// Ed "C++" Love
// Do not ask Prashant about this code.  He has no clue. 
//
// Constant Presssure/Volume Four Node Quadrilateral
// Plane Strain (NOT PLANE STRESS)
//

#include <iostream.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ConstantPressureVolumeQuad.h>
#include <Renderer.h>

//static data
Matrix ConstantPressureVolumeQuad :: stiff(8,8)   ;
Vector ConstantPressureVolumeQuad :: resid(8)     ;
Matrix ConstantPressureVolumeQuad :: mass(8,8)    ;
Matrix ConstantPressureVolumeQuad :: damping(8,8) ;
 
//volume-pressure constants
double ConstantPressureVolumeQuad :: one3  = 1.0 / 3.0 ;
double ConstantPressureVolumeQuad :: two3  = 2.0 / 3.0 ;
double ConstantPressureVolumeQuad :: four3 = 4.0 / 3.0 ;
double ConstantPressureVolumeQuad :: one9  = 1.0 / 9.0 ;
    
//quadrature data
double ConstantPressureVolumeQuad :: root3 = sqrt(3.0) ;
double ConstantPressureVolumeQuad :: one_over_root3 = 1.0 / root3 ;

double ConstantPressureVolumeQuad :: sg[] = { -one_over_root3,  
					       one_over_root3, 
					       one_over_root3, 
					      -one_over_root3 } ;

double ConstantPressureVolumeQuad :: tg[] = { -one_over_root3, 
                                              -one_over_root3, 
                                               one_over_root3,  
                                               one_over_root3 } ;

double ConstantPressureVolumeQuad :: wg[] = { 1.0, 1.0, 1.0, 1.0 } ;
  

//null constructor
ConstantPressureVolumeQuad :: ConstantPressureVolumeQuad( ) :
Element( 0, ELE_TAG_ConstantPressureVolumeQuad ),
connectedExternalNodes(4) 
{ 
  return ; 
}


//full constructor
ConstantPressureVolumeQuad :: ConstantPressureVolumeQuad( 
                            int tag, 
                  	    int node1,
			    int node2,
			    int node3,
			    int node4,
			    NDMaterial &theMaterial ) :
Element( tag, ELE_TAG_ConstantPressureVolumeQuad ),
connectedExternalNodes(4) 
{
  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

      materialPointers[i] = theMaterial.getCopy("AxiSymmetric2D") ;

      if (materialPointers[i] == 0) {

	  g3ErrorHandler->fatal("ConstantPressureVolumeQuad::constructor %s",
				"- failed to get a material of type: AxiSymmetric2D");
      } //end if
      
  } //end for i 

}


//destructor 
ConstantPressureVolumeQuad :: ~ConstantPressureVolumeQuad( )
{
  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ; 

    nodePointers[i] = 0 ;
  } //end for i

}


//set domain
void ConstantPressureVolumeQuad :: setDomain( Domain *theDomain ) 
{  
  int i ;
  for ( i = 0;  i < 4; i++ ) {

     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i)  ) ;

     if ( nodePointers[i] != 0 ) {
	 const Vector &coor = nodePointers[i]->getCrds( ) ;

	 xl[0][i] = coor(0) ;
	 xl[1][i] = coor(1) ; 
     } // end if 

  } //end for i 
  
  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int ConstantPressureVolumeQuad :: getNumExternalNodes( ) const
{
  return 4 ;
} 
 

//return connected external nodes
const ID& ConstantPressureVolumeQuad :: getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 


//return number of dofs
int ConstantPressureVolumeQuad :: getNumDOF( ) 
{
  return 8 ;
}


//commit state
int ConstantPressureVolumeQuad :: commitState( )
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int ConstantPressureVolumeQuad :: revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int ConstantPressureVolumeQuad :: revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void ConstantPressureVolumeQuad :: Print( ostream &s, int flag )
{
  s << '\n' ;
  s << "Four Node Quad -- Mixed Pressure/Volume -- Plane Strain \n" ;
  s << "Element Number " << this->getTag() << '\n' ;
  s << "Node 1 : " << connectedExternalNodes(0) << '\n' ;
  s << "Node 2 : " << connectedExternalNodes(1) << '\n' ;
  s << "Node 3 : " << connectedExternalNodes(2) << '\n' ;
  s << "Node 4 : " << connectedExternalNodes(3) << '\n' ;
  s << "Material Information : \n " ;

  materialPointers[0]->Print( s, flag ) ;

  s << endl ;
}

//return stiffness matrix 
const Matrix& ConstantPressureVolumeQuad :: getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    


//return secant matrix 
const Matrix& ConstantPressureVolumeQuad :: getSecantStiff( ) 
{
   int tang_flag = 1 ; //get the tangent

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  
    
  return stiff ;
}
    

//return damping matrix because frank is a dumb ass 
const Matrix& ConstantPressureVolumeQuad :: getDamp( ) 
{
  //not supported
  return damping ;
}    


//return mass matrix
const Matrix& ConstantPressureVolumeQuad :: getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 


//zero the load -- what load?
void ConstantPressureVolumeQuad :: zeroLoad( )
{
  return ;
}

//add load -- what load?
int  ConstantPressureVolumeQuad :: addLoad( const Vector &addP )
{
  return -1 ;
}


//get residual
const Vector& ConstantPressureVolumeQuad :: getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  return resid ;   
}


//get residual with inertia terms
const Vector& ConstantPressureVolumeQuad :: getResistingForceIncInertia( )
{
  int tang_flag = 0 ; //don't get the tangent

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  //inertia terms
  formInertiaTerms( tang_flag ) ;

  return resid ;
}

//*****************************************************************************
//form inertia terms

void   ConstantPressureVolumeQuad::formInertiaTerms( int tangFlag ) 
{

  static const int ndm = 2 ;

  static const int ndf = 2 ; 

  static const int numberNodes = 4 ;

  static const int numberGauss = 4 ;

  static const int nShape = 3 ;

  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix 

  double dvol ; //volume element

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static Vector momentum(ndf) ;

  static Matrix sx(ndm,ndm) ;

  int i, j, k, p ;
  int jj, kk ;

  double temp, rho, massJK ;


  //zero mass 
  mass.Zero( ) ;

  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj, sx ) ;
    
    //volume element
    dvol = wg[i] * xsj ;

    //node loop to compute acceleration
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
      momentum += shp[massIndex][j] * ( nodePointers[j]->getTrialAccel()  ) ; 

    //density
    rho = materialPointers[i]->getRho() ;

    //multiply acceleration by density to form momentum
    momentum *= rho ;


    //residual and tangent calculations node loops
    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      temp = shp[massIndex][j] * dvol ;

      for ( p = 0; p < ndf; p++ )
        resid( jj+p ) += ( temp * momentum(p) )  ;

      
      if ( tangFlag == 1 ) {

	 //multiply by density
	 temp *= rho ;

	 //node-node mass
         kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) {

	    massJK = temp * shp[massIndex][k] ;

            for ( p = 0; p < ndf; p++ )  
	      mass( jj+p, kk+p ) += massJK ;
            
            kk += ndf ;
          } // end for k loop

      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop


  } //end for i gauss loop 



}

//*********************************************************************

//form residual and tangent
void ConstantPressureVolumeQuad ::  formResidAndTangent( int tang_flag ) 
{
  // strains ordered  00, 11, 22, 01  
  //            i.e.  11, 22, 33, 12 
  //
  //            strain(0) =   eps_00
  //            strain(1) =   eps_11
  //            strain(2) =   eps_22
  //            strain(3) = 2*eps_01
  //
  //  same ordering for stresses but no 2 

  int i,  j,  k, l ;
  int jj, kk ;
  int node ;
  int success ;
  
  double tmp_shp[3][4] ; //shape functions

  double     shp[3][4][4] ; //shape functions at each gauss point

  double vol_avg_shp[3][4] ; // volume averaged shape functions

  double xsj ;  // determinant jacaobian matrix 

  static Matrix sx(2,2) ; // inverse jacobian matrix 

  double dvol[4] ; //volume elements

  double volume = 0.0 ; //volume of element

  double theta = 0.0 ; //average volume change (trace of strain) 

  double pressure = 0.0 ; //constitutive pressure  

  static Vector strain(4) ; //strain in vector form 

  static Vector sigBar(4) ; //stress in vector form
  static Vector sig(4) ; //mixed stress in vector form

  double trace = 0.0 ; //trace of the strain 

  static Matrix BJ(4,2) ; //strain-displacement matrices
  static Matrix BJtran(2,4) ; 
  static Matrix BK(4,2) ;

  static Matrix littleBJ(1,2) ; //volume strain-displacement matrices
  static Matrix littleBJtran(2,1) ;
  static Matrix littleBK(1,2) ; 

  static Matrix stiffJK(2,2) ; //nodeJ-nodeK 2x2 stiffness
  static Vector residJ(2) ; //nodeJ residual 
  
  static Vector one(4) ; //rank 2 identity as a vector
  static Matrix oneMatrix(4,1) ;
  static Matrix oneTran(1,4) ;
  
  static Matrix Pdev(4,4) ; //deviator projector

  static Matrix dd(4,4) ;  //material tangent

  static Matrix Pdev_dd_Pdev(4,4) ;
  static Matrix Pdev_dd_one(4,1) ; 
  static Matrix one_dd_Pdev(1,4) ;
  double bulk ;
  static Matrix BJtranD(2,4) ;
  static Matrix BJtranDone(2,1) ;
  static Matrix littleBJoneD(2,4) ;
  static Matrix littleBJtranBulk(2,1) ;
  
  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;

  //one vector
  one(0) = 1.0 ;
  one(1) = 1.0 ;
  one(2) = 1.0 ;
  one(3) = 0.0 ;

  //one matrix
  oneMatrix(0,0) = 1.0 ;
  oneMatrix(1,0) = 1.0 ;
  oneMatrix(2,0) = 1.0 ;
  oneMatrix(3,0) = 0.0 ;
  oneTran = this->transpose( 4, 1, oneMatrix ) ;

  //Pdev matrix
  Pdev.Zero( ) ;

  Pdev(0,0) =  two3 ;
  Pdev(0,1) = -one3 ;
  Pdev(0,2) = -one3 ;

  Pdev(1,0) = -one3 ;
  Pdev(1,1) =  two3 ;
  Pdev(1,2) = -one3 ;

  Pdev(2,0) = -one3 ;
  Pdev(2,1) = -one3 ;
  Pdev(2,2) =  two3 ;

  Pdev(3,3) = 1.0 ;


  //zero stuff
  volume = 0.0 ;

  for ( k = 0; k < 3; k++ ){
    for ( l = 0; l < 4; l++ ) 
        vol_avg_shp[k][l] = 0.0 ; 
  } //end for k


  //gauss loop to compute volume averaged shape functions

  for ( i = 0; i < 4; i++ ){
    
    shape2d( sg[i], tg[i], xl, tmp_shp, xsj, sx ) ;

    dvol[i] = wg[i] * xsj ;  // multiply by radius for axisymmetry 

    volume += dvol[i] ;

    for ( k = 0; k < 3; k++ ){
      for ( l = 0; l < 4; l++ ) {

	shp[k][l][i] = tmp_shp[k][l] ;

        vol_avg_shp[k][l] += tmp_shp[k][l] * dvol[i] ;

      } // end for l
    } //end for k

  } //end for i 


  //compute volume averaged shape functions
  for ( k = 0; k < 3; k++ ){
    for ( l = 0; l < 4; l++ ) 
        vol_avg_shp[k][l] /= volume ; 
  } //end for k


  //compute theta
  theta = 0.0 ;
  for ( i = 0; i < 4; i++ ) {

    strain.Zero( ) ;

    //node loop to compute strain
    for ( node = 0; node < 4; node++ ) {

      const Vector &ul = nodePointers[node]->getTrialDisp( ) ;

      strain(0) += shp[0][node][i] * ul(0) ;

      strain(1) += shp[1][node][i] * ul(1) ;

      strain(2) = 0.0 ;  // not zero for axisymmetry


    } // end for node

    trace  =  strain(0) + strain(1) + strain(2) ;

    theta +=  trace * dvol[i] ;

  } // end for i
  theta /= volume ;


  //compute pressure
  pressure = 0.0 ;
  for ( i = 0; i < 4; i++ ) {

    strain.Zero( ) ;

    //node loop to compute strain
    for ( node = 0; node < 4; node++ ) {

      const Vector &ul = nodePointers[node]->getTrialDisp( ) ;

      strain(0) += shp[0][node][i] * ul(0) ;

      strain(1) += shp[1][node][i] * ul(1) ;

      strain(2) = 0.0 ; // not zero for axisymmetry

      strain(3) +=  shp[1][node][i] * ul(0) 	
	          + shp[0][node][i] * ul(1) ; 

    } // end for node

    trace = strain(0) + strain(1) + strain(2) ;

    strain -= (one3*trace)*one ;
    strain += (one3*theta)*one ;
    
    success = materialPointers[i]->setTrialStrain( strain ) ;

    const Vector &sigBar = materialPointers[i]->getStress( ) ;

    pressure +=  one3 * ( sigBar(0) + sigBar(1) + sigBar(2) ) * dvol[i] ;

  } // end for i
  pressure /= volume ;


  //residual and tangent calculations gauss loop

  for ( i = 0; i < 4; i++ ) {

    const Vector &sigBar = materialPointers[i]->getStress( ) ;

    //stress for equilibrium

    trace = sigBar(0) + sigBar(1) + sigBar(2) ;

    sig  = sigBar ;
    sig -= (one3*trace)*one ;
    sig += pressure*one ;

    //multilply by volume elements and compute 
    // matrices for stiffness calculation

    sig *= dvol[i] ;

    if ( tang_flag == 1 ) {

	const Matrix &ddBar  = materialPointers[i]->getTangent( ) ;
	
	dd = ddBar * dvol[i] ;
        
	Pdev_dd_Pdev = Pdev * dd * Pdev ;

	Pdev_dd_one  = one3 * ( Pdev * dd * oneMatrix ) ;
	
	one_dd_Pdev  = one3 * ( oneTran * dd * Pdev ) ;

	bulk = one9 * (   dd(0,0) + dd(0,1) + dd(0,2)  
	                + dd(1,0) + dd(1,1) + dd(1,2) 
	                + dd(2,0) + dd(2,1) + dd(2,2) ) ;

    } //end if tang_flag
    

    //residual and tangent loop over nodes
    
    jj = 0 ;
    for ( j = 0; j < 4; j++ ) {
      
      BJ.Zero( );
      BJ(0,0) = shp[0][j][i] ;
      BJ(1,1) = shp[1][j][i] ; 

      // BJ(2,0) for axi-symmetry 

      BJ(3,0) = shp[1][j][i]  ;
      BJ(3,1) = shp[0][j][i]  ;


      littleBJ(0,0) = vol_avg_shp[0][j] ;
      littleBJ(0,1) = vol_avg_shp[1][j] ;

      BJtran = this->transpose( 4, 2, BJ ) ;

      littleBJtran = this->transpose( 1, 2, littleBJ ) ;

      //compute residual 

      residJ = BJtran * sig ; 

      resid( jj   ) += residJ(0) ;
      resid( jj+1 ) += residJ(1) ;

      if ( tang_flag == 1 ) { //stiffness matrix

	  BJtranD          =  BJtran * Pdev_dd_Pdev ;
	  BJtranDone       =  BJtran * Pdev_dd_one ;
	  littleBJoneD     =  littleBJtran * one_dd_Pdev ;
	  littleBJtranBulk =  bulk * littleBJtran ;
	  
 	  kk = 0 ;
	  for ( k = 0; k < 4; k++ ) {

	      BK.Zero( );
	      BK(0,0) = shp[0][k][i] ;
	      BK(1,1) = shp[1][k][i] ;

	      // BK(2,0) for axi-symmetry 

	      BK(3,0) = shp[1][k][i] ;
	      BK(3,1) = shp[0][k][i] ;


	      littleBK(0,0) = vol_avg_shp[0][k] ;
	      littleBK(0,1) = vol_avg_shp[1][k] ;

	      //compute stiffness matrix
        
	      stiffJK =  ( BJtranD + littleBJoneD ) * BK
                      +  ( BJtranDone + littleBJtranBulk ) * littleBK ; 

	      stiff( jj,   kk   ) += stiffJK(0,0) ;
	      stiff( jj+1, kk   ) += stiffJK(1,0) ;
	      stiff( jj,   kk+1 ) += stiffJK(0,1) ;
	      stiff( jj+1, kk+1 ) += stiffJK(1,1) ;

	      kk += 2 ;
	  } // end for k

      } // end if tang_flag 

      
      jj += 2 ;
    } // end for j 

  } //end for i

  
  return ;
}


//----------------------------------------------------------------------

Matrix ConstantPressureVolumeQuad :: transpose( int dim1, 
						int dim2, 
						const Matrix &M ) 
{
  int i ;
  int j ;

  Matrix Mtran( dim2, dim1 ) ;

  for ( i = 0; i < dim1; i++ ) {

    for ( j = 0; j < dim2; j++ ) 
      Mtran(j,i) = M(i,j) ;

  } // end for i

  return Mtran ;
}

//---------------------------------------------------------------------

//shape function routine for four node quads
void ConstantPressureVolumeQuad :: shape2d( double ss, double tt, 
		                            const double x[2][4], 
		                            double shp[3][4], 
		                            double &xsj, 
		                            Matrix &sx ) 
{ 

  int i, j, k ;

  double temp ;
     
  static const double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static const double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static Matrix xs(2,2) ;

  for ( i = 0; i < 4; i++ ) {
      shp[2][i] = ( 0.5 + s[i]*ss )*( 0.5 + t[i]*tt ) ;
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  } // end for i

  
  // Construct jacobian and its inverse
  
  xs.Zero( ) ;
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {

      for ( k = 0; k < 4; k++ )
	  xs(i,j) +=  x[i][k] * shp[j][k] ;

    } //end for j
  }  // end for i 

  xsj = xs(0,0)*xs(1,1) - xs(0,1)*xs(1,0) ;

  sx(0,0) =  xs(1,1) / xsj ;
  sx(1,1) =  xs(0,0) / xsj ;
  sx(0,1) = -xs(0,1) / xsj ; 
  sx(1,0) = -xs(1,0) / xsj ; 


  //form global derivatives 
  
  for ( i = 0; i < 4; i++ ) {
    temp      = shp[0][i]*sx(0,0) + shp[1][i]*sx(1,0) ;
    shp[1][i] = shp[0][i]*sx(0,1) + shp[1][i]*sx(1,1) ;
    shp[0][i] = temp ;
  } // end for i

  return ;
}
	
//***********************************************************************

int
ConstantPressureVolumeQuad::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = nodePointers[0]->getCrds();
    const Vector &end2Crd = nodePointers[1]->getCrds();	
    const Vector &end3Crd = nodePointers[2]->getCrds();	
    const Vector &end4Crd = nodePointers[3]->getCrds();	

    const Vector &end1Disp = nodePointers[0]->getDisp();
    const Vector &end2Disp = nodePointers[1]->getDisp();
    const Vector &end3Disp = nodePointers[2]->getDisp();
    const Vector &end4Disp = nodePointers[3]->getDisp();

    static Matrix coords(4,3) ;
    static Vector values(4) ;
    static Vector P(8) ;

    coords.Zero( ) ;

    values(0) = 1 ;
    values(1) = 1 ;
    values(2) = 1 ;
    values(3) = 1 ;

    if (displayMode < 3 && displayMode > 0)
      P = this->getResistingForce();

    for (int i = 0; i < 2; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
      coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
      coords(3,i) = end4Crd(i) + end4Disp(i)*fact;    
      /*      if (displayMode < 3 && displayMode > 0)
	values(i) = P(displayMode*2+i);
      else
      values(i) = 1;  */
    }

    //cerr << coords;
    int error = 0;

    error += theViewer.drawPolygon (coords, values);

    return error;
}
   
//-----------------------------------------------------------------------

int ConstantPressureVolumeQuad :: sendSelf (int commitTag, Channel &theChannel)
{
    return -1;
}
    
int ConstantPressureVolumeQuad :: recvSelf (int commitTag, 
					Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
    return -1;
}

