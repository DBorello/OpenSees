/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.3 $
// $Date: 2001-07-11 21:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/ShellMITC4.cpp,v $

// Ed "C++" Love
//
// B-bar four node shell element with membrane and drill
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
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ShellMITC4.h>
#include <R3vectors.h>
#include <Renderer.h>

#define min(a,b) ( (a)<(b) ? (a):(b) )

//static data
double  ShellMITC4::xl[2][4] ;
Matrix  ShellMITC4::stiff(24,24) ;
Vector  ShellMITC4::resid(24) ;
Matrix  ShellMITC4::mass(24,24) ;
Matrix  ShellMITC4::damping(24,24) ;

//static basis vectors
Vector  ShellMITC4::g1(3) ;
Vector  ShellMITC4::g2(3) ;
Vector  ShellMITC4::g3(3) ;

//intialize pointers to zero using intialization list 
Matrix  **ShellMITC4::GammaB1pointer = 0 ;
Matrix  **ShellMITC4::GammaD1pointer = 0 ;
Matrix  **ShellMITC4::GammaA2pointer = 0 ;
Matrix  **ShellMITC4::GammaC2pointer = 0 ; 
    
//quadrature data
const double  ShellMITC4::root3 = sqrt(3.0) ;
const double  ShellMITC4::one_over_root3 = 1.0 / root3 ;

const double  ShellMITC4::sg[] = { -one_over_root3,  
				one_over_root3, 
				one_over_root3, 
	                       -one_over_root3 } ;

const double  ShellMITC4::tg[] = { -one_over_root3, 
                               -one_over_root3, 
                                one_over_root3,  
                                one_over_root3 } ;

const double  ShellMITC4::wg[] = { 1.0, 1.0, 1.0, 1.0 } ;

  

//null constructor
ShellMITC4::ShellMITC4( ) :
Element( 0, ELE_TAG_ShellMITC4 ),
connectedExternalNodes(4) 
{ 

}


//*********************************************************************
//full constructor
ShellMITC4::ShellMITC4(  int tag, 
                         int node1,
                         int node2,
   	                 int node3,
                         int node4,
	                 SectionForceDeformation &theMaterial ) :
Element( tag, ELE_TAG_ShellMITC4 ),
connectedExternalNodes(4) 
{
  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

      materialPointers[i] = theMaterial.getCopy( ) ;

      if (materialPointers[i] == 0) {

	  g3ErrorHandler->fatal("ShellMITC4::constructor %s",
		"- failed to get a material of type: ShellSection");
      } //end if
      
  } //end for i 


  //shear matrix pointers

  if ( GammaB1pointer == 0 ) {
	GammaB1pointer = new Matrix*[4] ;    //four matrix pointers
	GammaB1pointer[0] = new Matrix(1,3) ;  //
	GammaB1pointer[1] = new Matrix(1,3) ;  //    four
	GammaB1pointer[2] = new Matrix(1,3) ;  //  1x3 matrices
	GammaB1pointer[3] = new Matrix(1,3) ;  //
  } //end if B1

  if ( GammaD1pointer == 0 ) {
	GammaD1pointer = new Matrix*[4] ;
	GammaD1pointer[0] = new Matrix(1,3) ;
	GammaD1pointer[1] = new Matrix(1,3) ;
	GammaD1pointer[2] = new Matrix(1,3) ;
	GammaD1pointer[3] = new Matrix(1,3) ;
  } //end if D1

  if ( GammaA2pointer == 0 ) {
	GammaA2pointer = new Matrix*[4] ;
	GammaA2pointer[0] = new Matrix(1,3) ;
	GammaA2pointer[1] = new Matrix(1,3) ;
	GammaA2pointer[2] = new Matrix(1,3) ;
	GammaA2pointer[3] = new Matrix(1,3) ;
  } //end if A2

  if ( GammaC2pointer == 0 ) {
	GammaC2pointer = new Matrix*[4] ;
	GammaC2pointer[0] = new Matrix(1,3) ;
	GammaC2pointer[1] = new Matrix(1,3) ;
	GammaC2pointer[2] = new Matrix(1,3) ;
	GammaC2pointer[3] = new Matrix(1,3) ;
  } //end if C2

}
//******************************************************************


//destructor 
ShellMITC4::~ShellMITC4( )
{
  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ; 

    nodePointers[i] = 0 ;

  } //end for i

}


//set domain
void  ShellMITC4::setDomain( Domain *theDomain ) 
{  

  static Matrix ddMembrane(3,3) ;

  static Vector eig(3) ;

  int i, j ;

  //node pointers
  for ( i = 0; i < 4; i++ ) 
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;


  //compute drilling stiffness penalty parameter

  static Vector strain(8) ;
  strain.Zero( ) ;

  //send the (zero) strain to the material 
  int success = materialPointers[0]->setTrialSectionDeformation( strain ) ;

  //compute the stress (should be zero)
  const Vector &stress = materialPointers[0]->getStressResultant( ) ;

  //compute the material tangent
  const Matrix &dd = materialPointers[0]->getSectionTangent( ) ;

  //assemble ddMembrane ;
  for ( i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j ++ )
         ddMembrane(i,j) = dd(i,j) ;
  } //end for i 
  
  //eigenvalues of ddMembrane
  eig = LovelyEig( ddMembrane ) ;

  //set ktt 
  //Ktt = dd(2,2) ;  //shear modulus 
  Ktt = min( eig(2), min( eig(0), eig(1) ) ) ;


  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int  ShellMITC4::getNumExternalNodes( ) const
{
  return 4 ;
} 
 

//return connected external nodes
const ID&  ShellMITC4::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 


//return number of dofs
int  ShellMITC4::getNumDOF( ) 
{
  return 24 ;
}


//commit state
int  ShellMITC4::commitState( )
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int  ShellMITC4::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  ShellMITC4::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void  ShellMITC4::Print( ostream &s, int flag )
{
  s << '\n' ;
  s << "MITC4 Bbar Non-Locking Four Node Shell \n" ;
  s << "Node 1 : " << connectedExternalNodes(0) << '\n' ;
  s << "Node 2 : " << connectedExternalNodes(1) << '\n' ;
  s << "Node 3 : " << connectedExternalNodes(2) << '\n' ;
  s << "Node 4 : " << connectedExternalNodes(3) << '\n' ;

  s << "Material Information : \n " ;
  materialPointers[0]->Print( s, flag ) ;

  s << endl ;
}

//return stiffness matrix 
const Matrix&  ShellMITC4::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    

//return secant matrix 
const Matrix&  ShellMITC4::getSecantStiff( ) 
{
   int tang_flag = 1 ; //get the tangent

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  
    
  return stiff ;
}
    

//return damping matrix because frank is a dumb ass 
const Matrix&  ShellMITC4::getDamp( ) 
{
  //not supported
  return damping ;
}    


//return mass matrix
const Matrix&  ShellMITC4::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 


//zero the load -- what load?
void  ShellMITC4::zeroLoad( )
{
  return ;
}

//add load -- what load?
int  ShellMITC4::addLoad( const Vector &addP )
{
  return -1 ;
}


//get residual
const Vector&  ShellMITC4::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  return resid ;   
}


//get residual with inertia terms
const Vector&  ShellMITC4::getResistingForceIncInertia( )
{
  int tang_flag = 0 ; //don't get the tangent

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  formInertiaTerms( tang_flag ) ;

  return resid ;
}


//*********************************************************************
//form inertia terms

void   ShellMITC4::formInertiaTerms( int tangFlag ) 
{

  //translational mass only
  //rotational inertia terms are neglected

  static const int ndm = 3 ;

  static const int ndf = 6 ; 

  static const int numberNodes = 4 ;

  static const int numberGauss = 4 ;

  static const int nShape = 3 ;

  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix 

  double dvol ; //volume element

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static Vector momentum(ndf) ;


  int i, j, k, p, q ;
  int jj, kk ;

  double temp, rhoH, massJK ;


  //zero mass 
  mass.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;


  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;

    //volume element to also be saved
    dvol = wg[i] * xsj ;  

    //node loop to compute accelerations
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
      momentum += ( shp[massIndex][j] * nodePointers[j]->getTrialAccel() ) ;


    //density
    rhoH = materialPointers[i]->getRho() ;

    //multiply acceleration by density to form momentum
    momentum *= rhoH ;


    //residual and tangent calculations node loops
    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      temp = shp[massIndex][j] * dvol ;

      for ( p = 0; p < 3; p++ )
        resid( jj+p ) += ( temp * momentum(p) ) ;

      
      if ( tangFlag == 1 ) {

	 //multiply by density
	 temp *= rhoH ;

	 //node-node translational mass
         kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) {

	   massJK = temp * shp[massIndex][k] ;

	   for ( p = 0; p < 3; p++ ) 
	      mass( jj+p, kk+p ) +=  massJK ;
            
            kk += ndf ;
          } // end for k loop

      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop


  } //end for i gauss loop 

}

//*********************************************************************

//form residual and tangent
void  ShellMITC4::formResidAndTangent( int tang_flag ) 
{
  //
  //  six(6) nodal dof's ordered :
  //
  //    -        - 
  //   |    u1    |   <---membrane
  //   |    u2    |
  //   |----------|
  //   |  w = u3  |   <---bending
  //   |  theta1  | 
  //   |  theta2  | 
  //   |----------|
  //   |  theta3  |   <---drill 
  //    -        -  
  //
  // membrane strains ordered :
  //
  //            strain(0) =   eps00     i.e.   (11)-strain
  //            strain(1) =   eps11     i.e.   (22)-strain
  //            strain(2) =   gamma01   i.e.   (12)-shear
  //
  // curvatures and shear strains ordered  :
  //
  //            strain(3) =     kappa00  i.e.   (11)-curvature
  //            strain(4) =     kappa11  i.e.   (22)-curvature
  //            strain(5) =   2*kappa01  i.e. 2*(12)-curvature 
  //
  //            strain(6) =     gamma02  i.e.   (13)-shear
  //            strain(7) =     gamma12  i.e.   (23)-shear
  //
  //  same ordering for moments/shears but no 2 
  //  
  //  Then, 
  //              epsilon00 = -z * kappa00      +    eps00_membrane
  //              epsilon11 = -z * kappa11      +    eps11_membrane
  //  gamma01 = 2*epsilon01 = -z * (2*kappa01)  +  gamma01_membrane 
  //
  //  Shear strains gamma02, gamma12 constant through cross section
  //

  static const int ndf = 6 ; //two membrane plus three bending plus one drill

  int i,  j,  k, p, q ;
  int jj, kk ;
  int node ;

  int success ;
  
  double volume = 0.0 ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[4] ; //volume element

  static Vector strain(8) ;  //strain

  static double shp[3][4] ;  //shape functions at a gauss point

  static double Shape[3][4][4] ; //all the shape functions

  static Vector residJ(ndf) ; //nodeJ residual 

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Vector stress(8) ;  //stress resultants

  static Matrix dd(8,8) ;  //material tangent

  static Matrix J0(2,2) ;  //Jacobian at center
 
  static Matrix J0inv(2,2) ; //inverse of Jacobian at center

  static double epsDrill ;  //drilling "strain"

  static double tauDrill ; //drilling "stress"

  //---------B-matrices------------------------------------

    static Matrix Bhat[4] ;      // tying integral 

    static Matrix BJ(8,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,8) ;

    static Matrix BK(8,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,8) ;


    static Matrix Bbend(3,3) ;  // bending B matrix

    static Matrix Bshear(2,3) ; // shear B matrix

    static Matrix Bmembrane(3,2) ; // membrane B matrix


    static Matrix BdrillJ(1,ndf) ; //drill B matrix

    static Matrix BdrillK(1,ndf) ; 

  //-------------------------------------------------------

  
  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;
  computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  //compute the gamma's
  computeGamma( xl, J0 ) ;

  //zero Bhat = \frac{1}{volume} \int{  B - \bar{B} } \diff A
  for ( node = 0;  node < 4;  node++ ) {
    Bhat[node].resize(2,3);
    Bhat[node].Zero( ) ;
  }

  //gauss loop to compute Bhat's 
  for ( i = 0; i < 4; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;

    //save shape functions
    for ( p = 0; p < 3; p++ ) {
       for ( q = 0; q < 4; q++ )
	  Shape[p][q][i] = shp[p][q] ;
    } // end for p

    //volume element to also be saved
    dvol[i] = wg[i] * xsj ;  

    volume += dvol[i] ;

    for ( node = 0; node < 4; node++ ) {

       //compute B shear matrix for this node
       Bhat[node] += (  dvol[i] * computeBshear(node,shp)  ) ;

       //compute B-bar shear matrix for this node
       Bhat[node] -= ( dvol[i] *
                  computeBbarShear( node, sg[i], tg[i], J0inv ) 
                      ) ;

    } //end for node   

  } // end for i gauss loop

  //compute Bhat 
    for ( node = 0;  node < 4;  node++ ) {
      Bhat[node] /= volume ;
      //      cerr << "Bhat[" << node << "]  = " << Bhat[node] ;
    }
    

  //gauss loop 
  for ( i = 0; i < 4; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < 3; p++ ) {
       for ( q = 0; q < 4; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p


    //zero the strains
    strain.Zero( ) ;
    epsDrill = 0.0 ;


    // j-node loop to compute strain 
    for ( j = 0; j < 4; j++ )  {

      //compute B matrix 

      Bmembrane = computeBmembrane( j, shp ) ;

      Bbend = computeBbend( j, shp ) ;

      Bshear = computeBbarShear( j, sg[i], tg[i], J0inv ) ;

      //add Bhat to shear terms 
      Bshear += Bhat[j] ;

      BJ = assembleB( Bmembrane, Bbend, Bshear ) ;
      
      //nodal "displacements" 
      const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

      //compute the strain
      strain += (BJ*ul) ; 


      //drilling B matrix
      BdrillJ = computeBdrill( j, shp ) ;

      //drilling "strain" 
      for ( p = 0; p < ndf; p++ )
	  epsDrill +=  BdrillJ(0,p)*ul(p) ;
 
    } // end for j
  

    //send the strain to the material 
    success = materialPointers[i]->setTrialSectionDeformation( strain ) ;

    //compute the stress
    stress = materialPointers[i]->getStressResultant( ) ;

    //drilling "stress" 
    tauDrill = Ktt * epsDrill ;

    //multiply by volume element
    stress   *= dvol[i] ;
    tauDrill *= dvol[i] ;

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getSectionTangent( ) ;
      dd *= dvol[i] ;
    } //end if tang_flag


    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j = 0; j < 4; j++ ) {

      Bmembrane = computeBmembrane( j, shp ) ;

      Bbend = computeBbend( j, shp ) ;

      //multiply bending terms by (-1.0) for correct statement
      // of equilibrium  
      Bbend *= (-1.0) ;

      Bshear = computeBbarShear( j, sg[i], tg[i], J0inv ) ;

      //add Bhat to shear terms
      Bshear += Bhat[j] ;

      BJ = assembleB( Bmembrane, Bbend, Bshear ) ;
   

      //transpose 
      BJtran = transpose( 8, ndf, BJ ) ;

      //residual
      residJ = BJtran * stress ;


      //drilling B matrix
      BdrillJ = computeBdrill( j, shp ) ;

      //residual 
      for ( p = 0; p < ndf; p++ )
        resid( jj + p ) += ( residJ(p) + BdrillJ(0,p)*tauDrill ) ;


      if ( tang_flag == 1 ) {

        BJtranD = BJtran * dd ;

        BdrillJ *= ( Ktt*dvol[i] ) ;

        kk = 0 ;
        for ( k = 0; k < 4; k++ ) {

          Bmembrane = computeBmembrane( k, shp ) ;

          Bbend = computeBbend( k, shp ) ;

          Bshear = computeBbarShear( k, sg[i], tg[i], J0inv ) ;

          Bshear += Bhat[k] ;


          BK = assembleB( Bmembrane, Bbend, Bshear ) ;
  
          BdrillK = computeBdrill( k, shp ) ;

 
          stiffJK =  BJtranD * BK  
                  +  transpose( 1,ndf,BdrillJ ) * BdrillK ; 

          for ( p = 0; p < ndf; p++ )  {
             for ( q = 0; q < ndf; q++ )
                stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
          } //end for p

          kk += ndf ;
        } // end for k loop

      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop


  } //end for i gauss loop 

  
  return ;
}


//************************************************************************
//compute local coordinates and basis

void   ShellMITC4::computeBasis( ) 
{
  //could compute derivatives \frac{ \partial {\bf x} }{ \partial L_1 } 
  //                     and  \frac{ \partial {\bf x} }{ \partial L_2 }
  //and use those as basis vectors but this is easier 
  //and the shell is flat anyway.

  
  //get two vectors (g1, g2) in plane of shell by 
  // nodal coordinate differences
  
  const Vector &coor0 = nodePointers[0]->getCrds( ) ;

  const Vector &coor1 = nodePointers[1]->getCrds( ) ;

  const Vector &coor2 = nodePointers[2]->getCrds( ) ;
  
  const Vector &coor3 = nodePointers[3]->getCrds( ) ;

  g1.Zero( ) ;

  g1 = 0.5 * ( coor2 + coor1 - coor3 - coor0 ) ;

  g2.Zero( ) ;

  g2 = 0.5 * ( coor3 + coor2 - coor1 - coor0 ) ;


 
  //normalize g1 

  double length = LovelyNorm( g1 ) ;

  g1 /= length ;


  //Gram-Schmidt process for g2 

  double alpha = LovelyInnerProduct( g2, g1 ) ;

  g2 = ( g2 - alpha*g1 ) ;


  //normalize g2 

  length = LovelyNorm( g2 ) ;

  g2 /= length ;


  //cross product for g3  

  g3 = LovelyCrossProduct( g1, g2 ) ;

  
  //local nodal coordinates in plane of shell

  int i ;
  for ( i = 0; i < 4; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( ) ;

       xl[0][i] = LovelyInnerProduct( coorI, g1 ) ;  
       xl[1][i] = LovelyInnerProduct( coorI, g2 ) ;

  }  //end for i 

}

//*************************************************************************
//compute Bdrill

Matrix   ShellMITC4::computeBdrill( int node, const double shp[3][4] )
{

  //  static Matrix Bdrill(1,6) ;

  static double B1 ;
  static double B2 ;
  static double B6 ;

  static Matrix BdrillShell(1,6) ;

//---Bdrill Matrix in standard {1,2,3} mechanics notation---------
//
//             -                                       -
//   Bdrill = | -0.5*N,2   +0.5*N,1    0    0    0   -N |   (1x6) 
//             -                                       -  
//
//----------------------------------------------------------------

  //  Bdrill.Zero( ) ;

  //Bdrill(0,0) = -0.5*shp[1][node] ;

  //Bdrill(0,1) = +0.5*shp[0][node] ;

  //Bdrill(0,5) =     -shp[2][node] ;


  B1 =  -0.5*shp[1][node] ; 

  B2 =  +0.5*shp[0][node] ;

  B6 =  -shp[2][node] ;

  
  BdrillShell.Zero( ) ;

  BdrillShell(0,0) = B1*g1(0) + B2*g2(0) ;
  BdrillShell(0,1) = B1*g1(1) + B2*g2(1) ; 
  BdrillShell(0,2) = B1*g1(2) + B2*g2(2) ;

  BdrillShell(0,3) = B6*g3(0) ;
  BdrillShell(0,4) = B6*g3(1) ; 
  BdrillShell(0,5) = B6*g3(2) ;
 
  return BdrillShell ;

}

//*************************************************************************
//compute the gamma's

void   ShellMITC4::computeGamma( const double xl[2][4], const Matrix &J )
{
  //static Matrix e1(1,2) ;
  //static Matrix e2(1,2) ;

  static Matrix Jtran(2,2) ;  // J transpose

  static Matrix e1Jtran(1,2) ; // e1*Jtran 

  static Matrix e2Jtran(1,2) ; // e2*Jtran

  static Matrix Bshear(2,3) ; // shear strain-displacement matrix

  static double shape[3][4] ; // shape functions 

  double xsj ; 

  double L1, L2 ;

  int node ;

 /*  
  e1(0,0) = 1.0 ;
  e1(0,1) = 0.0 ;

  e2(0,0) = 0.0 ;
  e2(0,1) = 1.0 ;
 */

  //Jtran = transpose( 2, 2, J ) ;
  Jtran(0,0) = J(0,0) ;
  Jtran(1,1) = J(1,1) ;
  Jtran(0,1) = J(1,0) ;
  Jtran(1,0) = J(0,1) ;

  //first row of Jtran
  //e1Jtran = e1*Jtran ; 
  e1Jtran(0,0) = Jtran(0,0) ;
  e1Jtran(0,1) = Jtran(0,1) ;

  //second row of Jtran 
  //e2Jtran = e2*Jtran ;
  e2Jtran(0,0) = Jtran(1,0) ;
  e2Jtran(0,1) = Jtran(1,1) ;


  //-------------GammaB1--------------------------------

    //set natural coordinate point location
    L1 =  0.0 ;
    L2 = -1.0 ;

    //get shape functions    
    shape2d( L1, L2, xl, shape, xsj ) ;

    for ( node = 0; node < 4; node++ ) {

        //compute shear B        
        Bshear = computeBshear( node, shape ) ;

        GammaB1pointer[node]->Zero( ) ;

        ( *GammaB1pointer[node] ) = e1Jtran*Bshear ;

    } //end for node
 
  //----------------------------------------------------


  //-------------GammaD1--------------------------------

    //set natural coordinate point location
    L1 =  0.0 ;
    L2 = +1.0 ;

    //get shape functions    
    shape2d( L1, L2, xl, shape, xsj ) ;

    for ( node = 0; node < 4; node++ ) {

        //compute shear B        
        Bshear = computeBshear( node, shape ) ;

        GammaD1pointer[node]->Zero( ) ;

        ( *GammaD1pointer[node] ) = e1Jtran*Bshear ;

    } //end for node
 
  //----------------------------------------------------    


  //-------------GammaA2--------------------------------

    //set natural coordinate point location
    L1 = -1.0 ;
    L2 =  0.0 ;

    //get shape functions    
    shape2d( L1, L2, xl, shape, xsj ) ;

    for ( node = 0; node < 4; node++ ) {

        //compute shear B        
        Bshear = computeBshear( node, shape ) ;

        GammaA2pointer[node]->Zero( ) ;

        ( *GammaA2pointer[node] ) = e2Jtran*Bshear ;

    } //end for node
 
  //----------------------------------------------------     


  //-------------GammaC2--------------------------------

    //set natural coordinate point location
    L1 = +1.0 ;
    L2 =  0.0 ;

    //get shape functions    
    shape2d( L1, L2, xl, shape, xsj ) ;

    for ( node = 0; node < 4; node++ ) {

        //compute shear B        
        Bshear = computeBshear( node, shape ) ;

        GammaC2pointer[node]->Zero( ) ;

        ( *GammaC2pointer[node] ) = e2Jtran*Bshear ;

    } //end for node
 
  //----------------------------------------------------     

 return ;

}

//********************************************************************
//assemble a B matrix

Matrix  
ShellMITC4::assembleB( const Matrix &Bmembrane,
                               const Matrix &Bbend, 
                               const Matrix &Bshear ) 
{

  //Matrix Bbend(3,3) ;  // plate bending B matrix

  //Matrix Bshear(2,3) ; // plate shear B matrix

  //Matrix Bmembrane(3,2) ; // plate membrane B matrix


    static Matrix B(8,6) ;

    static Matrix BmembraneShell(3,3) ; 
    
    static Matrix BbendShell(3,6) ; 

    static Matrix BshearShell(2,6) ;
 
    static Matrix Gmem(2,3) ;

    static Matrix Gbend(3,6) ;

    int p, q ;
    int pp ;

//    
// For Shell : 
//---B Matrices in standard {1,2,3} mechanics notation---------
//
//            -                     _          
//           | Bmembrane  |     0    |
//           | --------------------- |     
//    B =    |         Bbend         |   (8x6) 
//           | --------------------- |
//           |         Bshear        |
//            -           -         -
//
//-------------------------------------------------------------
//


    //shell modified membrane terms
    
    Gmem.Zero( ) ;

    Gmem(0,0) = g1(0) ;
    Gmem(0,1) = g1(1) ;
    Gmem(0,2) = g1(2) ;

    Gmem(1,0) = g2(0) ;
    Gmem(1,1) = g2(1) ;
    Gmem(1,2) = g2(2) ;

    BmembraneShell = Bmembrane * Gmem ;


    //shell modified bending terms 

    Gbend.Zero( ) ;

    Gbend(0,0) = g3(0) ;
    Gbend(0,1) = g3(1) ;
    Gbend(0,2) = g3(2) ;

    Gbend(1,3) = g1(0) ;
    Gbend(1,4) = g1(1) ;
    Gbend(1,5) = g1(2) ;

    Gbend(2,3) = g2(0) ;
    Gbend(2,4) = g2(1) ;
    Gbend(2,5) = g2(2) ;

    BbendShell = Bbend * Gbend ;
   

    //shell modified shear terms 

    BshearShell = Bshear * Gbend ;

   

  B.Zero( ) ;

  //assemble B from sub-matrices 

  //membrane terms 
  for ( p = 0; p < 3; p++ ) {

      for ( q = 0; q < 3; q++ ) 
          B(p,q) = BmembraneShell(p,q) ;

  } //end for p


  //bending terms
  for ( p = 0; p < 3; p++ ) {
      pp = p + 3 ;
      
      for ( q = 0; q < 6; q++ ) {
          B(pp,q) = BbendShell(p,q) ; 
      } // end for q
 
  } //end for p
    

  //shear terms 
  for ( p = 0; p < 2; p++ ) {
      pp = p + 6 ;
      
      for ( q = 0; q < 6; q++ ) {
          B(pp,q) = BshearShell(p,q) ; 
      } // end for q
 
  } //end for p
  
  return B ;

}

//***********************************************************************
//compute Bmembrane matrix

Matrix   
ShellMITC4::computeBmembrane( int node, const double shp[3][4] ) 
{

  static Matrix Bmembrane(3,2) ;

//---Bmembrane Matrix in standard {1,2,3} mechanics notation---------
//
//                -             -
//               | +N,1      0   | 
// Bmembrane =   |   0     +N,2  |    (3x2)
//               | +N,2    +N,1  |
//                -             -  
//
//-------------------------------------------------------------------

  Bmembrane.Zero( ) ;

  Bmembrane(0,0) = shp[0][node] ;
  Bmembrane(1,1) = shp[1][node] ;
  Bmembrane(2,0) = shp[1][node] ;
  Bmembrane(2,1) = shp[0][node] ;

  return Bmembrane ;

}

//***********************************************************************
//compute Bbend matrix

Matrix   ShellMITC4::computeBbend( int node, const double shp[3][4] )
{

    static Matrix Bbend(3,3) ;

//---Bbend Matrix in standard {1,2,3} mechanics notation---------
//
//            -                -
//   Bbend = | 0      0    -N,1 | 
//           | 0    +N,2     0  |    (3x3)
//           | 0    +N,1   -N,2 |
//            -                -  
//
//----------------------------------------------------------------

    Bbend.Zero( ) ;

    Bbend(0,2) = -shp[0][node] ;
    Bbend(1,1) =  shp[1][node] ;
    Bbend(2,1) =  shp[0][node] ;
    Bbend(2,2) = -shp[1][node] ; 

    return Bbend ;
}

//***********************************************************************
//compute Bbar shear matrix

Matrix  ShellMITC4::computeBbarShear( int node, double L1, double L2,
			              const Matrix &Jinv ) 
{
    static Matrix Bshear(2,3) ;

    static Matrix JinvTran(2,2) ;  // J-inverse-transpose

    static Matrix Gamma1(1,3) ;
    static Matrix Gamma2(1,3) ; 


    //JinvTran = transpose( 2, 2, Jinv ) ;
    JinvTran(0,0) = Jinv(0,0) ;
    JinvTran(1,1) = Jinv(1,1) ;
    JinvTran(0,1) = Jinv(1,0) ;
    JinvTran(1,0) = Jinv(0,1) ;

    // compute BShear from Bbar interpolation

    Gamma1 = ( 1.0 - L2 )*( *GammaB1pointer[node] )  +
             ( 1.0 + L2 )*( *GammaD1pointer[node] )   ;

    Gamma2 = ( 1.0 - L1 )*( *GammaA2pointer[node] )  +
             ( 1.0 + L1 )*( *GammaC2pointer[node] )   ;
 
    Gamma1 *= 0.50 ;
    Gamma2 *= 0.50 ;


    Bshear.Zero( ) ;

    Bshear(0,0) = Gamma1(0,0) ;
    Bshear(0,1) = Gamma1(0,1) ;
    Bshear(0,2) = Gamma1(0,2) ;

    Bshear(1,0) = Gamma2(0,0) ;
    Bshear(1,1) = Gamma2(0,1) ;
    Bshear(1,2) = Gamma2(0,2) ;

    //strain tensor push on Bshear
    Bshear = ( JinvTran * Bshear ) ;
    
    return Bshear ;
}

//***********************************************************************
//compute standard Bshear matrix

Matrix  ShellMITC4::computeBshear( int node, const double shp[3][4] )
{

  static Matrix Bshear(2,3) ;

//---Bshear Matrix in standard {1,2,3} mechanics notation------
//
//             -                -
//   Bshear = | +N,1      0    +N |  (2x3)
//            | +N,2     -N     0 |
//             -                -  
//
//-------------------------------------------------------------

  Bshear.Zero( ) ;

  Bshear(0,0) =  shp[0][node] ;
  Bshear(0,2) =  shp[2][node] ;
  Bshear(1,0) =  shp[1][node] ;
  Bshear(1,1) = -shp[2][node] ;

  return Bshear ;

}  

//***********************************************************************
//compute Jacobian matrix and inverse at point {L1,L2} 

void   ShellMITC4::computeJacobian( double L1, double L2,
				    const double x[2][4], 
                                    Matrix &JJ, 
                                    Matrix &JJinv )
{
  int i, j, k ;
     
  static const double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static const double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static double shp[2][4] ;

  double ss = L1 ;
  double tt = L2 ;

  for ( i = 0; i < 4; i++ ) {
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  } // end for i

  
  // Construct jacobian and its inverse
  
  JJ.Zero( ) ;
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {

      for ( k = 0; k < 4; k++ )
	  JJ(i,j) +=  x[i][k] * shp[j][k] ;

    } //end for j
  }  // end for i 

  double xsj = JJ(0,0)*JJ(1,1) - JJ(0,1)*JJ(1,0) ;

  //inverse jacobian
  JJinv(0,0) =  JJ(1,1) / xsj ;
  JJinv(1,1) =  JJ(0,0) / xsj ;
  JJinv(0,1) = -JJ(0,1) / xsj ; 
  JJinv(1,0) = -JJ(1,0) / xsj ; 

  return ;

}

//************************************************************************
//shape function routine for four node quads

void  ShellMITC4::shape2d( double ss, double tt, 
		           const double x[2][4], 
		           double shp[3][4], 
		           double &xsj            )

{ 

  int i, j, k ;

  double temp ;
     
  static const double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static const double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static Matrix xs(2,2) ;
  static Matrix sx(2,2) ;

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

  //inverse jacobian
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
	   
//**********************************************************************

Matrix  ShellMITC4::transpose( int dim1, 
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

//**********************************************************************

int  ShellMITC4::sendSelf (int commitTag, Channel &theChannel)
{
    return -1;
}
    
int  ShellMITC4::recvSelf (int commitTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
    return -1;
}
//**************************************************************************

int
ShellMITC4::displaySelf(Renderer &theViewer, int displayMode, float fact)
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

    static Matrix coords(4,3);
    static Vector values(4);
    static Vector P(24) ;

    values(0) = 1 ;
    values(1) = 1 ;
    values(2) = 1 ;
    values(3) = 1 ;

    if (displayMode < 3 && displayMode > 0)
      P = this->getResistingForce();

    for (int i = 0; i < 3; i++) {
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
