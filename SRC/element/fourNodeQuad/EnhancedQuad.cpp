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
                                                                        
// $Revision: 1.2 $
// $Date: 2001-08-07 21:05:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/EnhancedQuad.cpp,v $

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
#include <EnhancedQuad.h>
#include <Renderer.h>


//static data
double  EnhancedQuad::xl[2][4] ;
Matrix  EnhancedQuad::stiff(8,8) ;
Vector  EnhancedQuad::resid(8) ;
Matrix  EnhancedQuad::mass(8,8) ;
Matrix  EnhancedQuad::damping(8,8) ;

double  EnhancedQuad::stressData[3][4] ;
double  EnhancedQuad::tangentData[3][3][4] ;

    
//quadrature data
const double  EnhancedQuad::root3 = sqrt(3.0) ;
const double  EnhancedQuad::one_over_root3 = 1.0 / root3 ;

const double  EnhancedQuad::sg[] = { -one_over_root3,  
				     one_over_root3, 
				     one_over_root3, 
				     -one_over_root3 } ;

const double  EnhancedQuad::tg[] = { -one_over_root3, 
				     -one_over_root3, 
				     one_over_root3,  
				     one_over_root3 } ;

const double  EnhancedQuad::wg[] = { 1.0, 1.0, 1.0, 1.0 } ;

  

//null constructor
EnhancedQuad::EnhancedQuad( ) :
Element( 0, ELE_TAG_EnhancedQuad ),
connectedExternalNodes(4),
alpha(4) 
{ 

}


//full constructor
EnhancedQuad::EnhancedQuad(  int tag, 
                         int node1,
                         int node2,
   	                 int node3,
                         int node4,
	                 NDMaterial &theMaterial,
			 const char *type ) 
  :
Element( tag, ELE_TAG_EnhancedQuad ),
connectedExternalNodes(4),
alpha(4)
{

  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

      materialPointers[i] = theMaterial.getCopy(type) ;

      if (materialPointers[i] == 0) {

	  g3ErrorHandler->fatal("EnhancedQuad::constructor %s",
		"- failed to get a material of type: ShellSection");
      } //end if
      
  } //end for i 

}


//destructor 
EnhancedQuad::~EnhancedQuad( )
{
  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ; 

    nodePointers[i] = 0 ;

  } //end for i
}


//set domain
void  EnhancedQuad::setDomain( Domain *theDomain ) 
{  

  int i ;

  //zero enhanced parameters
  alpha.Zero( ) ;

  //node pointers
  for ( i = 0; i < 4; i++ ) 
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

  this->DomainComponent::setDomain(theDomain) ;
}


//get the number of external nodes
int  EnhancedQuad::getNumExternalNodes( ) const
{
  return 4 ;
} 
 

//return connected external nodes
const ID&  EnhancedQuad::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 


//return number of dofs
int  EnhancedQuad::getNumDOF( ) 
{
  return 8 ;
}


//commit state
int  EnhancedQuad::commitState( )
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int  EnhancedQuad::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  EnhancedQuad::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  //zero enhanced parameters
  this->alpha.Zero( ) ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void  EnhancedQuad::Print( ostream &s, int flag )
{
  s << '\n' ;
  s << "Enhanced Strain Four Node Quad \n" ;
  s << "Element Number: " << this->getTag() << '\n' ;
  s << "Node 1 : " << connectedExternalNodes(0) << '\n' ;
  s << "Node 2 : " << connectedExternalNodes(1) << '\n' ;
  s << "Node 3 : " << connectedExternalNodes(2) << '\n' ;
  s << "Node 4 : " << connectedExternalNodes(3) << '\n' ;

  s << "Material Information : \n " ;
  materialPointers[0]->Print( s, flag ) ;

  s << endl ;
}

//return stiffness matrix 
const Matrix&  EnhancedQuad::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    


//return secant matrix 
const Matrix&  EnhancedQuad::getSecantStiff( ) 
{
   int tang_flag = 1 ; //get the tangent

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  
    
  return stiff ;
}
    

//return damping matrix because frank is a dumb ass 
const Matrix&  EnhancedQuad::getDamp( ) 
{
  //not supported
  return damping ;
}    


//return mass matrix
const Matrix&  EnhancedQuad::getMass( ) 
{

  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;

} 


//zero the load -- what load?
void  EnhancedQuad::zeroLoad( )
{
  return ;
}

//add load -- what load?
int  EnhancedQuad::addLoad( const Vector &addP )
{
  return -1 ;
}


//get residual
const Vector&  EnhancedQuad::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  return resid ;   
}


//get residual with inertia terms
const Vector&  EnhancedQuad::getResistingForceIncInertia( )
{
  int tang_flag = 0 ; //don't get the tangent

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  //inertia terms
  formInertiaTerms( tang_flag ) ;

  return resid ;
}


//*********************************************************************
//form inertia terms

void   EnhancedQuad::formInertiaTerms( int tangFlag ) 
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

  int i, j, k, p ;
  int jj, kk ;

  double temp, rho, massJK ;


  //zero mass 
  mass.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;


  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;
    
    //volume element
    dvol = wg[i] * xsj ;

    //node loop to compute acceleration
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
      //momentum += shp[massIndex][j] * ( nodePointers[j]->getTrialAccel()  ) ; 
      momentum.addVector( 1.0,
			  nodePointers[j]->getTrialAccel(),
			  shp[massIndex][j] ) ;

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
void  EnhancedQuad::formResidAndTangent( int tang_flag ) 
{

  static const double tolerance = 1.0e-08 ;

  static const int nIterations = 10 ;

  static const int ndm = 2 ;

  static const int ndf = 2 ; 

  static const int nstress = 3 ;
 
  static const int numberNodes = 4 ;

  static const int numberGauss = 4 ;

  static const int nShape = 3 ;

  static const int nEnhanced = 4 ;
  
  static const int nModes = 2 ;

  static const int numberDOF = 8 ;


  int i, j, k, p, q ;
  int jj, kk ;

  int success ;
  

  static double xsj[numberGauss] ;  // determinant jacaobian matrix 

  static double dvol[numberGauss] ; //volume element

  static Vector strain(nstress) ;  //strain

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  static Vector residJ(ndf) ; //nodeJ residual 

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Matrix stiffKJ(ndf,ndf) ; //nodeKJ stiffness

  static Vector stress(nstress) ;  //stress

  static Matrix dd(nstress,nstress) ;  //material tangent

  static Matrix J0(ndm,ndm) ; //Jacobian matrix at center of element

  static Matrix J0inv(ndm,ndm) ; //inverse of above


  static Matrix Kee(nEnhanced,nEnhanced) ;

  static Vector residE(nEnhanced) ;

  static Vector Umode(ndf) ;

  static Vector dalpha(nEnhanced) ;

  static Matrix Kue(numberDOF,nEnhanced) ;

  static Matrix Keu(nEnhanced,numberDOF) ;

  static Matrix KeeInvKeu(nEnhanced,numberDOF) ;


  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BKtran(ndf,nstress) ;

    static Matrix BJtranD(ndf,nstress) ;

    static Matrix BKtranD(ndf,nstress) ;
  //-------------------------------------------------------

  
  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;

  Kee.Zero( ) ;
  residE.Zero( ) ;

  Kue.Zero( ) ;
  Keu.Zero( ) ;


  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;
  computeJacobian( L1, L2, xl, J0, J0inv ) ; 


  //gauss loop to compute and save shape functions 
  double det ;
  for ( i = 0; i < numberGauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, det ) ;

    //save shape functions
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
	  Shape[p][q][i] = shp[p][q] ;
    } // end for p

    //save jacobian determinant
    xsj[i] = det ;

    //volume element to also be saved
    dvol[i] = wg[i] * det ;  

  } // end for i gauss loop


  //-------------------------------------------------------------------
  //newton loop to solve for enhanced strain parameters

  int count = 0 ;
  do {

    residE.Zero( ) ;
    Kee.Zero( ) ;

    //gauss loop
    for ( i = 0; i < numberGauss; i++ ) {

      //extract shape functions from saved array
      for ( p = 0; p < nShape; p++ ) {
	for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
      } // end for p


      //zero the strain
      strain.Zero( ) ;

      // j-node loop to compute nodal strain contributions
      for ( j = 0; j < numberNodes; j++ )  {

        //compute B matrix 
        BJ = computeB( j, shp ) ;
      
        //nodal displacements 
        const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

        //compute the strain
        //strain += (BJ*ul) ; 
	strain.addMatrixVector(1.0, BJ,ul,1.0) ;

      } // end for j


      // j-node loop to compute enhanced strain contributions
      for ( j = 0; j < nModes; j++ )  {

        //compute B matrix 
        BJ = computeBenhanced( j, sg[i], tg[i], xsj[i], J0inv ) ; 
      
        //enhanced "displacements" 
        Umode(0) = this->alpha( 2*j     ) ;
        Umode(1) = this->alpha( 2*j + 1 ) ;

        //compute the strain
        //strain += (BJ*Umode) ; 
	strain.addMatrixVector(1.0, BJ,Umode,1.0) ;

      } // end for j


      //send the strain to the material 
      success = materialPointers[i]->setTrialStrain( strain ) ;

      //compute the stress
      stress = materialPointers[i]->getStress( ) ;

      //multiply by volume element
      stress  *= dvol[i] ;

      //tangent 
      dd = materialPointers[i]->getTangent( ) ;

      //multiply by volume element
      dd *= dvol[i] ;


      //save stress and tangent (already multiplied by volume element)
      saveData( i, stress, dd ) ; 


      //enhanced residual and tangent calculations loops

      jj = 0 ;
      for ( j = 0; j < nModes; j++ ) {

        //compute B matrix 
        BJ = computeBenhanced( j, sg[i], tg[i], xsj[i], J0inv ) ; 

	//transpose 
	BJtran = transpose( BJ ) ;


	//residual
	//residJ = (BJtran * stress) ;
	//residJ *= (-1.0) ;
	residJ.addMatrixVector(0.0, BJtran,stress,-1.0) ;

	//residual 
	for ( p = 0; p < ndf; p++ )
	  residE( jj+p ) += residJ(p)  ;


        //BJtranD = BJtran * dd ;
	BJtranD.addMatrixProduct(0.0, BJtran,dd,1.0) ;

	kk = 0 ;
	for ( k = 0; k < nModes; k++ ) {

          BK = computeBenhanced( k, sg[i], tg[i], xsj[i], J0inv ) ;
  
          //stiffJK =  BJtranD * BK  ;
	  stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0) ;

          for ( p = 0; p < ndf; p++ )  {
             for ( q = 0; q < ndf; q++ )
                Kee( jj+p, kk+q ) += stiffJK( p, q ) ;
          } //end for p

            kk += ndf ;
          } // end for k loop

	jj += ndf ;
      } // end for j loop


    } // end for i gauss loop 


    //solve for enhanced strain parameters

    dalpha.Zero( ) ;

    Kee.Solve( residE, dalpha ) ;

    this->alpha += dalpha ;

    count++ ;
    if ( count > nIterations ) {
      cerr << "Exceeded " << nIterations
	   << " solving for enhanced strain parameters " << endl ;
      break ;
    } //end if 


    //do at least 2 iterations so saved data is good
  } while ( residE.Norm() > tolerance  ||  count < 2 ) ;


  //end enhanced strain parameters newton loop
  //-------------------------------------------------------------------



  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p


    //recover stress and tangent from saved data
    getData( i, stress, dd ) ;


    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      BJ = computeB( j, shp ) ;
   
      //transpose 
      BJtran = transpose( BJ ) ;

      //residual
      //residJ = BJtran * stress ;
      residJ.addMatrixVector(0.0, BJtran,stress,1.0) ;

      for ( p = 0; p < ndf; p++ )
        resid( jj+p ) += residJ(p)  ;


      if ( tang_flag == 1 ) {

	//BJtranD = BJtran * dd ;
	BJtranD.addMatrixProduct(0.0, BJtran,dd,1.0) ;

	 //node-node stiffness
         kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) {

            BK = computeB( k, shp ) ;
  
            //stiffJK =  BJtranD * BK  ;
	    stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0) ;

            for ( p = 0; p < ndf; p++ )  {
               for ( q = 0; q < ndf; q++ )
                  stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
            } //end for p

            kk += ndf ;
          } // end for k loop


	 //node-enhanced stiffness Kue 
         kk = 0 ;
         for ( k = 0; k < nModes; k++ ) {

            BK = computeBenhanced( k, sg[i], tg[i], xsj[i], J0inv ) ;
  
            //stiffJK =  BJtranD * BK  ;
	    stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0) ;
	   
            for ( p = 0; p < ndf; p++ )  {
               for ( q = 0; q < ndf; q++ )
                  Kue( jj+p, kk+q ) += stiffJK( p, q ) ;
            } //end for p 

            kk += ndf ;
          } // end for k loop


	 //enhanced-node stiffness Keu 
         kk = 0 ;
         for ( k = 0; k < nModes; k++ ) {

            BK = computeBenhanced( k, sg[i], tg[i], xsj[i], J0inv ) ;

	    //transpose 
	    BKtran = transpose( BK ) ;
	    
	    //BKtranD = BKtran * dd ;
	    BKtranD.addMatrixProduct(0.0, BKtran,dd,1.0 ) ;

            //stiffKJ = (BKtran*dd)*BJ ;
	    stiffKJ.addMatrixProduct(0.0, BKtranD,BJ,1.0) ;

            for ( p = 0; p < ndf; p++ )  {
               for ( q = 0; q < ndf; q++ )
                  Keu( kk+p, jj+q ) += stiffKJ( p, q ) ;
            } //end for p  

            kk += ndf ;
          } // end for k loop



      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop


  } //end for i gauss loop 


  //static condensation of enhanced parameters

  if ( tang_flag == 1 ) {
    
     Kee.Solve( Keu, KeeInvKeu ) ;

     //stiff -= ( Kue * KeeInvKeu ) ;
     stiff.addMatrixProduct(1.0,  Kue,KeeInvKeu,-1.0 ) ;

  } //end if 

  return ;
}
  

//************************************************************************
void   EnhancedQuad::saveData( int gp, 
			       const Vector &stress,
			       const Matrix &tangent ) 
{

  int i, j ;

  //save stress
  for ( i=0; i<3; i++ )
    stressData[i][gp] = stress(i) ;

  //save tangent
  for ( i=0; i<3; i++ ) {
    for (j=0; j<3; j++ ) 
      tangentData[i][j][gp] = tangent(i,j) ;
  } //end for i

  return ;
}


//************************************************************************
//recover stress and tangent data

void  EnhancedQuad::getData( int gp,
			     Vector &stress,
			     Matrix &tangent ) 
{

  int i, j ;

  //get stress
  for ( i=0; i<3; i++ )
    stress(i) = stressData[i][gp] ;

  //get tangent
  for ( i=0; i<3; i++ ) {
    for ( j=0; j<3; j++ ) 
      tangent(i,j) = tangentData[i][j][gp] ;
  } //end for i

  return ;

}

//************************************************************************
//compute local coordinates and basis

void   EnhancedQuad::computeBasis( ) 
{
  //nodal coordinates 

  int i ;
  for ( i = 0; i < 4; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( ) ;

       xl[0][i] = coorI(0) ;
       xl[1][i] = coorI(1) ;

  }  //end for i 

}

//*************************************************************************
//compute B matrix

const Matrix&
EnhancedQuad::computeB( int node, const double shp[3][4] ) 
{

  static Matrix B(3,2) ;

//---B Matrix in standard {1,2,3} mechanics notation---------------
//
//                -             -
//               | +N,1      0   | 
// B    =        |   0     +N,2  |    (3x2)
//               | +N,2    +N,1  |
//                -             -  
//
//-------------------------------------------------------------------

  B.Zero( ) ;

  B(0,0) = shp[0][node] ;
  B(1,1) = shp[1][node] ;
  B(2,0) = shp[1][node] ;
  B(2,1) = shp[0][node] ;

  return B ;

}

//***********************************************************************
//compute enhanced strain B-matrices

const Matrix&   
EnhancedQuad::computeBenhanced( int node, 
					 double L1,
					 double L2,
					 double j, 
					 const Matrix &Jinv )
{
  static Matrix B(3,2) ;

  static double JinvTran[2][2] ;

  static double shape[2] ;

  static double parameter ;


  //compute JinvTran
  JinvTran[0][0] = Jinv(0,0) ;
  JinvTran[1][1] = Jinv(1,1) ;
  JinvTran[0][1] = Jinv(1,0) ;
  JinvTran[1][0] = Jinv(0,1) ;      //residual 



  if ( node == 0 ) {

    //first column of JinvTran 
    shape[0] = JinvTran[0][0] ;
    shape[1] = JinvTran[1][0] ;

    parameter = L1 / j ;

  }
  else if ( node == 1 ) {

    //second column of JinvTran 
    shape[0] = JinvTran[0][1] ;
    shape[1] = JinvTran[1][1] ;

    parameter = L2 / j ;

  } //end if 

  shape[0] *= parameter ;
  shape[1] *= parameter ; 


  B.Zero( ) ;

  B(0,0) = shape[0] ;
  B(1,1) = shape[1] ;
  B(2,0) = shape[1] ;
  B(2,1) = shape[0] ;
  
  return B ;

}

//***********************************************************************
//compute Jacobian matrix and inverse at point {L1,L2} 

void   EnhancedQuad::computeJacobian( double L1, double L2,
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
  double jinv = 1.0 / xsj ;
  JJinv(0,0) =  JJ(1,1) * jinv ;
  JJinv(1,1) =  JJ(0,0) * jinv ;
  JJinv(0,1) = -JJ(0,1) * jinv ;
  JJinv(1,0) = -JJ(1,0) * jinv ;

  return ;

}

//************************************************************************
//shape function routine for four node quads

void  EnhancedQuad::shape2d( double ss, double tt, 
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
  //inverse jacobian
  double jinv = 1.0 / xsj ;
  sx(0,0) =  xs(1,1) * jinv ;
  sx(1,1) =  xs(0,0) * jinv ;
  sx(0,1) = -xs(0,1) * jinv ;
  sx(1,0) = -xs(1,0) * jinv ;


  //form global derivatives 
  
  for ( i = 0; i < 4; i++ ) {
    temp      = shp[0][i]*sx(0,0) + shp[1][i]*sx(1,0) ;
    shp[1][i] = shp[0][i]*sx(0,1) + shp[1][i]*sx(1,1) ;
    shp[0][i] = temp ;
  } // end for i

  return ;
}
	   
//**********************************************************************

const Matrix&
EnhancedQuad::transpose( const Matrix &M ) 
{
  int i ;
  int j ;

  //we're always transposing 3x2 matrices for this element,
  //so always return a 2x3 .

  static int dim1 = 2 ;
  static int dim2 = 3 ;
  static Matrix Mtran(dim1,dim2) ;

  for ( i = 0; i < dim1; i++ ) {
     for ( j = 0; j < dim2; j++ ) 
         Mtran(i,j) = M(j,i) ;
  } // end for i

  return Mtran ;
}

//**********************************************************************

int  EnhancedQuad::sendSelf (int commitTag, Channel &theChannel)
{
    return -1;
}
    
int  EnhancedQuad::recvSelf (int commitTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
    return -1;
}
//**************************************************************************

int
EnhancedQuad::displaySelf(Renderer &theViewer, int displayMode, float fact)
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

//***************************************************************************
/* 

    //zero the strains
    strain.Zero( ) ;


    // j-node loop to compute strain 
    for ( j = 0; j < numberNodes; j++ )  {

      //compute B matrix 
      BJ = computeB( j, shp ) ;
      
      //nodal displacements 
      const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

      //compute the strain
      strain += (BJ*ul) ; 

    } // end for j


    // j-node loop to compute enhanced strain contributions
    for ( j = 0; j < nModes; j++ )  {

      //compute B matrix 
      BJ = computeBenhanced( j, sg[i], tg[i], xsj[i], J0inv ) ; 
      
      //enhanced "displacements" 
      Umode(0) = this->alpha( 2*j     ) ;
      Umode(1) = this->alpha( 2*j + 1 ) ;

      //compute the strain
      strain += (BJ*Umode) ; 

    } // end for j


    //send the strain to the material 
    success = materialPointers[i]->setTrialStrain( strain ) ;

    //compute the stress
    stress = materialPointers[i]->getStress( ) ;

    //multiply by volume element
    stress  *= dvol[i] ;

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getTangent( ) ;
      dd *= dvol[i] ;
    } //end if tang_flag
    
*/
