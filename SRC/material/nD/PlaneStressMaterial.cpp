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
                                                                        
// $Revision: 1.1 $
// $Date: 2001-08-07 20:59:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressMaterial.cpp,v $

//
// Ed "C++" Love
//
// Generic Plane Stress Material
//


#include <PlaneStressMaterial.h>


//static vector and matrices
Vector  PlaneStressMaterial::stress(3) ;
Matrix  PlaneStressMaterial::tangent(3,3) ;



//null constructor
PlaneStressMaterial::PlaneStressMaterial( ) : 
NDMaterial(0, ND_TAG_PlaneStressMaterial ), 
strain(3) 
{ }



//full constructor
PlaneStressMaterial::PlaneStressMaterial(    
				   int tag, 
                                   NDMaterial &the3DMaterial ) :
NDMaterial( tag, ND_TAG_PlaneStressMaterial ),
strain(3)
{
  theMaterial = the3DMaterial.getCopy("ThreeDimensional") ;

  strain22 = 0.0 ;
  gamma02 = 0.0 ;
  gamma12 = 0.0 ;
}



//destructor
PlaneStressMaterial::~PlaneStressMaterial( ) 
{ 
  delete theMaterial ;
} 



//make a clone of this material
NDMaterial*
PlaneStressMaterial::getCopy( ) 
{
  PlaneStressMaterial *clone ;   //new instance of this class

  clone = new PlaneStressMaterial( this->getTag(), 
                                   *theMaterial ) ; //make the copy

  clone->strain22 = this->strain22 ;
  clone->gamma02  = this->gamma02 ;
  clone->gamma12  = this->gamma12 ;

  return clone ;
}


//make a clone of this material
NDMaterial* 
PlaneStressMaterial::getCopy( const char *type ) 
{
  return this->getCopy( ) ;
}


//send back order of strain in vector form
int 
PlaneStressMaterial::getOrder( ) const
{
  return 3 ;
}


const char*
PlaneStressMaterial::getType( ) const 
{
  return "PlaneStress" ; 
}



//swap history variables
int 
PlaneStressMaterial::commitState( ) 
{
  return theMaterial->commitState( ) ;
}



//revert to last saved state
int 
PlaneStressMaterial::revertToLastCommit( )
{
  return theMaterial->revertToLastCommit( )  ;
}


//revert to start
int
PlaneStressMaterial::revertToStart( )
{
  this->strain22 = 0.0 ;
  this->gamma12  = 0.0 ;
  this->gamma02  = 0.0 ;

  return theMaterial->revertToStart( ) ;
}


//mass per unit volume
double
PlaneStressMaterial::getRho( )
{
  return theMaterial->getRho( ) ;
}


//receive the strain
int 
PlaneStressMaterial::setTrialStrain( const Vector &strainFromElement )
{

  this->strain(0) = strainFromElement(0) ;
  this->strain(1) = strainFromElement(1) ;
  this->strain(2) = strainFromElement(2) ;


  static Vector threeDstrain(6) ;

  threeDstrain(0) = this->strain(0) ;
  threeDstrain(1) = this->strain(1) ;

  threeDstrain(2) = this->strain22 ;
  
  threeDstrain(3) = this->strain(2) ; 

  threeDstrain(4) = this->gamma12 ;
  threeDstrain(5) = this->gamma02 ;

  return theMaterial->setTrialStrain( threeDstrain ) ;

}


//send back the strain
const Vector& 
PlaneStressMaterial::getStrain( )
{
  return this->strain ;
}


//send back the stress 
const Vector&  
PlaneStressMaterial::getStress( )
{
  static const double tolerance = 1.0e-08 ;

  int success ;

  double norm ;

  static Vector outOfPlaneStress(3) ;

  static Vector strainIncrement(3) ;

  static Vector threeDstrain(6) ;

  static Vector threeDstress(6) ;

  static Matrix threeDtangent(6,6) ;

  static Vector threeDstressCopy(6) ; 

  static Matrix threeDtangentCopy(6,6) ;

  static Matrix dd22(3,3) ;

  int i, j ;

  int ii, jj ;


  //newton loop to solve for out-of-plane strains
  do {

    //set three dimensional strain
    threeDstrain(0) = this->strain(0) ;
    threeDstrain(1) = this->strain(1) ;

    threeDstrain(2) = this->strain22 ;
  
    threeDstrain(3) = this->strain(2) ; 

    threeDstrain(4) = this->gamma12 ;
    threeDstrain(5) = this->gamma02 ;

    success = theMaterial->setTrialStrain( threeDstrain ) ;
   

    //three dimensional stress
    threeDstress = theMaterial->getStress( ) ;

    //three dimensional tangent 
    threeDtangent = theMaterial->getTangent( ) ;


    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
    //PlaneStressMaterial strain order = 11, 22, 12, 33, 23, 31 

    //swap matrix indices to sort out-of-plane components 
    for ( i=0; i<6; i++ ) {

      ii = this->indexMap(i) ;

      threeDstressCopy(ii) = threeDstress(i) ;

      for ( j=0; j<6; j++ ) {

	jj = this->indexMap(j) ;
	
	threeDtangentCopy(ii,jj) = threeDtangent(i,j) ;

      }//end for j
       
    }//end for i


    //partitioned stresses and tangent
    for ( i=0; i<3; i++ ) {

      this->stress(i)     = threeDstressCopy(i  ) ;
      outOfPlaneStress(i) = threeDstressCopy(i+3) ;

      for ( j=0; j<3; j++ ) 
	dd22(i,j) = threeDtangentCopy(i+3,j+3) ;

    }//end for i


    //set norm
    norm = outOfPlaneStress.Norm( ) ;

    //int Solve(const Vector &V, Vector &res) const;
    //int Solve(const Matrix &M, Matrix &res) const;
    //condensation 
    dd22.Solve( outOfPlaneStress, strainIncrement ) ;

    //update out of plane strains
    this->strain22 -= strainIncrement(0) ;
    this->gamma12  -= strainIncrement(1) ;
    this->gamma02  -= strainIncrement(2) ;


  } while ( norm > tolerance ) ;

  
  return this->stress ;
}


//send back the tangent 
const Matrix&  
PlaneStressMaterial::getTangent( )
{
  static const double tolerance = 1.0e-08 ;

  int success ;

  double norm ;

  static Vector outOfPlaneStress(3) ;

  static Vector strainIncrement(3) ;

  static Vector threeDstrain(6) ;

  static Vector threeDstress(6) ;

  static Matrix threeDtangent(6,6) ;

  static Vector threeDstressCopy(6) ; 

  static Matrix threeDtangentCopy(6,6) ;

  static Matrix dd11(3,3) ;
  static Matrix dd12(3,3) ;
  static Matrix dd21(3,3) ;
  static Matrix dd22(3,3) ;

  static Matrix dd22invdd21(3,3) ;

  int i, j ;

  int ii, jj ;


  //newton loop
  do {

    //set three dimensional strain
    threeDstrain(0) = this->strain(0) ;
    threeDstrain(1) = this->strain(1) ;

    threeDstrain(2) = this->strain22 ;
  
    threeDstrain(3) = this->strain(2) ; 

    threeDstrain(4) = this->gamma12 ;
    threeDstrain(5) = this->gamma02 ;

    success = theMaterial->setTrialStrain( threeDstrain ) ;
   

    //three dimensional stress
    threeDstress = theMaterial->getStress( ) ;

    //three dimensional tangent 
    threeDtangent = theMaterial->getTangent( ) ;


    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
    //PlaneStressMaterial strain order = 11, 22, 12, 33, 23, 31 

    //swap matrix indices to sort out-of-plane components 
    for ( i=0; i<6; i++ ) {

      ii = this->indexMap(i) ;

      threeDstressCopy(ii) = threeDstress(i) ;

      for ( j=0; j<6; j++ ) {

	jj = this->indexMap(j) ;
	
	threeDtangentCopy(ii,jj) = threeDtangent(i,j) ;

      }//end for j
       
    }//end for i


    //out of plane stress and tangents
    for ( i=0; i<3; i++ ) {

      outOfPlaneStress(i) = threeDstressCopy(i+3) ;

      for ( j=0; j<3; j++ ) {

	dd11(i,j) = threeDtangentCopy(i,  j  ) ;
	dd12(i,j) = threeDtangentCopy(i,  j+3) ;
	dd21(i,j) = threeDtangentCopy(i+3,j  ) ;
	dd22(i,j) = threeDtangentCopy(i+3,j+3) ;

      }//end for j

    }//end for i

    //set norm
    norm = outOfPlaneStress.Norm( ) ;

    //int Solve(const Vector &V, Vector &res) const;
    //int Solve(const Matrix &M, Matrix &res) const;
    //condensation 
    dd22.Solve( outOfPlaneStress, strainIncrement ) ;

    //update out of plane strains
    this->strain22 -= strainIncrement(0) ;
    this->gamma12  -= strainIncrement(1) ;
    this->gamma02  -= strainIncrement(2) ;
    

  } while ( norm > tolerance ) ;

    
  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  dd22.Solve( dd21, dd22invdd21 ) ;
  this->tangent   = dd11 ; 
  this->tangent  -= ( dd12*dd22invdd21 ) ;

  return this->tangent ;
}



int 
PlaneStressMaterial::indexMap( int i )
{
  int ii ;

  if ( i == 2 ) 
    ii = 3 ;
  else if ( i == 3 )
    ii = 2 ;
  else 
    ii = i ;

  return ii ;
}



//print out data
void  
PlaneStressMaterial::Print( ostream &s, int flag )
{
  s << "General Plane Stress Material \n" ;
  s << " Tag: " << this->getTag() << "\n" ; 
  s << "using the 3D material : \n" ;

  theMaterial->Print( s, flag ) ;

  return ;
}


int 
PlaneStressMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  return -1;
}

int 
PlaneStressMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}
 


