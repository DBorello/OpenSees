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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateFiberMaterial.cpp,v $

//
// Ed "C++" Love
//
// Generic Plate Fiber Material
//


#include <PlateFiberMaterial.h>


//static vector and matrices
Vector  PlateFiberMaterial::stress(5) ;
Matrix  PlateFiberMaterial::tangent(5,5) ;



//null constructor
PlateFiberMaterial::PlateFiberMaterial( ) : 
NDMaterial(0, ND_TAG_PlateFiberMaterial ), 
strain(5) 
{ }



//full constructor
PlateFiberMaterial::PlateFiberMaterial(    
				   int tag, 
                                   NDMaterial &the3DMaterial ) :
NDMaterial( tag, ND_TAG_PlateFiberMaterial ),
strain(5)
{
  theMaterial = the3DMaterial.getCopy("ThreeDimensional") ;

  strain22 = 0.0 ;
}



//destructor
PlateFiberMaterial::~PlateFiberMaterial( ) 
{ 
  delete theMaterial ;
} 



//make a clone of this material
NDMaterial*
PlateFiberMaterial::getCopy( ) 
{
  PlateFiberMaterial *clone ;   //new instance of this class

  clone = new PlateFiberMaterial( this->getTag(), 
                                   *theMaterial ) ; //make the copy

  clone->strain22 = this->strain22 ;

  return clone ;
}


//make a clone of this material
NDMaterial* 
PlateFiberMaterial::getCopy( const char *type ) 
{
  return this->getCopy( ) ;
}


//send back order of strain in vector form
int 
PlateFiberMaterial::getOrder( ) const
{
  return 5 ;
}


const char*
PlateFiberMaterial::getType( ) const 
{
  return "PlateFiber" ; 
}



//swap history variables
int 
PlateFiberMaterial::commitState( ) 
{
  return theMaterial->commitState( ) ;
}



//revert to last saved state
int 
PlateFiberMaterial::revertToLastCommit( )
{
  return theMaterial->revertToLastCommit( )  ;
}


//revert to start
int
PlateFiberMaterial::revertToStart( )
{
  this->strain22 = 0.0 ;

  return theMaterial->revertToStart( ) ;
}


//mass per unit volume
double
PlateFiberMaterial::getRho( )
{
  return theMaterial->getRho( ) ;
}


//receive the strain
int 
PlateFiberMaterial::setTrialStrain( const Vector &strainFromElement )
{

  this->strain(0) = strainFromElement(0) ; //11
  this->strain(1) = strainFromElement(1) ; //22
  this->strain(2) = strainFromElement(2) ; //12
  this->strain(3) = strainFromElement(3) ; //23
  this->strain(4) = strainFromElement(4) ; //31


  static Vector threeDstrain(6) ;

  threeDstrain(0) = this->strain(0) ;
  threeDstrain(1) = this->strain(1) ;

  threeDstrain(2) = this->strain22 ;
  
  threeDstrain(3) = this->strain(2) ; 
  threeDstrain(4) = this->strain(3) ;
  threeDstrain(5) = this->strain(4) ;

  return theMaterial->setTrialStrain( threeDstrain ) ;

}


//send back the strain
const Vector& 
PlateFiberMaterial::getStrain( )
{
  return this->strain ;
}


//send back the stress 
const Vector&  
PlateFiberMaterial::getStress( )
{
  static const double tolerance = 1.0e-08 ;

  int success ;

  double norm ;

  static Vector outOfPlaneStress(1) ;

  static Vector strainIncrement(1) ;

  static Vector threeDstrain(6) ;

  static Vector threeDstress(6) ;

  static Matrix threeDtangent(6,6) ;

  static Vector threeDstressCopy(6) ; 

  static Matrix threeDtangentCopy(6,6) ;

  static Matrix dd22(1,1) ;

  int i, j ;

  int ii, jj ;


  int count = 0 ;
  //newton loop to solve for out-of-plane strains
  do {

    //set three dimensional strain
    threeDstrain(0) = this->strain(0) ;
    threeDstrain(1) = this->strain(1) ;

    threeDstrain(2) = this->strain22 ;
  
    threeDstrain(3) = this->strain(2) ; 
    threeDstrain(4) = this->strain(3) ;
    threeDstrain(5) = this->strain(4);

    success = theMaterial->setTrialStrain( threeDstrain ) ;
   

    //three dimensional stress
    threeDstress = theMaterial->getStress( ) ;

    //three dimensional tangent 
    threeDtangent = theMaterial->getTangent( ) ;


    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
    //PlateFiberMaterial strain order =  11, 22, 12, 23, 31, 33 

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
    
    outOfPlaneStress(0) = threeDstress(2) ;

    dd22(0,0) = threeDtangentCopy(5,5) ;

    for ( i=0; i<5; i++ ) 
      this->stress(i)     = threeDstressCopy(i) ;

    //set norm
    norm = outOfPlaneStress.Norm( ) ;

    //int Solve(const Vector &V, Vector &res) const;
    //int Solve(const Matrix &M, Matrix &res) const;
    //condensation 
    dd22.Solve( outOfPlaneStress, strainIncrement ) ;

    //update out of plane strains
    this->strain22 -= strainIncrement(0) ;

    count++ ;
   } while ( norm > tolerance && count < 10 ) ;

  
  return this->stress ;
}

//send back the tangent 
const Matrix&  
PlateFiberMaterial::getTangent( )
{
  static const double tolerance = 1.0e-08 ;

  int success ;

  double norm ;

  static Vector outOfPlaneStress(1) ;

  static Vector strainIncrement(1) ;

  static Vector threeDstrain(6) ;

  static Vector threeDstress(6) ;

  static Matrix threeDtangent(6,6) ;

  static Vector threeDstressCopy(6) ; 

  static Matrix threeDtangentCopy(6,6) ;

  static Matrix dd11(5,5) ;
  static Matrix dd12(5,1) ;
  static Matrix dd21(1,5) ;
  static Matrix dd22(1,1) ;

  static Matrix dd22invdd21(1,5) ;

  int i, j ;

  int ii, jj ;


  int count = 0 ;
  //newton loop
  do {

    //set three dimensional strain
    threeDstrain(0) = this->strain(0) ;
    threeDstrain(1) = this->strain(1) ;

    threeDstrain(2) = this->strain22 ;
  
    threeDstrain(3) = this->strain(2) ; 
    threeDstrain(4) = this->strain(3) ;
    threeDstrain(5) = this->strain(4);

    success = theMaterial->setTrialStrain( threeDstrain ) ;
   

    //three dimensional stress
    threeDstress = theMaterial->getStress( ) ;

    //three dimensional tangent 
    threeDtangent = theMaterial->getTangent( ) ;


    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
    //PlateFiberMaterial strain order  = 11, 22, 12, 23, 31, 33 

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

    outOfPlaneStress(0) = threeDstressCopy(5) ;

    dd22(0,0) = threeDtangentCopy(5,5) ;

    for ( i=0; i<5; i++ ) {

	dd12(i,0) = threeDtangentCopy(i,5) ;
	dd21(0,i) = threeDtangentCopy(5,i) ;

      for ( j=0; j<5; j++) 
	dd11(i,j) = threeDtangentCopy(i,j) ;

    }//end for i


    //set norm
    norm = outOfPlaneStress.Norm( ) ;

    //int Solve(const Vector &V, Vector &res) const;
    //int Solve(const Matrix &M, Matrix &res) const;
    //condensation 
    dd22.Solve( outOfPlaneStress, strainIncrement ) ;

    //update out of plane strains
    this->strain22 -= strainIncrement(0) ;
    
    count++ ;
  } while ( norm > tolerance && count < 10 ) ;

    
  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  dd22.Solve( dd21, dd22invdd21 ) ;
  this->tangent   = dd11 ; 
  this->tangent  -= ( dd12*dd22invdd21 ) ;

  return this->tangent ;
}



int 
PlateFiberMaterial::indexMap( int i )
{
  int ii ;

  switch ( i+1 ) { //add 1 for standard vector indices

    case 1 :
      ii = 1 ; 
      break ;
 
    case 2 :
      ii = 2 ;
      break ;

    case 3 :
      ii = 6 ;
      break ;

    case 4 :
      ii = 3 ;
      break ;

    case 5 :
      ii = 4 ;
      break ;

    case 6 :
      ii = 5 ;
      break ;

    default :
      ii = 1 ;
      break ;

  } //end switch

  ii--;

  return ii ;
}



//print out data
void  
PlateFiberMaterial::Print( ostream &s, int flag )
{
  s << "General Plate Fiber Material \n" ;
  s << " Tag: " << this->getTag() << "\n" ; 
  s << "using the 3D material : \n" ;

  theMaterial->Print( s, flag ) ;

  return ;
}


int 
PlateFiberMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  return -1;
}

int 
PlateFiberMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}
 


