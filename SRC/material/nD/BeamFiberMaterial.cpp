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
// $Date: 2001-10-03 18:07:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BeamFiberMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition of BeamFiberMaterial.
// The BeamFiberMaterial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11, 12, and 13
// stress components which can then be integrated over an area to model a
// shear flexible 3D beam.

#include <BeamFiberMaterial.h>

Vector BeamFiberMaterial::stress(3);
Matrix BeamFiberMaterial::tangent(3,3);

BeamFiberMaterial::BeamFiberMaterial(void)
: NDMaterial(0, ND_TAG_BeamFiberMaterial),
strain22(0.0), strain33(0.0), gamma23(0.0),
theMaterial(0), strain(3)
{
	// Nothing to do
}

BeamFiberMaterial::BeamFiberMaterial(int tag, NDMaterial &theMat)
: NDMaterial(tag, ND_TAG_BeamFiberMaterial),
strain22(0.0), strain33(0.0), gamma23(0.0),
theMaterial(0), strain(3)

{
	// Get a copy of the material
	theMaterial = theMat.getCopy("ThreeDimensional");

	if (theMaterial == 0)
		g3ErrorHandler->fatal("%s -- failed to get copy of material",
			"BeamFiberMaterial::BeamFiberMaterial");
}

BeamFiberMaterial::~BeamFiberMaterial(void) 
{ 
	if (theMaterial != 0)
		delete theMaterial;
} 

NDMaterial*
BeamFiberMaterial::getCopy(void) 
{
	BeamFiberMaterial *theCopy =
		new BeamFiberMaterial(this->getTag(), *theMaterial);

	theCopy->strain22 = this->strain22;
	theCopy->strain33 = this->strain33;
	theCopy->gamma23  = this->gamma23;

	return theCopy;
}

NDMaterial* 
BeamFiberMaterial::getCopy(const char *type)
{
	if (strcmp(type, "BeamFiber") == 0)
		return this->getCopy();
	else
		return 0;
}

int 
BeamFiberMaterial::getOrder(void) const
{
	return 3;
}

const char*
BeamFiberMaterial::getType(void) const 
{
	return "BeamFiber";
}

int 
BeamFiberMaterial::commitState(void)
{
	return theMaterial->commitState();
}

int 
BeamFiberMaterial::revertToLastCommit(void)
{
	return theMaterial->revertToLastCommit();
}

int
BeamFiberMaterial::revertToStart()
{
	this->strain22 = 0.0;
	this->strain33 = 0.0;
	this->gamma23  = 0.0;

	return theMaterial->revertToStart( );
}

double
BeamFiberMaterial::getRho(void)
{
	return theMaterial->getRho();
}


//receive the strain
//NDmaterial strain order        = 11, 22, 33, 12, 23, 31  
//BeamFiberMaterial strain order = 11, 12, 31, 22, 33, 23
int 
BeamFiberMaterial::setTrialStrain(const Vector &strainFromElement)
{
  this->strain(0) = strainFromElement(0);
  this->strain(1) = strainFromElement(1);
  this->strain(2) = strainFromElement(2);

  static Vector threeDstrain(6);

  threeDstrain(0) = this->strain(0);
  threeDstrain(1) = this->strain22;
  threeDstrain(2) = this->strain33;
  threeDstrain(3) = this->strain(1);
  threeDstrain(4) = this->gamma23;
  threeDstrain(5) = this->strain(2);

  return theMaterial->setTrialStrain(threeDstrain);
}

const Vector& 
BeamFiberMaterial::getStrain(void)
{
	return this->strain;
}

const Vector&  
BeamFiberMaterial::getStress( )
{
  static const double tolerance = 1.0e-08 ;

  int success ;

  double norm ;

  static Vector condensedStress(3) ;

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
    threeDstrain(1) = this->strain22 ;
    threeDstrain(2) = this->strain33 ;
    threeDstrain(3) = this->strain(1) ; 
    threeDstrain(4) = this->gamma23 ;
    threeDstrain(5) = this->strain(2) ;

    success = theMaterial->setTrialStrain( threeDstrain ) ;
   

    //three dimensional stress
    threeDstress = theMaterial->getStress( ) ;

    //three dimensional tangent 
    threeDtangent = theMaterial->getTangent( ) ;


    //NDmaterial strain order        = 11, 22, 33, 12, 23, 31  
	//BeamFiberMaterial strain order = 11, 12, 31, 22, 33, 23

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
      condensedStress(i) = threeDstressCopy(i+3) ;

      for ( j=0; j<3; j++ ) 
	dd22(i,j) = threeDtangentCopy(i+3,j+3) ;

    }//end for i


    //set norm
    norm = condensedStress.Norm( ) ;

    //int Solve(const Vector &V, Vector &res) const;
    //int Solve(const Matrix &M, Matrix &res) const;
    //condensation 
    dd22.Solve( condensedStress, strainIncrement ) ;

    //update out of plane strains
    this->strain22 -= strainIncrement(0) ;
    this->strain33 -= strainIncrement(1) ;
    this->gamma23  -= strainIncrement(2) ;

  } while ( norm > tolerance ) ;

  
  return this->stress ;
}

const Matrix&  
BeamFiberMaterial::getTangent( )
{
  static const double tolerance = 1.0e-08 ;

  int success ;

  double norm ;

  static Vector condensedStress(3) ;

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
    threeDstrain(1) = this->strain22 ;
    threeDstrain(2) = this->strain33 ;
    threeDstrain(3) = this->strain(1) ; 
    threeDstrain(4) = this->gamma23 ;
    threeDstrain(5) = this->strain(2) ;

    success = theMaterial->setTrialStrain( threeDstrain ) ;
   

    //three dimensional stress
    threeDstress = theMaterial->getStress( ) ;

    //three dimensional tangent 
    threeDtangent = theMaterial->getTangent( ) ;


    //NDmaterial strain order        = 11, 22, 33, 12, 23, 31 
	//BeamFiberMaterial strain order = 11, 12, 31, 22, 33, 23

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

      condensedStress(i) = threeDstressCopy(i+3) ;

      for ( j=0; j<3; j++ ) {

	dd11(i,j) = threeDtangentCopy(i,  j  ) ;
	dd12(i,j) = threeDtangentCopy(i,  j+3) ;
	dd21(i,j) = threeDtangentCopy(i+3,j  ) ;
	dd22(i,j) = threeDtangentCopy(i+3,j+3) ;

      }//end for j

    }//end for i

    //set norm
    norm = condensedStress.Norm( ) ;

    //int Solve(const Vector &V, Vector &res) const;
    //int Solve(const Matrix &M, Matrix &res) const;
    //condensation 
    dd22.Solve( condensedStress, strainIncrement ) ;

    //update out of plane strains
    this->strain22 -= strainIncrement(0) ;
    this->strain33 -= strainIncrement(1) ;
    this->gamma23  -= strainIncrement(2) ;
    

  } while ( norm > tolerance ) ;

    
  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  dd22.Solve( dd21, dd22invdd21 ) ;
  this->tangent   = dd11 ; 
  this->tangent  -= ( dd12*dd22invdd21 ) ;

  return this->tangent ;
}

//NDmaterial strain order        = 11, 22, 33, 12, 23, 31 
//BeamFiberMaterial strain order = 11, 12, 31, 22, 33, 23
int 
BeamFiberMaterial::indexMap( int i )
{
  int ii ;

  if ( i == 3 ) 
	  ii = 1 ;
  else if ( i == 5 )
	  ii = 2 ;
  else if ( i == 1 )
	  ii = 3;
  else if ( i == 2 )
	  ii = 4;
  else if ( i == 4 )
	  ii = 5;
  else 
	  ii = i ;

  return ii ;
}

void  
BeamFiberMaterial::Print( ostream &s, int flag )
{
	s << "BeamFiberMaterial, tag: " << this->getTag() << endl;
	s << "\tWrapped material: "<< theMaterial->getTag() << endl;

	theMaterial->Print( s, flag );
}

int 
BeamFiberMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
	return -1;
}

int 
BeamFiberMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return -1;
}
