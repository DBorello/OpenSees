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
                                                                        
// $Revision: 1.19 $
// $Date: 2003-03-11 02:56:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/Truss.cpp,v $
                                                                        
                                                                        
// File: ~/element/truss/Truss.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the implementation for the Truss class.
//
// What: "@(#) Truss.C, revA"

#include "Truss.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

//#include <fstream>

// initialise the class wide variables
Matrix Truss::trussM2(2,2);
Matrix Truss::trussM4(4,4);
Matrix Truss::trussM6(6,6);
Matrix Truss::trussM12(12,12);
Vector Truss::trussV2(2);
Vector Truss::trussV4(4);
Vector Truss::trussV6(6);
Vector Truss::trussV12(12);

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the truss end nodes.
Truss::Truss(int tag, 
	     int dim,
	     int Nd1, int Nd2, 
	     UniaxialMaterial &theMat,
	     double a, double rho)
 :Element(tag,ELE_TAG_Truss),     
  theMaterial(0), connectedExternalNodes(2),
  dimension(dim), numDOF(0), theLoad(0),
  theMatrix(0), theVector(0), t(0),
  L(0.0), A(a), M(rho)
{
    // get a copy of the material and check we obtained a valid copy
    theMaterial = theMat.getCopy();
    if (theMaterial == 0) {
      opserr << "FATAL Truss::Truss - " << tag <<
	"failed to get a copy of material with tag " << theMat.getTag() << endln;
      exit(-1);
    }
    
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL Truss::Truss - " <<  tag << "failed to create an ID of size 2\n";
      exit(-1);
    }

    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;        

    // set node pointers to NULL
    for (int i=0; i<2; i++)
      theNodes[i] = 0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	theLoadSens = 0;
// AddingSensitivity:END //////////////////////////////////////
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
Truss::Truss()
:Element(0,ELE_TAG_Truss),     
 theMaterial(0),connectedExternalNodes(2),
 dimension(0), numDOF(0),
 theMatrix(0), theVector(0), t(0), 
  L(0.0), A(0.0), M(0.0)
{
    // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL Truss::Truss - failed to create an ID of size 2\n";
      exit(-1);
  }

  for (int i=0; i<2; i++)
    theNodes[i] = 0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	theLoadSens = 0;
// AddingSensitivity:END //////////////////////////////////////
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
Truss::~Truss()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theMaterial != 0)
	delete theMaterial;
    if (t != 0)
	delete t;
    if (theLoad != 0)
	delete theLoad;
    if (theLoadSens != 0)
	delete theLoadSens;
}


int
Truss::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
Truss::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
Truss::getNodePtrs(void) 
{
  return theNodes;
}

int
Truss::getNumDOF(void) 
{
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the truss element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
Truss::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	L = 0;
	return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	
    
    // if can't find both - send a warning message
    if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
      if (theNodes[0] == 0)
	opserr <<"Truss::setDomain() - truss" << this->getTag() << " node " << Nd1 <<
	  "does not exist in the model\n";
      else
	opserr <<"Truss::setDomain() - truss" << this->getTag() << " node " << Nd2 <<
	  "does not exist in the model\n";

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	

      return;
    }

    // now determine the number of dof and the dimesnion    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) {
      opserr <<"WARNING Truss::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
	"have differing dof at ends for truss " << this->getTag() << endln;

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
	
      return;
    }	

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // now set the number of dof for element and set matrix and vector pointer
    if (dimension == 1 && dofNd1 == 1) {
	numDOF = 2;    
	theMatrix = &trussM2;
	theVector = &trussV2;
    }
    else if (dimension == 2 && dofNd1 == 2) {
	numDOF = 4;
	theMatrix = &trussM4;
	theVector = &trussV4;	
    }
    else if (dimension == 2 && dofNd1 == 3) {
	numDOF = 6;	
	theMatrix = &trussM6;
	theVector = &trussV6;		
    }
    else if (dimension == 3 && dofNd1 == 3) {
	numDOF = 6;	
	theMatrix = &trussM6;
	theVector = &trussV6;			
    }
    else if (dimension == 3 && dofNd1 == 6) {
	numDOF = 12;	    
	theMatrix = &trussM12;
	theVector = &trussV12;			
    }
    else {
      opserr <<"WARNING Truss::setDomain cannot handle " << dimension << " dofs at nodes in " << 
	dofNd1  << " problem\n";

      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
      return;
    }

    // create a transformation matrix for the element
    t = new Matrix(1,numDOF);
    theLoad = new Vector(numDOF);
    
    if (t == 0 || (t->noCols() != numDOF)) {
      opserr << "Truss::setDomain - truss " << this->getTag() <<
	"out of memory creating T matrix of size 1 x" << numDOF << endln;
      exit(-1);
      return;
    }      
    
    if (theLoad == 0) {
      opserr << "Truss::setDomain - truss " << this->getTag() << 
	"out of memory creating vector of size" << numDOF << endln;
      exit(-1);
      return;
    }          
    
    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    if (dimension == 1) {
	Matrix &trans = *t;      
	trans(0,0) = -1;
	trans(0,1) = 1;

	double dx = end2Crd(0)-end1Crd(0);	
	L = sqrt(dx*dx);
	

	if (L == 0.0) {
	  opserr <<"WARNING Truss::setDomain() - truss " << this->getTag() << " has zero length\n";
	  return;
	}	

    } else if (dimension == 2) {
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
    
	L = sqrt(dx*dx + dy*dy);
    
	if (L == 0.0) {
	  opserr <<"WARNING Truss::setDomain() - truss " << this->getTag() << " has zero length\n";
	  return;
	}
	
	double cs = dx/L;
	double sn = dy/L;

	Matrix &trans = *t;
	if (numDOF == 4) {
	    trans(0,0) = -cs;
	    trans(0,1) = -sn;
	    trans(0,2) = cs;
	    trans(0,3) = sn;
	} else { // it must be 6
	    trans(0,0) = -cs;
	    trans(0,1) = -sn;
	    trans(0,2) = 0.0;
	    trans(0,3) = cs;
	    trans(0,4) = sn;	
	    trans(0,5) = 0.0;
	}     

    } else {
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
	double dz = end2Crd(2)-end1Crd(2);		
    
	L = sqrt(dx*dx + dy*dy + dz*dz);
    
	if (L == 0.0) {
	  opserr <<"WARNING Truss::setDomain() - truss " << this->getTag() << " has zero length\n";
	  return;
	}
	
	double cx = dx/L;
	double cy = dy/L;
	double cz = dz/L;	

	Matrix &trans = *t;	
	if (numDOF == 6) {
	    trans(0,0) = -cx;
	    trans(0,1) = -cy;
	    trans(0,2) = -cz;
	    trans(0,3) = cx;
	    trans(0,4) = cy;
	    trans(0,5) = cz;	    
	} else { // it must be 12
	    trans(0,0) = -cx;
	    trans(0,1) = -cy;
	    trans(0,2) = -cz;
	    trans(0,3) = 0;
	    trans(0,4) = 0;
	    trans(0,5) = 0;	    
	    trans(0,6) = cx;
	    trans(0,7) = cy;
	    trans(0,8) = cz;
	    trans(0,9) = 0;
	    trans(0,10) = 0;
	    trans(0,11) = 0;	    	    
	}     
    }

    // determine the nodal mass for lumped mass approach
    M = M * A * L/2;
}   	 


int
Truss::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "Truss::commitState () - failed in base class";
  }    
  retVal = theMaterial->commitState();
  return retVal;
}

int
Truss::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
Truss::revertToStart()
{
    return theMaterial->revertToStart();
}

int
Truss::update(void)
{
    // determine the current strain given trial displacements at nodes
    double strain = this->computeCurrentStrain();
    double rate = this->computeCurrentStrainRate();
    return theMaterial->setTrialStrain(strain, rate);
}


const Matrix &
Truss::getTangentStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    double E = theMaterial->getTangent();

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;
    Matrix &trans = *t;

    stiff = trans^trans;

    stiff *= A*E/L;  

    return *theMatrix;
}


const Matrix &
Truss::getInitialStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    double E = theMaterial->getInitialTangent();

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;
    Matrix &trans = *t;

    stiff = trans^trans;

    stiff *= A*E/L;  

    return *theMatrix;
}
/*
const Matrix &
Truss::getSecantStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    // get the current E from the material for this strain
    double strain = this->computeCurrentStrain();

    double stress = theMaterial->getStress();    
    double E = stress/strain;

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;
    Matrix &trans = *t;
    stiff = trans^trans;
    stiff *= A*E/L;  
    
    return *theMatrix;
}
*/  
const Matrix &
Truss::getDamp(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    double eta = theMaterial->getDampTangent();

    // come back later and redo this if too slow
    Matrix &damp = *theMatrix;
    Matrix &trans = *t;

    damp = trans^trans;

    damp *= A*eta/L;

    return *theMatrix;
}


const Matrix &
Truss::getMass(void)
{   // zero the matrix
    theMatrix->Zero();    
  
    // check for quick return
    if (L == 0.0 || M == 0.0) { // - problem in setDomain() no further warnings
	return *theMatrix;
    }    

    Matrix &mass = *theMatrix;
    if (dimension == 1 && numDOF == 2) {
	mass(0,0) = M; 
	mass(1,1) = M;
    }
    else if (dimension == 2 && numDOF == 4) {
	mass(0,0) = M; 
	mass(1,1) = M;
	mass(2,2) = M; 
	mass(3,3) = M;	
    }
    else if (dimension == 2 && numDOF == 6) {
	mass(0,0) = M; 
	mass(1,1) = M;
	mass(3,3) = M;
	mass(4,4) = M; 
    }
    else if (dimension == 3 && numDOF == 6) {
	mass(0,0) = M; 
	mass(1,1) = M;
	mass(2,2) = M; 
	mass(3,3) = M;
	mass(4,4) = M; 
	mass(5,5) = M;		
    }
    else if (dimension == 3 && numDOF == 12) {
	mass(0,0) = M; 
	mass(1,1) = M;
	mass(2,2) = M; 
	mass(6,6) = M; 
	mass(7,7) = M;
	mass(8,8) = M; 
    }

    return *theMatrix; // so it will compile
}



void 
Truss::zeroLoad(void)
{
    theLoad->Zero();
}

int 
Truss::addLoad(ElementalLoad *theLoad, double loadFactor)

{  
  opserr <<"Truss::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  
  return -1;
}

int 
Truss::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for a quick return
    if (L == 0.0 || M == 0.0) 
	return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);    

  int nodalDOF = numDOF/2;
    
#ifdef _G3DEBUG    
  if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
    opserr <<"Truss::addInertiaLoadToUnbalance " <<
      "matrix and vector sizes are incompatable\n";
    return -1;
  }
#endif
    
    // want to add ( - fact * M R * accel ) to unbalance
    for (int i=0; i<dimension; i++) {
	double val1 = Raccel1(i);
	double val2 = Raccel2(i);	
	
	// perform - fact * M*(R * accel) // remember M a diagonal matrix
	val1 *= -M;
	val2 *= -M;
	
	(*theLoad)(i) += val1;
	(*theLoad)(i+nodalDOF) += val2;
    }	

    return 0;
}


int 
Truss::addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool somethingRandomInMotions)
{

  if (theLoadSens == 0) {
    theLoadSens = new Vector(numDOF);
  }
  else {
    theLoadSens->Zero();
  }
  
  
  if (somethingRandomInMotions) {
    
    
    // check for a quick return
    if (L == 0.0 || M == 0.0) 
      return 0;
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    
    
    int nodalDOF = numDOF/2;
    
#ifdef _G3DEBUG    
    if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
      opserr << "Truss::addInertiaLoadToUnbalance " <<
	"matrix and vector sizes are incompatable\n";
      return -1;
    }
#endif
    
    // want to add ( - fact * M R * accel ) to unbalance
    for (int i=0; i<dimension; i++) {
      double val1 = Raccel1(i);
      double val2 = Raccel2(i);	
      
      // perform - fact * M*(R * accel) // remember M a diagonal matrix
      val1 *= M;
      val2 *= M;
      
      (*theLoadSens)(i) = val1;
      (*theLoadSens)(i+nodalDOF) = val2;
    }	
  }
  else {
    
    // check for a quick return
    if (L == 0.0 || M == 0.0) 
      return 0;
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    
    
    int nodalDOF = numDOF/2;
    
#ifdef _G3DEBUG    
    if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
      opserr << "Truss::addInertiaLoadToUnbalance " <<
	"matrix and vector sizes are incompatable\n";
      return -1;
    }
#endif
    
    // want to add ( - fact * M R * accel ) to unbalance
    for (int i=0; i<dimension; i++) {
      double val1 = Raccel1(i);
      double val2 = Raccel2(i);	
      
      // perform - fact * M*(R * accel) // remember M a diagonal matrix
      
      double massDerivative = 0.0;
      if (parameterID == 1) {
	double rho = M*2.0/(A*L);
	massDerivative = rho*1.0*L/2.0;
      }
      else if (parameterID == 2) {
	massDerivative = 1.0*A*L/2.0;
      }
      
      val1 *= massDerivative;
      val2 *= massDerivative;
      
      (*theLoadSens)(i) = val1;
      (*theLoadSens)(i+nodalDOF) = val2;
    }	
  }
  return 0;
}

const Vector &
Truss::getResistingForce()
{	
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theVector->Zero();
	return *theVector;
    }
    
    // R = Ku - Pext
    // Ku = F * transformation
    double force = A*theMaterial->getStress();
    for (int i=0; i<numDOF; i++)
	(*theVector)(i) = (*t)(0,i)*force;

    // subtract external load:  Ku - P
    (*theVector) -= *theLoad;


/*
ofstream outputFile1( "tmatrix_p.out", ios::out );
char mystring[100] = "";
for (int iii=0; iii<t->noCols(); iii++) {
	sprintf(mystring,"%25.16e  ",(*t)(0,iii));
	outputFile1 << mystring << endl;
}
outputFile1.close();


ofstream outputFile2( "thevector_p.out", ios::out );
for (iii=0; iii<theVector->Size(); iii++) {
	sprintf(mystring,"%25.16e  ",(*theVector)(iii));
	outputFile2 << mystring << endl;
}
outputFile2.close();

double streess = theMaterial->getStress();
ofstream outputFile3( "stress_p.out", ios::out );
sprintf(mystring,"%25.16e  ",streess);
outputFile3 << mystring << endl;
outputFile3.close();

*/

    return *theVector;
}


const Vector &
Truss::getResistingForceIncInertia()
{	
    this->getResistingForce();
    
    // now include the mass portion
    if (L != 0.0 && M != 0.0) {
	
	// remember we set M = M*A*L/2 in setDoamin()
	const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();	
	
	int dof = dimension;
	int start = numDOF/2;
	for (int i=0; i<dof; i++) {
	    (*theVector)(i) = (*theVector)(i) + M*accel1(i);
	    (*theVector)(i+start) = (*theVector)(i+start) + M*accel2(i);
	}
    }    
    
    return *theVector;
}


int
Truss::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // truss packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(7);
  data(0) = this->getTag();
  data(1) = dimension;
  data(2) = numDOF;
  data(3) = A;
  if (L != 0)
      data(6) = M * 2 / (L*A);
  else
      data(6) = M;
  
  data(4) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();

  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(5) = matDbTag;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING Truss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }	      

  // truss then sends the tags of it's two end nodes

  res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr <<"WARNING Truss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  // finally truss asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr <<"WARNING Truss::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }

  return 0;
}

int
Truss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{

  int res;
  int dataTag = this->getDbTag();

  // truss creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(7);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING Truss::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  dimension = data(1);
  numDOF = data(2);
  A = data(3);
  M = data(6);
  
  // truss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr <<"WARNING Truss::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally truss creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int matClass = data(4);
  int matDb = data(5);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create a new material object
    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr <<"WARNING Truss::recvSelf() - " << this->getTag() 
	<< " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING Truss::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }

  return 0;
}


int
Truss::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // ensure setDomain() worked
    if (L == 0.0)
       return 0;

    // first determine the two end points of the truss based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();    

    if (displayMode == 1 || displayMode == 2) {
	Vector v1(3);
	Vector v2(3);
	for (int i=0; i<dimension; i++) {
	    v1(i) = end1Crd(i)+end1Disp(i)*fact;
	    v2(i) = end2Crd(i)+end2Disp(i)*fact;    
	}
	
	// compute the strain and axial force in the member
	double strain, force;
	if (L == 0.0) {
	    strain = 0.0;
	    force = 0.0;
	} else {
	    strain = this->computeCurrentStrain();
	    theMaterial->setTrialStrain(strain);
	    force = A*theMaterial->getStress();    
	}
    
	if (displayMode == 2) // use the strain as the drawing measure
	    return theViewer.drawLine(v1, v2, strain, strain);	
	else { // otherwise use the axial force as measure
	    return theViewer.drawLine(v1,v2, force, force);
	}
    }
    return 0;
}


void
Truss::Print(OPS_Stream &s, int flag)
{
    // compute the strain and axial force in the member
    double strain, force;
    strain = theMaterial->getStrain();
    force = A * theMaterial->getStress();
    
    if (flag == 0) { // print everything
	s << "Element: " << this->getTag(); 
	s << " type: Truss  iNode: " << connectedExternalNodes(0);
	s << " jNode: " << connectedExternalNodes(1);
	s << " Area: " << A << " Total Mass: " << M*2;
	
	s << " \n\t strain: " << strain;
	s << " axial load: " <<  force;
	if (L != 0.0) {
	  for (int i=0; i<numDOF; i++)
	    (*theVector)(i) = (*t)(0,i)*force;
	  s << " \n\t unbalanced load: " << *theVector;	
	}

	s << " \t Material: " << *theMaterial;
	s << endln;
    } else if (flag == 1) {
	s << this->getTag() << "  " << strain << "  ";
	s << force << endln;
    }
}

double
Truss::computeCurrentStrain(void) const
{
    // NOTE method will not be called if L == 0

    // determine the strain
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();	

    double dLength = 0.0;
    for (int i=0; i<dimension; i++){
	dLength -= (disp2(i)-disp1(i))* (*t)(0,i);
    }

//cerr << "disp1: " << disp1 << endl;
//cerr << "disp2: " << disp2 << endl;
//cerr << "strain from element: " << dLength/L << endl;
//cerr << "***********************" << endl;
  
	// this method should never be called with L == 0
    return dLength/L;
}

double
Truss::computeCurrentStrainRate(void) const
{
    // NOTE method will not be called if L == 0

    // determine the strain
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();	

    double dLength = 0.0;
    for (int i=0; i<dimension; i++){
	dLength -= (vel2(i)-vel1(i))* (*t)(0,i);
    }

    // this method should never be called with L == 0
    return dLength/L;
}

Response*
Truss::setResponse(const char **argv, int argc, Information &eleInfo)
{
  
  //
  // we compare argv[0] for known response types for the Truss
  //

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"axialForce") == 0)
    return new ElementResponse(this, 1, 0.0);

  else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	   strcmp(argv[0],"deformation") == 0)
    return new ElementResponse(this, 2, 0.0);

  // tangent stiffness matrix
  else if (strcmp(argv[0],"stiff") == 0)
    return new ElementResponse(this, 3, *theMatrix);

  // a material quantity    
  else if (strcmp(argv[0],"-material") == 0)
    return theMaterial->setResponse(&argv[1], argc-1, eleInfo);
  
  else
    return 0;
}

int 
Truss::getResponse(int responseID, Information &eleInfo)
{
  switch (responseID) {
    case 1:
      return eleInfo.setDouble(A * theMaterial->getStress());
      
    case 2:
      return eleInfo.setDouble(L * theMaterial->getStrain());
      
    case 3:
      return eleInfo.setMatrix(this->getTangentStiff());

    default:
      return 0;
  }
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
Truss::setParameter (const char **argv, int argc, Information &info)
{
    if (argc < 1)
        return -1;

    // Cross sectional area of the truss itself
    if (strcmp(argv[0],"A") == 0) {
        info.theType = DoubleType;
        return 1;
    }

    // Mass densitity (per unit volume) of the truss itself
    if (strcmp(argv[0],"rho") == 0) {
        info.theType = DoubleType;
        return 2;
    }

    // a material parameter
    if (strcmp(argv[0],"-material") == 0 || strcmp(argv[0],"material") == 0) {
      int ok = theMaterial->setParameter(&argv[1], argc-1, info);
      if (ok < 0)
	return -1;
      else
	return ok + 100;
    } 
    
    // otherwise parameter is unknown for the Truss class
    else
      return -1;
}

int
Truss::updateParameter (int parameterID, Information &info)
{
  switch (parameterID) {
    case -1:
      return -1;
      
    case 1:
		this->M = (M*2.0/(A*L)) * info.theDouble * L/2.0;
        this->A = info.theDouble;
        return 0;

    case 2:
        this->M = info.theDouble * A * L/2.0;
        return 0;

    default:
      if (parameterID >= 100)
	  return theMaterial->updateParameter(parameterID-100, info);
      else
	  return -1;
  }
}
int
Truss::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	// The identifier needs to be passed "downwards" also when it's zero
	if (passedParameterID == 0 || passedParameterID == 1 || passedParameterID == 2) {
		theMaterial->activateParameter(0);
	}

	// If the identifier is non-zero and the parameter belongs to the material
	else if ( passedParameterID > 100) {
		theMaterial->activateParameter(passedParameterID-100);
	}

	return 0;
}
/*
const Vector &
Truss::getResistingForceDerivativeIncInertia(int gradNumber, const Vector &vel,const Vector &acc,double a2,double a3,double a4,double a6,double a7,double a8,double alphaM,double betaK,double dAlphaMdh,double dBetaKdh)
{

	// The term -dPint/dh|u fixed
	// Note that this term is delivered by another method in this class
	// which is also using "theVector".  There the minus sign is implemented.
	Vector dPintdh = this->getResistingForceDerivative(gradNumber);
	theVector->Zero();
	theVector->addVector(1.0,dPintdh,1.0);


	// The term -dM/dh*acc
	double area = 0.0;
	double rho = 0.0;
	if (parameterID == 1) {
		area = 1.0;
		rho = M*2.0/(A*L);
	}
	else if (parameterID == 2) {
		rho = 1.0;
		area = A;
	}
	Matrix *dMdh = 0;
    if (dimension == 1 && numDOF == 2) {
		dMdh = new Matrix(2,2);
		(*dMdh)(0,0) = rho*area*L/2.0; 
		(*dMdh)(1,1) = rho*area*L/2.0;
    }
    else if (dimension == 2 && numDOF == 4) {
		dMdh = new Matrix(4,4);
		(*dMdh)(0,0) = rho*area*L/2.0; 
		(*dMdh)(1,1) = rho*area*L/2.0;
		(*dMdh)(2,2) = rho*area*L/2.0; 
		(*dMdh)(3,3) = rho*area*L/2.0;	
    }
    else if (dimension == 2 && numDOF == 6) {
		dMdh = new Matrix(6,6);
		(*dMdh)(0,0) = rho*area*L/2.0; 
		(*dMdh)(1,1) = rho*area*L/2.0;
		(*dMdh)(3,3) = rho*area*L/2.0;
		(*dMdh)(4,4) = rho*area*L/2.0; 
    }
    else if (dimension == 3 && numDOF == 6) {
		dMdh = new Matrix(6,6);
		(*dMdh)(0,0) = rho*area*L/2.0; 
		(*dMdh)(1,1) = rho*area*L/2.0;
		(*dMdh)(2,2) = rho*area*L/2.0; 
		(*dMdh)(3,3) = rho*area*L/2.0;
		(*dMdh)(4,4) = rho*area*L/2.0; 
		(*dMdh)(5,5) = rho*area*L/2.0;		
    }
    else if (dimension == 3 && numDOF == 12) {
		dMdh = new Matrix(12,12);
		(*dMdh)(0,0) = rho*area*L/2.0; 
		(*dMdh)(1,1) = rho*area*L/2.0;
		(*dMdh)(2,2) = rho*area*L/2.0; 
		(*dMdh)(6,6) = rho*area*L/2.0; 
		(*dMdh)(7,7) = rho*area*L/2.0;
		(*dMdh)(8,8) = rho*area*L/2.0; 
    }
	theVector->addVector(1.0,((*dMdh)^acc),-1.0);


	// The term -(dAlphaMdh*M)*vel
	theVector->addVector(1.0,(dAlphaMdh*((this->getMass())^vel)),-1.0);


	// The term -(alphaM*dMdh)*vel
	theVector->addVector(1.0,(alphaM*((*dMdh)^vel)),-1.0);

	
	// dMdh is used for the last time; delete it
	if (dMdh != 0) {
		delete dMdh;
	}

	// The term -(dBetaKdh*K)*vel
	theVector->addVector(1.0,(dBetaKdh*((this->getTangentStiff())^vel)),-1.0);


	// The term -(betaK*dKdh)*vel
	Matrix &dKdh = *theMatrix;
	if (parameterID == 1) {
		double E = theMaterial->getTangent();
		Matrix &trans = *t;
		dKdh = trans^trans;
		dKdh *= A*E/L;  
	}
	else {
//		double dEdh = theMaterial->getTangentDerivative();
double dEdh = 0.0;
		Matrix &trans = *t;
		dKdh = trans^trans;
		dKdh *= A*dEdh/L;  
	}
	theVector->addVector(1.0,(betaK*(dKdh^vel)),-1.0);


	// Recover displacement sensitivity results
	Vector v(theVector->Size());
	int counter = 0;
	int i;
	for (i=0; i<dimension; i++){
		v(counter) = end1Ptr->getDisplSensitivity(i+1, gradNumber);
		counter++;
	}
	for (i=0; i<dimension; i++){
		v(counter) = end2Ptr->getDisplSensitivity(i+1, gradNumber);
		counter++;
	}


	// The term -M*(a2*v)
	theVector->addVector(1.0,(((this->getMass())^v)*a2),-1.0);


	// The term -C*(a6*v)
	theVector->addVector(1.0,(((this->getDamp())^v)*a6),-1.0);


	// Recover velocity sensitivity results
	Vector vdot(theVector->Size());
	counter = 0;
	for (i=0; i<dimension; i++){
		vdot(counter) = end1Ptr->getVelSensitivity(i+1, gradNumber);
		counter++;
	}
	for (i=0; i<dimension; i++){
		vdot(counter) = end2Ptr->getVelSensitivity(i+1, gradNumber);
		counter++;
	}


	// The term -M*(a3*vdot)
	theVector->addVector(1.0,(((this->getMass())^vdot)*a3),-1.0);


	// The term -C*(a7*vdot)
	theVector->addVector(1.0,(((this->getDamp())^vdot)*a7),-1.0);


	// Recover acceleration sensitivity results
	Vector vdotdot(theVector->Size());
	counter = 0;
	for (i=0; i<dimension; i++){
		vdotdot(counter) = end1Ptr->getAccSensitivity(i+1, gradNumber);
		counter++;
	}
	for (i=0; i<dimension; i++){
		vdotdot(counter) = end2Ptr->getAccSensitivity(i+1, gradNumber);
		counter++;
	}


	// The term -M*(a4*vdotdot)
	theVector->addVector(1.0,(((this->getMass())^vdotdot)*a4),-1.0);


	// The term -C*(a8*vdotdot)
	theVector->addVector(1.0,(((this->getDamp())^vdotdot)*a8),-1.0);

	return *theVector;
}
*/

const Matrix &
Truss::getKiSensitivity(int gradNumber)
{
	theMatrix->Zero();

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;
    Matrix &trans = *t;

	if (parameterID == 0) {
	}
    else if (parameterID == 1) {
		// If cross sectional area is random
	    double E = theMaterial->getInitialTangent();
		stiff = trans^trans;
		stiff *= 1.0*E/L;  
	}
	else if (parameterID == 2) {
		// Nothing here when 'rho' is random
	}
	else {
		double Esens = theMaterial->getInitialTangentSensitivity(gradNumber);
		stiff = trans^trans;
	    stiff *= A*Esens/L;  
	}

	return *theMatrix;
}

const Matrix &
Truss::getMassSensitivity(int gradNumber)
{
	theMatrix->Zero();
	

	double massDerivative = 0.0;
	if (parameterID == 1) {
		double rho = M*2.0/(A*L);
		massDerivative = rho*1.0*L/2.0;
	}
	else if (parameterID == 2) {
		massDerivative = 1.0*A*L/2.0;
	}
    if (dimension == 1 && numDOF == 2) {
		(*theMatrix)(0,0) = massDerivative; 
		(*theMatrix)(1,1) = massDerivative;
    }
    else if (dimension == 2 && numDOF == 4) {
		(*theMatrix)(0,0) = massDerivative; 
		(*theMatrix)(1,1) = massDerivative;
		(*theMatrix)(2,2) = massDerivative; 
		(*theMatrix)(3,3) = massDerivative;	
    }
    else if (dimension == 2 && numDOF == 6) {
		(*theMatrix)(0,0) = massDerivative; 
		(*theMatrix)(1,1) = massDerivative;
		(*theMatrix)(3,3) = massDerivative;
		(*theMatrix)(4,4) = massDerivative; 
    }
    else if (dimension == 3 && numDOF == 6) {
		(*theMatrix)(0,0) = massDerivative; 
		(*theMatrix)(1,1) = massDerivative;
		(*theMatrix)(2,2) = massDerivative; 
		(*theMatrix)(3,3) = massDerivative;
		(*theMatrix)(4,4) = massDerivative; 
		(*theMatrix)(5,5) = massDerivative;		
    }
    else if (dimension == 3 && numDOF == 12) {
		(*theMatrix)(0,0) = massDerivative; 
		(*theMatrix)(1,1) = massDerivative;
		(*theMatrix)(2,2) = massDerivative; 
		(*theMatrix)(6,6) = massDerivative; 
		(*theMatrix)(7,7) = massDerivative;
		(*theMatrix)(8,8) = massDerivative; 
    }

	return *theMatrix;
}

const Vector &
Truss::getResistingForceSensitivity(int gradNumber)
{
	theVector->Zero();


	// Initial declarations
	int i;
	double stressSensitivity, strain, temp1, temp2;
	double coordStressSensitivity = 0.0;

	// Check if a nodal coordinate is random
	Vector nodeParameterID(2);
	nodeParameterID(0) = (double)theNodes[0]->getCrdsSensitivity();
	nodeParameterID(1) = (double)theNodes[1]->getCrdsSensitivity();
	bool nodeCoordIsRandom = false;
	if (nodeParameterID.Norm() != 0.0) {
		nodeCoordIsRandom = true;
	}

	// Determine original geometry and elemet stretch
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();	
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();	

	double x1 = end1Crd(0);
	double x2 = end2Crd(0);
	double y1 = end1Crd(1);	
	double y2 = end2Crd(1);
	double dx = (x2-x1);
	double dy = (y2-y1);
	double L = sqrt(dx*dx + dy*dy);

	double dLength = 0.0;
	for (i=0; i<dimension; i++){
		dLength -= (disp2(i)-disp1(i))* (*t)(0,i);
	}
	double rate = this->computeCurrentStrainRate();

	// Make sure the material is up to date
	strain = dLength/L;
	theMaterial->setTrialStrain(strain);


	// Compute derivative of transformation matrix (assume 4 dofs)
	Vector dtdh(4);
	if (nodeCoordIsRandom) {

		const Vector &end1Crd = theNodes[0]->getCrds();
		const Vector &end2Crd = theNodes[1]->getCrds();	
		double dx = end2Crd(0)-end1Crd(0);
		double dy = end2Crd(1)-end1Crd(1);	
		double L = sqrt(dx*dx + dy*dy);

		for (i=0; i<2; i++) {


			if (i==0) {

				if ( ((int)nodeParameterID(i))==1 ) { // here x1 is random
					temp1 = (-L+dx*dx/L)/(L*L);
					temp2 = dx*dy/(L*L*L);
					dtdh(0) = -temp1;
					dtdh(1) = -temp2;
					dtdh(2) = temp1;
					dtdh(3) = temp2;
				}
				else if ( ((int)nodeParameterID(i))==2 ) { // here y1 is random
					temp1 = (-L+dy*dy/L)/(L*L);
					temp2 = dx*dy/(L*L*L);
					dtdh(0) = -temp2;
					dtdh(1) = -temp1;
					dtdh(2) = temp2;
					dtdh(3) = temp1;
				}
			}
			else if (i==1) {

				if ( ((int)nodeParameterID(i))==1 ) { // here x2 is random
					temp1 = (L-dx*dx/L)/(L*L);
					temp2 = -dx*dy/(L*L*L);
					dtdh(0) = -temp1;
					dtdh(1) = -temp2;
					dtdh(2) = temp1;
					dtdh(3) = temp2;
				}
				else if ( ((int)nodeParameterID(i))==2 ) { // here y2 is random
					temp1 = (L-dy*dy/L)/(L*L);
					temp2 = -dx*dy/(L*L*L);
					dtdh(0) = -temp2;
					dtdh(1) = -temp1;
					dtdh(2) = temp2;
					dtdh(3) = temp1;
				}

			}
		}
	}



	// Determine stress sensitivity 
	if (nodeCoordIsRandom) {

		double strainSensitivity, materialTangent;
		double dLengthDerivative = 0.0;
		for (i=0; i<dimension; i++){
			dLengthDerivative -= (disp2(i)-disp1(i))* dtdh(i);
		}

		for (i=0; i<2; i++) {
			if (i==0) {
				if ( ((int)nodeParameterID(i))==1 ) {		// here x1 is random
					materialTangent = theMaterial->getTangent();
					strainSensitivity = (dLengthDerivative*L+dLength/L*(x2-x1))/(L*L);
					coordStressSensitivity = materialTangent * strainSensitivity;
				}
				else if ( ((int)nodeParameterID(i))==2 ) {	// here y1 is random
					materialTangent = theMaterial->getTangent();
					strainSensitivity = (dLengthDerivative*L+dLength/L*(y2-y1))/(L*L);
					coordStressSensitivity = materialTangent * strainSensitivity;
				}
			}
			else {
				if ( ((int)nodeParameterID(i))==1 ) {		// here x2 is random
					materialTangent = theMaterial->getTangent();
					strainSensitivity = (dLengthDerivative*L-dLength/L*(x2-x1))/(L*L);
					coordStressSensitivity = materialTangent * strainSensitivity;
				}
				else if ( ((int)nodeParameterID(i))==2 ) {	// here y2 is random
					materialTangent = theMaterial->getTangent();
					strainSensitivity = (dLengthDerivative*L-dLength/L*(y2-y1))/(L*L);
					coordStressSensitivity = materialTangent * strainSensitivity;
				}
			}
		}
	}

	stressSensitivity = coordStressSensitivity + theMaterial->getStressSensitivity(gradNumber,true);



	// Compute sensitivity depending on 'parameter'
	double stress = theMaterial->getStress();
	if( parameterID == 1 ) {			// Cross-sectional area
		for (i=0; i<numDOF; i++)
			(*theVector)(i) = (-1)*(*t)(0,i)*stress +
			(-1)*(*t)(0,i)*A*stressSensitivity;
	}
	else {		// Density, material parameter or nodal coordinate
		for (i=0; i<numDOF; i++) {	
			(*theVector)(i) = (-1)*(*t)(0,i)*A*stressSensitivity
				            + (-1)* dtdh(i) *A*stress;
		}
	}



	// subtract external load sensitivity
	if (theLoadSens == 0) {
		theLoadSens = new Vector(numDOF);
	}
	(*theVector) -= *theLoadSens;

	return *theVector;
}








int
Truss::commitSensitivity(int gradNumber, int numGrads)
{
	// Initial declarations
	double dLdh = 0.0; 
	int i; 
	double strainSensitivity, temp1, temp2;


	// Check if a nodal coordinate is random
	Vector nodeParameterID(2);
	nodeParameterID(0) = (double)theNodes[0]->getCrdsSensitivity();
	nodeParameterID(1) = (double)theNodes[1]->getCrdsSensitivity();
	bool nodeCoordIsRandom = false;
	if (nodeParameterID.Norm() != 0.0) {
		nodeCoordIsRandom = true;
	}


	// Displacement difference between the two ends
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();	
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();	
	double dLength = 0.0;
	for (i=0; i<dimension; i++){
		dLength -= (disp2(i)-disp1(i))* (*t)(0,i);
	}

	// Displacement sensitivity difference between the two ends
	double sens1;
	double sens2;
	double dSensitivity = 0.0;
	for (i=0; i<dimension; i++){
		sens1 = theNodes[0]->getDispSensitivity(i+1, gradNumber);
		sens2 = theNodes[1]->getDispSensitivity(i+1, gradNumber);
		dSensitivity -= (sens2-sens1)* (*t)(0,i);
	}

	double x1 = end1Crd(0);
	double x2 = end2Crd(0);
	double y1 = end1Crd(1);	
	double y2 = end2Crd(1);
	double dx = (x2-x1);
	double dy = (y2-y1);
	double L = sqrt(dx*dx + dy*dy);


	// Compute derivative of transformation matrix (assume 4 dofs)
	Vector dtdh(4);
	if (nodeCoordIsRandom) {

		const Vector &end1Crd = theNodes[0]->getCrds();
		const Vector &end2Crd = theNodes[1]->getCrds();	
		double dx = end2Crd(0)-end1Crd(0);
		double dy = end2Crd(1)-end1Crd(1);	
		double L = sqrt(dx*dx + dy*dy);

		for (i=0; i<2; i++) {


			if (i==0) {

				if ( ((int)nodeParameterID(i))==1 ) { // here x1 is random
					temp1 = (-L+dx*dx/L)/(L*L);
					temp2 = dx*dy/(L*L*L);
					dtdh(0) = -temp1;
					dtdh(1) = -temp2;
					dtdh(2) = temp1;
					dtdh(3) = temp2;
				}
				else if ( ((int)nodeParameterID(i))==2 ) { // here y1 is random
					temp1 = (-L+dy*dy/L)/(L*L);
					temp2 = dx*dy/(L*L*L);
					dtdh(0) = -temp2;
					dtdh(1) = -temp1;
					dtdh(2) = temp2;
					dtdh(3) = temp1;
				}
			}
			else if (i==1) {

				if ( ((int)nodeParameterID(i))==1 ) { // here x2 is random
					temp1 = (L-dx*dx/L)/(L*L);
					temp2 = -dx*dy/(L*L*L);
					dtdh(0) = -temp1;
					dtdh(1) = -temp2;
					dtdh(2) = temp1;
					dtdh(3) = temp2;
				}
				else if ( ((int)nodeParameterID(i))==2 ) { // here y2 is random
					temp1 = (L-dy*dy/L)/(L*L);
					temp2 = -dx*dy/(L*L*L);
					dtdh(0) = -temp2;
					dtdh(1) = -temp1;
					dtdh(2) = temp2;
					dtdh(3) = temp1;
				}

			}
		}
	}











	// Determine strain sensitivity
	if (nodeCoordIsRandom) {

		double dLengthDerivative = 0.0;
		for (i=0; i<dimension; i++){
			dLengthDerivative -= (disp2(i)-disp1(i))* dtdh(i);
		}

		for (i=0; i<2; i++) {
			if (i==0) {
				if ( ((int)nodeParameterID(i))==1 ) {		// here x1 is random
					strainSensitivity = dLength/(L*L*L)*(x2-x1)
									  + dSensitivity/L
									  + dLengthDerivative/L;
				}
				else if ( ((int)nodeParameterID(i))==2 ) {	// here y1 is random
					strainSensitivity = dLength/(L*L*L)*(y2-y1)
									  + dSensitivity/L
									  + dLengthDerivative/L;
				}
			}
			else {
				if ( ((int)nodeParameterID(i))==1 ) {		// here x2 is random
					strainSensitivity = -dLength/(L*L*L)*(x2-x1)
									  + dSensitivity/L
									  + dLengthDerivative/L;
				}
				else if ( ((int)nodeParameterID(i))==2 ) {	// here y2 is random
					strainSensitivity = -dLength/(L*L*L)*(y2-y1)
									  + dSensitivity/L
									  + dLengthDerivative/L;
				}
			}
		}
	}
	else {
		strainSensitivity = dSensitivity/L;
	}

	
	// Pass it down to the material
	theMaterial->commitSensitivity(strainSensitivity, gradNumber, numGrads);

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////
