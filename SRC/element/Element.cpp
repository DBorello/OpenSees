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
                                                                        
// $Revision: 1.8 $
// $Date: 2002-12-06 20:26:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Element.cpp,v $
                                                                        
                                                                        
// File: ~/model/Element.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for Element.
// Element is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// The interface:
//

#include "Element.h"
#include <Renderer.h>
#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Domain.h>

Matrix **Element::theMatrices; 
Vector **Element::theVectors1; 
Vector **Element::theVectors2; 
int  Element::numMatrices(0);

// Element(int tag, int noExtNodes);
// 	constructor that takes the element's unique tag and the number
//	of external nodes for the element.

Element::Element(int tag, int cTag) 
  :DomainComponent(tag, cTag), alphaM(0), betaK(0), betaK0(0), index(-1)

{
    // does nothing
}


Element::~Element() 
{

}

int
Element::update(void)
{
    return 0;
}

int
Element::revertToStart(void)
{
    return 0;
}


int
Element::setRayleighDampingFactors(double alpham, double betak, double betak0)
{
  alphaM = alpham;
  betaK  = betak;
  betaK0 = betak0;

  // check that memory has been allocated to store compute/return
  // damping matrix & residual force calculations
  if (index == -1) {
    int numDOF = this->getNumDOF();
    for (int i=0; i<numMatrices; i++) {
      Matrix *aMatrix = theMatrices[i];
      if (aMatrix->noRows() == numDOF) {
	index = i;
	i = numMatrices;
      }
    }
    if (index == -1) {
      Matrix **nextMatrices = new Matrix *[numMatrices+1];
      if (nextMatrices == 0) {
	g3ErrorHandler->fatal("Element::getTheMatrix - out of memory");
      }
	  int j;
      for (j=0; j<numMatrices; j++)
	nextMatrices[j] = theMatrices[j];
      Matrix *theMatrix = new Matrix(numDOF, numDOF);
      if (theMatrix == 0) {
	g3ErrorHandler->fatal("Element::getTheMatrix - out of memory");
      }
      nextMatrices[numMatrices] = theMatrix;

      Vector **nextVectors1 = new Vector *[numMatrices+1];
      Vector **nextVectors2 = new Vector *[numMatrices+1];
      if (nextVectors1 == 0 || nextVectors2 == 0) {
	g3ErrorHandler->fatal("Element::getTheVector - out of memory");
      }

      for (j=0; j<numMatrices; j++) {
	nextVectors1[j] = theVectors1[j];
	nextVectors2[j] = theVectors2[j];
      }
	
      Vector *theVector1 = new Vector(numDOF);
      Vector *theVector2 = new Vector(numDOF);
      if (theVector1 == 0 || theVector2 == 0) {
	g3ErrorHandler->fatal("Element::getTheVector - out of memory");
      }
      nextVectors1[numMatrices] = theVector1;
      nextVectors2[numMatrices] = theVector2;

      if (numMatrices != 0) {
	delete [] theMatrices;
	delete [] theVectors1;
	delete [] theVectors2;
      }
      index = numMatrices;
      numMatrices++;
      theMatrices = nextMatrices;
      theVectors1 = nextVectors1;
      theVectors2 = nextVectors2;
    }
  }

  return 0;
}

const Matrix &
Element::getDamp(void) 
{
  if (index  == -1) {
    this->setRayleighDampingFactors(0.0, 0.0, 0.0);
  }

  // now compute the damping matrix
  Matrix *theMatrix = theMatrices[index]; 
  theMatrix->Zero();
  if (alphaM != 0.0)
    theMatrix->addMatrix(0.0, this->getMass(), alphaM);
  if (betaK != 0.0)
    theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);      
  if (betaK0 != 0.0)
    theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);      

  // return the computed matrix
  return *theMatrix;
}



const Matrix &
Element::getMass(void)
{
  if (index  == -1) {
    this->setRayleighDampingFactors(0.0, 0.0, 0.0);
  }

  // zero the matrix & return it
  Matrix *theMatrix = theMatrices[index]; 
  theMatrix->Zero();
  return *theMatrix;
}

const Vector &
Element::getResistingForceIncInertia(void) 
{
  if (index == -1) {
    this->setRayleighDampingFactors(0.0, 0.0, 0.0);
  }

  Matrix *theMatrix = theMatrices[index]; 
  Vector *theVector = theVectors2[index];
  Vector *theVector2 = theVectors1[index];

  //
  // perform: R = P(U) - Pext(t);
  //

  (*theVector) = this->getResistingForce();

  //
  // perform: R = R - M * a
  //

  int loc = 0;
  Node **theNodes = this->getNodePtrs();
  int numNodes = this->getNumExternalNodes();

  int i;
  for (i=0; i<numNodes; i++) {
    const Vector &acc = theNodes[i]->getAccel();
    for (int i=0; i<acc.Size(); i++) {
      (*theVector2)(loc++) = acc(i);
    }
  }
  theVector->addMatrixVector(1.0, this->getMass(), *theVector2, +1.0);

  //
  // perform: R = R + (alphaM * M + betaK0 * K0 + betaK * K) * v
  //            = R + D * v
  //

  // determine the vel vector from ele nodes
  loc = 0;
  for (i=0; i<numNodes; i++) {
    const Vector &vel = theNodes[i]->getTrialVel();
    for (int i=0; i<vel.Size(); i++) {
      (*theVector2)(loc++) = vel[i];
    }
  }

  // now compute the damping matrix
  theMatrix->Zero();
  if (alphaM != 0.0)
    theMatrix->addMatrix(0.0, this->getMass(), alphaM);
  if (betaK != 0.0)
    theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);      
  if (betaK0 != 0.0)
    theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);      

  // finally the D * v
  theVector->addMatrixVector(1.0, *theMatrix, *theVector2, 1.0);

  return *theVector;
}


const Vector &
Element::getRayleighDampingForces(void) 
{

  if (index == -1) {
    this->setRayleighDampingFactors(0.0, 0.0, 0.0);
  }

  Matrix *theMatrix = theMatrices[index]; 
  Vector *theVector = theVectors2[index];
  Vector *theVector2 = theVectors1[index];

  //
  // perform: R = (alphaM * M + betaK0 * K0 + betaK * K) * v
  //            = D * v
  //

  // determine the vel vector from ele nodes
  Node **theNodes = this->getNodePtrs();
  int numNodes = this->getNumExternalNodes();
  int loc = 0;
  for (int i=0; i<numNodes; i++) {
    const Vector &vel = theNodes[i]->getTrialVel();
    for (int i=0; i<vel.Size(); i++) {
      (*theVector2)(loc++) = vel[i];
    }
  }

  // now compute the damping matrix
  theMatrix->Zero();
  if (alphaM != 0.0)
    theMatrix->addMatrix(0.0, this->getMass(), alphaM);
  if (betaK != 0.0)
    theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);      
  if (betaK0 != 0.0)
    theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);      

  // finally the D * v
  theVector->addMatrixVector(0.0, *theMatrix, *theVector2, 1.0);

  return *theVector;
}

int 
Element::addLoad(ElementalLoad *theLoad, double loadFactor) {
  return 0;
}

/*
int 
Element::addInertiaLoadToUnbalance(const Vector &accel)
{
  // some vectors to hold the load increment and RV
  int ndof = this->getNumDOF();
  Vector load(ndof);
  Vector RV(ndof);

  // 
  // for each node we will add it's R*accel contribution to RV
  //

  const ID &theNodes = this->getExternalNodes();
  int numNodes = theNodes.Size();
  int loc = 0;
  Domain *theDomain = this->getDomain();
  for (int i=0; i<numNodes; i++) {
    Node *theNode = theDomain->getNode(theNodes(i));
    if (theNode == 0)
      return -1;
    else {
      int numNodeDOF = theNode->getNumberDOF();
      const Vector &nodeRV = theNode->getRV(accel);
      for (int j=0; j<numNodeDOF; j++)
#ifdef _G3DEBUG
 	if (loc<ndof)
#endif
	  RV(loc++) = nodeRV(j);
    }
  }

  //
  // now we determine - M * R * accel
  //
  const Matrix &mass = this->getMass();
  load = mass * RV;
  load *= -1.0;

  return this->addLoad(load);
}
*/

bool
Element::isSubdomain(void)
{
    return false;
}

Response*
Element::setResponse(char **argv, int argc, Information &eleInfo)
{
	return 0;
}

int
Element::getResponse(int responseID, Information &eleInformation)
{
    return -1;
}
