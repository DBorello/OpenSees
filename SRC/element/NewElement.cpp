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
// $Date: 2001-08-16 21:37:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/NewElement.cpp,v $
                                                                        
// Written: fmk 
// Created: 08/01
//
// Description: This file contains the implementation for the NewElement class.
//
// What: "@(#) NewElement.cpp, revA"

#include "NewElement.h"
#include <Information.h>
#include <ElementResponse.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>


static int NUM_NODE =2;
static int NUM_DOF  =4;

// constructors:
NewElement::NewElement(int tag)
 :Element(tag,ELE_TAG_NewElement),     
  connectedExternalNodes(2), theMatrix(NUM_DOF, NUM_DOF), theVector(NUM_DOF)
{

}

NewElement::NewElement()
 :Element(0,ELE_TAG_NewElement),     
  connectedExternalNodes(2), theMatrix(NUM_DOF, NUM_DOF), theVector(NUM_DOF)
{

}

//  destructor:
NewElement::~NewElement()
{

}


int
NewElement::getNumExternalNodes(void) const
{
    return NUM_NODE;
}

const ID &
NewElement::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

int
NewElement::getNumDOF(void) 
{
    return NUM_DOF;
}


void
NewElement::setDomain(Domain *theDomain)
{
    // call the base class method
    this->DomainComponent::setDomain(theDomain);

}   	 


int
NewElement::commitState()
{
  return 0;
}

int
NewElement::revertToLastCommit()
{
  return 0;
}

int
NewElement::revertToStart()
{
  return 0;
}

int
NewElement::update(void)
{
  return 0;
}


const Matrix &
NewElement::getTangentStiff(void)
{
  return theMatrix;
}

const Matrix &
NewElement::getSecantStiff(void)
{
  return theMatrix;
}
    
const Matrix &
NewElement::getDamp(void)
{
  return theMatrix;
}


const Matrix &
NewElement::getMass(void)
{ 
  return theMatrix;
}



void 
NewElement::zeroLoad(void)
{
  return;
}

int 
NewElement::addLoad(const Vector &addP)
{
  return 0;
}

int 
NewElement::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
NewElement::getResistingForce()
{	
  return theVector;
}


const Vector &
NewElement::getResistingForceIncInertia()
{	
  return theVector;
}


int
NewElement::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
NewElement::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}


int
NewElement::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  return 0;
}


void
NewElement::Print(ostream &s, int flag)
{
  return;
}


Response*
NewElement::setResponse(char **argv, int argc, Information &eleInfo)
{
  return 0;
}


int 
NewElement::getResponse(int responseID, Information &eleInfo)
{
  return -1;
}


int
NewElement::setParameter (char **argv, int argc, Information &info)
{
  return -1;
}
    

int
NewElement::updateParameter (int parameterID, Information &info)
{
  return -1;
}
