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
                                                                        
// $Revision: 1.5 $
// $Date: 2001-07-11 23:09:00 $
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

// Element(int tag, int noExtNodes);
// 	constructor that takes the element's unique tag and the number
//	of external nodes for the element.

Element::Element(int tag, int cTag) 
  :DomainComponent(tag, cTag), Ki(0)
{
    // does nothing
}


Element::~Element() 
{
  if (Ki != 0)
    delete Ki;
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
Element::setKi(void) {
  if (Ki != 0)
    delete Ki;
  Ki = new Matrix(this->getTangentStiff());
  if (Ki == 0) {
    cerr << "Element::setKi(void) - out of memory for element " << this->getTag() << endl;
    return -1;
  }
  return 0;
}
      
  
const Matrix &
Element::getKi(void) {
  if (Ki != 0)
    return *Ki;
  else {
    this->setKi();
    if (Ki == 0)
      return this->getTangentStiff();
    else
      return *Ki;
  }
}



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
