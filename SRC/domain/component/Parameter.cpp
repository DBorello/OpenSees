/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
// $Date: 2006-09-05 20:21:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/Parameter.cpp,v $

#include <classTags.h>
#include <Parameter.h>
#include <DomainComponent.h>

Parameter::Parameter(int passedTag,
		     DomainComponent *parentObject,
		     const char **argv, int argc)
  :TaggedObject(passedTag),
   theObjects(0), parameterID(0),
   numObjects(0), maxNumObjects(0)
{
  int ok = -1;

  maxNumObjects = initialSize;

  theObjects = new MovableObject *[maxNumObjects];
  parameterID = new int[maxNumObjects];

  if (parentObject != 0)
    ok = parentObject->setParameter(argv, argc, *this);
  
  if (ok < 0)
    opserr << "Parameter::Parameter "<< this->getTag() <<" -- unable to set parameter" << endln;
}

Parameter::Parameter(const Parameter &param):
  TaggedObject(param.getTag())
{
  theInfo = param.theInfo;
  numObjects = param.numObjects;
  maxNumObjects = param.maxNumObjects;

  theObjects = new MovableObject *[maxNumObjects];
  parameterID = new int[maxNumObjects];

  for (int i = 0; i < numObjects; i++) {
    theObjects[i] = param.theObjects[i];
    parameterID[i] = param.parameterID[i];
  }
}

Parameter::~Parameter()
{
  if (theObjects != 0)
    delete [] theObjects;

  if (parameterID != 0)
    delete [] parameterID;
}

int
Parameter::addObject(DomainComponent *parentObject,
		     const char **argv, int argc)
{
  return (parentObject != 0) ?
    parentObject->setParameter(argv, argc, *this) : -1;
}

int
Parameter::update(int newValue)
{
  theInfo.theInt = newValue;

  int ok = 0;

  for (int i = 0; i < numObjects; i++)
    ok += theObjects[i]->updateParameter(parameterID[i], theInfo);

  return ok;
}

int
Parameter::update(double newValue)
{
  theInfo.theDouble = newValue;

  int ok = 0;

  for (int i = 0; i < numObjects; i++)
    ok += theObjects[i]->updateParameter(parameterID[i], theInfo);

  return ok;
}

int
Parameter::activate(bool active)
{
  int ok = 0;
  for (int i = 0; i < numObjects; i++) {
    ok += theObjects[i]->activateParameter(active ? parameterID[i] : 0);
  }

  return ok;
}

void
Parameter::Print(OPS_Stream &s, int flag)  
{
  s << "Parameter, tag = " << this->getTag() << endln;
  //s << "\tparameterID = " << parameterID << endln;
}

int
Parameter::addObject(int paramID, MovableObject *object)
{
  //opserr << "Parameter::addObject " << numObjects << endln;

  if (numObjects == maxNumObjects) {
    maxNumObjects += expandSize;
    MovableObject **newObjects = new MovableObject*[maxNumObjects];
    int *newParameterID = new int[maxNumObjects];

    for (int i = 0; i < numObjects; i++) {
      newObjects[i] = theObjects[i];
      newParameterID[i] = parameterID[i];
    }

    delete [] theObjects;
    delete [] parameterID;

    theObjects = newObjects;
    parameterID = newParameterID;
  }

  parameterID[numObjects] = paramID;
  theObjects[numObjects] = object;
  numObjects++;

  return 0;
}
