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
// $Date: 2001-06-14 06:21:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/NodalLoad.cpp,v $
                                                                        
                                                                        
// File: ~/domain/node/NodalLoad.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of NodalLoad

#include <NodalLoad.h>
#include <stdlib.h>
#include <Domain.h>
#include <Channel.h>
#include <Information.h>

NodalLoad::NodalLoad(int theClasTag)
:Load(0,theClasTag), 
 myNode(0), myNodePtr(0), load(0), konstant(false)
{
    // constructor for a FEM_ObjectBroker
}

NodalLoad::NodalLoad(int tag, int node, int theClassTag)
:Load(tag,theClassTag), 
 myNode(node), myNodePtr(0), load(0), konstant(false)
{
}

NodalLoad::NodalLoad(int tag, int node, const Vector &theLoad, bool isLoadConstant)
:Load(tag, LOAD_TAG_NodalLoad), 
 myNode(node), myNodePtr(0), load(0), konstant(isLoadConstant)
{
    load = new Vector(theLoad);    

    if (load == 0) {
	cerr << "FATAL NodalLoad::NodalLoad(int node, const Vector &theLoad) -";
	cerr << " ran out of memory for load on Node " << node << endl;
	exit(-1);
    }
}


// ~NodalLoad()
// 	destructor

NodalLoad::~NodalLoad()
{
    if (load != 0)
	delete load;
}

void 
NodalLoad::setDomain(Domain *newDomain)
{
    // first get myNodePtr
  if (newDomain == 0)
    return;

  // invoke the ancestor class method
  this->DomainComponent::setDomain(newDomain);    
}

int 
NodalLoad::getNodeTag(void) const
{
    return myNode;
}


void
NodalLoad::applyLoad(double loadFactor)
{
    if (myNodePtr == 0) {
      Domain *theDomain=this->getDomain();
      if ((theDomain == 0) || 
	  (myNodePtr = theDomain->getNode(myNode)) == 0) {
	cerr << "WARNING NodalLoad::applyLoad() - No associated Node node " ;
	cerr << " for NodalLoad " << *this;
	return;
      }
    }

    // add the load times the loadfactor to nodal unbalanced load
    if (konstant == false)
	myNodePtr->addUnbalancedLoad(*load,loadFactor);
    else
	myNodePtr->addUnbalancedLoad(*load,1.0);	
} 



int 
NodalLoad::sendSelf(int cTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    ID data(5);
    data(0) = this->getTag();    
    data(1) = myNode;
    if (load != 0)
	data(2) = load->Size();
    else
	data(2) = 0;
    data(3) = konstant;
    data(4) = this->getLoadPatternTag();
    
    int result = theChannel.sendID(dataTag, cTag, data);
    if (result < 0) {
	cerr << "NodalLoad::sendSelf - failed to send data\n";
	return result;
    }

    if (load != 0){
	int result = theChannel.sendVector(dataTag, cTag, *load);
	if (result < 0) {
	    cerr << "NodalLoad::sendSelf - failed to Load data\n";
	    return result;
	}
    }    

    return 0;
}

int 
NodalLoad::recvSelf(int cTag, Channel &theChannel, 
		    FEM_ObjectBroker &theBroker)
{	
    int result;
    int dataTag = this->getDbTag();
    ID data(5);
    result = theChannel.recvID(dataTag, cTag, data);
    if (result < 0) {
      cerr << "NodalLoad::recvSelf() - failed to recv data\n";
      return result;
    }    
    this->setTag(data(0));
    myNode = data(1);
    int loadSize = data(2);
    konstant = data(3);
    this->setLoadPatternTag(data(4));
    if (loadSize != 0) {
	load = new Vector(data(2));
	result = theChannel.recvVector(dataTag, cTag, *load);
	if (result < 0) {
	  cerr << "NodalLoad::recvSelf() - failed to recv load\n";
	  return result;
	}    
    }

    return 0;
}


void
NodalLoad::Print(ostream &s, int flag)
{
     s << "Nodal Load: " << myNode;
     if (load != 0)
	 s << " load : " << *load;
}


int
NodalLoad::setParameter(char **argv, int argc, Information &info)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"direction1") == 0) {
		info.theType = DoubleType;
		return 1;
	}
	if (strcmp(argv[0],"direction2") == 0) {
		info.theType = DoubleType;
		return 2;
	}
	if (strcmp(argv[0],"direction3") == 0) {
		info.theType = DoubleType;
		return 3;
	}
	else
		return -1;
}

int
NodalLoad::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		(*load)(0) = info.theDouble;
		return 0;
	case 2:
		(*load)(1) = info.theDouble;
		return 0;
	case 3:
		(*load)(2) = info.theDouble;
		return 0;
	default:
		return -1;
	}
}
