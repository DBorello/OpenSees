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
                                                                        
// $Revision: 1.4 $
// $Date: 2002-06-07 22:05:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/LoadPattern.cpp,v $
                                                                        
                                                                        
// File: ~/domain/pattern/LoadPattern.C
//
// Written: fmk 07/99
// Revised:
//
// Purpose: This file contains the method definitions for class LoadPattern.
// LoadPattern is a container class.
//
// The LoadPattern interface:


#include <LoadPattern.h>
#include <stdlib.h>
#include <ID.h>
#include <TimeSeries.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <SP_Constraint.h>
#include <ArrayOfTaggedObjects.h>
#include <ElementalLoadIter.h>
#include <NodalLoadIter.h>
#include <SingleDomSP_Iter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <G3Globals.h>

LoadPattern::LoadPattern(int tag, int clasTag)
:DomainComponent(tag,clasTag),
 loadFactor(0), isConstant(1), 
 theSeries(0), 
 currentGeoTag(0), lastGeoSendTag(-1),
 theNodalLoads(0), theElementalLoads(0), theSPs(0),
 theNodIter(0), theEleIter(0), theSpIter(0)
{
    // constructor for subclass
    theNodalLoads = new ArrayOfTaggedObjects(32);
    theElementalLoads = new ArrayOfTaggedObjects(32);
    theSPs = new ArrayOfTaggedObjects(32);

    if (theNodalLoads == 0 || theElementalLoads == 0 || theSPs == 0) {
	cerr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }    

    theEleIter = new ElementalLoadIter(theElementalLoads);    
    theNodIter = new NodalLoadIter(theNodalLoads);
    theSpIter = new SingleDomSP_Iter(theSPs);
    
    if (theEleIter == 0 || theNodIter == 0 || theSpIter == 0) {
	cerr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }  
// AddingSensitivity:BEGIN /////////////////////////////
	randomLoads = 0;
// AddingSensitivity:END ///////////////////////////////
}


LoadPattern::LoadPattern()
:DomainComponent(0,PATTERN_TAG_LoadPattern),
 loadFactor(0), isConstant(1), 
 theSeries(0), 
 currentGeoTag(0), lastGeoSendTag(-1),
 dbSPs(0), dbNod(0), dbEle(0), 
 theNodalLoads(0), theElementalLoads(0), theSPs(0),
 theNodIter(0), theEleIter(0), theSpIter(0)
{
    theNodalLoads = new ArrayOfTaggedObjects(32);
    theElementalLoads = new ArrayOfTaggedObjects(32);
    theSPs = new ArrayOfTaggedObjects(32);

    if (theNodalLoads == 0 || theElementalLoads == 0 || theSPs == 0) {
	cerr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }    

    theEleIter = new ElementalLoadIter(theElementalLoads);    
    theNodIter = new NodalLoadIter(theNodalLoads);
    theSpIter = new SingleDomSP_Iter(theSPs);
    
    if (theEleIter == 0 || theNodIter == 0 || theSpIter == 0) {
	cerr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }
// AddingSensitivity:BEGIN /////////////////////////////
	randomLoads = 0;
// AddingSensitivity:END ///////////////////////////////
}


LoadPattern::LoadPattern(int tag)
:DomainComponent(tag,PATTERN_TAG_LoadPattern),
 loadFactor(0), isConstant(1), 
 theSeries(0), 
 currentGeoTag(0), lastGeoSendTag(-1),
 dbSPs(0), dbNod(0), dbEle(0), 
 theNodalLoads(0), theElementalLoads(0), theSPs(0),
 theNodIter(0), theEleIter(0), theSpIter(0)
{
    theNodalLoads = new ArrayOfTaggedObjects(32);
    theElementalLoads = new ArrayOfTaggedObjects(32);
    theSPs = new ArrayOfTaggedObjects(32);

    if (theNodalLoads == 0 || theElementalLoads == 0 || theSPs == 0) {
	cerr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }    

    theEleIter = new ElementalLoadIter(theElementalLoads);    
    theNodIter = new NodalLoadIter(theNodalLoads);
    theSpIter = new SingleDomSP_Iter(theSPs);
    
    if (theEleIter == 0 || theNodIter == 0 || theSpIter == 0) {
	cerr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }
// AddingSensitivity:BEGIN /////////////////////////////
	randomLoads = 0;
// AddingSensitivity:END ///////////////////////////////
}

    
// ~LoadPattern()
//	destructor

LoadPattern::~LoadPattern()
{
    if (theSeries != 0)
	delete theSeries;
    
    if (theNodalLoads != 0)
      delete theNodalLoads;

    if (theElementalLoads != 0)
      delete theElementalLoads;

    if (theSPs != 0)
      delete theSPs;

    if (theEleIter != 0)
      delete theEleIter;

    if (theNodIter != 0)
      delete theNodIter;

    if (theSpIter != 0)
      delete theSpIter;

// AddingSensitivity:BEGIN /////////////////////////////
	if (randomLoads != 0)
		delete randomLoads;
// AddingSensitivity:END ///////////////////////////////
}


void
LoadPattern::setTimeSeries(TimeSeries *theTimeSeries)
{
    // invoke the destructor on the old TimeSeries
    if (theSeries != 0)
	delete theSeries;

    // set the pointer to the new series object
    theSeries = theTimeSeries;
}


void
LoadPattern::setDomain(Domain *theDomain)
{
    // if subclass does not implement .. check for 0 pointer
    if (theNodalLoads != 0) {
	NodalLoad *nodLoad;
	NodalLoadIter &theNodalIter = this->getNodalLoads();
	while ((nodLoad = theNodalIter()) != 0)
	    nodLoad->setDomain(theDomain);
    
	ElementalLoad *eleLoad;
	ElementalLoadIter &theElementalIter = this->getElementalLoads();
	while ((eleLoad = theElementalIter()) != 0)
	    eleLoad->setDomain(theDomain);    

	SP_Constraint *theSP;
	SP_ConstraintIter &theSpConstraints = this->getSPs();
	while ((theSP = theSpConstraints()) != 0)
	    theSP->setDomain(theDomain);
    }

    // now we set this load patterns domain
    this->DomainComponent::setDomain(theDomain);
}



bool
LoadPattern::addNodalLoad(NodalLoad *load)
{
    Domain *theDomain = this->getDomain();    

    bool result = theNodalLoads->addComponent(load);
    if (result == true) {
	if (theDomain != 0)
	    load->setDomain(theDomain);
	load->setLoadPatternTag(this->getTag());
	currentGeoTag++;
    } else  
	cerr << "WARNING: LoadPattern::addNodalLoad() - load could not be added\n";

    return result;
}

bool
LoadPattern::addElementalLoad(ElementalLoad *load)
{
    Domain *theDomain = this->getDomain();
    
    bool result = theElementalLoads->addComponent(load);
    if (result == true) {
	if (theDomain != 0)
	    load->setDomain(theDomain);
	load->setLoadPatternTag(this->getTag());
	currentGeoTag++;
    } else
	cerr << "WARNING: LoadPattern::addElementalLoad() - load could not be added\n";	
    
    return result;
}

bool
LoadPattern::addSP_Constraint(SP_Constraint *theSp)
{
    Domain *theDomain = this->getDomain();
    
    bool result = theSPs->addComponent(theSp);
    if (result == true) {
	if (theDomain != 0)
	    theSp->setDomain(theDomain);
	theSp->setLoadPatternTag(this->getTag());
	currentGeoTag++;
    } else
	cerr << "WARNING: LoadPattern::addSP_Constraint() - load could not be added\n";
    return result;
}


NodalLoadIter  &
LoadPattern::getNodalLoads(void)
{
    theNodIter->reset();
    return *theNodIter;
}
    
ElementalLoadIter &
LoadPattern::getElementalLoads(void)  
{
    theEleIter->reset();
    return *theEleIter;    
}

SP_ConstraintIter &
LoadPattern::getSPs(void)  
{
    theSpIter->reset();
    return *theSpIter;    
}

void
LoadPattern::clearAll(void)
{
    theElementalLoads->clearAll();
    theNodalLoads->clearAll();
    theSPs->clearAll();
    currentGeoTag++;
}

NodalLoad *
LoadPattern::removeNodalLoad(int tag)
{
    TaggedObject *obj = theNodalLoads->removeComponent(tag);
    if (obj == 0)
	return 0;
    NodalLoad *result = (NodalLoad *)obj;
    result->setDomain(0);
    currentGeoTag++;
    return result;
}

ElementalLoad *
LoadPattern::removeElementalLoad(int tag)
{
    TaggedObject *obj = theElementalLoads->removeComponent(tag);
    if (obj == 0)
	return 0;
    ElementalLoad *result = (ElementalLoad *)obj;
    result->setDomain(0);
    currentGeoTag++;
    return result;    
}    

SP_Constraint *
LoadPattern::removeSP_Constraint(int tag)
{
    TaggedObject *obj = theSPs->removeComponent(tag);
    if (obj == 0)
	return 0;
    SP_Constraint *result = (SP_Constraint *)obj;
    result->setDomain(0);
    currentGeoTag++;
    return result;    
}    


void
LoadPattern::applyLoad(double pseudoTime)
{
  // first determine the load factor
  if (theSeries != 0 && isConstant != 0)
    loadFactor = theSeries->getFactor(pseudoTime);

  NodalLoad *nodLoad;
  NodalLoadIter &theNodalIter = this->getNodalLoads();
  while ((nodLoad = theNodalIter()) != 0)
    nodLoad->applyLoad(loadFactor);
    
  ElementalLoad *eleLoad;
  ElementalLoadIter &theElementalIter = this->getElementalLoads();
  while ((eleLoad = theElementalIter()) != 0)
    eleLoad->applyLoad(loadFactor);

  SP_Constraint *sp;
  SP_ConstraintIter &theIter = this->getSPs();
  while ((sp = theIter()) != 0)
    sp->applyConstraint(loadFactor);
}

void
LoadPattern::setLoadConstant(void) 
{
  isConstant = 0;
}


int
LoadPattern::sendSelf(int cTag, Channel &theChannel)
{
  // get my current database tag
  // NOTE - dbTag equals 0 if not sending to a database OR has not yet been sent
  int myDbTag = this->getDbTag();

  // into an ID we place all info needed to determine state of LoadPattern
  int numNodLd, numEleLd, numSPs;
  ID lpData(12);

  numNodLd = theNodalLoads->getNumComponents();
  numEleLd = theElementalLoads->getNumComponents();
  numSPs = theSPs->getNumComponents();

  lpData(11) = this->getTag();
  lpData(0) = currentGeoTag;
  lpData(1) = numNodLd;
  lpData(2) = numEleLd;
  lpData(3) = numSPs;

  if (dbNod == 0) {
    dbNod = theChannel.getDbTag();
    dbEle = theChannel.getDbTag();
    dbSPs = theChannel.getDbTag();
  } 

  lpData(4) = dbNod;
  lpData(5) = dbEle;
  lpData(6) = dbSPs;

  lpData(7) = isConstant;

  if (theSeries != 0) {
    int dbtag = theSeries->getDbTag();
    int classtag = theSeries->getClassTag();
    if (dbtag == 0) {
      dbtag = theChannel.getDbTag();
      theSeries->setDbTag(dbtag);
    }
    lpData(8) = classtag;
    lpData(9) = dbtag;
  } else
    lpData(8) = -1;

  // see if we can save sending the vector containing just the load factor
  // will happen in parallel if sending the loadPattern .. not in database
  if (loadFactor == 0.0)
    lpData(10) = 0;
  else
    lpData(10) = 1;

  if (theChannel.sendID(myDbTag, cTag, lpData) < 0) {
    g3ErrorHandler->warning("LoadPattern::sendSelf - channel failed to send the initial ID");
    return -1;
  }    
  
  if (loadFactor != 0.0) {
    Vector data(1);
    data(0) = loadFactor;
    if (theChannel.sendVector(myDbTag, cTag, data) < 0) {
      g3ErrorHandler->warning("LoadPattern::sendSelf - channel failed to send the Vector");
      return -2;
    }
  }

  if (theSeries != 0)
    if (theSeries->sendSelf(cTag, theChannel) < 0) {
      g3ErrorHandler->warning("LoadPattern::sendSelf - the TimeSeries failed to send");
      return -3;
    }

  // now check if data defining the objects in the LoadPAttern needs to be sent 
  // NOTE THIS APPROACH MAY NEED TO CHANGE FOR VERY LARGE PROBLEMS IF CHANNEL CANNOT
  // HANDLE VERY LARGE ID OBJECTS.
  if (lastGeoSendTag != currentGeoTag) {


    //
    // into an ID we are gonna place the class and db tags for each node so can rebuild
    // this ID we then send to the channel
    //

    // create the ID and get the node iter
    if (numNodLd != 0) {
      ID nodeData(numNodLd*2);
      NodalLoad *theNode;
      NodalLoadIter &theNodes = this->getNodalLoads();
      int loc =0;

      // loop over nodes in domain adding their classTag and dbTag to the ID
      while ((theNode = theNodes()) != 0) {
	nodeData(loc) = theNode->getClassTag();
	int dbTag = theNode->getDbTag();
	
	// if dbTag still 0 get one from Channel; 
	// if this tag != 0 set the dbTag in node
	if (dbTag == 0 && myDbTag != 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theNode->setDbTag(dbTag);
	}
	
	nodeData(loc+1) = dbTag;
	loc+=2;
      }    

      // now send the ID
      if (theChannel.sendID(dbNod, currentGeoTag, nodeData) < 0) {
	g3ErrorHandler->warning("LoadPattern::sendSelf - channel failed to send the NodalLoads ID");
	return -4;
      }
    }

    // we do the same for elemental loads as we did for nodal loads above .. see comments above!

    if (numEleLd != 0) {
      ID elementData(numEleLd*2);
      ElementalLoad *theEle;
      ElementalLoadIter &theElements = this->getElementalLoads();
      int loc = 0;
    
      while ((theEle = theElements()) != 0) {
	elementData(loc) = theEle->getClassTag();
	int dbTag = theEle->getDbTag();

	if (dbTag == 0 && myDbTag != 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theEle->setDbTag(dbTag);
	}
      
	elementData(loc+1) = dbTag;
	loc+=2;
      }

      // now send the ID
      if (theChannel.sendID(dbEle, currentGeoTag, elementData) < 0) {
	g3ErrorHandler->warning("Domain::send - channel failed to send the element ID");
	return -5;
      }
    }

    // we do the same for SP_Constraints as for NodalLoads above .. see comments above!
    
    if (numSPs != 0) {
      ID spData(numSPs*2);
      SP_Constraint *theSP;
      SP_ConstraintIter &theSPs = this->getSPs();
      int loc = 0;
    
      while ((theSP = theSPs()) != 0) {
	spData(loc) = theSP->getClassTag();
	int dbTag = theSP->getDbTag();

	if (dbTag == 0 && myDbTag != 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theSP->setDbTag(dbTag);
	}
	
	spData(loc+1) = dbTag;
	loc+=2;
      }    

      if (theChannel.sendID(dbSPs, currentGeoTag, spData) < 0) {
	g3ErrorHandler->warning("LoadPAttern::sendSelf - channel failed sending SP_Constraint ID");
	return -6;
      }
    }

    // set the lst send db tag so we don't have to do all that again
    lastGeoSendTag = currentGeoTag;
  }
  

  // now we invoke sendSelf on all the NodalLoads, ElementalLoads and SP_Constraints
  // which have been added to the LoadCase
  NodalLoad *theNode;
  NodalLoadIter &theNodes = this->getNodalLoads();
  while ((theNode = theNodes()) != 0) {
    if (theNode->sendSelf(cTag, theChannel) < 0) {
      g3ErrorHandler->warning("LoadPattern::sendSelf - node with tag %d failed in sendSelf",
			      theNode->getTag());
      return -7;
    }
  }

  ElementalLoad *theEle;
  ElementalLoadIter &theElements = this->getElementalLoads();
  while ((theEle = theElements()) != 0) {
    if (theEle->sendSelf(cTag, theChannel) < 0) {
      g3ErrorHandler->warning("LoadPattern::sendSelf - element with tag %d failed in sendSelf",
			      theEle->getTag());
      return -8;
    }
  }

  SP_Constraint *theSP;
  SP_ConstraintIter &theSPs = this->getSPs();
  while ((theSP = theSPs()) != 0) {
    if (theSP->sendSelf(cTag, theChannel) < 0) {
      g3ErrorHandler->warning("LoadPattern::sendSelf - SP_Constraint tagged %d failed sendSelf",
			      theSP->getTag());
      return -9;
    }
  }    

  // if we get here we are successfull
  return 0;
}



int
LoadPattern::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{

  // get my current database tag
  // NOTE - dbTag equals 0 if not sending to a database OR has not yet been sent
  int myDbTag = this->getDbTag();

  // into an ID we place all info needed to determine state of LoadPattern
  int numNod, numEle, numSPs;
  ID lpData(12);

  if (theChannel.recvID(myDbTag, cTag, lpData) < 0) {
    g3ErrorHandler->warning("LoadPattern::recvSelf - channel failed to recv the initial ID");
    return -1;
  }    

  this->setTag(lpData(11));

  isConstant = lpData(7);

  if (lpData(10) == 1) { // we must recv the load factor in a Vector
    Vector data(1);
    data(0) = loadFactor;
    if (theChannel.recvVector(myDbTag, cTag, data) < 0) {
      g3ErrorHandler->warning("LoadPattern::recvSelf - channel failed to recv the Vector");
      return -2;
    }
  }
  
  // read data about the time series
  if (lpData(8) != -1) {
    if (theSeries == 0) {
      theSeries = theBroker.getNewTimeSeries(lpData(8));
    } else if (theSeries->getClassTag() != lpData(8)) {
      delete theSeries;    
      theSeries = theBroker.getNewTimeSeries(lpData(8));
    }
    if (theSeries == 0) {
      g3ErrorHandler->warning("LoadPattern::recvSelf - failed to create TimeSeries");
      return -3;
    }
  
    theSeries->setDbTag(lpData(9));

    if (theSeries->recvSelf(cTag, theChannel, theBroker) < 0) {
      g3ErrorHandler->warning("LoadPattern::recvSelf - the TimeSeries failed to recv");
      return -3;
    }
  }

  if (currentGeoTag != lpData(0)) {

    // clear out the all the components in the current load pattern
    this->clearAll();

    currentGeoTag = lpData(0);

    numNod = lpData(1);
    numEle = lpData(2);
    numSPs = lpData(3);
    dbNod = lpData(4);
    dbEle = lpData(5);
    dbSPs = lpData(6);    

    // 
    // now we rebuild the nodal loads
    //
    
    // first get the information from the domainData about the nodes
    if (numNod != 0) {
      ID nodeData(2*numNod);

      // now receive the ID about the nodes, class tag and dbTags
      if (theChannel.recvID(dbNod, currentGeoTag, nodeData) < 0) {
	g3ErrorHandler->warning("LoadPAttern::recvSelf - channel failed to recv the NodalLoad ID");
	return -2;
      }

      // now for each NodalLoad we 1) get a new node of the correct type from the ObjectBroker
      // 2) ensure the node exists and set it's dbTag, 3) we invoke recvSelf on this new 
      // blank node and 4) add this node to the domain

      int loc = 0;

      for (int i=0; i<numNod; i++) {
	int classTag = nodeData(loc);
	int dbTag = nodeData(loc+1);
	
	NodalLoad *theNode = theBroker.getNewNodalLoad(classTag);

	if (theNode == 0) {
	  g3ErrorHandler->warning("LoadPattern::recv - cannot create NodalLoad with classTag %d ",
				  classTag);
	  return -2;
	}			
	
	theNode->setDbTag(dbTag);
	
	if (theNode->recvSelf(cTag, theChannel, theBroker) < 0) {
	  g3ErrorHandler->warning("LoadPattern::recvSelf - NodalLoad with dbTag %d failed in recvSelf",
				  dbTag);
	  return -2;
	}			

	if (this->addNodalLoad(theNode) == false) {
	  g3ErrorHandler->warning("LoadPattern::recvSelf - failed adding NodalLoad tagged %d into LP!",
				  theNode->getTag());
	  return -3;
	}			
	  
	loc+=2;
      }   
    }

    // 
    // now we rebuild the ElementalLoads .. same as NodalLoads above .. see comments above
    //
    
    if (numEle != 0) {
      ID eleData(2*numEle);
      
      if (theChannel.recvID(dbEle, currentGeoTag, eleData) < 0) {
	g3ErrorHandler->warning("LoadPattern::recvSelf - channel failed to recv the EleLoad ID");
	return -2;
      }

      int loc = 0;
      for (int i=0; i<numEle; i++) {
	int classTag = eleData(loc);
	int dbTag = eleData(loc+1);
      
	ElementalLoad *theEle = theBroker.getNewElementalLoad(classTag);
	if (theEle == 0) {
	  g3ErrorHandler->warning("LoadPattern::recv - cannot create ElementalLoad with classTag %d ",
				  classTag);
	  return -2;
	}			

	theEle->setDbTag(dbTag);
	
	if (theEle->recvSelf(cTag, theChannel, theBroker) < 0) {
	  g3ErrorHandler->warning("LoadPattern::recvSelf - Ele with dbTag %d failed in recvSelf",
				  dbTag);
	  return -2;
	}			
	
	if (this->addElementalLoad(theEle) == false) {
	  g3ErrorHandler->warning("LoadPattern::recvSelf - could not add Ele with tag %d into LP!",
				  theEle->getTag());
	  return -3;
	}			
	
	loc+=2;
      }
    }

    // 
    // now we rebuild the SP_Constraints .. same as nodes above .. see above if can't understand!!
    //
    
    if (numSPs != 0) {
      ID spData(2*numSPs);

      if (theChannel.recvID(dbSPs, currentGeoTag, spData) < 0) {
	g3ErrorHandler->warning("LoadPattern::recvSelf - channel failed to recv the SP_Constraints ID");
	return -2;
      }

      int loc = 0;
      for (int i=0; i<numSPs; i++) {
	int classTag = spData(loc);
	int dbTag = spData(loc+1);
      
	SP_Constraint *theSP = theBroker.getNewSP(classTag);
	if (theSP == 0) {
	  g3ErrorHandler->warning("LoadPattern::recv - cannot create SP_Constraint with classTag %d ",
				  classTag);
	  return -2;
	}			
	theSP->setDbTag(dbTag);
      
	if (theSP->recvSelf(cTag, theChannel, theBroker) < 0) {
	  g3ErrorHandler->warning("LoadPattern::recvSelf - SP_Constraint with dbTag %d failed in recvSelf",
				  dbTag);
	  return -2;
	}			
	
	if (this->addSP_Constraint(theSP) == false) {
	  g3ErrorHandler->warning("LoadPattern::recvSelf - could not add SP_Constraint with tag %d into LP!",
				  theSP->getTag());
	  return -3;
	}			
	
	loc+=2;
      }
    }

    // now set the load pattern db count
    currentGeoTag = lpData(0);
    lastGeoSendTag  = currentGeoTag;

  } else {
    if (theSeries != 0)
      if (theSeries->recvSelf(cTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("LoadPattern::recvSelf - the TimeSeries failed to recv");
	return -3;
      }

    
    NodalLoad *theNode;
    NodalLoadIter &theNodes = this->getNodalLoads();
    while ((theNode = theNodes()) != 0) {
      if (theNode->recvSelf(cTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("LoadPattern::recvSelf - node with tag %d failed in recvSelf",
				theNode->getTag());
	return -7;
      }
    }

    ElementalLoad *theEle;
    ElementalLoadIter &theElements = this->getElementalLoads();
    while ((theEle = theElements()) != 0) {
      if (theEle->recvSelf(cTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("LoadPattern::recvSelf - element with tag %d failed in recvSelf",
				theEle->getTag());
	return -8;
      }
    }

    SP_Constraint *theSP;
    SP_ConstraintIter &theSPs = this->getSPs();
    while ((theSP = theSPs()) != 0) {
      if (theSP->recvSelf(cTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("LoadPattern::recvSelf - SP_Constraint tagged %d failed recvSelf",
				theSP->getTag());
	return -9;
      }
    }    
  }

  // if we get here we are successfull
  return 0;
}

void
LoadPattern::Print(ostream &s, int flag)
{
    s << "Load Pattern: " << this->getTag() << "\n";
    if (theSeries != 0)
      theSeries->Print(s,flag);
    cerr << "  Nodal Loads: \n";
    theNodalLoads->Print(s,flag);
    cerr << "\n  Elemental Loads: \n";
    theElementalLoads->Print(s, flag);
    cerr << "\n  Single Point Constraints: \n";
    theSPs->Print(s, flag);
}


LoadPattern *
LoadPattern::getCopy(void)
{
  LoadPattern *theCopy = new LoadPattern(this->getTag());
  if (theCopy == 0) {
    g3ErrorHandler->fatal("LoadPattern::getCopy() - ran out of memory\n");
    return theCopy; // in case fatal() does not exit
  }
  theCopy->loadFactor = loadFactor;
  theCopy->isConstant = isConstant;
  theCopy->theSeries = theSeries;
  return theCopy;
}
    
    
int
LoadPattern::setParameter(char **argv, int argc, Information &info)
{
	// Nodal load
	if (strcmp(argv[0],"loadAtNode") == 0) {

		int nodeNumber = atoi(argv[1]);
		NodalLoad *thePossibleNodalLoad;
		NodalLoad *theNodalLoad = 0;
		NodalLoadIter &theNodalIter = this->getNodalLoads();

		while ((thePossibleNodalLoad = theNodalIter()) != 0) {
			if ( nodeNumber == thePossibleNodalLoad->getNodeTag() ) {
				theNodalLoad = thePossibleNodalLoad;
			}
		}

		int ok = -1;
		ok = theNodalLoad->setParameter (&argv[2], argc, info);

		if (ok > 0 )
			return ok*1000 + nodeNumber;
		else
			return -1;
	}

	// Unknown parameter
	else
		return -1;
}

int
LoadPattern::updateParameter(int parameterID, Information &info)
{
	NodalLoad *thePossibleNodalLoad = 0;
	NodalLoad *theNodalLoad = 0;
	NodalLoadIter &theNodalIter = this->getNodalLoads();

	switch (parameterID) {
	case 1: case -1:  // Not implemented.
		return -1;
    default:
		if (parameterID > 1000  &&  parameterID < 2000)  {
			int nodeNumber = parameterID-1000;
			while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
				if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
					theNodalLoad = thePossibleNodalLoad;
				}
			}
			return theNodalLoad->updateParameter(1, info);
		}
        else if (parameterID > 2000  &&  parameterID < 3000)  {
            int nodeNumber = parameterID-2000;
            while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
                if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
                    theNodalLoad = thePossibleNodalLoad;
                }
            }
            return theNodalLoad->updateParameter(2, info);
        }
        else if (parameterID > 3000  &&  parameterID < 4000)  {
            int nodeNumber = parameterID-3000;
            while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
                if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
                        theNodalLoad = thePossibleNodalLoad;
                }
			}
			return theNodalLoad->updateParameter(3, info);
        }
        else
			return -1;
    }
}

// AddingSensitivity:BEGIN ////////////////////////////////////
const Vector &
LoadPattern::gradient(bool compute, int identifier)
{
	// If compute==true: return the randomLoads vector, 
	// if not: just set the gradientIdentifier of
	// the nodal loads of this pattern. 

	// For the purpose of setting flags, the identifier
	// will either be 0 (nothing random) or parameterID. 
	// For the purpose of computing gradients the
	// identifier has no use. 


	// Initial declarations
	Vector tempRandomLoads(1);
	int sizeRandomLoads;

	// Start with a fresh return vector
	if (randomLoads == 0) {
		randomLoads = new Vector(1);
	}
	else {
		delete randomLoads;
		randomLoads = new Vector(1);
	}
	
	// Prepare the vector identifying which loads are random. 
	if (compute) { 

		NodalLoad *theNodalLoad = 0;
		NodalLoadIter &theNodalIter = this->getNodalLoads();
		int i;

		// Loop through the nodal loads to pick up possible contributions
		int nodeNumber;
		int dofNumber;
		while ((theNodalLoad = theNodalIter()) != 0)  {
			const Vector &gradientVector = theNodalLoad->gradient(true,0);
			if (gradientVector(0) != 0.0 ) {

				// Found a random load! Get nodeNumber and dofNumber
				nodeNumber = theNodalLoad->getNodeTag();
				dofNumber = (int)gradientVector(0);
				
				// Update the randomLoads vector
				sizeRandomLoads = randomLoads->Size();
				if (sizeRandomLoads == 1) {
					delete randomLoads;
					randomLoads = new Vector(2);
					(*randomLoads)(0) = (double)nodeNumber;
					(*randomLoads)(1) = (double)dofNumber;
				}
				else {
					tempRandomLoads = (*randomLoads);
					delete randomLoads;
					randomLoads = new Vector(sizeRandomLoads+2);
					for (i=0; i<sizeRandomLoads; i++) {
						(*randomLoads)(i) = tempRandomLoads(i);
					}
					(*randomLoads)(sizeRandomLoads) = nodeNumber;
					(*randomLoads)(sizeRandomLoads+1) = dofNumber;
				}
			}	
		}
	}
	else {
		
		// Don't set flag here in the load pattern itself.
		// (Assume there always may be random loads)

		NodalLoad *theNodalLoad = 0;
		NodalLoadIter &theNodalIter = this->getNodalLoads();

		if (identifier == 0) {

			// Go through all nodal loads and zero out gradientIdentifier
			// (Remember: the identifier is only zero if we are in 
			// the process of zeroing out all sensitivity flags).
			while ((theNodalLoad = theNodalIter()) != 0)  {
				theNodalLoad->gradient(false,0);
			}

		}
		else {

			// Find the right nodal load and set the flag
			if (identifier > 1000  &&  identifier < 2000)  {
				int nodeNumber = identifier-1000;
				while ((theNodalLoad = theNodalIter()) != 0)  {
					if ( nodeNumber == theNodalLoad->getNodeTag() )  {
						theNodalLoad->gradient(false,1);
					}
				}
			}
			else if (identifier > 2000  &&  identifier < 3000)  {
				int nodeNumber = identifier-2000;
				while ((theNodalLoad = theNodalIter()) != 0)  {
					if ( nodeNumber == theNodalLoad->getNodeTag() )  {
						theNodalLoad->gradient(false,2);
					}
				}
			}
			else if (identifier > 3000  &&  identifier < 4000)  {
				int nodeNumber = identifier-3000;
				while ((theNodalLoad = theNodalIter()) != 0)  {
					if ( nodeNumber == theNodalLoad->getNodeTag() )  {
						theNodalLoad->gradient(false,3);
					}
				}
			}
			else {
				cerr << "LoadPattern::gradient() -- error in identifier. " << endl;
			}
		}
	}
	return (*randomLoads);
}
// AddingSensitivity:END //////////////////////////////////////
