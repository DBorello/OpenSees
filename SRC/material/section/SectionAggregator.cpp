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
// $Date: 2002-05-17 01:24:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/SectionAggregator.cpp,v $
                                                                        
                                                                        
// File: ~/section/SectionAggregator.C
//
// Written: MHS
// Created: June 2000
// Revision: A
//
// Purpose: This file contains the implementation for the SectionAggregator class. 

#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <classTags.h>
#include <SectionAggregator.h>
#include <MaterialResponse.h>
#include <ID.h>

#include <classTags.h>

#define maxOrder 10

// Assumes section order is less than or equal to maxOrder.
// Can increase if needed!!!
double SectionAggregator::workArea[(maxOrder+1)*maxOrder];
int    SectionAggregator::codeArea[maxOrder];

// constructors:
SectionAggregator::SectionAggregator (int tag, SectionForceDeformation &theSec,
				      int numAdds, UniaxialMaterial **theAdds,
				      const ID &addCodes): 
  SectionForceDeformation(tag, SEC_TAG_Aggregator), 
  theSection(0), theAdditions(0), matCodes(0), numMats(numAdds),
  theVector(0),
  theMatrix(0),
  theCode(0),
  otherDbTag(0)
{
    theSection = theSec.getCopy();
    
    if (!theSection)
	g3ErrorHandler->fatal("%s -- failed to get copy of section",
			      "SectionAggregator::SectionAggregator");

    if (!theAdds)
	g3ErrorHandler->fatal("%s -- null uniaxial material array passed",
			      "SectionAggregator::SectionAggregator");

    theAdditions = new UniaxialMaterial *[numMats];

    if (!theAdditions)
	g3ErrorHandler->fatal("%s -- failed to allocate pointers",
			      "SectionAggregator::SectionAggregator");
    
    int i;
    
    for (i = 0; i < numMats; i++) {
	if (!theAdds[i])
	    g3ErrorHandler->fatal("%s -- null uniaxial material pointer passed",
				  "SectionAggregator::SectionAggregator");
	
	theAdditions[i] = theAdds[i]->getCopy(this);
	
	if (!theAdditions[i])
	    g3ErrorHandler->fatal("%s -- failed to copy uniaxial material",
				  "SectionAggregator::SectionAggregator");
    }

    if (theSec.getOrder()+numAdds > maxOrder) {
      g3ErrorHandler->fatal("%s -- order too big, need to modify the #define in SectionAggregator.cpp to %d",
			    "SectionAggregator::SectionAggregator", theSec.getOrder()+numAdds);
    }

    theCode = new ID(codeArea, theSec.getOrder()+numAdds);
    theMatrix = new Matrix(&workArea[maxOrder], theSec.getOrder()+numAdds, theSec.getOrder()+numAdds);
    theVector = new Vector(workArea, theSec.getOrder()+numAdds);
    matCodes = new ID(addCodes);

    if (theCode == 0 || theMatrix == 0 || theVector == 0 || matCodes == 0)
	g3ErrorHandler->fatal("%s -- out of memory",
			      "SectionAggregator::SectionAggregator");
}

SectionAggregator::SectionAggregator (int tag, int numAdds,
				      UniaxialMaterial **theAdds,
				      const ID &addCodes): 
  SectionForceDeformation(tag, SEC_TAG_Aggregator), 
  theSection(0), theAdditions(0), matCodes(0), numMats(numAdds),
  theVector(0), theMatrix(0), theCode(0),
  otherDbTag(0)
{
    if (!theAdds)
	g3ErrorHandler->fatal("%s -- null uniaxial material array passed",
			      "SectionAggregator::SectionAggregator");

    theAdditions = new UniaxialMaterial *[numMats];

    if (!theAdditions)
	g3ErrorHandler->fatal("%s -- failed to allocate pointers",
			      "SectionAggregator::SectionAggregator");
    
    int i;
    
    for (i = 0; i < numMats; i++) {
	if (!theAdds[i])
	    g3ErrorHandler->fatal("%s -- null uniaxial material pointer passed",
				  "SectionAggregator::SectionAggregator");
	
	theAdditions[i] = theAdds[i]->getCopy(this);
	
	if (!theAdditions[i])
	    g3ErrorHandler->fatal("%s -- failed to copy uniaxial material",
				  "SectionAggregator::SectionAggregator");
    }

    if (numAdds > maxOrder) {
      	    g3ErrorHandler->fatal("%s -- order too big, need to modify the #define in SectionAggregator.cpp to %d",
				  "SectionAggregator::SectionAggregator", numAdds);
    }

    theCode = new ID(codeArea, numAdds);
    theMatrix = new Matrix(&workArea[maxOrder], numAdds, numAdds);
    theVector = new Vector(workArea, numAdds);
    matCodes = new ID(addCodes);

    if (theCode == 0 || theMatrix == 0 || theVector == 0 || matCodes == 0)
      g3ErrorHandler->fatal("%s -- out of memory",
			    "SectionAggregator::SectionAggregator");
}

SectionAggregator::SectionAggregator (int tag, SectionForceDeformation &theSec,
				      UniaxialMaterial &theAddition, int c) :
  SectionForceDeformation(tag, SEC_TAG_Aggregator),
  theSection(0), theAdditions(0), matCodes(0), numMats(1),
  theVector(0), theMatrix(0), theCode(0),
  otherDbTag(0)
{
    theSection = theSec.getCopy();
    
    if (!theSection)
	g3ErrorHandler->fatal("%s -- failed to get copy of section",
			      "SectionAggregator::SectionAggregator");

    theAdditions = new UniaxialMaterial *[1];

    theAdditions[0] = theAddition.getCopy(this);
    
    if (!theAdditions[0])
	g3ErrorHandler->fatal("%s -- failed to copy uniaxial material",
			      "SectionAggregator::SectionAggregator");
    
    matCodes = new ID(1);
    (*matCodes)(0) = c;

    if (theSec.getOrder()+1 > maxOrder) {
      g3ErrorHandler->fatal("%s -- order too big, need to modify the #define in SectionAggregator.cpp to %d",
			    "SectionAggregator::SectionAggregator", theSec.getOrder()+1);
    }

    theVector = new Vector(workArea, theSec.getOrder()+1);
    theMatrix = new Matrix(&workArea[maxOrder], theSec.getOrder()+1, theSec.getOrder()+1);
    theCode   = new ID(codeArea, theSec.getOrder()+1);
    
    if (theCode == 0 || theMatrix == 0 || theVector == 0 || matCodes == 0)
	g3ErrorHandler->fatal("%s -- out of memory",
			      "SectionAggregator::SectionAggregator");
}

// constructor for blank object that recvSelf needs to be invoked upon
SectionAggregator::SectionAggregator():
  SectionForceDeformation(0, SEC_TAG_Aggregator),
  theSection(0), theAdditions(0), matCodes(0), numMats(0), 
  theVector(0), theMatrix(0), theCode(0),
  otherDbTag(0)
{

}

// destructor:
SectionAggregator::~SectionAggregator()
{
   int i;

   if (theSection)
       delete theSection;

   for (i = 0; i < numMats; i++)
       if (theAdditions[i])
	   delete theAdditions[i];

   if (theAdditions)
       delete [] theAdditions;
   
   if (theVector != 0)
     delete theVector;

   if (theMatrix != 0)
     delete theMatrix;

   if (theCode != 0)
     delete theCode;

   if (matCodes != 0)
     delete matCodes;
}

int SectionAggregator::setTrialSectionDeformation (const Vector &def)
{
  int ret = 0;
  int i = 0;

  int theSectionOrder = 0;

  if (theSection) {
    theSectionOrder = theSection->getOrder();
    Vector v(workArea, theSectionOrder);
    
    for (i = 0; i < theSectionOrder; i++)
      v(i) = def(i);
    
    ret = theSection->setTrialSectionDeformation(v);
  }

  int order = theSectionOrder + numMats;
  
  for ( ; i < order; i++)
    ret += theAdditions[i-theSectionOrder]->setTrialStrain(def(i));
  
  return ret;
}

const Vector &
SectionAggregator::getSectionDeformation(void)
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const Vector &v = theSection->getSectionDeformation();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*theVector)(i) = v(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*theVector)(i) = theAdditions[i-theSectionOrder]->getStrain();

  return *theVector;
}

const Matrix &
SectionAggregator::getSectionTangent(void)
{
  int i = 0;

  int theSectionOrder = 0;

  // Zero before assembly
  theMatrix->Zero();

  if (theSection) {
    const Matrix &ks = theSection->getSectionTangent();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*theMatrix)(i,j) = ks(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*theMatrix)(i,i) = theAdditions[i-theSectionOrder]->getTangent();
  
  return *theMatrix;
}

const Matrix &
SectionAggregator::getSectionFlexibility(void)
{
  int i = 0;
    
  int theSectionOrder = 0;

  // Zero before assembly
  theMatrix->Zero();

  if (theSection) {
    const Matrix &fs = theSection->getSectionFlexibility();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*theMatrix)(i,j) = fs(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++) {
    double k = theAdditions[i-theSectionOrder]->getTangent();
    if (k == 0.0) {
      g3ErrorHandler->warning("%s -- singular section stiffness",
			       "SectionAggregator::getSectionFlexibility");
      (*theMatrix)(i,i) = 1.e14;
    }
    else
      (*theMatrix)(i,i) = 1/k;
  }	
  
  return *theMatrix;
}

const Vector &
SectionAggregator::getStressResultant(void)
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const Vector &s = theSection->getStressResultant();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*theVector)(i) = s(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*theVector)(i) = theAdditions[i-theSectionOrder]->getStress();
  
  return *theVector;
}

SectionForceDeformation *
SectionAggregator::getCopy(void)
{
  SectionAggregator *theCopy = 0;
    
  if (theSection)
    theCopy = new SectionAggregator(this->getTag(), *theSection,
				    numMats, theAdditions, *matCodes);
  else
    theCopy = new SectionAggregator(this->getTag(), numMats,
				    theAdditions, *matCodes);
  
  if (theCopy == 0)
    g3ErrorHandler->fatal("%s -- failed to allocate copy",
			  "SectionAggregator::getCopy");
  
  return theCopy;
}

const ID&
SectionAggregator::getType ()
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const ID &secType = theSection->getType();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*theCode)(i) = secType(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*theCode)(i) = (*matCodes)(i-theSectionOrder);

  return *theCode;
}

int
SectionAggregator::getOrder () const
{
  int order = numMats;

  if (theSection != 0)
    order += theSection->getOrder();

  return order;
}

int
SectionAggregator::commitState(void)
{
  int err = 0;
    
  if (theSection)
    err += theSection->commitState();
  
  for (int i = 0; i < numMats; i++)
    err += theAdditions[i]->commitState();
  
  return err;
}

int
SectionAggregator::revertToLastCommit(void)
{
  int err = 0;
  
  int i = 0;
  
  // Revert the section
  if (theSection)
    err += theSection->revertToLastCommit();
  
  // Do the same for the uniaxial materials
  for (i = 0; i < numMats; i++)
    err += theAdditions[i]->revertToLastCommit();
  
  return err;
}	

int
SectionAggregator::revertToStart(void)
{
  int err = 0;
  
  // Revert the section
  if (theSection)
    err += theSection->revertToStart();
  
  // Do the same for the uniaxial materials
  for (int i = 0; i < numMats; i++)
    err += theAdditions[i]->revertToStart();
  
  return err;
}

int
SectionAggregator::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  // Need otherDbTag since classTags ID and data ID may be the same size
  if (otherDbTag == 0) 
    otherDbTag = theChannel.getDbTag();
  
  // Create ID for tag and section order data
  static ID data(5);
  
  int order = this->getOrder();
  
  data(0) = this->getTag();
  data(1) = otherDbTag;
  data(2) = order;
  data(3) = (theSection != 0) ? theSection->getOrder() : 0;
  data(4) = numMats;

  // Send the tag and section order data
  res += theChannel.sendID(this->getDbTag(), cTag, data);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not send data ID",
			    "SectionAggregator::sendSelf");
    return res;
  }
  
  // Determine how many classTags there are and allocate ID vector
  // for the tags and section code
  int numTags = (theSection == 0) ? numMats : numMats + 1;
  ID classTags(2*numTags + numMats);
  
  // Loop over the UniaxialMaterials filling in class and db tags
  int i, dbTag;
  for (i = 0; i < numMats; i++) {
    classTags(i) = theAdditions[i]->getClassTag();
    
    dbTag = theAdditions[i]->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theAdditions[i]->setDbTag(dbTag);
    }
    
    classTags(i+numTags) = dbTag;
  }
  
  // Put the Section class and db tags into the ID vector
  if (theSection != 0) {
    classTags(numTags-1) = theSection->getClassTag();
    
    dbTag = theSection->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theSection->setDbTag(dbTag);
    }
    
    classTags(2*numTags-1) = dbTag;
  }
  
  // Put the UniaxialMaterial codes into the ID vector
  int j = 2*numTags;
  for (i = 0; i < numMats; i++, j++)
    classTags(j) = (*matCodes)(i);
  
  // Send the material class and db tags and section code
  res += theChannel.sendID(otherDbTag, cTag, classTags);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not send classTags ID",
			    "SectionAggregator::sendSelf");
    return res;
  }

  // Ask the UniaxialMaterials to send themselves
  for (i = 0; i < numMats; i++) {
    res += theAdditions[i]->sendSelf(cTag, theChannel);
    if (res < 0) {
      g3ErrorHandler->warning("%s -- could not send UniaxialMaterial, i = %d",
			      "SectionAggregator::sendSelf", i);
      return res;
    }
  }
  
  // Ask the Section to send itself
  if (theSection != 0) {
    res += theSection->sendSelf(cTag, theChannel);
    if (res < 0) {
      g3ErrorHandler->warning("%s -- could not send SectionForceDeformation",
			      "SectionAggregator::sendSelf");
      return res;
    }
  }
  
  return res;
}


int
SectionAggregator::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // Create an ID and receive tag and section order
  static ID data(5);
  res += theChannel.recvID(this->getDbTag(), cTag, data);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not receive data ID",
			    "SectionAggregator::recvSelf");
    return res;
  }
  
  this->setTag(data(0));
  otherDbTag = data(1);
  int order = data(2);
  int theSectionOrder = data(3);
  numMats = data(4);

  if (order > 0) {
    if (theVector == 0 || theVector->Size() != order) {
      if (theVector != 0) {
	delete theVector;
	delete theMatrix;
	delete theCode;
      }
      theMatrix = new Matrix(&workArea[maxOrder], order,order);
      theVector = new Vector(workArea, order);
      theCode = new ID(codeArea, order);
    }
  }

  if (numMats > 0) {
    if (matCodes == 0 || matCodes->Size() != numMats) {
      if (matCodes != 0)
	delete matCodes;

      matCodes = new ID(numMats);
    }
  }

  // Determine how many classTags there are and allocate ID vector
  int numTags = (theSectionOrder == 0) ? numMats : numMats + 1;
  ID classTags(numTags*2 + numMats);
  
  // Receive the material class and db tags
  res += theChannel.recvID(otherDbTag, cTag, classTags);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not receive classTags ID",
			    "SectionAggregator::recvSelf");
    return res;
  }

  // Check if null pointer, allocate if so
  if (theAdditions == 0) {
    theAdditions = new UniaxialMaterial *[numMats];
    if (theAdditions == 0) {
      g3ErrorHandler->warning("%s -- could not allocate UniaxialMaterial array",
			      "SectionAggregator::recvSelf");
      return -1;
    }
    // Set pointers to null ... will get allocated by theBroker
    for (int j = 0; j < numMats; j++)
      theAdditions[j] = 0;
  }
  
  // Loop over the UniaxialMaterials
  int i, classTag;
  for (i = 0; i < numMats; i++) {
    classTag = classTags(i);
    
    // Check if the UniaxialMaterial is null; if so, get a new one
    if (theAdditions[i] == 0)
      theAdditions[i] = theBroker.getNewUniaxialMaterial(classTag);
    
    // Check that the UniaxialMaterial is of the right type; if not, delete
    // the current one and get a new one of the right type
    else if (theAdditions[i]->getClassTag() != classTag) {
      delete theAdditions[i];
      theAdditions[i] = theBroker.getNewUniaxialMaterial(classTag);
    }
    
    // Check if either allocation failed
    if (theAdditions[i] == 0) {
      g3ErrorHandler->warning("%s -- could not get UniaxialMaterial, i = %d",
			      "SectionAggregator::recvSelf", i);
      return -1;
    }
    
    // Now, receive the UniaxialMaterial
    theAdditions[i]->setDbTag(classTags(i+numTags));
    res += theAdditions[i]->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      g3ErrorHandler->warning("%s -- could not receive UniaxialMaterial, i = %d",
			      "SectionAggregator::recvSelf", i);
      return res;
    }
  }

  // If there is no Section to receive, return
  if (theSectionOrder == 0)
    return res;
  
  classTag = classTags(numTags-1);
  
  // Check if the Section is null; if so, get a new one
  if (theSection == 0)
    theSection = theBroker.getNewSection(classTag);
  
  // Check that the Section is of the right type; if not, delete
  // the current one and get a new one of the right type
  else if (theSection->getClassTag() != classTag) {
    delete theSection;
    theSection = theBroker.getNewSection(classTag);
  }
  
  // Check if either allocation failed
  if (theSection == 0) {
    g3ErrorHandler->warning("%s -- could not get a SectionForceDeformation",
			    "SectionAggregator::recvSelf");
    return -1;
  }

  // Now, receive the Section
  theSection->setDbTag(classTags(2*numTags-1));
  res += theSection->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not receive SectionForceDeformation",
			    "SectionAggregator::recvSelf");
    return res;
  }
  
  // Fill in the section code
  int j = 2*numTags;
  for (i = 0; i < numMats; i++, j++)
    (*matCodes)(i) = classTags(j);

  return res;
}

void
SectionAggregator::Print(ostream &s, int flag)
{
  s << "\nSection Aggregator, tag: " << this->getTag() << endl;
  if (theSection) {
    s << "\tSection, tag: " << theSection->getTag() << endl;
    theSection->Print(s, flag);
  }
  s << "\tUniaxial Additions" << endl;
  for (int i = 0; i < numMats; i++)
    s << "\t\tUniaxial Material, tag: " << theAdditions[i]->getTag() << endl;
  s << "\tUniaxial codes " << *matCodes << endl;
}

Response*
SectionAggregator::setResponse(char **argv, int argc, Information &info)
{
	// See if the response is one of the defaults
	Response *res = SectionForceDeformation::setResponse(argv, argc, info);
	if (res != 0)
		return res;

	// If not, forward the request to the section (need to do this to get fiber response)
	// CURRENTLY NOT SENDING ANYTHING OFF TO THE UniaxialMaterials ... Probably
	// don't need anything more from them than stress, strain, and stiffness, 
	// which are covered in base class method ... can change if need arises
	else if (theSection != 0)
		return theSection->setResponse(argv, argc, info);

	else
		return 0;
}

int
SectionAggregator::getResponse(int responseID, Information &info)
{
	// Just call the base class method ... don't need to define
	// this function, but keeping it here just for clarity
	return SectionForceDeformation::getResponse(responseID, info);
}

int
SectionAggregator::setVariable(const char *argv)
{
	// Axial strain
	if (strcmp(argv,"axialStrain") == 0)
		return 1;
	// Curvature about the section z-axis
	else if (strcmp(argv,"curvatureZ") == 0)
		return 2;
	// Curvature about the section y-axis
	else if (strcmp(argv,"curvatureY") == 0)
		return 3;
	else
		return -1;
}

int
SectionAggregator::getVariable(int variableID, double &info)
{
  int i;

  info = 0.0;

  int order = numMats;
  if (theSection != 0)
    order += theSection->getOrder();

  const Vector &e = this->getSectionDeformation();
  const ID &code  = this->getType();

  switch (variableID) {
  case 1:	// Axial strain
    // Series model ... add all sources of deformation
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_P)
	info += e(i);
    return 0;
  case 2:	// Curvature about the section z-axis
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_MZ)
	info += e(i);
    return 0;
  case 3:	// Curvature about the section y-axis
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_MY)
	info += e(i);
    return 0;
  default:
    return -1;
  }
}
