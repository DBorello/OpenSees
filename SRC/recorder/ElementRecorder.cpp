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
                                                                        
// $Revision: 1.6 $
// $Date: 2001-10-25 21:31:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ElementRecorder.cpp,v $
                                                                        
                                                                        
// File: ~/recorder/ElementRecorder.C
//
// Written: fmk 
// Created: 09/99
// Revision: A
//
// Description: This file contains the class implementatation of ElementRecorder.
//
// What: "@(#) ElementRecorder.C, revA"

#include <ElementRecorder.h>
#include <Domain.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>

ElementRecorder::ElementRecorder(const ID &eleID, Domain &theDom, 
				 char **argv, int argc,
				 bool echoTime, double dT, char *fileName)
:numEle(eleID.Size()), responseID(eleID.Size()), theDomain(&theDom),
 echoTimeFlag(echoTime), deltaT(dT), nextTimeStampToRecord(0.0)
{
  theElements = new Element *[numEle];
  for (int ii=0; ii<numEle; ii++)
    theElements[ii] = 0;

  theResponses = new Response *[numEle];
  for (int j=0; j<numEle; j++)
    theResponses[j] = 0;

  eleInfoObjects = new Information[numEle];
  for (int i=0; i<numEle; i++) {
    Element *theEle = theDom.getElement(eleID(i));
    if (theEle == 0) {
      cerr << "WARNING ElementRecorder::ElementRecorder() -";
      cerr << " no element with tag: " << eleID(i) << " exists in Domain\n";
      numEle = 0;
      return;
    } else {
      theResponses[i] = theEle->setResponse(argv, argc, eleInfoObjects[i]);
      theElements[i] = theEle;
    }
  }

  // if file is specified, copy name and open the file
  if (fileName != 0) {
    if (strlen(fileName) > MAX_FILENAMELENGTH) 
      g3ErrorHandler->warning("WARNING - ElementRecorder::ElementRecorder() - fileName %s too long, max %d\n",
			      fileName, MAX_FILENAMELENGTH);
    else {
      strcpy(theFileName, fileName);    
      theFile.open(fileName, ios::out);
      if (theFile.bad()) {
	cerr << "WARNING - ElementRecorder::ElementRecorder()";
	cerr << " - could not open file " << fileName << endl;
      }    
    } 
  }
}

ElementRecorder::~ElementRecorder()
{
    // close the file
    if (!theFile.bad())
	theFile.close();    

    if (theElements != 0)
      delete [] theElements;

    if (theResponses != 0) {
      for (int i = 0; i < numEle; i++)
	delete theResponses[i];
      delete [] theResponses;
    }

    if (eleInfoObjects != 0)
      delete [] eleInfoObjects;
}


int 
ElementRecorder::record(int commitTag, double timeStamp)
{
  int result = 0;
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
      
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    // print out the pseudo time if requested
    if (echoTimeFlag == true) {
      if (!theFile) 
	cerr << timeStamp << " ";			
      else 
	theFile << timeStamp << " ";	
    }

    // for each element do a getResponse() & print the result
    for (int i=0; i< numEle; i++) {
      int theID = responseID(i);
      if (theResponses[i] != 0) {
      
	// ask the element for the reponse
	int res;
	Information &eleInfo = eleInfoObjects[i];
	if (( res = theResponses[i]->getResponse()) < 0)
	  result = res;
	else {
	  // print results to file or stderr depending on whether
	  // a file was opened
	
	  if (theFile.bad())
	    theResponses[i]->Print(cerr);	    
	  else {
	    theResponses[i]->Print(theFile);
	    theFile << "  ";  // added for OSP
	  }
	}
      } 
    }

    if (theFile.bad()) 
      cerr << endl;
    else {
      theFile << " \n";
      theFile.flush();
    }
  }

  // succesfull completion - return 0
  return result;
}


int 
ElementRecorder::playback(int commitTag)
{
    return 0;
}

void 
ElementRecorder::restart(void)
{
  theFile.close();
  theFile.open(theFileName, ios::out);
  if (!theFile) {
    cerr << "WARNING - FileNodeDispRecorder::restart() - could not open file ";
    cerr << theFileName << endl;
  }    
}
