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
                                                                        
// $Revision: 1.7 $
// $Date: 2002-04-02 18:49:50 $
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
#include <FE_Datastore.h>

ElementRecorder::ElementRecorder(const ID &eleID, Domain &theDom, 
				 char **argv, int argc,
				 bool echoTime, double dT, char *theFileName)
:numEle(eleID.Size()), responseID(eleID.Size()), theDomain(&theDom),
 echoTimeFlag(echoTime), deltaT(dT), nextTimeStampToRecord(0.0), 
 db(0), dbColumns(0), numDbColumns(0), data(0)
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
  if (theFileName != 0) {

    // create char array to store file name
    int fileNameLength = strlen(theFileName) + 1;
    fileName = new char[fileNameLength];
    if (fileName == 0) {
      g3ErrorHandler->fatal("FileNodeDispRecorder::FileNodeDispRecorder - out of memory creating string %d long\n",
			    fileNameLength);
    }

    // copy file name string
    strcpy(fileName, theFileName);    

    // open the file
    theFile.open(fileName, ios::out);
    if (theFile.bad()) {
      cerr << "WARNING - ElementRecorder::ElementRecorder()";
      cerr << " - could not open file " << fileName << endl;
    }    
  } 
}


ElementRecorder::ElementRecorder(const ID &eleID, Domain &theDom, 
				 char **argv, int argc,
				 bool echoTime, 
				 FE_Datastore *database, 
				 char *tableName, 
				 double dT)
:numEle(eleID.Size()), responseID(eleID.Size()), theDomain(&theDom),
 echoTimeFlag(echoTime), deltaT(dT), nextTimeStampToRecord(0.0), 
 db(database), dbColumns(0), numDbColumns(0), data(0)
{

  int numDbColumns = 1;  // 1 for the pseudo-time

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
      if (theResponses[i] != 0) {
      
	// from the response type determine no of cols for each
	Information &eleInfo = theResponses[i]->myInfo;
	if (eleInfo.theType == IntType || eleInfo.theType == DoubleType)
	  numDbColumns += 1;
	else if (eleInfo.theType == VectorType) 
	  numDbColumns += eleInfo.theVector->Size();
	else if (eleInfo.theType == IdType) 
	  numDbColumns += eleInfo.theID->Size();
      }
      theElements[i] = theEle;
    }
  }

  // create char array to store table name
  int fileNameLength = strlen(tableName) + 1;
  fileName = new char[fileNameLength];
  if (fileName == 0) {
    g3ErrorHandler->fatal("FileNodeDispRecorder::FileNodeDispRecorder - out of memory creating string %d long\n",
			  fileNameLength);
  }

  // copy the strings
  strcpy(fileName, tableName);    

  // now create the columns strings for the database
  // for each element do a getResponse() & print the result
  dbColumns = new (char *)[numDbColumns];

  char aColumn[256]; // assumes a column name will not be longer than 256 characters
  
  char *newColumn = new char[10];
  sprintf(newColumn, "%s","time");  
  dbColumns[0] = newColumn;

  char *dataToStore = argv[argc-1];

  int counter = 1;
  for (int i=0; i<eleID.Size(); i++) {
    int eleTag = eleID(i);
    int numVariables = 0;
    if (theResponses[i]!= 0) {
      Information &eleInfo = theResponses[i]->myInfo;
      
      if (eleInfo.theType == IntType || eleInfo.theType == DoubleType) {
	// create column heading for single data item for element
	numVariables = 0;
	sprintf(aColumn, "%s_%d",dataToStore,eleTag);
	int lenColumn = strlen(aColumn+1);
	char *newColumn = new char[lenColumn];
	sprintf(newColumn, "%s",aColumn);
	dbColumns[counter] = newColumn;
	counter++;
      }
      else if (eleInfo.theType == VectorType) 
	numVariables = eleInfo.theVector->Size();
      else if (eleInfo.theType == IdType) 
	numVariables = eleInfo.theVector->Size();

      // create the column headings for multiple data for the element
      for (int j=1; j<=numVariables; j++) {
	sprintf(aColumn, "%s_%d_%d",dataToStore,eleTag,j);
	int lenColumn = strlen(aColumn+1);
	char *newColumn = new char[lenColumn];
	sprintf(newColumn, "%s",aColumn);
	dbColumns[counter] = newColumn;
	counter++;
      }
    }
  }
  
  // create the table in the database
  db->createTable(tableName, numDbColumns, dbColumns);
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

    if (fileName != 0)
      delete [] fileName;
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
  theFile.open(fileName, ios::out);
  if (!theFile) {
    cerr << "WARNING - FileNodeDispRecorder::restart() - could not open file ";
    cerr << fileName << endl;
  }    
}
