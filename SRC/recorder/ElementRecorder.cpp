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
                                                                        
// $Revision: 1.17 $
// $Date: 2004-11-13 00:57:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ElementRecorder.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 09/99
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
#include <DataOutputHandler.h>
#include <OPS_Globals.h>


ElementRecorder::ElementRecorder(const ID &eleID, 
				 const char **argv, 
				 int argc,
				 bool echoTime, 
				 Domain &theDom, 
				 DataOutputHandler &theOutputHandler,
				 double dT)
:numEle(eleID.Size()), responseID(eleID.Size()), theResponses(0), 
 theDomain(&theDom), theHandler(&theOutputHandler),
 echoTimeFlag(echoTime), deltaT(dT), nextTimeStampToRecord(0.0), data(0)
{
  // 
  // Set the response objects:
  //   1. create an array of pointers for them
  //   2. iterate over the elements invoking setResponse() to get the new objects & determine size of data
  //

  int numDbColumns = 0;
  if (echoTimeFlag == true) 
    numDbColumns = 1;  // 1 for the pseudo-time

  theResponses = new Response *[numEle];
  for (int j=0; j<numEle; j++)
    theResponses[j] = 0;

  Information eleInfo(1.0);
  int i;
  for (i=0; i<numEle; i++) {
    Element *theEle = theDom.getElement(eleID(i));
    if (theEle == 0) {
      opserr << "WARNING ElementRecorder::ElementRecorder() -";
      opserr << " no element with tag: " << eleID(i) << " exists in Domain\n";
      theResponses[i] = 0;
    } else {
      theResponses[i] = theEle->setResponse(argv, argc, eleInfo);
      if (theResponses[i] != 0) {
	// from the response type determine no of cols for each
	Information &eleInfo = theResponses[i]->getInformation();
	const Vector &eleData = eleInfo.getData();
	numDbColumns += eleData.Size();
      }
    }
  }

  //
  // now create the columns strings for the data description
  // for each element do a getResponse() 
  //

  char **dbColumns = new char *[numDbColumns];
  static char aColumn[1012]; // assumes a column name will not be longer than 256 characters
  char *newColumn = new char[5];

  sprintf(newColumn, "%s","time");  
  dbColumns[0] = newColumn;
  
  int lengthString = 0;
  for (i=0; i<argc; i++)
    lengthString += strlen(argv[i])+1;
  char *dataToStore = new char[lengthString];
  lengthString = 0;
  for (int j=0; j<argc; j++) {
    int argLength = strlen(argv[j]);
    strcpy(&dataToStore[lengthString], argv[j]);
    if (j<(argc-1)) {
      lengthString += argLength;
      dataToStore[lengthString] = ' ';
      lengthString ++;
    } else
      lengthString += argLength+1;
  }
  
  int counter = 1;
  for (i=0; i<eleID.Size(); i++) {
    int eleTag = eleID(i);
    int numVariables = 0;
    if (theResponses[i]!= 0) {
      const Information &eleInfo = theResponses[i]->getInformation();
      
      if (eleInfo.theType == IntType || eleInfo.theType == DoubleType) {
	// create column heading for single data item for element
	numVariables = 0;
	sprintf(aColumn, "Element%d_%s", eleTag, dataToStore);
	int lenColumn = strlen(aColumn);
	char *newColumn = new char[lenColumn+1];
	strcpy(newColumn, aColumn);
	dbColumns[counter] = newColumn;
	counter++;
      }
      
      else if (eleInfo.theType == VectorType) 
	numVariables = eleInfo.theVector->Size();
      else if (eleInfo.theType == IdType) 
	numVariables = eleInfo.theID->Size();
      
      // create the column headings for multiple data for the element
      for (int j=1; j<=numVariables; j++) {
	sprintf(aColumn, "Element%d_%s_%d",eleTag, dataToStore, j);
	int lenColumn = strlen(aColumn);
	char *newColumn = new char[lenColumn+1];
	strcpy(newColumn, aColumn);
	dbColumns[counter] = newColumn;
	counter++;
      }
    }
  }

  // replace spaces with undescore for tables
  for (i=0; i<numDbColumns; i++) {
    char *data = dbColumns[i];
    int length = strlen(data);
    for (int j=0; j<length; j++)
      if (data[j] == ' ') data[j]='_';
  }

  //
  // call open in the handler with the data description
  //

  theHandler->open(dbColumns, numDbColumns);

  //
  // clean up the data description
  //

  if (dbColumns != 0) {

    for (int i=0; i<numDbColumns; i++) 
      delete [] dbColumns[i];

      delete [] dbColumns;
  }
  
  // create the vector to hold the data
  data = new Vector(numDbColumns);
}


ElementRecorder::~ElementRecorder()
{
  //
  // invoke the destructor on the respons eobjects
  //
  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }
}


int 
ElementRecorder::record(int commitTag, double timeStamp)
{
  int result = 0;

  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {

    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;
    

    int loc = 0;
    if (echoTimeFlag == true) 
      (*data)(loc++) = timeStamp;
    
    // for each element do a getResponse() & print the result
    for (int i=0; i< numEle; i++) {
      if (theResponses[i] != 0) {
	
	// ask the element for the reponse
	int res;
	if (( res = theResponses[i]->getResponse()) < 0)
	  result = res;
	else {
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  for (int j=0; j<eleData.Size(); j++)
	    (*data)(loc++) = eleData(j);
	}
      } 
    }
    
    theHandler->write(*data);
  }
  
  // succesfull completion - return 0
  return result;
}

void 
ElementRecorder::restart(void)
{
  if (data != 0)
    data->Zero();
}
