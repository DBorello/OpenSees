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
// $Date: 2004-11-13 00:57:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeElementRecorder.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class implementatation of 
// EnvelopeElementRecorder.
//
// What: "@(#) EnvelopeElementRecorder.C, revA"

#include <EnvelopeElementRecorder.h>
#include <Domain.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <FE_Datastore.h>
#include <Information.h>
#include <DataOutputHandler.h>

EnvelopeElementRecorder::EnvelopeElementRecorder(const ID &eleID, 
						 const char **argv, 
						 int argc,
						 Domain &theDom, 
						 DataOutputHandler &theOutputHandler,
						 double dT)
:numEle(eleID.Size()), responseID(eleID.Size()), theDomain(&theDom),
 theHandler(&theOutputHandler), deltaT(dT), nextTimeStampToRecord(0.0), 
 data(0), currentData(0), first(true)
{
  theResponses = new Response *[numEle];
  for (int j=0; j<numEle; j++)
    theResponses[j] = 0;

  Information eleInfo(1.0);
  int numDbColumns = 0;

  for (int ii=0; ii<numEle; ii++) {
    Element *theEle = theDom.getElement(eleID(ii));
    if (theEle == 0) {
      opserr << "WARNING EnvelopeElementRecorder::EnvelopeElementRecorder() -";
      opserr << " no element with tag: " << eleID(ii) << " exists in Domain\n";
      theResponses[ii] = 0;
    } else {
      theResponses[ii] = theEle->setResponse(argv, argc, eleInfo);
      Information &eleInfo = theResponses[ii]->getInformation();
      const Vector &eleData = eleInfo.getData();
      numDbColumns += eleData.Size();
    }
  }
  
  // create the matrix & vector that holds the data
  data = new Matrix(3, numDbColumns);
  currentData = new Vector(numDbColumns);
  if (data == 0 || currentData == 0) {
    opserr << "EnvelopeElementRecorder::EnvelopeElementRecorder() - out of memory\n";
    exit(-1);
  }
  
  // now create the columns strings for the database
  // for each element do a getResponse() & print the result
  char **dbColumns = new char *[numDbColumns];
  static char aColumn[1012]; // assumes a column name will not be longer than 256 characters
  

  int lengthString = 0;
  for (int l=0; l<argc; l++)
    lengthString += strlen(argv[l])+1;
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

  int counter = 0;
  for (int i=0; i<eleID.Size(); i++) {
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
  for (int k=0; k<numDbColumns; k++) {
    char *data = dbColumns[k];
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

}
  
EnvelopeElementRecorder::~EnvelopeElementRecorder()
{
  //
  // write the data
  //

  for (int i=0; i<3; i++) {
    int size = currentData->Size();
    for (int j=0; j<size; j++)
      (*currentData)(j) = (*data)(i,j);
    theHandler->write(*currentData);
  }

  //
  // clean up the memory
  //

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }
  
  if (data != 0)
    delete data;
  
  if (currentData != 0)
    delete currentData;
}


int 
EnvelopeElementRecorder::record(int commitTag, double timeStamp)
{
  int result = 0;
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
      
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;
    
    int loc = 0;
    // for each element do a getResponse() & put the result in current data
    for (int i=0; i< numEle; i++) {
      if (theResponses[i] != 0) {
	// ask the element for the reponse
	int res;
	if (( res = theResponses[i]->getResponse()) < 0)
	  result += res;
	else {
	  // print results to file or stderr depending on whether
	  // a file was opened
	  
	  // from the response determine no of cols for each
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  for (int j=0; j<eleData.Size(); j++) 
	    (*currentData)(loc++) = eleData(j);
	}
      }
    }

    // check if max or min
    // check if currentData modifies the saved data
    bool writeIt = false;
    int size = currentData->Size();
    if (first == true) {
      for (int i=0; i<size; i++) {
	(*data)(0,i) = (*currentData)(i);
	(*data)(1,i) = (*currentData)(i);
	(*data)(2,i) = fabs((*currentData)(i));
	first = false;
	writeIt = true;
      } 
    } else {
      for (int i=0; i<size; i++) {
	double value = (*currentData)(i);
	if ((*data)(0,i) > value) {
	  (*data)(0,i) = value;
	  double absValue = fabs(value);
	  if ((*data)(2,i) < absValue) 
	    (*data)(2,i) = absValue;
	  writeIt = true;
	} else if ((*data)(1,i) < value) {
	  (*data)(1,i) = value;
	  double absValue = fabs(value);
	  if ((*data)(2,i) < absValue) 
	    (*data)(2,i) = absValue;
	  writeIt = true;
	}
      }
    }
  }    

  // succesfull completion - return 0
  return result;
}

void 
EnvelopeElementRecorder::restart(void)
{
  data->Zero();
  first = true;
}
