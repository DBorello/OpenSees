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
// $Date: 2004-11-13 00:57:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeNodeRecorder.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class definition for EnvelopeNodeRecorder.
// A EnvelopeNodeRecorder is used to record the envelop of specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) EnvelopeNodeRecorder.C, revA"

#include <math.h>

#include <EnvelopeNodeRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <FE_Datastore.h>

#include <DataOutputHandler.h>

#include <string.h>
#include <iomanip>
using std::ios;

EnvelopeNodeRecorder::EnvelopeNodeRecorder(const ID &dofs, 
					   const ID &nodes, 
					   const char *dataToStore,
					   Domain &theDom,
					   DataOutputHandler &theOutputHandler,
					   double dT)
:theDofs(0), theNodes(0), 
 currentData(0), data(0), 
 theDomain(&theDom), theHandler(&theOutputHandler),
 deltaT(dT), nextTimeStampToRecord(0.0), 
 first(true)
{
  // verify dof are valid 
  int numDOF = dofs.Size();
  theDofs = new ID(0, numDOF);

  int count = 0;
  int i;
  for (i=0; i<numDOF; i++) {
    int dof = dofs(i);
    if (dof >= 0) {
      (*theDofs)[count] = dof;
      count++;
    } else {
      opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - invalid dof  " << dof;
      opserr << " will be ignored\n";
    }
  }

  //
  // verify the nodes exist 
  //

  count = 0;
  int numNode = nodes.Size();
  theNodes = new ID(1, numNode);
  for (i=0; i<numNode; i++) {
    int nodeTag = nodes(i);
    Node *theNode = theDomain->getNode(nodeTag);
    if (theNode == 0) {
      opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - invalid node  " << nodeTag;
      opserr << " does not exist in domain - will be ignored\n";
    } else {
      (*theNodes)[count++] = nodeTag;
    }
  }

  currentData = new Vector(theNodes->Size()*theDofs->Size());
  data = new Matrix(3, theNodes->Size()*theDofs->Size());
  data->Zero();

  //
  // set the data flag used as aswitch to get the response in a record
  //

  if ((strcmp(dataToStore, "disp") == 0)) {
    dataFlag = 0;
  } else if ((strcmp(dataToStore, "vel") == 0)) {
    dataFlag = 1;
  } else if ((strcmp(dataToStore, "accel") == 0)) {
    dataFlag = 2;
  } else if ((strcmp(dataToStore, "incrDisp") == 0)) {
    dataFlag = 3;
  } else if ((strcmp(dataToStore, "incrDeltaDisp") == 0)) {
    dataFlag = 4;
  } else if ((strncmp(dataToStore, "eigen",5) == 0)) {
    int mode = atoi(&(dataToStore[5]));
    if (mode > 0)
      dataFlag = 10 + mode;
    else
      dataFlag = 6;
  } else {
    dataFlag = 6;
    opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - dataToStore " << dataToStore;
    opserr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }

  //
  // need to create the data description, i.e. what each column of data is
  //

  int numDbColumns = nodes.Size()*dofs.Size();
  char **dbColumns = new char *[numDbColumns];
  static char aColumn[128]; // assumes a column name will not be longer than 256 characters

  int counter = 0;
  for (i=0; i<theNodes->Size(); i++) {
    int nodeTag = (*theNodes)(i);
    for (int j=0; j<theDofs->Size(); j++) {
      int dof = (*theDofs)(j);
      sprintf(aColumn, "Node%d_%s_%d", nodeTag, dataToStore, dof);
      int lenColumn = strlen(aColumn);
      char *newColumn = new char[lenColumn+1];
      strcpy(newColumn, aColumn);
      dbColumns[counter] = newColumn;
      counter++;
    }
  }
  
  //
  // call open in the handler with the data description
  //

  if (theHandler != 0)
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


EnvelopeNodeRecorder::~EnvelopeNodeRecorder()
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

  if (theDofs != 0)
    delete theDofs;
  
  if (theNodes != 0)
    delete theNodes;
  
  if (theHandler != 0)
    delete theHandler;

  if (currentData != 0)
    delete currentData;

  if (data != 0)
    delete data;
}

int 
EnvelopeNodeRecorder::record(int commitTag, double timeStamp)
{
  // now we go get the displacements from the nodes
  int numDOF = theDofs->Size();
  int numNodes = theNodes->Size();
  
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
    
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;
    
    for (int i=0; i<numNodes; i++) {
      int cnt = i*numDOF;
      Node *theNode = theDomain->getNode((*theNodes)(i));
      if (theNode != 0) {
	if (dataFlag == 0) {
	  const Vector &response = theNode->getTrialDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (response.Size() > dof) {
	      (*currentData)(cnt) = response(dof);
	    }else 
	      (*currentData)(cnt) = 0.0;
	    
	    cnt++;
	  }
	} else if (dataFlag == 1) {
	  const Vector &response = theNode->getTrialVel();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (response.Size() > dof) {
	      (*currentData)(cnt) = response(dof);
	    } else 
	      (*currentData)(cnt) = 0.0;
	    
	    cnt++;
	  }
	} else if (dataFlag == 2) {
	  const Vector &response = theNode->getTrialAccel();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (response.Size() > dof) {
	      (*currentData)(cnt) = response(dof);
	    } else 
	      (*currentData)(cnt) = 0.0;
	    
	    cnt++;
	  }
	} else if (dataFlag == 3) {
	  const Vector &response = theNode->getIncrDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (response.Size() > dof) {
	      (*currentData)(cnt) = response(dof);
	    } else 
	      (*currentData)(cnt) = 0.0;
	    
	    cnt++;
	  }
	} else if (dataFlag == 4) {
	  const Vector &response = theNode->getIncrDeltaDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (response.Size() > dof) {
	      (*currentData)(cnt) = response(dof);
	    } else 
	      (*currentData)(cnt) = 0.0;
	    
	    cnt++;
	  }
	} else if (dataFlag > 10) {
	  int mode = dataFlag - 10;
	  int column = mode - 1;
	  const Matrix &theEigenvectors = theNode->getEigenvectors();
	  if (theEigenvectors.noCols() > column) {
	    int noRows = theEigenvectors.noRows();
	    for (int j=0; j<numDOF; j++) {
	      int dof = (*theDofs)(j);
	      if (noRows > dof) {
		(*currentData)(cnt) = theEigenvectors(dof,column);
	      } else 
		(*currentData)(cnt) = 0.0;
	      cnt++;		
	    }
	  } else {
	    for (int j=0; j<numDOF; j++) {
	      (*currentData)(cnt) = 0.0;
	      cnt++;		
	    }
	  }
	}
      }
    }
    
    // check if currentData modifies the saved data
    bool writeIt = false;
    if (first == true) {
      for (int i=0; i<numNodes*numDOF; i++) {
	(*data)(0,i) = (*currentData)(i);
	(*data)(1,i) = (*currentData)(i);
	(*data)(2,i) = fabs((*currentData)(i));
	first = false;
	writeIt = true;
      } 
    } else {
      for (int i=0; i<numNodes*numDOF; i++) {
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

  return 0;
}


void
EnvelopeNodeRecorder::restart(void)
{
  data->Zero();
  first = true;
}





