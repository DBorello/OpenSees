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
                                                                        
// $Revision: 1.18 $
// $Date: 2004-11-13 00:57:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/NodeRecorder.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class definition for NodeRecorder.
// A NodeRecorder is used to record the specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) NodeRecorder.C, revA"

#include <NodeRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <FE_Datastore.h>
#include <DataOutputHandler.h>

#include <string.h>

NodeRecorder::NodeRecorder(const ID &dofs, 
			   const ID &nodes, 
			   int psensitivity,
			   const char *dataToStore,
			   Domain &theDom,
			   DataOutputHandler &theOutputHandler,
			   double dT,
			   int startFlag)
:theDofs(0), theNodes(0), response(1 + nodes.Size()*dofs.Size()), 
 theDomain(&theDom), theHandler(&theOutputHandler),
 flag(startFlag), dataFlag(0), 
 deltaT(dT), nextTimeStampToRecord(0.0), 
 sensitivity(psensitivity)
{
  //
  // store copy of dof's to be recorder, verifying dof are valid, i.e. >= 0
  //

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
      opserr << "NodeRecorder::NodeRecorder - invalid dof  " << dof;
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
      opserr << "NodeRecorder::NodeRecorder - invalid node  " << nodeTag;
      opserr << " does not exist in domain - will be ignored\n";
    } else {
      (*theNodes)[count++] = nodeTag;
    }
  }

  response.Zero();

  //
  // set the data flag used as aswitch to get the response in a record
  //

  if (dataToStore == 0 || (strcmp(dataToStore, "disp") == 0)) {
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
    opserr << "NodeRecorder::NodeRecorder - dataToStore " << dataToStore;
    opserr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }

  //
  // need to create the data description, i.e. what each column of data is
  //

  int numDbColumns = 1 + nodes.Size()*dofs.Size();
  char **dbColumns = new char *[numDbColumns];
    
  static char aColumn[128]; // assumes a column name will not be longer than 256 characters
    
  char *newColumn = new char[5];
  sprintf(newColumn, "%s","time");  
  dbColumns[0] = newColumn;
  
  int counter = 1;
  for (i=0; i<theNodes->Size(); i++) {
    int nodeTag = (*theNodes)(i);
    for (int j=0; j<theDofs->Size(); j++) {
      int dof = (*theDofs)(j);
      sprintf(aColumn, "Node%d_%s_%d", nodeTag, dataToStore, dof+1);
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




NodeRecorder::~NodeRecorder()
{
  if (theHandler != 0)
    delete theHandler;

  if (theDofs != 0)
    delete theDofs;
  
  if (theNodes != 0)
    delete theNodes;
}

int 
NodeRecorder::record(int commitTag, double timeStamp)
{
  if (theHandler == 0) {
    opserr << "NodeRecorder::record() - no DataOutputHandler has been set\n";
    return -1;
  }


  int numDOF = theDofs->Size();
  int numNodes = theNodes->Size();
  
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {

    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    //
    // now we go get the responses from the nodes & place them in disp vector
    //

    for (int i=0; i<numNodes; i++) {
      int cnt = i*numDOF + 1; 
      Node *theNode = theDomain->getNode((*theNodes)(i));
      if (theNode != 0) {
	if (dataFlag == 0) {
	  // AddingSensitivity:BEGIN ///////////////////////////////////
	  if (sensitivity==0) {
	    const Vector &theResponse = theNode->getTrialDisp();
	    for (int j=0; j<numDOF; j++) {
	      int dof = (*theDofs)(j);
	      if (theResponse.Size() > dof) {
		response(cnt) = theResponse(dof);
	      }
	      else {
		response(cnt) = 0.0;
	      }
	      cnt++;
	    }
	  }
	  else {
	    for (int j=0; j<numDOF; j++) {
	      int dof = (*theDofs)(j);
	      response(cnt) = theNode->getDispSensitivity(dof+1, sensitivity);
	    }
	  }
	  // AddingSensitivity:END /////////////////////////////////////
	} else if (dataFlag == 1) {
	  const Vector &theResponse = theNode->getTrialVel();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else 
	      response(cnt) = 0.0;
	    
	    cnt++;
	  }
	} else if (dataFlag == 2) {
	  const Vector &theResponse = theNode->getTrialAccel();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else 
	      response(cnt) = 0.0;
	    
	    cnt++;
	  }
	} else if (dataFlag == 3) {
	  const Vector &theResponse = theNode->getIncrDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else 
	      response(cnt) = 0.0;
	    
	    cnt++;
	  }
	} else if (dataFlag == 4) {
	  const Vector &theResponse = theNode->getIncrDeltaDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else 
	      response(cnt) = 0.0;
	    
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
		response(cnt) = theEigenvectors(dof,column);
	      } else 
		response(cnt) = 0.0;
	      cnt++;		
	    }
	  } else {
	    for (int j=0; j<numDOF; j++) {
	      response(cnt) = 0.0;
	    }
	  }
	}
      }
    }
    
    // insert the data into the database
    response(0) = timeStamp;
    theHandler->write(response);
  }
    
  return 0;
}

void
NodeRecorder::restart(void)
{
  opserr << "WARNING NodeRecorder::restart() - does nothing, OutputHandlers will contain old & new data\n";
  opserr << "   - typically most users want to use \"remove -recorders\" before this as want to save all data\n";
}





