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
                                                                        


// $Revision: 1.8 $
// $Date: 2002-04-03 00:06:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/NodeRecorder.cpp,v $
                                                                        

// File: ~/recorder/NodeRecorder.C
//
// Written: fmk 
// Created: 11/98
// Revision: A
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

#include <iostream.h>
#include <fstream.h>
#include <string.h>

NodeRecorder::NodeRecorder(const ID &dofs, 
			   const ID &nodes, 
			   Domain &theDom,
			   char *theFileName,
			   char *dataToStore,
			   double dT,
			   int startFlag)
:theDofs(dofs), theNodes(nodes), disp(1 + nodes.Size()*dofs.Size()), 
  theDomain(&theDom), flag(startFlag), deltaT(dT), nextTimeStampToRecord(0.0), db(0), dbColumns(0)
{
  // create char array to store file name
  int fileNameLength = strlen(theFileName) + 1;
  fileName = new char[fileNameLength];
  if (fileName == 0) {
    g3ErrorHandler->fatal("FileNodeDispRecorder::FileNodeDispRecorder - out of memory creating string %d long\n",
			  fileNameLength);
  }

  // copy the strings
  strcpy(fileName, theFileName);    

  // open the file
  theFile.open(fileName, ios::out);
  if (theFile.bad()) {
    cerr << "WARNING - FileNodeDispRecorder::FileNodeDispRecorder()";
    cerr << " - could not open file " << fileName << endl;
  }    

  disp.Zero();

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
    cerr << "NodeRecorder::NodeRecorder - dataToStore " << dataToStore;
    cerr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }
}


NodeRecorder::NodeRecorder(const ID &dofs, 
			   const ID &nodes, 
			   Domain &theDom,
			   FE_Datastore *database,
			   char *dbTable,
			   char *dataToStore,
			   double dT,
			   int startFlag)
:theDofs(dofs), theNodes(nodes), disp(1 + nodes.Size()*dofs.Size()), 
  theDomain(&theDom), flag(startFlag), deltaT(dT), nextTimeStampToRecord(0.0), 
  db(database), dbColumns(0)
{

  // create char array to store table name
  int fileNameLength = strlen(dbTable) + 1;
  fileName = new char[fileNameLength];
  if (fileName == 0) {
    g3ErrorHandler->fatal("FileNodeDispRecorder::FileNodeDispRecorder - out of memory creating string %d long\n",
			  fileNameLength);
  }

  // copy the strings
  strcpy(fileName, dbTable);

  disp.Zero();

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
    cerr << "NodeRecorder::NodeRecorder - dataToStore " << dataToStore;
    cerr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }

  // now create the columns strings for the database
  numDbColumns = 1 + nodes.Size()*dofs.Size();
  dbColumns = new char *[numDbColumns];

  char aColumn[256]; // assumes a column name will not be longer than 256 characters
  
  char *newColumn = new char[5];
  sprintf(newColumn, "%s","time");  
  dbColumns[0] = newColumn;
  
  int counter = 1;
  for (int i=0; i<nodes.Size(); i++) {
    int nodeTag = nodes(i);
    for (int j=0; j<dofs.Size(); j++) {
      int dof = dofs(j);
      sprintf(aColumn, "%s_%d_%d",dataToStore,nodeTag,dof);
      int lenColumn = strlen(aColumn+1);
      char *newColumn = new char[lenColumn];
      sprintf(newColumn, "%s",aColumn);
      dbColumns[counter] = newColumn;
      counter++;
    }
  }

  // create the table in the database
  db->createTable(dbTable, numDbColumns, dbColumns);
}


NodeRecorder::~NodeRecorder()
{
    if (!theFile)
	theFile.close();

    if (fileName != 0)
      delete fileName;

    if (dbColumns != 0) {
      for (int i=0; i<numDbColumns; i++)
	delete [] dbColumns[i];
      
      delete [] dbColumns;
    }
}

int 
NodeRecorder::record(int commitTag, double timeStamp)
{
    // now we go get the displacements from the nodes
    int numDOF = theDofs.Size();
    int numNodes = theNodes.Size();
    
    if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
      
      if (deltaT != 0.0) 
	nextTimeStampToRecord = timeStamp + deltaT;

      for (int i=0; i<numNodes; i++) {
	int cnt = i*numDOF + 1;
	Node *theNode = theDomain->getNode(theNodes(i));
	if (theNode != 0) {
	  if (dataFlag == 0) {
	    const Vector &theDisp = theNode->getTrialDisp();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		}else 
		  disp(cnt) = 0.0;

		cnt++;
	    }
	  } else if (dataFlag == 1) {
	    const Vector &theDisp = theNode->getTrialVel();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		} else 
		  disp(cnt) = 0.0;

		cnt++;
	    }
	  } else if (dataFlag == 2) {
	    const Vector &theDisp = theNode->getTrialAccel();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		} else 
		  disp(cnt) = 0.0;

		cnt++;
	    }
	  } else if (dataFlag == 3) {
	    const Vector &theDisp = theNode->getIncrDisp();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		} else 
		  disp(cnt) = 0.0;

		cnt++;
	    }
	  } else if (dataFlag == 4) {
	    const Vector &theDisp = theNode->getIncrDeltaDisp();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		} else 
		  disp(cnt) = 0.0;

		cnt++;
	    }
	  } else if (dataFlag > 10) {
	    int mode = dataFlag - 10;
	    int column = mode - 1;
	    const Matrix &theEigenvectors = theNode->getEigenvectors();
	    if (theEigenvectors.noCols() > column) {
	      int noRows = theEigenvectors.noRows();
	      for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (noRows > dof) {
		  disp(cnt) = theEigenvectors(dof,column);
		} else 
		  disp(cnt) = 0.0;
		cnt++;		
	      }
	    } else {
	      for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		disp(cnt) = 0.0;
	      }
	    }
	  }
	}
      }

      if (db == 0) {      

	// write them to the file
	if (flag == 1)
	  theFile << timeStamp << " ";
	else if (flag == 2)
	  theFile << timeStamp << " ";
	
	for (int j=1; j<=numNodes*numDOF; j++)
	  theFile << disp(j) << " ";
      
	theFile << endl;
	theFile.flush();

      } else {

	// insert the data into the database
	disp(0) = timeStamp;
	db->insertData(fileName, dbColumns, commitTag, disp);
      }
    }    

    return 0;
}

int 
NodeRecorder::playback(int commitTag)
{
  if (theFile.bad())
    return 0;
    
  // close o/p file to ensure all buffered data gets written to file
  theFile.close(); 

  int numDOF = theDofs.Size();
  int numNodes = theNodes.Size();  
  
  // open a stream for reading from the file
  ifstream inputFile;
  inputFile.open(fileName, ios::in);
  if (inputFile.bad()) {
    cerr << "WARNING - FileNodeDispRecorder::playback() - could not open file ";
    cerr << fileName << endl;
    return -1;
  }   

  double data;
  // read file up until line we want
  for (int i=0; i<(commitTag-1); i++)
    // now read in a line
    if (flag == 1 || flag == 2) {
	inputFile >> data;
	for (int i=0; i<numNodes*numDOF; i++)
	    inputFile >> data;
    }

  // now read in our line and print out
  if (flag == 1 || flag == 2) {
      inputFile >> data;
      cerr << data << " ";
      for (int i=0; i<numNodes*numDOF; i++) {
	  inputFile >> data;
	  cerr << data << " ";
      }	
      cerr << endl;
  }
  inputFile.close();
      
    // open file again for writing
    theFile.open(fileName, ios::app);
    if (theFile.bad()) {
      cerr << "WARNING - FileNodeDispRecorder::playback() - could not open file ";
      cerr << fileName << endl;
      return -1;
    }    
  
  // does nothing
  return 0;
}

void
NodeRecorder::restart(void)
{
  theFile.close();
  theFile.open(fileName, ios::out);
  if (theFile.bad()) {
    cerr << "WARNING - FileNodeDispRecorder::restart() - could not open file ";
    cerr << fileName << endl;
  }
}





