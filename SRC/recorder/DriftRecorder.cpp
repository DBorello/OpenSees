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
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DriftRecorder.cpp,v $

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for DriftRecorder.

#include <math.h>

#include <DriftRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <DataOutputHandler.h>
#include <string.h>

DriftRecorder::DriftRecorder(int ni, int nj, int df, int dirn,
			     Domain &theDom, DataOutputHandler &theDataOutputHandler)
  :ndI(1), ndJ(1), dof(df), perpDirn(dirn), oneOverL(1), data(2),
   theDomain(&theDom), theHandler(&theDataOutputHandler)
{
  ndI(0) = ni;
  ndJ(0) = nj;
  
  Node *nodeI = theDomain->getNode(ni);
  Node *nodeJ = theDomain->getNode(nj);

  if (nodeI == 0 || nodeJ == 0) {
    opserr << "DriftRecorder::DriftRecorder-- node " << nodeI << " or " << nodeJ << " not in domain\n";    
    oneOverL(0) = 0.0;
  }

  int numDbColumns = 2;
  char **dbColumns = new char *[numDbColumns];
  static char aColumn[128]; // assumes a column name will not be longer than 256 characters
  
  char *timeColumn = new char[5];
  sprintf(timeColumn, "%s","time");  
  dbColumns[0] = timeColumn;
  
  const Vector &crdI = nodeI->getCrds();
  const Vector &crdJ = nodeJ->getCrds();

  if (crdI(dirn) == crdJ(dirn)) {
    opserr << "DriftRecorder::DriftRecorder-- Nodal projection has zero component along chosen direction\n";
    oneOverL(0) = 0.0;
  }
  else 
    oneOverL(0) = 1.0/fabs(crdJ(dirn) - crdI(dirn));

  sprintf(aColumn, "Drift%d_%d_%d", ni, nj, dirn);
  int lenColumn = strlen(aColumn);
  char *newColumn = new char[lenColumn+1];
  strcpy(newColumn, aColumn);
  dbColumns[1] = newColumn;
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


DriftRecorder::DriftRecorder(const ID &nI, const ID &nJ, int df, int dirn,
			     Domain &theDom, DataOutputHandler &theDataOutputHandler)
  :ndI(nI), ndJ(nJ), dof(df), perpDirn(dirn), oneOverL(nI.Size()), data(nI.Size()+1),
   theDomain(&theDom), theHandler(&theDataOutputHandler)
 {
   int ndIsize = ndI.Size();
   int ndJsize = ndJ.Size();
   if (ndIsize != ndJsize) {
     opserr << "FATAL - DriftRecorder::DriftRecorder() - error node arrays differ in size\n";
     exit(-1);
   }

   //
   // create the data description for the OutputHandler & compute length between nodes
   //

   int numDbColumns = 1 + ndIsize;
   char **dbColumns = new char *[numDbColumns];
   
   static char aColumn[128]; // assumes a column name will not be longer than 256 characters
   
   char *newColumn = new char[5];
   sprintf(newColumn, "%s","time");  
   dbColumns[0] = newColumn;

   for (int i=0; i<ndIsize; i++) {
     int ni = ndI(i);
     int nj = ndJ(i);
     
     Node *nodeI = theDomain->getNode(ni);
     Node *nodeJ = theDomain->getNode(nj);
     
     if (nodeI == 0 || nodeJ == 0) {
       opserr << "DriftRecorder::DriftRecorder-- node " << nodeI << " or " << nodeJ << " not in domain\n";    
       oneOverL(0) = 0.0;
     }
     
     const Vector &crdI = nodeI->getCrds();
     const Vector &crdJ = nodeJ->getCrds();
     
     if (crdI(dirn) == crdJ(dirn)) {
       opserr << "DriftRecorder::DriftRecorder-- Nodal projection has zero component along chosen direction\n";
       oneOverL(i) = 0.0;
     }
     else 
       oneOverL(i) = 1.0/fabs(crdJ(dirn) - crdI(dirn));

     sprintf(aColumn, "Drift%d_%d_%d", ni, nj, dirn);
     int lenColumn = strlen(aColumn);
     char *newColumn = new char[lenColumn+1];
     strcpy(newColumn, aColumn);
     dbColumns[i+1] = newColumn;
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

DriftRecorder::~DriftRecorder()
{
  if (theHandler != 0)
    delete theHandler;
}

int 
DriftRecorder::record(int commitTag, double timeStamp)
{
  data(0) = theDomain->getCurrentTime();
  
  for (int i=0; i<ndI.Size(); i++) {
    if (oneOverL(i) != 0.0) {
      Node *nodeI = theDomain->getNode(ndI(i));
      Node *nodeJ = theDomain->getNode(ndJ(i));
      
      const Vector &dispI = nodeI->getTrialDisp();
      const Vector &dispJ = nodeJ->getTrialDisp();
      
      double dx = dispJ(dof)-dispI(dof);
      
      
      data(i+1) =  dx*oneOverL(i);
    }
    else
      data(i+1) = 0.0;
  }

  theHandler->write(data);
  return 0;
}

void
DriftRecorder::restart(void)
{
  opserr << "DriftRecorder::restart() - should not be called\n";
  //  theHandler->restart();
}
