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
                                                                        
// $Revision: 1.10 $
// $Date: 2004-11-13 00:57:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/NodeRecorder.h,v $
                                                                        
#ifndef NodeRecorder_h
#define NodeRecorder_h

// Written: fmk 
// Created: 09/00
// Revision: A
//
// Description: This file contains the class definition for 
// NodeRecorder. A NodeRecorder is used to store the specified nodal dof responses
// for the specified nodes in a file.
//
// What: "@(#) NodeRecorder.h, revA"


#include <Recorder.h>
#include <ID.h>
#include <Vector.h>

class Domain;
class FE_Datastore;
class DataOutputHandler;

class NodeRecorder: public Recorder
{
  public:
  NodeRecorder(const ID &theDof, 
	       const ID &theNodes, 
	       int sensitivity,
	       const char *dataToStore,
	       Domain &theDomain,
	       DataOutputHandler &theOutputHandler,
	       double deltaT = 0.0,
	       int startFlag = 0); 
    
    ~NodeRecorder();

    int record(int commitTag, double timeStamp);
    void restart(void);    
    
  protected:
    
  private:	
    ID *theDofs;
    ID *theNodes;
    Vector response;

    Domain *theDomain;
    DataOutputHandler *theHandler;

    int flag;   // flag indicating whether time, load factor or nothing printed
	        // at start of each line in file

    int dataFlag; // flag indicating what it is to be stored in recorder

    double deltaT;
    double nextTimeStampToRecord;

    // AddingSensitivity:BEGIN //////////////////////////////
    int sensitivity;
    // AddingSensitivity:END ////////////////////////////////
};

#endif
