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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-06-18 23:18:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/PatternRecorder.h,v $
                                                                        
#ifndef PatternRecorder_h
#define PatternRecorder_h

// File: ~/recorder/PatternRecorder.h
//
// Written: fmk 
// Created: 09/00
// Revision: A
//
// Description: This file contains the class definition for 
// PatternRecorder. A PatternRecorder is used to store the specified nodal dof responses
// for the specified nodes in a file.
//
// What: "@(#) PatternRecorder.h, revA"


#include <Recorder.h>
#include <ID.h>
#include <Vector.h>

#include <fstream>
using std::ofstream;
class Domain;
class FE_Datastore;

class PatternRecorder: public Recorder
{
  public:
    PatternRecorder(int thePattern,
		 Domain &theDomain,
		 char *fileName,
		 double deltaT = 0.0,
		 int startFlag = 0); 

    ~PatternRecorder();
    int record(int commitTag, double timeStamp);
    int playback(int commitTag);
    void restart(void);    
    
  protected:
    
  private:	
	int thePattern;
    Domain *theDomain;

    int flag;   // flag indicating whether time, load factor or nothing printed
	        // at start of each line in file
    char *fileName;
    ofstream theFile;     

    double deltaT;
    double nextTimeStampToRecord;
};

#endif
