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
                                                                        
// $Revision: 1.4 $
// $Date: 2002-04-02 18:49:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ElementRecorder.h,v $
                                                                        
                                                                        
#ifndef ElementRecorder_h
#define ElementRecorder_h

// File: ~/recorder/ElementRecorder.h
//
// Written: fmk 
// Created: 09/99
// Revision: A
//
// Description: This file contains the class definition for ElementRecorder.
// A ElementRecorder is used to obtain a response from an element during 
// the analysis.
//
// What: "@(#) ElementRecorder.h, revA"

#include <Recorder.h>
#include <fstream.h>

#include <Information.h>
#include <G3Globals.h>
#include <ID.h>


class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;

class ElementRecorder: public Recorder
{
  public:
    ElementRecorder(const ID &eleID, Domain &theDomain, char **argv, int argc,
		    bool echoTime, double deltaT = 0.0, char *fileName =0);

    ElementRecorder(const ID &eleID, Domain &theDomain, char **argv, int argc,
		    bool echoTime, FE_Datastore *db, char *tableName, double deltaT = 0.0);

    ~ElementRecorder();
    int record(int commitTag, double timeStamp);
    int playback(int commitTag);

    void restart(void);    
    
  protected:
    
  private:	
    int numEle;
    ID responseID;                 // integer element returns in setResponse
    Information	*eleInfoObjects;   // object into whcih element places the response
    Element **theElements;         // pointer to the elements

    Response **theResponses;

    Domain *theDomain;
    bool echoTimeFlag;             // flag indicating if pseudo time also printed
    char *fileName;                // file name  
    ofstream theFile; 	           // output stream

    double deltaT;
    double nextTimeStampToRecord;

    FE_Datastore *db;
    char **dbColumns;
    int numDbColumns;
    Vector *data;
};


#endif
