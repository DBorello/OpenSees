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
                                                                        
// $Revision: 1.5 $
// $Date: 2004-11-13 00:57:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeElementRecorder.h,v $
                                                                        
#ifndef EnvelopeElementRecorder_h
#define EnvelopeElementRecorder_h

// Written: fmk 
//
// What: "@(#) EnvelopeElementRecorder.h, revA"

#include <Recorder.h>
#include <Information.h>
#include <OPS_Globals.h>
#include <ID.h>


class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;
class DataOutputHandler;

class EnvelopeElementRecorder: public Recorder
{
  public:
    EnvelopeElementRecorder(const ID &eleID, 
			    const char **argv, 
			    int argc,
			    Domain &theDomain, 
			    DataOutputHandler &theOutputHandler,
			    double deltaT = 0.0);


    ~EnvelopeElementRecorder();

    int record(int commitTag, double timeStamp);
    void restart(void);    
    
  protected:
    
  private:	
    int numEle;
    ID responseID;                 // integer element returns in setResponse
    Response **theResponses;

    Domain *theDomain;
    DataOutputHandler *theHandler;

    double deltaT;
    double nextTimeStampToRecord;

    Matrix *data;
    Vector *currentData;
    bool first;
};


#endif
