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
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DriftRecorder.h,v $
                                                                        
#ifndef DriftRecorder_h
#define DriftRecorder_h

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for 
// DriftRecorder. 

#include <Recorder.h>
#include <ID.h>
#include <Vector.h>

class Domain;
class DataOutputHandler;

class DriftRecorder: public Recorder
{
 public:
  DriftRecorder(int ndI, int ndJ, int dof, int perpDirn,
		Domain &theDomain, 
		DataOutputHandler &theHandler);

  DriftRecorder(const ID &ndI, const ID &ndJ, int dof, int perpDirn,
		Domain &theDomain, 
		DataOutputHandler &theHandler);
  
  ~DriftRecorder();

  int record(int commitTag, double timeStamp);
  void restart(void);    
  
 protected:
  
 private:	
  ID ndI;
  ID ndJ;
  int dof;
  int perpDirn;
  Vector oneOverL;
  Vector data;

  Domain *theDomain;
  DataOutputHandler *theHandler;
};

#endif
