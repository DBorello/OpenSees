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
// $Date: 2001-05-16 04:19:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/AlgorithmIncrements.h,v $
                                                                        
                                                                        
// File: ~/recorder/tcl/AlgorithmIncrements.h.h
// 
// Written: fmk 
// Created: 01/01
// Revision: A
//
// Description: This file contains the class definition for AlgorithmIncrements.
// A AlgorithmIncrements will display the X and B in the SOE associated with the
// algorithm on a record.

//
// What: "@(#) ModelBuilder.h, revA"

#ifndef AlgorithmIncrements_h
#define AlgorithmIncrements_h

#include <Recorder.h>
#include <G3Globals.h>


#include <fstream.h>

class EquiSolnAlgo;
class Renderer;
class ColorMap;
class ID;
class Vector;

class AlgorithmIncrements : public Recorder
{
  public:
    AlgorithmIncrements(EquiSolnAlgo *theAlgo,
			char *windowTitle, 
			int xLoc, int yLoc, int width, int height,
			bool displayRecord = false,
			char *fileName = 0);
    
    ~AlgorithmIncrements();    

    int plotData(const Vector &X, const Vector &B);

    int record(int cTag);
    int playback(int cTag);
    void restart(void);    

  protected:

  private:
    ColorMap *theMap;
    Renderer *theRenderer;
    EquiSolnAlgo *theAlgo;

    int numRecord;
    bool displayRecord;
    char *fileName;
    ofstream theFile;     
};

#endif







