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
                                                                        
// $Revision: 1.2 $
// $Date: 2002-01-09 23:13:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/20nbrick/TclTwentyNodeBrickCommand.cpp,v $
                                                                        
                                                                        
// File: ~/element/TclEightNodeBrickCommand.C
// 
// Written: fmk 
// Created: 07/99
// Modified: 08/01  Zhaohui Yang, Boris Jeremic @ucdavis
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_addEightNodeBrick() 
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <ErrorHandler.h>
#include <TwentyNodeBrick.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addTwentyNodeBrick(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  char **argv, Domain*theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      g3ErrorHandler->warning("command: element Brick20N - no modelbuilder");
      return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 27) {
      g3ErrorHandler->warning("command: element Brick20N - insufficient args - want %s",
          "element Brick20N eleTag? node1? node2? .. node20?  matTag? bforce1? bforce2? bforce3? massDensity?\n");
      return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleID, matID;
  int nodes[20];
  double bodyforces[3], massDensity;
  //char *type;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N - invalid integer tag %s",      
			      argv[1+eleArgStart]);

      return TCL_ERROR;
  }
  
  // read the 20 node tags
  int i;
  for (i=0; i<20; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	  g3ErrorHandler->warning("command: element Brick20N %d - invalid integer tag %s",      
				  eleID, argv[2+i+eleArgStart]);
	  return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[22+eleArgStart], &matID) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N %d - invalid matID tag %s",      
			      eleID, argv[22+eleArgStart]);      
      return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
      g3ErrorHandler->warning("command: element Brick20N %d - no NDMaterial with tag %s exists",
			      eleID, argv[22+eleArgStart]);      
      return TCL_ERROR;      
  }
  
  //type = argv[11+eleArgStart];

  // read the 3 bodyforce accel's 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[23+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	  g3ErrorHandler->warning("command: element Brick20N %d - invalid bodyforces tag %s",      
				  eleID, argv[23+i+eleArgStart]);
	  return TCL_ERROR;
      }
  }

  // now get the massDensity
  if (Tcl_GetDouble(interp, argv[26+eleArgStart], &massDensity) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N %d - invalid massDensity %s",      
    			  eleID, argv[26+eleArgStart]);      
      return TCL_ERROR;
  }  
  
  // now create the EightNodeBrick and add it to the Domain
  TwentyNodeBrick *theEle = new TwentyNodeBrick(eleID, nodes[ 0], nodes [1], nodes[ 2], nodes[ 3], nodes[ 4],
                                                      nodes[ 5], nodes [6], nodes[ 7], nodes[ 8], nodes[ 9],
                                                      nodes[10], nodes[11], nodes[12], nodes[13], nodes[14],
                                                      nodes[15], nodes[16], nodes[17], nodes[18], nodes[19],
					       theMaterial, bodyforces[0], bodyforces[1], bodyforces[2], massDensity, 0.0);
					      
  if (theEle == 0) {
      g3ErrorHandler->warning("command: element Brick20N %d - out of memory", eleID);      
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
      g3ErrorHandler->warning("command: element Brick20N %d - could not add ele to domain", 
			      eleID);      
      delete theEle;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



