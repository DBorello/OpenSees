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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-01-30 04:45:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/8nbrick/TclEightNodeBrickCommand.cpp,v $
                                                                        
                                                                        
// File: ~/element/TclEightNodeBrickCommand.C
// 
// Written: fmk 
// Created: 07/99
// Modified: 01/01  Zhaohui Yang, Boris Jeremic @ucdavis
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
#include <EightNodeBrick.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addEightNodeBrick(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  char **argv, Domain*theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      g3ErrorHandler->warning("command: element brick - no modelbuilder");
      return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 15) {
      g3ErrorHandler->warning("command: element brick - insufficient args - want %s",
          "element brick eleTag? node1? node2? .. node8? matTag? bforce1? bforce2? bforce3? massDensity?\n");
      return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleID, matID;
  int nodes[8];
  double bodyforces[3], massDensity;
  //char *type;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
      g3ErrorHandler->warning("command: element brick - invalid integer tag %s",      
			      argv[1+eleArgStart]);

      return TCL_ERROR;
  }
  
  // read the 8 node tags
  int i;
  for (i=0; i<8; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	  g3ErrorHandler->warning("command: element brick %d - invalid integer tag %s",      
				  eleID, argv[2+i+eleArgStart]);
	  return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[10+eleArgStart], &matID) != TCL_OK) {
      g3ErrorHandler->warning("command: element brick %d - invalid matID tag %s",      
			      eleID, argv[10+eleArgStart]);      
      return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
      g3ErrorHandler->warning("command: element brick %d - no NDMaterial with tag %s exists",
			      eleID, argv[10+eleArgStart]);      
      return TCL_ERROR;      
  }
  
  //type = argv[11+eleArgStart];

  // read the 3 bodyforce accel's 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[11+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	  g3ErrorHandler->warning("command: element brick %d - invalid bodyforces tag %s",      
				  eleID, argv[11+i+eleArgStart]);
	  return TCL_ERROR;
      }
  }

  // now get the massDensity
  if (Tcl_GetDouble(interp, argv[14+eleArgStart], &massDensity) != TCL_OK) {
      g3ErrorHandler->warning("command: element brick %d - invalid massDensity %s",      
    			  eleID, argv[14+eleArgStart]);      
      return TCL_ERROR;
  }  
  
  // now create the EightNodeBrick and add it to the Domain
  EightNodeBrick *theEle = new EightNodeBrick(eleID,
	                                             nodes[0], 
	                                             nodes[1], 
																																													 nodes[2], 
																																													 nodes[3], 
																																													 nodes[4],
                                              nodes[5],
																																													 nodes[6], 
																																													 nodes[7], 
																																													 theMaterial, 
                                   					      bodyforces[0], 
																																													 bodyforces[1],
																																													 bodyforces[2], 
																																													 massDensity, 
																																													 0.0);
					      
  if (theEle == 0) {
      g3ErrorHandler->warning("command: element brick %d - out of memory", eleID);      
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
      g3ErrorHandler->warning("command: element brick %d - could not add ele to domain", 
			      eleID);      
      delete theEle;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}



