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
// $Date: 2002-04-12 22:48:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/upU/TclTwentyNodeBrick_u_p_U.cpp,v $
                                                                        
                                                                        
// File: ~/element/upU/TclTwentyNodeBrick_u_p_U.cpp
// 
// Written: fmk 
// Created: 07/99
// Modified: 01/02  Xiaoyan Wu, Boris Jeremic @ucdavis
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_addTwentyNodeBrick_u_p_U() 
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <ErrorHandler.h>
#include <TwentyNodeBrick.h>
#include <TwentyNodeBrick_u_p_U.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addTwentyNodeBrick_u_p_U(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  char **argv, Domain*theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U - no modelbuilder");
      return TCL_ERROR;
  }

  // check the number of arguments is correct
  // should have 36 arguments.  
  //  element  Brick20N_u_p_U eleTag?  node1? node2? .. node20?  matTag?   bforce1?  bforce2?  bforce3?\n
  //  argv[0]	argv[1]	      argv[2] argv[3]			 argv[23]  argv[24]  argv[25]  argv[26]
  //  porosity?  alpha?  solidDensity? fluidDensity? x_permeability? y_permeability? z_permeability?
  //  argv[27]  argv[28]   argv[29]	 argv[30]       argv[31]          argv[32]      argv[33]
  //  solid_bulk_modulus? fluid_bulk_modulus?   pressure?
  //       argv[34]            argv[35]     	 argv[36] 
  // Xiaoyan added this comments. 01/07/2002
  if ((argc-eleArgStart) < 35) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U - insufficient args - want %s",
          "element Brick20N_u_p_U eleTag? node1? node2? .. node20?  matTag? bforce1? bforce2? bforce3?\n" 
	  "porosity? alpha?  solidDensity? fluidDensity? \n"
 	  "permeability_in_x_dir? permeability_in_y_dir? permeability_in_z_dir?"
	  "solid_bulk_modulus? fluid_bulk_modulus? pressure?\n");
      return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleID, matID;
  int nodes[20];
  double bodyforces[3], porosity, alpha, solidDensity, fluidDensity;
  double perm_x, perm_y, perm_z, kks,kkf;
  //char *type;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U - invalid integer tag %s",      
			      argv[1+eleArgStart]);

      return TCL_ERROR;
  }
  
  // read the 20 node tags
  int i;
  for (i=0; i<20; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	  g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid integer tag %s",      
				  eleID, argv[2+i+eleArgStart]);
	  return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[22+eleArgStart], &matID) != TCL_OK) {	     // argv[23]  Xiaoyan 01/07/2002
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid matID tag %s",      
			      eleID, argv[22+eleArgStart]);      
      return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - no NDMaterial with tag %s exists",
			      eleID, argv[22+eleArgStart]);      
      return TCL_ERROR;      
  }
  
  //type = argv[11+eleArgStart];

  // read the 3 bodyforce accel's 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[23+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	  g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid bodyforces tag %s",      
				  eleID, argv[23+i+eleArgStart]);
	  return TCL_ERROR;
      }
  }

  // now get the porosity
  if (Tcl_GetDouble(interp, argv[26+eleArgStart], &porosity) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid porosity %s",      
    			  eleID, argv[26+eleArgStart]);      
      return TCL_ERROR;
  } 
 
  // now get the alpha for solid alpha=1.0
   if (Tcl_GetDouble(interp, argv[27+eleArgStart], &alpha) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid alpha %s",      
    			  eleID, argv[27+eleArgStart]);      
      return TCL_ERROR;
  }  
 
   // now get the solidDensity
   if (Tcl_GetDouble(interp, argv[28+eleArgStart], &solidDensity) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid solidDensity %s",      
    			  eleID, argv[28+eleArgStart]);      
      return TCL_ERROR;
  }  
 
   // now get the fluidDensity
   if (Tcl_GetDouble(interp, argv[29+eleArgStart], &fluidDensity) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid fluidDensity %s",      
    			  eleID, argv[29+eleArgStart]);      
      return TCL_ERROR;
  }  
  // permeability in x direction
  if (Tcl_GetDouble(interp, argv[30+eleArgStart], &perm_x) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid permeability in x direction %s",      
    			  eleID, argv[30+eleArgStart]);      
      return TCL_ERROR;
  }  
  // permeability in y direction
  if (Tcl_GetDouble(interp, argv[31+eleArgStart], &perm_y) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid permeability in y direction %s",      
    			  eleID, argv[31+eleArgStart]);      
      return TCL_ERROR;
  }  
  // permeability in z direction
  if (Tcl_GetDouble(interp, argv[32+eleArgStart], &perm_z) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - invalid permeability in z direction %s",      
    			  eleID, argv[32+eleArgStart]);      
      return TCL_ERROR;
  }  
  // now get the bulk modulus of solid
     if (Tcl_GetDouble(interp, argv[33+eleArgStart], &kks) != TCL_OK) {        // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid bulk modulus of solid %s",      
    			  eleID, argv[33+eleArgStart]);      
      return TCL_ERROR;
  }  
  // now get the bulk modulus of fluid
     if (Tcl_GetDouble(interp, argv[34+eleArgStart], &kkf) != TCL_OK) {        // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid bulk modulus of fluid %s",      
    			  eleID, argv[34+eleArgStart]);      
      return TCL_ERROR;
  }  
  
  // now create the EightNodeBrick and add it to the Domain
  TwentyNodeBrick_u_p_U *theEle = new TwentyNodeBrick_u_p_U(eleID, nodes[ 0], nodes [1], nodes[ 2], nodes[ 3], nodes[ 4],
                                                      nodes[ 5], nodes [6], nodes[ 7], nodes[ 8], nodes[ 9],
                                                      nodes[10], nodes[11], nodes[12], nodes[13], nodes[14],
                                                      nodes[15], nodes[16], nodes[17], nodes[18], nodes[19],
					              theMaterial, bodyforces[0], bodyforces[1], bodyforces[2], 
						      porosity, alpha, solidDensity, fluidDensity,
						      perm_x, perm_y, perm_z, kks, kkf, 0.0);
					      
  if (theEle == 0) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - out of memory", eleID);      
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
      g3ErrorHandler->warning("command: element Brick20N_u_p_U %d - could not add ele to domain", 
			      eleID);      
      delete theEle;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



