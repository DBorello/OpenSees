///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:             
// CLASS:            
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class for coupled system
// RETURN:
// VERSION:
// LANGUAGE:          C++.ver >= 3.0
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Xiaoyan Wu
// DATE:              Sept. 2001
//
//  "Coupled system": Solid and fluid coexist.
//                    u-- Solid displacement
//                    p-- Pore pressure
//                    U-- Absolute fluid displacement
//
//
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <ErrorHandler.h>
#include <EightNodeBrick_u_p_U.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addEightNodeBrick_u_p_U(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  char **argv, Domain*theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      g3ErrorHandler->warning("command: element Brick8_u_p_U - no modelbuilder");
      return TCL_ERROR;
  }

  // check the number of arguments is correct
 
  //  element  Brick20N_u_p_U eleTag?  node1? node2? .. node8?  matTag?   bforce1?  bforce2?  bforce3?\n
  //  argv[0]	argv[1]	      argv[2] argv[3]		       argv[11]  argv[12]  argv[13]  argv[14]
  //  porosity?  alpha?  solidDensity? fluidDensity?  x_permeability? y_permeability? z_permeability? 
  //  argv[15]  argv[16]   argv[17]	 argv[18]        argv[19]          argv[20]      argv[21]     
  //  solid_bulk_modulus? fluid_bulk_modulus?   pressure?
  //       argv[22]            argv[23]     	 argv[24] 
  // Xiaoyan added this comments. 01/07/2002

  if ((argc-eleArgStart) < 23) {
      g3ErrorHandler->warning("command: element Brick8_u_p_U - insufficient args - want %s",
          "element Brick8_u_p_U eleTag? node1? node2? ... node8? matTag? bforce1? bforce2? bforce3? \n"
	  " porosity? alpha? solidDensity? fluidDensity? \n"
 	  "permeability_in_x_dir? permeability_in_y_dir? permeability_in_z_dir?"
	  "solid_bulk_modulus? fluid_bulk_modulus? pressure?\n");
      return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleID, matID;
  int nodes[8];					  
  double bodyforces[3], porosity, alpha, solidDensity, fluidDensity;
  double perm_x, perm_y, perm_z, kks, kkf;
  //char *type;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick8_u_p_U - invalid integer tag %s",      
			      argv[1+eleArgStart]);

      return TCL_ERROR;
  }
  
  // read the 8 node tags
  int i;
  for (i=0; i<8; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	  g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid integer tag %s",      
				  eleID, argv[2+i+eleArgStart]);
	  return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[10+eleArgStart], &matID) != TCL_OK) {
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid matID tag %s",      
			      eleID, argv[10+eleArgStart]);      
      return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - no NDMaterial with tag %s exists",
			      eleID, argv[10+eleArgStart]);      
      return TCL_ERROR;      
  }
  
  //type = argv[11+eleArgStart];

  // read the 3 bodyforce accel's 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[11+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	  g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid bodyforces tag %s",      
				  eleID, argv[11+i+eleArgStart]);
	  return TCL_ERROR;
      }
  }

  // now get the porosity
  if (Tcl_GetDouble(interp, argv[14+eleArgStart], &porosity) != TCL_OK) {	       // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid porosity %s",      
    			  eleID, argv[14+eleArgStart]);      
      return TCL_ERROR;
  }  

  // now get the alpha
    if (Tcl_GetDouble(interp, argv[15+eleArgStart], &alpha) != TCL_OK) {	      // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid alpha %s",      
    			  eleID, argv[15+eleArgStart]);      
      return TCL_ERROR;
  }  

  // now get the soldDensity

  if (Tcl_GetDouble(interp, argv[16+eleArgStart], &solidDensity) != TCL_OK) {	      // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid solidDensity %s",      
    			  eleID, argv[16+eleArgStart]);      
      return TCL_ERROR;
  }    
  
  // now get the fludiDensity
     if (Tcl_GetDouble(interp, argv[17+eleArgStart], &fluidDensity) != TCL_OK) {        // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid fluidDensity %s",      
    			  eleID, argv[17+eleArgStart]);      
      return TCL_ERROR;
  }  

  
  // now get the permeability in x direction
     if (Tcl_GetDouble(interp, argv[18+eleArgStart], &perm_x) != TCL_OK) {        // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid permeability in x direction %s",      
    			  eleID, argv[18+eleArgStart]);      
      return TCL_ERROR;
  }  
  // now get the permeability in y direction
     if (Tcl_GetDouble(interp, argv[19+eleArgStart], &perm_y) != TCL_OK) {        // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid permeability in y direction %s",      
    			  eleID, argv[19+eleArgStart]);      
      return TCL_ERROR;
  }  
  // now get the permeability in z direction
     if (Tcl_GetDouble(interp, argv[20+eleArgStart], &perm_z) != TCL_OK) {        // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid permeability in z direction %s",      
    			  eleID, argv[20+eleArgStart]);      
      return TCL_ERROR;
  }  
 
  // now get the bulk modulus of solid
     if (Tcl_GetDouble(interp, argv[21+eleArgStart], &kks) != TCL_OK) {        // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid bulk modulus of solid %s",      
    			  eleID, argv[21+eleArgStart]);      
      return TCL_ERROR;
  }  
  // now get the bulk modulus of fluid
     if (Tcl_GetDouble(interp, argv[22+eleArgStart], &kkf) != TCL_OK) {        // wxy added 01/07/2002
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - invalid bulk modulus of fluid %s",      
    			  eleID, argv[22+eleArgStart]);      
      return TCL_ERROR;
  }  
  // now create the EightNodeBrick and add it to the Domain
  EightNodeBrick_u_p_U *theEle = new EightNodeBrick_u_p_U(eleID,nodes[0], nodes[1], nodes[2], nodes[3], nodes[4],
                                              nodes[5],nodes[6], nodes[7], theMaterial, 
					      bodyforces[0], bodyforces[1], bodyforces[2], 
					      porosity, alpha, solidDensity, fluidDensity, 
					      perm_x, perm_y, perm_z, kks, kkf,0.0);
					      
  if (theEle == 0) {
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - out of memory", eleID);      
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
      g3ErrorHandler->warning("command: element Brick8_u_p_U %d - could not add ele to domain", 
			      eleID);      
      delete theEle;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



