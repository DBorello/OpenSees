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
// $Date: 2005-08-22 20:50:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TclReinforcingSteel.cpp,v $
                                                                        

#include <TclModelBuilder.h>
#include <ReinforcingSteel.h>   // Jon Mohle

#include <Vector.h>
#include <string.h>
#include <tcl.h>

int
TclCommand_ReinforcingSteel(ClientData clientData, Tcl_Interp *interp, int argc, 
			    TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
  UniaxialMaterial *theMaterial = 0;

  if (argc < 10) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial ReinforcingSteel tag? fy? fu? Es? Esh? esh? eult? slenderness ratio? alpha? r? gama? Fatigue1? Fatigue2? Degrade1?" << endln;
    return TCL_ERROR;
  }
  
  int tag;
  double fy, fu, Es, Esh, esh, eult;
  double slen = 0.0;
  double Fat1 = 0.114;
  double Fat2 = -4.46;
  double Deg1 = 0.00679;
  double alpha = 1.0;
  double r = 1.6;
  double gama = 0.4;
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial ReinforcingSteel tag" << endln;
    return TCL_ERROR;		
  }
  
  if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
    opserr << "WARNING invalid fy\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[4], &fu) != TCL_OK) {
    opserr << "WARNING invalid fu\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetDouble(interp, argv[5], &Es) != TCL_OK) {
    opserr << "WARNING invalid Es\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[6], &Esh) != TCL_OK) {
    opserr << "WARNING invalid Esh\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[7], &esh) != TCL_OK) {
    opserr << "WARNING invalid esh\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[8], &eult) != TCL_OK) {
    opserr << "WARNING invalid eult\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  if (argc > 9) {
    if (Tcl_GetDouble(interp, argv[9], &slen) != TCL_OK) {
      opserr << "WARNING invalid slenderness ratio\n";
      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
      return TCL_ERROR;	
    }
  }
  if (argc > 10) {
    if (Tcl_GetDouble(interp, argv[10], &alpha) != TCL_OK) {
      opserr << "WARNING invalid alpha factor\n";
      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[11], &r) != TCL_OK) {
      opserr << "WARNING invalid r factor\n";
      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[12], &gama) != TCL_OK) {
      opserr << "WARNING invalid gama factor\n";
      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
      return TCL_ERROR;	
    }
  }
  if (argc > 13) {
    if (Tcl_GetDouble(interp, argv[13], &Fat1) != TCL_OK) {
      opserr << "WARNING invalid Fatigue1\n";
      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
      return TCL_ERROR;	
    }
    if (Tcl_GetDouble(interp, argv[14], &Fat2) != TCL_OK) {
      opserr << "WARNING invalid Fatigue2\n";
      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
      return TCL_ERROR;	
    }
    if (Tcl_GetDouble(interp, argv[15], &Deg1) != TCL_OK) {
      opserr << "WARNING invalid Degrade1\n";
      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
      return TCL_ERROR;	
    }
  } 
  // Parsing was successful, allocate the material
  theMaterial = new ReinforcingSteel(tag, fy, fu, Es, Esh, esh, eult, slen, alpha, r, gama, Fat1, Fat2, Deg1);

  if (theMaterial != 0) 
    return theTclBuilder->addUniaxialMaterial(*theMaterial);
  else
    return -1;
}
