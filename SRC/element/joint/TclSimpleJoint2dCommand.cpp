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
                                                                        
// Written: Arash	Created: 04/01
// Revision: 
//				AAA		05/01
//
// Description: This file contains the implementation of the TclModelBuilder_addSimpleJoint2D()
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <SimpleJoint2D.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addSimpleJoint2D(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				char **argv, 
				Domain *theTclDomain,
				TclModelBuilder *theTclBuilder)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 3) {
		cerr << "WARNING -- model dimensions and/or nodal DOF not compatible with BCconnect element\n";
		return TCL_ERROR;
	}

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc-argStart) < 9) {
    cerr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    cerr << "Want: element SimpleJoint2D eleTag? iNode? jNode? kNode? lNode? matTag? jntNodeTag? rotEleTag?\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int SimpleJoint2DId, iNode, jNode, kNode, lNode;
  SimpleJoint2D *theSimpleJoint2D;

  if (Tcl_GetInt(interp, argv[argStart], &SimpleJoint2DId) != TCL_OK) {
    cerr << "WARNING invalid SimpleJoint2D eleTag" << endl;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[1+argStart], &iNode) != TCL_OK) {
    cerr << "WARNING invalid iNode\n";
    cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2+argStart], &jNode) != TCL_OK) {
    cerr << "WARNING invalid jNode\n";
    cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[3+argStart], &kNode) != TCL_OK) {
    cerr << "WARNING invalid kNode\n";
    cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }  
  
  if (Tcl_GetInt(interp, argv[4+argStart], &lNode) != TCL_OK) {
    cerr << "WARNING invalid lNode\n";
    cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }  
  
  int PanelMatId;

  if (Tcl_GetInt(interp, argv[5+argStart], &PanelMatId) != TCL_OK) {
    cerr << "WARNING invalid matID\n";
    cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }
  
  int LargeDisp;

  if (Tcl_GetInt(interp, argv[6+argStart], &LargeDisp) != TCL_OK) {
    LargeDisp = 0;
  }

  int newNodeTag, newRotEleTag;
  if (Tcl_GetInt(interp, argv[7+argStart], &newNodeTag) != TCL_OK) {
    cerr << "WARNING invalid NodeTag \n";
    cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[8+argStart], &newRotEleTag) != TCL_OK) {
    cerr << "WARNING invalid RotEleTag\n";
    cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }

  UniaxialMaterial *PanelMaterial = theTclBuilder->getUniaxialMaterial(PanelMatId);
  
  if (PanelMaterial == 0) {
    cerr << "WARNING material not found\n";
    cerr << "Material: " << PanelMatId;
    cerr << "\nSimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }

  theSimpleJoint2D = new SimpleJoint2D( SimpleJoint2DId,
					iNode,jNode,kNode,lNode, 
					*PanelMaterial , 
					theTclDomain, 
					newNodeTag,
					newRotEleTag,
					LargeDisp);

  if (theSimpleJoint2D == 0) {
    cerr << "WARNING ran out of memory creating element\n";
    cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theSimpleJoint2D) == false) {
      cerr << "WARNING could not add element to the domain\n";
      cerr << "SimpleJoint2D element: " << SimpleJoint2DId << endl;
      delete theSimpleJoint2D;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}



