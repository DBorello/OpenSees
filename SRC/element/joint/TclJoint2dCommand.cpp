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
// Description: This file contains the implementation of the TclModelBuilder_addJoint2D()
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <SimpleJoint2D.h>
#include <Joint2D.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addJoint2D(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				char **argv, 
				Domain *theTclDomain,
				TclModelBuilder *theTclBuilder)
{
	// ensure the destructor has not been called
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
		cerr << "Want:\n";
		cerr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatC? RotSpringTag? LrgDsp?\n";
		cerr << "or:\n";
		cerr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatI? MatJ? MatK? MatL? MatC? LrgDsp?\n";
		return TCL_ERROR;
	}
	
	// get the id and end nodes
	int Joint2DId, iNode, jNode, kNode, lNode;
	if (Tcl_GetInt(interp, argv[argStart], &Joint2DId) != TCL_OK) {
		cerr << "WARNING invalid Joint2D eleTag" << endl;
		return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp, argv[1+argStart], &iNode) != TCL_OK) {
		cerr << "WARNING invalid iNode\n";
		cerr << "Joint2D element: " << Joint2DId << endl;
		return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp, argv[2+argStart], &jNode) != TCL_OK) {
		cerr << "WARNING invalid jNode\n";
		cerr << "Joint2D element: " << Joint2DId << endl;
		return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp, argv[3+argStart], &kNode) != TCL_OK) {
		cerr << "WARNING invalid kNode\n";
		cerr << "Joint2D element: " << Joint2DId << endl;
		return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp, argv[4+argStart], &lNode) != TCL_OK) {
		cerr << "WARNING invalid lNode\n";
		cerr << "Joint2D element: " << Joint2DId << endl;
		return TCL_ERROR;
	}

	// Get the center node
	int CenterNodeTag;
	if (Tcl_GetInt(interp, argv[5+argStart], &CenterNodeTag) != TCL_OK) {
		cerr << "WARNING invalid tag for center node\n";
		cerr << "Joint2D element: " << Joint2DId << endl;
		return TCL_ERROR;
	}
		
	// check domain for existence of internal node tag
	Node *CenterNode = theTclDomain->getNode(CenterNodeTag);
	if (CenterNode != 0) {
		cerr << "WARNING node tag specified for the center node already exists.\n";
		cerr << "Use a new node tag.\n";
		cerr << "Joint2D element: " << Joint2DId << endl;
		return TCL_ERROR;
	}


	// Decide to use SimpleJoint2D or Joint2D constructor, based on the number of arguments
	if ((argc-argStart) == 9 ) {
		// Using SimpleJoint2D constructor
		
		int PanelMatId;
		if (Tcl_GetInt(interp, argv[6+argStart], &PanelMatId) != TCL_OK) {
			cerr << "WARNING invalid matID\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}

		int RotEleTag;
		if (Tcl_GetInt(interp, argv[7+argStart], &RotEleTag) != TCL_OK) {
			cerr << "WARNING invalid tag for rotational spring.\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}

		// check domain for existence of spring element tag
		Element *RotEle = theTclDomain->getElement(RotEleTag);
		if (RotEle != 0) {
			cerr << "WARNING element tag specified for the rotational spring already exists.\n";
			cerr << "Use a new element tag.\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		
		int LargeDisp;
		if (Tcl_GetInt(interp, argv[8+argStart], &LargeDisp) != TCL_OK) {
			// use 0 as default
			LargeDisp = 0;
		}

		UniaxialMaterial *PanelMaterial = theTclBuilder->getUniaxialMaterial(PanelMatId);
		
		if ( PanelMaterial == 0 ) {
			cerr << "WARNING material not found\n";
			cerr << "Material: " << PanelMatId;
			cerr << "\nJoint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		
		SimpleJoint2D *theSimpleJoint2D;
		theSimpleJoint2D = new SimpleJoint2D( Joint2DId,
					iNode,jNode,kNode,lNode, 
					*PanelMaterial , 
					theTclDomain, 
					CenterNodeTag,
					RotEleTag,
					LargeDisp);
		
		if (theSimpleJoint2D == 0) {
			cerr << "WARNING ran out of memory creating element\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}

		if (theTclDomain->addElement(theSimpleJoint2D) == false) {
			cerr << "WARNING could not add element to the domain\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			delete theSimpleJoint2D;
			return TCL_ERROR;
		}

		// if get here we have sucessfully created the element and added it to the domain
		return TCL_OK;

		}
	else if ((argc-argStart) == 12 ){
		// Using Joint2D constructor
		int MatIid;
		if (Tcl_GetInt(interp, argv[6+argStart], &MatIid) != TCL_OK) {
			cerr << "WARNING invalid material ID for spring I\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		UniaxialMaterial *MatI = theTclBuilder->getUniaxialMaterial(MatIid);
		
		if ( MatI == 0 ) {
			cerr << "WARNING material not found\n";
			cerr << "Material: " << MatIid;
			cerr << "\nJoint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		
		int MatJid;
		if (Tcl_GetInt(interp, argv[7+argStart], &MatJid) != TCL_OK) {
			cerr << "WARNING invalid material ID for spring J\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		UniaxialMaterial *MatJ = theTclBuilder->getUniaxialMaterial(MatJid);
		
		if ( MatJ == 0 ) {
			cerr << "WARNING material not found\n";
			cerr << "Material: " << MatJid;
			cerr << "\nJoint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}

		int MatKid;
		if (Tcl_GetInt(interp, argv[8+argStart], &MatKid) != TCL_OK) {
			cerr << "WARNING invalid material ID for spring K\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		UniaxialMaterial *MatK = theTclBuilder->getUniaxialMaterial(MatKid);
		
		if ( MatK == 0 ) {
			cerr << "WARNING material not found\n";
			cerr << "Material: " << MatKid;
			cerr << "\nJoint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}

		int MatLid;
		if (Tcl_GetInt(interp, argv[9+argStart], &MatLid) != TCL_OK) {
			cerr << "WARNING invalid material ID for spring L\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		UniaxialMaterial *MatL = theTclBuilder->getUniaxialMaterial(MatLid);
		
		if ( MatL == 0 ) {
			cerr << "WARNING material not found\n";
			cerr << "Material: " << MatLid;
			cerr << "\nJoint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		
		int PanelMatId;
		if (Tcl_GetInt(interp, argv[10+argStart], &PanelMatId) != TCL_OK) {
			cerr << "WARNING invalid matID\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}
		UniaxialMaterial *PanelMaterial = theTclBuilder->getUniaxialMaterial(PanelMatId);
		
		if ( PanelMaterial == 0 ) {
			cerr << "WARNING material not found\n";
			cerr << "Material: " << PanelMatId;
			cerr << "\nJoint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}

		int LargeDisp;
		if (Tcl_GetInt(interp, argv[11+argStart], &LargeDisp) != TCL_OK) {
			// use 0 as default
			LargeDisp = 0;
		}

		Joint2D *theJoint2D;
		theJoint2D = new Joint2D( Joint2DId,
					iNode,jNode,kNode,lNode,CenterNodeTag,
					*MatI,*MatJ,*MatK,*MatL,*PanelMaterial, 
					theTclDomain, 
					LargeDisp);
		
		if (theJoint2D == 0) {
			cerr << "WARNING ran out of memory creating element\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			return TCL_ERROR;
		}

		if (theTclDomain->addElement(theJoint2D) == false) {
			cerr << "WARNING could not add element to the domain\n";
			cerr << "Joint2D element: " << Joint2DId << endl;
			delete theJoint2D;
			return TCL_ERROR;
		}

		// if get here we have sucessfully created the element and added it to the domain
		return TCL_OK;
		}
	else {
		return TCL_ERROR;
	}
}



