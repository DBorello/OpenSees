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
// $Date: 2001-07-13 21:11:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/TclDispBeamColumnCommand.cpp,v $
                                                                        
// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the implementation of the 
// TclModelBuilder_addDispBeamColumn() command. 

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <DispBeamColumn2d.h>
#include <DispBeamColumn3d.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addDispBeamColumn(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				char **argv, 
				Domain*theTclDomain,
				TclModelBuilder *theTclBuilder)
{
	// ensure the destructor has not been called - 
	if (theTclBuilder == 0) {
		cerr << "WARNING builder has been destroyed\n";    
		return TCL_ERROR;
	}

	int ndm = theTclBuilder->getNDM();
	int ndf = theTclBuilder->getNDF();

	int ok = 0;
	if (ndm == 2 && ndf == 3)
		ok = 1;
	if (ndm == 3 && ndf == 6)
		ok = 1;

	if (ok == 0) {
		cerr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
			<< " not compatible with dispBeamColumn element" << endl;
		return TCL_ERROR;
	}

	if (argc < 8) {
		cerr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		cerr << "Want: element dispBeamColumn eleTag? iNode? jNode? numSec? secTag? transfTag?\n";
		return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, numSec, secTag, transfTag;

	if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
		cerr << "WARNING invalid dispBeamColumn eleTag" << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
		cerr << "WARNING invalid iNode\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
		cerr << "WARNING invalid jNode\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}
  
	if (Tcl_GetInt(interp, argv[5], &numSec) != TCL_OK) {
		cerr << "WARNING invalid numSec\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}  
  
	if (Tcl_GetInt(interp, argv[6], &secTag) != TCL_OK) {
		cerr << "WARNING invalid secTag\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}  

	if (Tcl_GetInt(interp, argv[7], &transfTag) != TCL_OK) {
		cerr << "WARNING invalid transfTag\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}

	SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
      
	if (theSection == 0) {
		cerr << "WARNING section not found\n";
		cerr << "Section: " << secTag;
		cerr << "\ndispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}

	Element *theElement = 0;

	if (ndm == 2) {

		CrdTransf2d *theTransf = theTclBuilder->getCrdTransf2d(transfTag);
      
		if (theTransf == 0) {
			cerr << "WARNING transformation not found\n";
			cerr << "transformation: " << transfTag;
			cerr << "\ndispBeamColumn element: " << eleTag << endl;
			return TCL_ERROR;
		}

		// now create the DispBeamColumn and add it to the Domain
		theElement = new DispBeamColumn2d(eleTag,iNode,jNode,numSec,*theSection,*theTransf);
	}

	if (ndm == 3) {

		CrdTransf3d *theTransf = theTclBuilder->getCrdTransf3d(transfTag);
      
		if (theTransf == 0) {
			cerr << "WARNING transformation not found\n";
			cerr << "transformation: " << transfTag;
			cerr << "\ndispBeamColumn element: " << eleTag << endl;
			return TCL_ERROR;
		}

		// now create the DispBeamColumn and add it to the Domain
		theElement = new DispBeamColumn3d(eleTag,iNode,jNode,numSec,*theSection,*theTransf);
	}

	if (theElement == 0) {
		cerr << "WARNING ran out of memory creating element\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
		cerr << "WARNING could not add element to the domain\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		delete theElement;
		return TCL_ERROR;
	}

	// if get here we have sucessfully created the element and added it to the domain
	return TCL_OK;
}



