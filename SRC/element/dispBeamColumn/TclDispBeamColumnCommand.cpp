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
// $Date: 2002-12-10 03:04:26 $
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
		cerr << "Want: element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
		return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, transfTag;
	int secTag[10]; // Max size of integration rule ... can change if needed
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
		cerr << "WARNING invalid dispBeamColumn eleTag" << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
		cerr << "WARNING invalid iNode\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
		cerr << "WARNING invalid jNode\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}
  
	if (Tcl_GetInt(interp, argv[argi++], &nIP) != TCL_OK) {
		cerr << "WARNING invalid nIP\n";
		cerr << "dispBeamColumn element: " << eleTag << endl;
		return TCL_ERROR;
	}  
  
	if (strcmp(argv[argi], "-sections") == 0) {
	  argi++;
	  if (argi+nIP > argc) {
	    interp->result = "WARNING insufficient number of section tags - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?";
	    return TCL_ERROR;
	  }
	  int section;
	  for (int i = 0; i < nIP; i++) {
	    if (Tcl_GetInt(interp, argv[argi+i], &section) != TCL_OK) {
	      interp->result = "WARNING invalid secTag - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?";
	      return TCL_ERROR;
	    }
	    secTag[i] = section;
	  }
	  argi += nIP;
	}
	
	else {
	  int section;
	  if (Tcl_GetInt(interp, argv[argi++], &section) != TCL_OK) {
	    interp->result = "WARNING invalid secTag - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?";
	    return TCL_ERROR;
	  }
	  for (int i = 0; i < nIP; i++)
	    secTag[i] = section;
	}
	
	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
	  interp->result = "WARNING invalid transfTag? - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?";
	  return TCL_ERROR;
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [nIP];
	
	if (!sections) {
	  interp->result = "WARNING TclElmtBuilder - addFrameElement - Insufficient memory to create sections";
	  return TCL_ERROR;
	}
	
	for (int j=0; j<nIP; j++) {
	  SectionForceDeformation *theSection = theTclBuilder->getSection(secTag[j]);
	  
	  if (theSection == 0) {
	    cerr << "WARNING TclElmtBuilder - frameElement - no Section found with tag ";
	    cerr << secTag[j] << endl;
	    delete [] sections;
	    return TCL_ERROR;
	  }

	  sections[j] = theSection;
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
		theElement = new DispBeamColumn2d(eleTag,iNode,jNode,nIP,sections,*theTransf);

		delete [] sections;
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
		theElement = new DispBeamColumn3d(eleTag,iNode,jNode,nIP,sections,*theTransf);

		delete [] sections;
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



