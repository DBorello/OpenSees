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
// $Date: 2001-08-18 21:39:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/drain/TclDrainMaterialCommand.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the implementation of the
// TclModelBuilder_addDrainMaterial() function. 

#include <DrainHardeningMaterial.h>
#include <DrainBilinearMaterial.h>
#include <DrainClough1Material.h>
#include <DrainClough2Material.h>
#include <DrainPinch1Material.h>

#include <TclModelBuilder.h>
#include <Vector.h>
#include <string.h>

static void printCommand(int argc, char **argv)
{
    cerr << "Input command: ";
    for (int i=0; i<argc; i++)
	cerr << argv[i] << " ";
    cerr << endl;
} 

UniaxialMaterial *
TclModelBuilder_addDrainMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				char **argv, TclModelBuilder *theTclBuilder)
{
	if (argc < 3) {
		cerr << "WARNING insufficient number of arguments\n";
		printCommand(argc, argv);
		return 0;
	}

	int tag;
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial tag\n";
		printCommand(argc, argv);
	    return 0;
	}

	UniaxialMaterial *theMaterial = 0;

	if (strcmp(argv[1],"Hardening2") == 0 || strcmp(argv[1],"Hardening02") == 0) {
		if (argc < 7) {
			cerr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Hardening02 tag? E? sigY? Hiso? Hkin?" << endl;
			return 0;
		}
		
		double E, sigY, Hiso, Hkin;

		if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
			cerr << "WARNING invalid E\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &sigY) != TCL_OK) {
			cerr << "WARNING invalid sigY\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &Hiso) != TCL_OK) {
			cerr << "WARNING invalid Hiso\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &Hkin) != TCL_OK) {
			cerr << "WARNING invalid Hkin\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new DrainHardeningMaterial(tag, E, sigY, Hiso, Hkin);
	}
       
	else if (strcmp(argv[1],"BiLinear") == 0) {
		if (argc < 19) {
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial BiLinear tag? ..." << endl;
			return 0;
		}
		
		Vector input(16);
		double temp;

		for (int i = 3, j = 0; j < 16; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				cerr << "WARNING invalid input, data " << i << '\n';
				printCommand(argc, argv);
				return 0;
			}
			input(j) = temp;
		}

		theMaterial = new DrainBilinearMaterial(tag, input);
	}

	else if (strcmp(argv[1],"Clough1") == 0) {
		if (argc < 19) {
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Clough1 tag? ..." << endl;
			return 0;
		}
		
		Vector input(16);
		double temp;

		for (int i = 3, j = 0; j < 16; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				cerr << "WARNING invalid input, data " << i << '\n';
				printCommand(argc, argv);
				return 0;
			}
			input(j) = temp;
		}

		theMaterial = new DrainClough1Material(tag, input);
	}

	else if (strcmp(argv[1],"Clough2") == 0) {
		if (argc < 19) {
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Clough2 tag? ..." << endl;
			return 0;
		}
		
		Vector input(16);
		double temp;

		for (int i = 3, j = 0; j < 16; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				cerr << "WARNING invalid input, data " << i << '\n';
				printCommand(argc, argv);
				return 0;
			}
			input(j) = temp;
		}

		theMaterial = new DrainClough2Material(tag, input);
	}

	else if (strcmp(argv[1],"Pinch1") == 0) {
		if (argc < 22) {
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Pinch1 tag? ..." << endl;
			return 0;
		}
		
		Vector input(19);
		double temp;

		for (int i = 3, j = 0; j < 19; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				cerr << "WARNING invalid input, data " << i << '\n';
				printCommand(argc, argv);
				return 0;
			}
			input(j) = temp;
		}

		theMaterial = new DrainPinch1Material(tag, input);
	}

	return theMaterial;
}
