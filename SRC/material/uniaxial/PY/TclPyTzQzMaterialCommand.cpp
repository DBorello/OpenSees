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
// $Date: 2002-04-29 00:04:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PY/TclPyTzQzMaterialCommand.cpp,v $

//temp out BJ #include <Domain.h>     // RWB bringing in Domain for PyLiq
//temp out BJ // !!!!!!!!!!!!! IFF to be discussed and changed, this Domain, ROss needs
//temp out BJ // that to make PYliq working...

#include <TclModelBuilder.h>

//PY Springs: RWBoulanger and BJeremic
#include <PySimple1.h>// RWB
#include <TzSimple1.h>// RWB
#include <QzSimple1.h>// RWB
// out for the moment#include <PyLiq1.h>// RWB
// out for the moment#include <TzLiq1.h>// RWB

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
TclModelBuilder_addPyTzQzMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				char **argv, TclModelBuilder *theTclBuilder)
//temp out BJ					char **argv, Domain *theDomain, TclModelBuilder *theTclBuilder)
//temp out BJ // !!!!!!!!!!!!! IFF to be discussed and changed, this theDomain argument, ROss needs
//temp out BJ // that to make PYliq working...
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

//  INSERTING THE EXTRA LINES FOR PySimple1 //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
	if (strcmp(argv[1],"PySimple1") == 0) {
	if (argc < 7) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial PySimple1 tag? soilType? pult? y50? drag? dashpot? " << endl;
	    return 0;
	}

	int tag, soilType;
	double pult, y50, drag, dashpot;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial PySimple1 tag" << endl;
	    return 0;		
	}

	if (Tcl_GetInt(interp, argv[3], &soilType) != TCL_OK) {
	    cerr << "WARNING invalid soilType\n";
	    cerr << "uniaxialMaterial PySimple1: " << tag << endl;
	    return 0;	
	}

	if (Tcl_GetDouble(interp, argv[4], &pult) != TCL_OK) {
	    cerr << "WARNING invalid pult\n";
	    cerr << "uniaxialMaterial PySimple1: " << tag << endl;
	    return 0;
	}

	if (Tcl_GetDouble(interp, argv[5], &y50) != TCL_OK) {
	    cerr << "WARNING invalid y50\n";
	    cerr << "uniaxialMaterial PySimple1: " << tag << endl;
	    return 0;	
	}

	if (Tcl_GetDouble(interp, argv[6], &drag) != TCL_OK) {
	    cerr << "WARNING invalid drag\n";
	    cerr << "uniaxialMaterial PySimple1: " << tag << endl;
	    return 0;	
	}

	if (argc == 7) dashpot = 0.0;

	if (argc > 7) {
		if (Tcl_GetDouble(interp, argv[7], &dashpot) != TCL_OK) {
			cerr << "WARNING invalid dashpot\n";
			cerr << "uniaxialMaterial PySimple1: " << tag << endl;
			return 0;	
		}
	}

	// Parsing was successful, allocate the material
	theMaterial = new PySimple1(tag,MAT_TAG_PySimple1,soilType, pult, y50, drag, dashpot);
    }

//  INSERTING THE EXTRA LINES FOR PyLiq1 //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
//temp out BJ 	else if (strcmp(argv[1],"PyLiq1") == 0) {
//temp out BJ 	if (argc < 10) {
//temp out BJ 	    cerr << "WARNING insufficient arguments\n";
//temp out BJ 	    printCommand(argc,argv);
//temp out BJ 	    cerr << "Want: uniaxialMaterial PyLiq1 tag? soilType? pult? y50? drag? dashpot? solidElem1? solidElem2?" << endl;
//temp out BJ 	    return 0;
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	int tag, soilType, solidElem1, solidElem2;
//temp out BJ 	double pult, y50, drag, dashpot;
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid uniaxialMaterial PyLiq1 tag" << endl;
//temp out BJ 	    return 0;		
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[3], &soilType) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid soilType\n";
//temp out BJ 	    cerr << "uniaxialMaterial PyLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[4], &pult) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid pult\n";
//temp out BJ 	    cerr << "uniaxialMaterial PyLiq1: " << tag << endl;
//temp out BJ 	    return 0;
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[5], &y50) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid y50\n";
//temp out BJ 	    cerr << "uniaxialMaterial PyLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[6], &drag) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid drag\n";
//temp out BJ 	    cerr << "uniaxialMaterial PyLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[7], &dashpot) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid dashpot\n";
//temp out BJ 	    cerr << "uniaxialMaterial PyLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[8], &solidElem1) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid solidElem\n";
//temp out BJ 	    cerr << "uniaxialMaterial PyLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 	
//temp out BJ 	if (Tcl_GetInt(interp, argv[9], &solidElem2) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid solidElem\n";
//temp out BJ 	    cerr << "uniaxialMaterial PyLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	// Parsing was successful, allocate the material
//temp out BJ 	theMaterial = new PyLiq1(tag, MAT_TAG_PyLiq1,soilType, pult, y50, drag, dashpot,
//temp out BJ 							solidElem1, solidElem2, theDomain);
//temp out BJ     }

//  INSERTING THE EXTRA LINES FOR QzSimple1 //////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

	else if (strcmp(argv[1],"QzSimple1") == 0) {
	if (argc < 6) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial QzSimple1 tag? QzType? Qult? z50? suction? dashpot? " << endl;
	    return 0;
	}

	int tag, QzType;
	double Qult, z50, suction, dashpot;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial QzSimple1 tag" << endl;
	    return 0;		
	}

	if (Tcl_GetInt(interp, argv[3], &QzType) != TCL_OK) {
	    cerr << "WARNING invalid QzType\n";
	    cerr << "uniaxialMaterial QzSimple1: " << tag << endl;
	    return 0;	
	}

	if (Tcl_GetDouble(interp, argv[4], &Qult) != TCL_OK) {
	    cerr << "WARNING invalid Qult\n";
	    cerr << "uniaxialMaterial QzSimple1: " << tag << endl;
	    return 0;
	}

	if (Tcl_GetDouble(interp, argv[5], &z50) != TCL_OK) {
	    cerr << "WARNING invalid z50\n";
	    cerr << "uniaxialMaterial QzSimple1: " << tag << endl;
	    return 0;	
	}

	if (argc == 6) {
		suction = 0.0;
		dashpot = 0.0;
	}

	if (argc > 6) {
		if (Tcl_GetDouble(interp, argv[6], &suction) != TCL_OK) {
		    cerr << "WARNING invalid suction\n";
			cerr << "uniaxialMaterial QzSimple1: " << tag << endl;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[7], &dashpot) != TCL_OK) {
			cerr << "WARNING invalid dashpot\n";
			cerr << "uniaxialMaterial QzSimple1: " << tag << endl;
			return 0;	
		}
	}

	// Parsing was successful, allocate the material
	theMaterial = new QzSimple1(tag, QzType, Qult, z50, suction, dashpot);
    }


//  INSERTING THE EXTRA LINES FOR TzSimple1 //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
	else if (strcmp(argv[1],"TzSimple1") == 0) {
	if (argc < 6) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial TzSimple1 tag? tzType? tult? z50? dashpot? " << endl;
	    return 0;
	}

	int tag, tzType;
	double tult, z50, dashpot;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial TzSimple1 tag" << endl;
	    return 0;		
	}

	if (Tcl_GetInt(interp, argv[3], &tzType) != TCL_OK) {
	    cerr << "WARNING invalid tzType\n";
	    cerr << "uniaxialMaterial TzSimple1: " << tag << endl;
	    return 0;	
	}

	if (Tcl_GetDouble(interp, argv[4], &tult) != TCL_OK) {
	    cerr << "WARNING invalid tult\n";
	    cerr << "uniaxialMaterial TzSimple1: " << tag << endl;
	    return 0;
	}

	if (Tcl_GetDouble(interp, argv[5], &z50) != TCL_OK) {
	    cerr << "WARNING invalid z50\n";
	    cerr << "uniaxialMaterial TzSimple1: " << tag << endl;
	    return 0;	
	}

	if (argc == 6) dashpot = 0.0;

	if (argc > 6) {
		if (Tcl_GetDouble(interp, argv[6], &dashpot) != TCL_OK) {
			cerr << "WARNING invalid dashpot\n";
			cerr << "uniaxialMaterial TzSimple1: " << tag << endl;
			return 0;	
		}
	}

	// Parsing was successful, allocate the material
	theMaterial = new TzSimple1(tag, MAT_TAG_TzSimple1, tzType, tult, z50, dashpot);
    }

	//  INSERTING THE EXTRA LINES FOR TzLiq1 //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
//temp out BJ 	else if (strcmp(argv[1],"TzLiq1") == 0) {
//temp out BJ 	if (argc < 9) {
//temp out BJ 	    cerr << "WARNING insufficient arguments\n";
//temp out BJ 	    printCommand(argc,argv);
//temp out BJ 	    cerr << "Want: uniaxialMaterial TzLiq1 tag? tzType? tult? z50? dashpot? solidElem1? solidElem2?" << endl;
//temp out BJ 	    return 0;
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	int tag, tzType, solidElem1, solidElem2;
//temp out BJ 	double tult, z50, dashpot;
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid uniaxialMaterial TzLiq1 tag" << endl;
//temp out BJ 	    return 0;		
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[3], &tzType) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid tzType\n";
//temp out BJ 	    cerr << "uniaxialMaterial TzLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[4], &tult) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid tult\n";
//temp out BJ 	    cerr << "uniaxialMaterial TzLiq1: " << tag << endl;
//temp out BJ 	    return 0;
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[5], &z50) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid z50\n";
//temp out BJ 	    cerr << "uniaxialMaterial TzLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[6], &dashpot) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid dashpot\n";
//temp out BJ 	    cerr << "uniaxialMaterial TzLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[7], &solidElem1) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid solidElem\n";
//temp out BJ 	    cerr << "uniaxialMaterial TzLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 	
//temp out BJ 	if (Tcl_GetInt(interp, argv[8], &solidElem2) != TCL_OK) {
//temp out BJ 	    cerr << "WARNING invalid solidElem\n";
//temp out BJ 	    cerr << "uniaxialMaterial TzLiq1: " << tag << endl;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	// Parsing was successful, allocate the material
//temp out BJ 	theMaterial = new TzLiq1(tag, MAT_TAG_TzLiq1,tzType, tult, z50, dashpot,
//temp out BJ 							solidElem1, solidElem2, theDomain);
//temp out BJ     }
//temp out BJ 
//////////////////////////////////////////////////////////////////////

	return theMaterial;
}
