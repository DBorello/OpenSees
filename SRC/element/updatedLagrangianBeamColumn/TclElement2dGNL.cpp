#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <Elastic2DGNL.h>

#include <TclModelBuilder.h>

#define  tcl_debug 1

// Elastic2DGNL(int tag, double A, double E, double I, int Nd1, int Nd2,
//             double rho = 0.0, bool islinear = false);

int
TclModelBuilder_addElastic2dGNL (ClientData clientData, Tcl_Interp *interp,
				 int argc, char **argv,
				 Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

    if(tcl_debug)
        cout << " TclModelBuilder_addElastic2dGNL \n";

	if (argc < 8)
	{
		cerr << "WARNING insufficient arguments\n";
		cerr << "element element2dGNL int tag, int Nd1, int Nd2, double A, double E, double Iz, <int linear>\n";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
	double massDens = 0.0;
	bool   linear = false;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		cerr << "WARNING invalid Elastic2dGNL tag" << endl;
		return TCL_ERROR;
	}
    if(tcl_debug) cout << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		cerr << "WARNING invalid node I\n";
		cerr << "Elastic2dGNL: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		cerr << "WARNING invalid node J\n";
		cerr << "Elastic2dGNL: " << tag << endl;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		cerr << "WARNING invalid A\n";
		cerr << "Elastic2dGNL: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		cerr << "WARNING invalid E\n";
		cerr << "Elastic2dGNL: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		cerr << "WARNING invalid I\n";
		cerr << "Elastic2dGNL: " << tag << endl;
		return TCL_ERROR;
	}
	
	if(argc == 9)
	{
		int lin = 0;
		if(Tcl_GetInt(interp, argv[8], &lin) != TCL_OK)
		{
			cerr << "WARNING invalid Linear Flag\n";
			cerr << "Elastic2dGNL: " << tag << endl;
			return TCL_ERROR;
		}
		
		if(lin == 1)
			linear = true;
		
		if(tcl_debug)
			cout << " 9 arguments - " << lin << endl;
	}
	
	// if(tcl_debug) cout << "\tAdded upto mass - input parameters\n";


	Element *theElement = new Elastic2dGNL(tag, A, E, I, ndI, ndJ, linear);//, false, massDens);

	if(tcl_debug) cout << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		cerr << "WARNING ran out of memory creating element\n";
		cerr << "Elastic2dGNL: " << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false)
	{
		cerr << "WARNING TclElmtBuilder - addElastic2dGNL - could not add element to domain ";
		cerr << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if(tcl_debug) cout << "\tElement number " << tag << " added to domain - returning\n";

	return TCL_OK;
}
