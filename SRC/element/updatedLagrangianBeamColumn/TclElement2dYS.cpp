#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <Inelastic2DYS01.h>
#include <Inelastic2DYS02.h>
#include <Inelastic2DYS03.h>
//#include <Inelastic2DYS04.h>
//#include <Inelastic2DYS05.h>

#include <YieldSurface_BC.h>
#include <TclModelBuilder.h>

#define  tcl_debug 0

// Element2dGNL(int tag, double A, double E, double I, int Nd1, int Nd2,
//             double rho = 0.0, bool islinear = false);

int
TclModelBuilder_addElement2dYS01 (ClientData clientData, Tcl_Interp *interp,
								   int argc, char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

    if(tcl_debug)
        cout << " TclModelBuilder_addElement2dGNL \n";

	if (argc < 11)
	{
		cerr << "WARNING insufficient arguments\n";
		cerr << "element element2dYS tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
//	double massDens = 0.0;
	int ysID1, ysID2;
	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		cerr << "WARNING invalid element2dYS tag" << endl;
		return TCL_ERROR;
	}
    if(tcl_debug) cout << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		cerr << "WARNING invalid node I\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		cerr << "WARNING invalid node J\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		cerr << "WARNING invalid A\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		cerr << "WARNING invalid E\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		cerr << "WARNING invalid I\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[8], &ysID1) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[9], &ysID2) != TCL_OK)
	{
		cerr << "WARNING invalid ysID2\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[10], &rf_algo) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		cerr << "WARNING element2dYS: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID1 << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		cerr << "WARNING element2dYS: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID2 << endl;
		return TCL_ERROR;
	}

// 		Inelastic2DYS(	int tag, double A, double E, double I, int Nd1, int Nd2,
// 						YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
// 						int rf_algo = -1, // updated
	Element *theElement = new Inelastic2DYS01(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

	if(tcl_debug) cout << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		cerr << "WARNING ran out of memory creating element\n";
		cerr << "element2dYS: " << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false)
	{
		cerr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		cerr << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if(tcl_debug) cout << "\tElement number " << tag << " added to domain - returning\n";

	return TCL_OK;
}


int
TclModelBuilder_addElement2dYS02 (ClientData clientData, Tcl_Interp *interp,
								   int argc, char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

    if(tcl_debug)
        cout << " TclModelBuilder_addElement2dGNL \n";

	if (argc < 12)
	{
		cerr << "WARNING insufficient arguments\n";
		cerr << "element element2dYS tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? power? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
//	double massDens = 0.0;
	int ysID1, ysID2;

	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		cerr << "WARNING invalid element2dYS tag" << endl;
		return TCL_ERROR;
	}
    if(tcl_debug) cout << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		cerr << "WARNING invalid node I\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		cerr << "WARNING invalid node J\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		cerr << "WARNING invalid A\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		cerr << "WARNING invalid E\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		cerr << "WARNING invalid I\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[8], &ysID1) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[9], &ysID2) != TCL_OK)
	{
		cerr << "WARNING invalid ysID2\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	double pow_min, pow_Kun;
	if (Tcl_GetDouble (interp, argv[10], &pow_min) != TCL_OK)
	{
		cerr << "WARNING invalid power\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}
	
	if (Tcl_GetDouble (interp, argv[11], &pow_Kun) != TCL_OK)
	{
		cerr << "WARNING invalid power\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[12], &rf_algo) != TCL_OK)
	{
		cerr << "WARNING invalid rfalgo\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		cerr << "WARNING element2dYS: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID1 << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		cerr << "WARNING element2dYS: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID2 << endl;
		return TCL_ERROR;
	}

//Inelastic2DYS02(int tag, double a, double e, double i, int Nd1, int Nd2,
//				YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
//				int rf_algo, bool islinear, double rho)


Element *theElement = new Inelastic2DYS02(tag, A, E, I, ndI, ndJ, theYS1, theYS2, pow_min, pow_Kun, rf_algo);

    cout << "Inelastic2DYS02 created\n";

	if(tcl_debug) cout << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		cerr << "WARNING ran out of memory creating element\n";
		cerr << "element2dYS: " << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	cout << "Inelastic2DYS02 adding to domain\n";

	if (theDomain->addElement(theElement) == false)
	{
		cerr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		cerr << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	//if(tcl_debug)
		cout << "Inelastic2DYS02 #" << tag << " added to domain - returning\n";

	return TCL_OK;
}


int
TclModelBuilder_addElement2dYS03 (ClientData clientData, Tcl_Interp *interp,
								   int argc, char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

    if(tcl_debug)
        cout << " TclModelBuilder_addElement2dGNL \n";

	if (argc < 11)
	{
		cerr << "WARNING insufficient arguments\n";
		cerr << "element element2dYS03 tag? Nd1? Nd2? A_ten? A_com? E? IzPos? IzNeg? ysID1? ysID2? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, aTens, aComp, Ipos, Ineg;
//	double massDens = 0.0;
	int ysID1, ysID2;

	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		cerr << "WARNING invalid element2dYS tag" << endl;
		return TCL_ERROR;
	}
    if(tcl_debug) cout << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		cerr << "WARNING invalid node I\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		cerr << "WARNING invalid node J\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &aTens) != TCL_OK)
	{
		cerr << "WARNING invalid A\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &aComp) != TCL_OK)
	{
		cerr << "WARNING invalid A\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &E) != TCL_OK)
	{
		cerr << "WARNING invalid E\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[8], &Ipos) != TCL_OK)
	{
		cerr << "WARNING invalid I\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[9], &Ineg) != TCL_OK)
	{
		cerr << "WARNING invalid I\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[10], &ysID1) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[11], &ysID2) != TCL_OK)
	{
		cerr << "WARNING invalid ysID2\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[12], &rf_algo) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS: " << tag << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		cerr << "WARNING element2dYS: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID1 << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		cerr << "WARNING element2dYS: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID2 << endl;
		return TCL_ERROR;
	}

//	Inelastic2DYS03(int tag, double a_ten, double a_com, double e,
//	                double iz_pos, double iz_neg, int Nd1, int Nd2,
//                    YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
//                    int rf_algo, bool islinear, double rho);

Element *theElement = new Inelastic2DYS03(tag, aTens, aComp, E,
                                          Ipos, Ineg, ndI, ndJ,
                                          theYS1, theYS2, rf_algo);

    cout << "Inelastic2DYS03 created\n";

	if(tcl_debug) cout << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		cerr << "WARNING ran out of memory creating element\n";
		cerr << "element2dYS: " << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	cout << "Inelastic2DYS03 adding to domain\n";

	if (theDomain->addElement(theElement) == false)
	{
		cerr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		cerr << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if(tcl_debug)
		cout << "Inelastic2DYS03 #" << tag << " added to domain - returning\n";

	return TCL_OK;
}


/*
int
TclModelBuilder_addElement2dYS04 (ClientData clientData, Tcl_Interp *interp,
								   int argc, char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

	if (argc < 11)
	{
		cerr << "WARNING insufficient arguments\n";
		cerr << "element element2dYS04 tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
//	double massDens = 0.0;
	int ysID1, ysID2;
	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		cerr << "WARNING invalid element2dYS04 tag" << endl;
		return TCL_ERROR;
	}
    if(tcl_debug) cout << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		cerr << "WARNING invalid node I\n";
		cerr << "element2dYS04: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		cerr << "WARNING invalid node J\n";
		cerr << "element2dYS04: " << tag << endl;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		cerr << "WARNING invalid A\n";
		cerr << "element2dYS04: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		cerr << "WARNING invalid E\n";
		cerr << "element2dYS04: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		cerr << "WARNING invalid I\n";
		cerr << "element2dYS04: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[8], &ysID1) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS04: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[9], &ysID2) != TCL_OK)
	{
		cerr << "WARNING invalid ysID2\n";
		cerr << "element2dYS04: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[10], &rf_algo) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS04: " << tag << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		cerr << "WARNING element2dYS04: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID1 << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		cerr << "WARNING element2dYS04: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID2 << endl;
		return TCL_ERROR;
	}

// 		Inelastic2DYS(	int tag, double A, double E, double I, int Nd1, int Nd2,
// 						YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
// 						int rf_algo = -1, // updated

	Element *theElement = new Inelastic2DYS04(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

	if(tcl_debug) cout << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		cerr << "WARNING ran out of memory creating element\n";
		cerr << "element2dYS04: " << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false)
	{
		cerr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		cerr << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if(tcl_debug) cout << "\tElement number " << tag << " added to domain - returning\n";

	return TCL_OK;
}


int
TclModelBuilder_addElement2dYS05 (ClientData clientData, Tcl_Interp *interp,
								   int argc, char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

	if (argc < 11)
	{
		cerr << "WARNING insufficient arguments\n";
		cerr << "element element2dYS04 tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
//	double massDens = 0.0;
	int ysID1, ysID2;
	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		cerr << "WARNING invalid element2dYS05 tag" << endl;
		return TCL_ERROR;
	}
    if(tcl_debug) cout << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		cerr << "WARNING invalid node I\n";
		cerr << "element2dYS05: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		cerr << "WARNING invalid node J\n";
		cerr << "element2dYS05: " << tag << endl;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		cerr << "WARNING invalid A\n";
		cerr << "element2dYS05: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		cerr << "WARNING invalid E\n";
		cerr << "element2dYS05: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		cerr << "WARNING invalid I\n";
		cerr << "element2dYS05: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[8], &ysID1) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS05: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[9], &ysID2) != TCL_OK)
	{
		cerr << "WARNING invalid ysID2\n";
		cerr << "element2dYS05: " << tag << endl;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[10], &rf_algo) != TCL_OK)
	{
		cerr << "WARNING invalid ysID1\n";
		cerr << "element2dYS05: " << tag << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		cerr << "WARNING element2dYS05: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID1 << endl;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		cerr << "WARNING element2dYS05: " << tag << "\n";
		cerr <<  " no yield surface exists with tag: " << ysID2 << endl;
		return TCL_ERROR;
	}

// 		Inelastic2DYS(	int tag, double A, double E, double I, int Nd1, int Nd2,
// 						YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
// 						int rf_algo = -1, // updated

	Element *theElement = new Inelastic2DYS05(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

	if(tcl_debug) cout << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		cerr << "WARNING ran out of memory creating element\n";
		cerr << "element2dYS05: " << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false)
	{
		cerr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		cerr << tag << endl;
		cin.get();
		return TCL_ERROR;
	}

	if(tcl_debug) cout << "\tElement number " << tag << " added to domain - returning\n";

	return TCL_OK;
}
*/

/*******************************************************************************************/
int
TclModelBuilder_addElement2dYS (ClientData clientData, Tcl_Interp *interp,
								   int argc, char **argv,
								   Domain *theTclDomain, TclModelBuilder *theTclBuilder)
{

  if (strcmp(argv[1],"inelastic2dYS01") == 0) {
	  int result = TclModelBuilder_addElement2dYS01(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }
  else if (strcmp(argv[1],"inelastic2dYS02") == 0) {
	  int result = TclModelBuilder_addElement2dYS02(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }
  
   else if (strcmp(argv[1],"inelastic2dYS03") == 0) {
	  int result = TclModelBuilder_addElement2dYS03(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }
  
  /*else if (strcmp(argv[1],"inelastic2dYS04") == 0) {
	  int result = TclModelBuilder_addElement2dYS04
	  (clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }
  else if (strcmp(argv[1],"inelastic2dYS05") == 0) {
	  int result = TclModelBuilder_addElement2dYS05
	  (clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }*/
  else

	return TCL_ERROR;

}
