//
// Description: This file contains the implementation of the TclModelBuilder_addFourNodeQuadUP()
// command. 
//

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <FourNodeQuadUP.h>

#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

/*  *****************************************************************************
    
    Q U A D  U_P

    ***************************************************************************** */

int
TclModelBuilder_addFourNodeQuadUP(ClientData clientData, Tcl_Interp *interp,  
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

	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 3) {
		cerr << "WARNING -- model dimensions and/or nodal DOF not compatible with QuadUP element\n";
		return TCL_ERROR;
	}

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc-argStart) < 12) {
    cerr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    cerr << "Want: element FourNodeQuadUP eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? bulk? rho? perm_x? perm_y? <b1? b2? pressure? dM? dK?>\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int FourNodeQuadUPId, iNode, jNode, kNode, lNode, matID;
  double thickness, bk, r, perm1, perm2;
	double p = 0.0;		// uniform normal traction (pressure)
	double b1 = 0.0;
	double b2 = 0.0;
	double dM = 0.0;
	double dK = 0.0;

  char *type;
  if (Tcl_GetInt(interp, argv[argStart], &FourNodeQuadUPId) != TCL_OK) {
    cerr << "WARNING invalid FourNodeQuadUP eleTag" << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1+argStart], &iNode) != TCL_OK) {
    cerr << "WARNING invalid iNode\n";
    cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2+argStart], &jNode) != TCL_OK) {
     cerr << "WARNING invalid jNode\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[3+argStart], &kNode) != TCL_OK) {
     cerr << "WARNING invalid kNode\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }  
  
  if (Tcl_GetInt(interp, argv[4+argStart], &lNode) != TCL_OK) {
     cerr << "WARNING invalid lNode\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[5+argStart], &thickness) != TCL_OK) {
     cerr << "WARNING invalid thickness\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }  
  
  type = argv[6+argStart];
  
  if (Tcl_GetInt(interp, argv[7+argStart], &matID) != TCL_OK) {
     cerr << "WARNING invalid matID\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }

	if (Tcl_GetDouble(interp, argv[8+argStart], &bk) != TCL_OK) {
     cerr << "WARNING invalid fluid bulk modulus\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[9+argStart], &r) != TCL_OK) {
     cerr << "WARNING invalid fluid mass density\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[10+argStart], &perm1) != TCL_OK) {
     cerr << "WARNING invalid lateral permeability\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[11+argStart], &perm2) != TCL_OK) {
     cerr << "WARNING invalid vertical permeability\n";
     cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
     return TCL_ERROR;
  }  

	if ((argc-argStart) >= 13) {
		if (Tcl_GetDouble(interp, argv[12+argStart], &b1) != TCL_OK) {
			cerr << "WARNING invalid b1\n";
			cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 14) {
		if (Tcl_GetDouble(interp, argv[13+argStart], &b2) != TCL_OK) {
			cerr << "WARNING invalid b2\n";
			cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 15) {
		if (Tcl_GetDouble(interp, argv[14+argStart], &p) != TCL_OK) {
			cerr << "WARNING invalid pressure\n";
			cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 16) {
		if (Tcl_GetDouble(interp, argv[15+argStart], &dM) != TCL_OK) {
			cerr << "WARNING invalid dM\n";
			cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 17) {
		if (Tcl_GetDouble(interp, argv[16+argStart], &dK) != TCL_OK) {
			cerr << "WARNING invalid dK\n";
			cerr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endl;
			return TCL_ERROR;
		}
	}

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
      
  if (theMaterial == 0) {
      cerr << "WARNING material not found\n";
      cerr << "Material: " << matID;
      cerr << "\nFourNodeQuadUP element: " << FourNodeQuadUPId << endl;
      return TCL_ERROR;
  }
  
  // now create the FourNodeQuadUP and add it to the Domain
  FourNodeQuadUP *theFourNodeQuadUP = 
      new FourNodeQuadUP(FourNodeQuadUPId,iNode,jNode,kNode,lNode,
		       *theMaterial, type, thickness, bk, r, perm1, perm2, b1, b2, p, dM, dK);
  if (theFourNodeQuadUP == 0) {
      cerr << "WARNING ran out of memory creating element\n";
      cerr << "FourNodeQuad element: " << FourNodeQuadUPId << endl;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theFourNodeQuadUP) == false) {
      cerr << "WARNING could not add element to the domain\n";
      cerr << "FourNodeQuad element: " << FourNodeQuadUPId << endl;
      delete theFourNodeQuadUP;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}


