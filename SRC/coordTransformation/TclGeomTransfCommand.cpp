/* ****************************************************************** **
**    Opensees - Open System for Earthquake Engineering Simulation    **
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
// $Date: 2002-09-25 20:40:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/TclGeomTransfCommand.cpp,v $

#include <TclModelBuilder.h>

#include <LinearCrdTransf2d.h>
#include <LinearCrdTransf3d.h>
#include <PDeltaCrdTransf2d.h>
#include <PDeltaCrdTransf3d.h>
#include <CorotCrdTransf2d.h>
#include <CorotCrdTransf3d.h>

static Domain *theTclModelBuilderDomain = 0;
static TclModelBuilder *theTclModelBuilder = 0;

// 
// to create a coordinate transformation 
//
int
TclModelBuilder_addGeomTransf(ClientData clientData, Tcl_Interp *interp,
			      int argc, char **argv,
			      Domain *theDomain,
			      TclModelBuilder *theBuilder)
  
{
  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    cerr << "WARNING insufficient number of geomTransf arguments\n";
    cerr << "Want: geomTransf type? tag? <specific transf args>" << endl;
    return TCL_ERROR;
  }   
	
  theTclModelBuilderDomain = theDomain;
  theTclModelBuilder = theBuilder;
    
  int NDM, NDF;
     
  NDM = theTclModelBuilder->getNDM();   // dimension of the structure (1d, 2d, or 3d)
  NDF = theTclModelBuilder->getNDF();   // number of degrees of freedom per node

  // create 2d coordinate transformation
  if (NDM == 2 && NDF == 3) {

    int crdTransfTag;
    Vector jntOffsetI(2), jntOffsetJ(2);
	 
    if (argc < 3) {
      interp->result = "WARNING insufficient arguments - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>"; 
      return TCL_ERROR;
    }
	    
    int argi = 2;  
    if (Tcl_GetInt(interp, argv[argi++], &crdTransfTag) != TCL_OK) {	
      interp->result = "WARNING invalid tag - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>";
      return  TCL_ERROR;
    }

    // allow additional options at end of command
    int i;
    
    while (argi != argc) {
      if (strcmp(argv[argi],"-jntOffset") == 0) {
	argi++;
	for (i = 0; i < 2; i++) {
	  if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) {
	    interp->result = "WARNING invalid jntOffset value - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>";
	    return TCL_ERROR;
	  }
	}
 
	for (i = 0; i < 2; i++) {
	  if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) {
	    interp->result = "WARNING invalid jntOffset value - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>";
	    return TCL_ERROR;
	  }
	}
      }
      
      else {
	interp->result = "WARNING bad command - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>";
	cerr << "invalid: " << argv[argi] << endl;
	return TCL_ERROR;
      }
    }

    // construct the transformation object
    
    CrdTransf2d *crdTransf2d;
    
    if (strcmp(argv[1],"Linear") == 0)
      crdTransf2d = new LinearCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ);
    
    else if (strcmp(argv[1],"PDelta") == 0 || strcmp(argv[1],"LinearWithPDelta") == 0)
      crdTransf2d = new PDeltaCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ);
    
#ifdef _COROTATIONAL
    else if (strcmp(argv[1],"Corotational") == 0)
      crdTransf2d = new CorotCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ);
#endif
    else {
      interp->result = "WARNING TclElmtBuilder - addGeomTransf - invalid Type";
      cerr << argv[1] << endl;
      return TCL_ERROR;
    }
    
    if (crdTransf2d == 0) {
      interp->result = "WARNING TclElmtBuilder - addGeomTransf - ran out of memory to create geometric transformation object";
      return TCL_ERROR;
    }
    
    // add the transformation to the modelBuilder
    if (theTclModelBuilder->addCrdTransf2d(*crdTransf2d)) {
      interp->result = "WARNING TclElmtBuilder - addGeomTransf  - could not add geometric transformation to model Builder";
      return TCL_ERROR;
    }
  }

  else if  (NDM == 3 && NDF == 6) {
    int crdTransfTag;
    Vector vecxzPlane(3);                  // vector that defines local xz plane
    Vector jntOffsetI(3), jntOffsetJ(3);   // joint offsets in global coordinates
    
    if (argc < 6) {
      interp->result = "WARNING insufficient arguments - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
      return TCL_ERROR;
    }
    
    int argi = 2;  
    if (Tcl_GetInt(interp, argv[argi++], &crdTransfTag) != TCL_OK) {	
      interp->result = "WARNING invalid tag - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
      return  TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(0)) != TCL_OK) {
      interp->result = "WARNING invalid vecxzPlaneX - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(1)) != TCL_OK) {
      interp->result = "WARNING invalid vecxzPlaneY - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(2)) != TCL_OK) {
      interp->result = "WARNING invalid vecxzPlaneZ - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
      return TCL_ERROR;
    }
    
    // allow additional options at end of command
    int i;
    
    while (argi != argc) {
      if (strcmp(argv[argi],"-jntOffset") == 0) {
	argi++;
	for (i = 0; i < 3; i++) {
	  if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) {
	    interp->result = "WARNING invalid jntOffset value - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >"; 
	    return TCL_ERROR;
	  }
	}
	
	for (i = 0; i < 3; i++) {
	  if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) {
	    interp->result = "WARNING invalid jntOffset value - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";  
	    return TCL_ERROR;
	  }
	}
      }
      else {
	interp->result = "WARNING bad command - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >"; 
	cerr << "invalid: " << argv[argi] << endl;
	return TCL_ERROR;
      }
    }
    
    // construct the transformation object
    
    CrdTransf3d *crdTransf3d;
    
    if (strcmp(argv[1],"Linear") == 0)
      crdTransf3d = new LinearCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);
    
    else if (strcmp(argv[1],"PDelta") == 0 || strcmp(argv[1],"LinearWithPDelta") == 0)
      crdTransf3d = new PDeltaCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);
    
#ifdef _COROTATIONAL
    else if (strcmp(argv[1],"Corotational") == 0)
      crdTransf3d = new CorotCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);
#endif
    
    else {
      interp->result = "WARNING TclElmtBuilder - addGeomTransf - invalid Type";
      return TCL_ERROR;
    }
    
    if (crdTransf3d == 0) {
      interp->result = "WARNING TclElmtBuilder - addGeomTransf - ran out of memory to create geometric transformation object";
      return TCL_ERROR;
    }
    
    // add the transformation to the modelBuilder
    if (theTclModelBuilder->addCrdTransf3d(*crdTransf3d)) {
      interp->result = "WARNING TclElmtBuilder - addGeomTransf  - could not add geometric transformation to model Builder";
      return TCL_ERROR;
    }
  }
  else {
    cerr << "WARNING NDM = " << NDM << " and NDF = " << NDF << "is imcompatible with available frame elements";
    return TCL_ERROR;
  }
  
  //  Tcl_Free ((char *)argv);
  
  // if get here we have sucessfully created the element and added it to the domain
  
  return TCL_OK;
}
