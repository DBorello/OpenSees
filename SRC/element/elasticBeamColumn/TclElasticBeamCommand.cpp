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
// $Date: 2001-08-31 17:18:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/elasticBeamColumn/TclElasticBeamCommand.cpp,v $
                                                                        
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the elasticBeamColumn element.

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <ElasticBeam2d.h>
#include <ElasticBeam3d.h>

#include <CrdTransf2d.h>
#include <CrdTransf3d.h>

#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addElasticBeam(ClientData clientData, Tcl_Interp *interp, int argc, 
			       char **argv, Domain *theTclDomain, TclModelBuilder *theTclBuilder,
			       int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - elasticBeamColumn \n";    
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  Element *theBeam = 0;

  if (ndm == 2) {
    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      cerr << "WARNING invalid ndf: " << ndf;
      cerr << ", for plane problem need 3 - elasticBeamColumn \n";    
      return TCL_ERROR;
    } 

    // check the number of arguments
    if ((argc-eleArgStart) < 8) {
      cerr << "WARNING bad command - want: elasticBeamColumn beamId iNode jNode A E I transTag\n";
      printCommand(argc, argv);
      return TCL_ERROR;
    }    

    // get the id, end nodes, and section properties
    int beamId, iNode, jNode, transTag;
    double A,E,I;
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &beamId) != TCL_OK) {
      cerr << "WARNING invalid beamId: " << argv[1+eleArgStart];
      cerr << " - elasticBeamColumn beamId iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
      cerr << "WARNING invalid iNode - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
      cerr << "WARNING invalid jNode - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
      cerr << "WARNING invalid A - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5+eleArgStart], &E) != TCL_OK) {
      cerr << "WARNING invalid E - elasticBeam " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[6+eleArgStart], &I) != TCL_OK) {
      cerr << "WARNING invalid I - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetInt(interp, argv[7+eleArgStart], &transTag) != TCL_OK) {
      cerr << "WARNING invalid transTag - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    
    CrdTransf2d *theTrans = theTclBuilder->getCrdTransf2d(transTag);
    
    if (theTrans == 0) {
	cerr << "WARNING transformation object not found - elasticBeamColumn " << beamId;
	return TCL_ERROR;
    }
    
    // now create the beam and add it to the Domain
    theBeam = new ElasticBeam2d (beamId,A,E,I,iNode,jNode, *theTrans);
    
    if (theBeam == 0) {
      cerr << "WARNING ran out of memory creating beam - elasticBeamColumn ";	
      cerr << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      cerr << "WARNING invalid ndof: " << ndf;
      cerr << ", for 3d problem  need 6 - elasticBeamColumn \n";    
      return TCL_ERROR;
    } 

    // check the number of arguments
    if ((argc-eleArgStart) < 11) {
      cerr << "WARNING bad command - want: elasticBeamColumn beamId iNode jNode";
      cerr << " A E G Jx Iy Iz transTag" << endl;
      printCommand(argc, argv);
      return TCL_ERROR;
    }    

    // get the id, end nodes, and section properties
    int beamId, iNode, jNode, transTag;
    double A,E,G,Jx,Iy,Iz;
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &beamId) != TCL_OK) {
      cerr << "WARNING invalid beamId: " << argv[1+eleArgStart];
      cerr << " - elasticBeamColumn beamId iNode jNode A E G Jx Iy Iz\n ";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
      cerr << "WARNING invalid iNode - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
      cerr << "WARNING invalid jNode - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
      cerr << "WARNING invalid A - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5+eleArgStart], &E) != TCL_OK) {
      cerr << "WARNING invalid E - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6+eleArgStart], &G) != TCL_OK) {
      cerr << "WARNING invalid G - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetDouble(interp, argv[7+eleArgStart], &Jx) != TCL_OK) {
      cerr << "WARNING invalid Jx - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetDouble(interp, argv[8+eleArgStart], &Iy) != TCL_OK) {
      cerr << "WARNING invalid Iy - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetDouble(interp, argv[9+eleArgStart], &Iz) != TCL_OK) {
      cerr << "WARNING invalid Iz - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetInt(interp, argv[10+eleArgStart], &transTag) != TCL_OK) {
      cerr << "WARNING invalid transTag - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }      

    CrdTransf3d *theTrans = theTclBuilder->getCrdTransf3d(transTag);
    
    if (theTrans == 0) {
	cerr << "WARNING transformation object not found - elasticBeamColumn " << beamId;
	return TCL_ERROR;
    }
    
    // now create the beam and add it to the Domain
    theBeam = new ElasticBeam3d (beamId,A,E,G,Jx,Iy,Iz,iNode,jNode, *theTrans);
    
    if (theBeam == 0) {
      cerr << "WARNING ran out of memory creating beam - elasticBeamColumn ";	
      cerr << beamId << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    
  }

  else {
    cerr << "WARNING elasticBeamColumn command only works when ndm is 2 or 3, ndm: ";
    cerr << ndm << endl;
    return TCL_ERROR;
  }

  // now add the beam to the domain
  if (theTclDomain->addElement(theBeam) == false) {
    cerr << "WARNING TclModelBuilder - addBeam - could not add beam to domain ";
    cerr << *theBeam;
    delete theBeam; // clean up the memory to avoid leaks
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}
