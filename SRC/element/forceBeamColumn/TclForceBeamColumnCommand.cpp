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
                                                                        
// $Revision: 1.3 $
// $Date: 2002-12-19 21:06:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/TclForceBeamColumnCommand.cpp,v $
                                                                        
// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the implementation of the 
// TclModelBuilder_addDispBeamColumn() command. 

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <TclModelBuilder.h>

#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn3d.h>

#include <LobattoBeamIntegration.h>
#include <UserDefinedBeamIntegration.h>

#include <HingeMidpointBeamIntegration2d.h>
#include <HingeMidpointBeamIntegration3d.h>
#include <HingeRadauBeamIntegration2d.h>
#include <HingeRadauBeamIntegration3d.h>
#include <UserDefinedHingeIntegration2d.h>
#include <UserDefinedHingeIntegration3d.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addForceBeamColumn(ClientData clientData, Tcl_Interp *interp,  
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
	 << " not compatible with forceBeamColumn element" << endl;
    return TCL_ERROR;
  }
  
  if (argc < 6) {
    cerr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    cerr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? ...\n";
    return TCL_ERROR;
  }
  
  // get the id and end nodes 
  int eleTag, iNode, jNode, transfTag;
  
  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    cerr << "WARNING invalid forceBeamColumn eleTag" << endl;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    cerr << "WARNING invalid iNode\n";
    cerr << "forceBeamColumn element: " << eleTag << endl;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    cerr << "WARNING invalid jNode\n";
    cerr << "forceBeamColumn element: " << eleTag << endl;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[5], &transfTag) != TCL_OK) {
    cerr << "WARNING invalid transfTag\n";
    cerr << "forceBeamColumn element: " << eleTag << endl;
    return TCL_ERROR;
  }
  
  CrdTransf2d *theTransf2d = 0;
  CrdTransf3d *theTransf3d = 0;
  
  if (ndm == 2) {
    
    theTransf2d = theTclBuilder->getCrdTransf2d(transfTag);
    
    if (theTransf2d == 0) {
      cerr << "WARNING transformation not found\n";
      cerr << "transformation: " << transfTag;
      cerr << "\nforceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
  }
  
  if (ndm == 3) {
    
    theTransf3d = theTclBuilder->getCrdTransf3d(transfTag);
    
    if (theTransf3d == 0) {
      cerr << "WARNING transformation not found\n";
      cerr << "transformation: " << transfTag;
      cerr << "\nforceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
  }
  
  Element *theElement = 0;
  
  if (strcmp(argv[6],"Lobatto") == 0) {
    int secTag, nIP;
    
    if (argc < 9) {
      cerr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      cerr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? Lobatto secTag? nIP?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTag) != TCL_OK) {
      cerr << "WARNING invalid secTag\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8], &nIP) != TCL_OK) {
      cerr << "WARNING invalid nIP\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    
    SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
    if (theSection == 0) {
      cerr << "WARNING section not found\n";
      cerr << "Section: " << secTag;
      cerr << "\nforceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (int i = 0; i < nIP; i++)
      sections[i] = theSection;
    
    LobattoBeamIntegration beamIntegr;

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);

    delete [] sections;
  }

  else if (strcmp(argv[6],"HingeMidpoint") == 0 ||
	   strcmp(argv[6],"HingeRadau") == 0) {
    
    if (argc < 14) {
      cerr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      cerr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? type secTagI? lpI? secTagJ? lpJ? E? A? Iz? <Iy? G? J?>\n";
      return TCL_ERROR;
    }

    int secTagI, secTagJ;
    double lpI, lpJ;
    double E, A, Iz, Iy, G, J;
    
    if (Tcl_GetInt(interp, argv[7], &secTagI) != TCL_OK) {
      cerr << "WARNING invalid secTagI\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &lpI) != TCL_OK) {
      cerr << "WARNING invalid lpI\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &secTagJ) != TCL_OK) {
      cerr << "WARNING invalid secTagJ\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &lpJ) != TCL_OK) {
      cerr << "WARNING invalid lpJ\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11], &E) != TCL_OK) {
      cerr << "WARNING invalid E\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[12], &A) != TCL_OK) {
      cerr << "WARNING invalid A\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[13], &Iz) != TCL_OK) {
      cerr << "WARNING invalid I\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    
    if (ndm == 3 && argc > 16) {
      if (Tcl_GetDouble(interp, argv[14], &Iy) != TCL_OK) {
	cerr << "WARNING invalid Iy\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[15], &G) != TCL_OK) {
	cerr << "WARNING invalid G\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[16], &J) != TCL_OK) {
	cerr << "WARNING invalid J\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
    }

    SectionForceDeformation *sectionI = theTclBuilder->getSection(secTagI);
    if (sectionI == 0) {
      cerr << "WARNING section not found\n";
      cerr << "Section: " << secTagI;
      cerr << "\nforceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    SectionForceDeformation *sectionJ = theTclBuilder->getSection(secTagJ);
    if (sectionJ == 0) {
      cerr << "WARNING section not found\n";
      cerr << "Section: " << secTagJ;
      cerr << "\nforceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    
    SectionForceDeformation *sections[2];
    sections[0] = sectionI;
    sections[1] = sectionJ;

    BeamIntegration *beamIntegr = 0;

    if (ndm == 2) {
      if (strcmp(argv[6],"HingeMidpoint") == 0)
	beamIntegr = new HingeMidpointBeamIntegration2d(E, A, Iz, lpI, lpJ);
      else
	beamIntegr = new HingeRadauBeamIntegration2d(E, A, Iz, lpI, lpJ);
      
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, 2, sections,
					 *beamIntegr, *theTransf2d);
    }
    else {
      if (strcmp(argv[6],"HingeMidpoint") == 0)
	beamIntegr =
	  new HingeMidpointBeamIntegration3d(E, A, Iz, Iy, G, J, lpI, lpJ);
      else
	beamIntegr =
	  new HingeRadauBeamIntegration3d(E, A, Iz, Iy, G, J, lpI, lpJ);
      
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, 2, sections,
					 *beamIntegr, *theTransf3d);
    }

    delete beamIntegr;
  }

  else if (strcmp(argv[6],"UserDefined") == 0) {

    if (argc < 9) {
      cerr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      cerr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? UserDefined nIP? secTag1? ... pt1? ... wt1? ...\n";
      return TCL_ERROR;
    }

    int nIP;
    
    if (Tcl_GetInt(interp, argv[7], &nIP) != TCL_OK) {
      cerr << "WARNING invalid nIP\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    
    ID secs(nIP);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    for (i = 0, j = 8; i < nIP; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	cerr << "WARNING invalid sec\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+nIP], &pt) != TCL_OK) {
	cerr << "WARNING invalid pt\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*nIP], &wt) != TCL_OK) {
	cerr << "WARNING invalid wt\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      secs(i) = sec;
      pts(i)  = pt;
      wts(i)  = wt;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	cerr << "WARNING section not found\n";
	cerr << "Section: " << secs(i);
	cerr << "\nforceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }
    
    UserDefinedBeamIntegration beamIntegr(nIP, pts, wts);

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    
    delete [] sections;
  }

  else if (strcmp(argv[6],"UserHinge") == 0) {

    if (argc < 9) {
      cerr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      cerr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? UserHinge E? A? Iz? <Iy? G? J?> npL? secTagL1? ... ptL1? ... wtL1? ... npR? secTagR1? ... ptR1? ... wtR1? ...\n";
      return TCL_ERROR;
    }

    double E, A, Iz, Iy, G, J;
    
    if (Tcl_GetDouble(interp, argv[7], &E) != TCL_OK) {
      cerr << "WARNING invalid E\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &A) != TCL_OK) {
      cerr << "WARNING invalid A\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9], &Iz) != TCL_OK) {
      cerr << "WARNING invalid I\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }

    int argStart = 10;
    if (ndm == 3) {
      if (Tcl_GetDouble(interp, argv[10], &Iy) != TCL_OK) {
	cerr << "WARNING invalid I\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[11], &G) != TCL_OK) {
	cerr << "WARNING invalid G\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[12], &J) != TCL_OK) {
	cerr << "WARNING invalid J\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      argStart = 13;
    }

    int npL, npR;
      
    if (Tcl_GetInt(interp, argv[argStart], &npL) != TCL_OK) {
      cerr << "WARNING invalid npL\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[argStart+3*npL+1], &npR) != TCL_OK) {
      cerr << "WARNING invalid npR\n";
      cerr << "forceBeamColumn element: " << eleTag << endl;
      return TCL_ERROR;
    }

    int nIP = npL+npR;

    ID secs(nIP);
    Vector ptsL(npL);
    Vector wtsL(npL);
    Vector ptsR(npR);
    Vector wtsR(npR);

    int i, j;
    for (i = 0, j = argStart+1; i < npL; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	cerr << "WARNING invalid sec\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+npL], &pt) != TCL_OK) {
	cerr << "WARNING invalid pt\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*npL], &wt) != TCL_OK) {
	cerr << "WARNING invalid wt\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      secs(i) = sec;
      ptsL(i) = pt;
      wtsL(i) = wt;
    }
    for (i = 0, j = 1+(argStart+1)+3*npL; i < npR; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	cerr << "WARNING invalid sec\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+npR], &pt) != TCL_OK) {
	cerr << "WARNING invalid pt\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*npR], &wt) != TCL_OK) {
	cerr << "WARNING invalid wt\n";
	cerr << "forceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      secs(i+npL) = sec;
      ptsR(i)     = pt;
      wtsR(i)     = wt;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	cerr << "WARNING section not found\n";
	cerr << "Section: " << secs(i);
	cerr << "\nforceBeamColumn element: " << eleTag << endl;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }
    
    if (ndm == 2) {
      UserDefinedHingeIntegration2d beamIntegr(npL, ptsL, wtsL,
					       npR, ptsR, wtsR,
					       E, A, Iz);
    
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    }
    else {
      UserDefinedHingeIntegration3d beamIntegr(npL, ptsL, wtsL,
					       npR, ptsR, wtsR,
					       E, A, Iz, Iy, G, J);
    
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    }
    
    delete [] sections;
  }

  else {
    cerr << "Unknown integration type: " << argv[6] << endl;
    cerr << "forceBeamColumn element: " << eleTag << endl;
    return TCL_ERROR;
  }
  
  if (theElement == 0) {
    cerr << "WARNING ran out of memory creating element\n";
    cerr << "forceBeamColumn element: " << eleTag << endl;
    return TCL_ERROR;
  }
  
  if (theTclDomain->addElement(theElement) == false) {
    cerr << "WARNING could not add element to the domain\n";
    cerr << "forceBeamColumn element: " << eleTag << endl;
    delete theElement;
    return TCL_ERROR;
  }
  
  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}
