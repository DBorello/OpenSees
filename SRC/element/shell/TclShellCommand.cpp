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
// $Date: 2001-05-19 06:00:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/TclShellCommand.cpp,v $
                                                                        
// File: ~/element/TclMembranePlateCommand.C
// 
// Written: fmk 
// Created: 03/01
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_addMembranePlate()
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <ShellMITC4.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addShellMITC4(ClientData clientData, Tcl_Interp *interp,  int argc, 
			 char **argv, Domain*theTclDomain,
			 TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 7) {
    cerr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    cerr << "Want: element ShellMITC4 eleTag? iNode? jNode? kNode? lNode? secTag?\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int ShellMITC4Id, iNode, jNode, kNode, lNode, matID;
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &ShellMITC4Id) != TCL_OK) {
    cerr << "WARNING invalid ShellMITC4 eleTag" << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    cerr << "WARNING invalid iNode\n";
    cerr << "ShellMITC4 element: " << ShellMITC4Id << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     cerr << "WARNING invalid jNode\n";
     cerr << "ShellMITC4 element: " << ShellMITC4Id << endl;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4+eleArgStart], &kNode) != TCL_OK) {
     cerr << "WARNING invalid jNode\n";
     cerr << "ShellMITC4 element: " << ShellMITC4Id << endl;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5+eleArgStart], &lNode) != TCL_OK) {
     cerr << "WARNING invalid jNode\n";
     cerr << "ShellMITC4 element: " << ShellMITC4Id << endl;
     return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6+eleArgStart], &matID) != TCL_OK) {
    cerr << "WARNING invalid matTag\n";
    cerr << "ShellMITC4 element: " << ShellMITC4Id << endl;
    return TCL_ERROR;
  }

  SectionForceDeformation *theSection = theTclBuilder->getSection(matID);      
      
  if (theSection == 0) {
    cerr << "WARNING section not found\n";
    cerr << "section tag: " << matID;
    cerr << "\nShellMITC4 element: " << ShellMITC4Id << endl;
    return TCL_ERROR;
  }
  
  // now create the ShellMITC4 and add it to the Domain
  ShellMITC4 *theShellMITC4 = new ShellMITC4(ShellMITC4Id,iNode,jNode,kNode,lNode,*theSection);
  if (theShellMITC4 == 0) {
    cerr << "WARNING ran out of memory creating element\n";
    cerr << "ShellMITC4 element: " << ShellMITC4Id << endl;
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theShellMITC4) == false) {
    cerr << "WARNING could not add element to the domain\n";
    cerr << "ShellMITC4 element: " << ShellMITC4Id << endl;
    delete theShellMITC4;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

