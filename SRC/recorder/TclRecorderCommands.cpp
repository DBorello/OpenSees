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
// $Date: 2001-05-16 04:19:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/TclRecorderCommands.cpp,v $
                                                                        
                                                                        
// File: ~/recorders/TclRecordersCommand.C
// 
// Written: fmk 
// Created: 04/98
// Revision: A
//
// Description: This file contains the function that is invoked
// by the interpreter when the comand 'record' is invoked by the 
// user.
//
// What: "@(#) commands.C, revA"


#include <tcl.h>
#include <tk.h>

#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <Domain.h>
#include <EquiSolnAlgo.h>

// recorders
#include <MaxNodeDispRecorder.h>
#include <NodeRecorder.h>
#include <ElementRecorder.h>
#include <TclFeViewer.h>
#include <FilePlotter.h>
#include <AlgorithmIncrements.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <Node.h>
#include <Element.h>

#include <EquiSolnAlgo.h>

static EquiSolnAlgo *theAlgorithm =0;

int 
TclCreateRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
		  char **argv, Domain &theDomain, Recorder **theRecorder)
{
    // make sure at least one other argument to contain integrator
    if (argc < 2) {
	interp->result = "WARNING need to specify a Recorder type ";
	return TCL_ERROR;
    }    

    //
    // check argv[1] for type of Recorder, parse in rest of arguments
    // needed for the type of Recorder, create the object and add to Domain
    //

    (*theRecorder) = 0;

    // an Element Recorder
    if (strcmp(argv[1],"Element") == 0) {
      int eleID;
	if (argc < 4) {
	    cerr << "WARNING recorder Element eleID1? eleID2? ...  <-time> "
		<< "<-file fileName?> parameters";
	    return TCL_ERROR;
	}    

	int endEleIDs = 2;
	int allFlag = 0;
	while (Tcl_GetInt(interp, argv[endEleIDs], &eleID) == TCL_OK) {
	  endEleIDs++;
	}

	// determine the number of elements
	int numEle = 0;
	if (strcmp(argv[endEleIDs],"all") == 0) {
	  endEleIDs += 1;
	  allFlag = 1;
	  numEle = theDomain.getNumElements();
	} else
	  numEle = endEleIDs-2;
	
        if ((endEleIDs-2) == 0) {
	    cerr << "WARNING recorder Element eleID1? eleID2? .. <-time> "
		<< "<-file fileName?> parameters";	    
	    return TCL_ERROR;
	}	    	    

	// create an ID to hold ele tags
        ID eleIDs(numEle); 

	// read in the ele tags to the ID
	if (allFlag == 1) {
	  int loc = 0;
	  ElementIter &theEleIter = theDomain.getElements();
	  Element *theEle;
	  while ((theEle = theEleIter()) != 0)
	    eleIDs[loc++] = theEle->getTag();
	} else {
	  for (int i=2; i<endEleIDs; i++) {
	    if (Tcl_GetInt(interp, argv[i], &eleID) != TCL_OK)	
	      return TCL_ERROR;	
	    eleIDs[i-2] = eleID;	  
	  }
	}

	bool echoTime = false;
	char *fileName = 0;
	int endMarker = endEleIDs;
	int flags = 0;
	while (flags == 0) {
	  if (strcmp(argv[endMarker],"-time") == 0) {
	    // allow user to specify const load
	    echoTime = true;
	    endMarker++;
	  } 
	  else if (strcmp(argv[endMarker],"-file") == 0) {
	    // allow user to specify load pattern other than current
	    endMarker++;
	    fileName = argv[endMarker];
	    endMarker++;
	  }
	  else
	    flags = 1;
	}
	
	(*theRecorder) = new ElementRecorder(eleIDs, theDomain, &argv[endMarker], 
					  argc-endMarker, echoTime, fileName);
    }
    
    // a MaxNodeDisp Recorder
    else if (strcmp(argv[1],"MaxNodeDisp") == 0) {
	int dof;
	if (argc < 5) {
	    cerr << "WARNING recorder MaxNodeDisp <dof> node <list nodes>";
	    return TCL_ERROR;
	}    
	if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK)	
	    return TCL_ERROR;	
	
	dof--; // subtract 1 for interpreter to C indexing

	int numNodes = argc -4;
	ID theNodes(numNodes);
	for (int i=0; i<numNodes; i++) {
	    int node;
	    if (Tcl_GetInt(interp, argv[i+4], &node) != TCL_OK)	
		return TCL_ERROR;		  
	    theNodes(i) = node;
	}
	
	(*theRecorder) = new MaxNodeDispRecorder(dof, theNodes, theDomain);
    }

    // create a recorder to write nodal displacement quantities to a file
    else if (strcmp(argv[1],"Node") == 0) {
	if (argc < 7) {
	    cerr << "WARNING recorder Node filename? response? <startFlag> ";
	    cerr << "-node <list nodes> -dof <doflist>";
	    return TCL_ERROR;
	}    


	char *responseID = argv[3];

	int flag = 0;
	if (strcmp(argv[4],"-time") == 0) 
	    flag = 1;
	if (strcmp(argv[4],"-load") == 0)
	    flag = 2;      
      
	int pos = 5;
	if (flag != 0) pos = 6;

	int numNodes = 0;
	int allFlag = 0;
	// get the number of nodes;
	if (strcmp(argv[pos],"all") == 0) {
	  numNodes = theDomain.getNumNodes();
	  allFlag = 1;
	} else {
	  numNodes = 0;
	  for (int j=pos; j< argc; j++) 
	    if (strcmp(argv[j],"-dof") != 0)
	      numNodes++;
	    else
	      j = argc+2;
	  if (numNodes >= (argc-pos)){
	    cerr << "WARNING no -dof to indicate end of node list ";
	    cerr << "- recorder Node filename? response? <startFlag> ";
	    cerr << "-node <list nodes> -dof <doflist>";
	    return TCL_ERROR;
	  }
	}

	// create an id containing node tags
	ID theNodes(numNodes);

	// read in the node tags
	if (allFlag == 1) {
	  NodeIter &theNodeIter = theDomain.getNodes();
	  Node *theNode;
	  int loc=0;
	  while ((theNode= theNodeIter()) != 0) {
	    int tag = theNode->getTag();
	    theNodes(loc++) = tag;
	  }
	  pos += 2;
	} else {
	  for (int i=0; i<numNodes; i++) {
	    int node;
	    if (Tcl_GetInt(interp, argv[i+pos], &node) != TCL_OK)	
	      return TCL_ERROR;		  
	    theNodes(i) = node;
	  }
	  pos += numNodes+1;
	}

	int numDof = argc - pos;
	ID theDofs(numDof);
	for (int k=0; k<numDof; k++) {
	    int dof;
	    if (Tcl_GetInt(interp, argv[k+pos], &dof) != TCL_OK)	
		return TCL_ERROR;		  
	    dof--;  // subtract 1 for interpreter to C indexing
	    theDofs(k) = dof;
	}	  

	(*theRecorder) = new NodeRecorder(theDofs, theNodes, 
					  theDomain,
					  argv[2], responseID, flag);
    } 


    // a recorder for the graphical display of the domain
    else if (strcmp(argv[1],"display") == 0) {

	{
	int xLoc, yLoc, width, height;
	if (argc < 7) {
	    cerr << "WARNING recorder display title xLoc yLoc pixelsX pixelsY <-file fileName?>";
	    return TCL_ERROR;
	}    
	if (Tcl_GetInt(interp, argv[3], &xLoc) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[4], &yLoc) != TCL_OK)	
	    return TCL_ERROR;	      
	if (Tcl_GetInt(interp, argv[5], &width) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[6], &height) != TCL_OK)	
	    return TCL_ERROR;
	
	// check if we are to wipe image on each redraw
	int wipeFlag = 0;
	if (argc == 8) 
	  if (strcmp(argv[7],"-wipe") == 0) 
	    wipeFlag = 1;

	if (argc == 7 || argc == 8)
	  (*theRecorder) = new TclFeViewer(argv[2], xLoc, yLoc, width, height, theDomain, wipeFlag, interp);
	else if (argc == 9)
	  (*theRecorder) = new TclFeViewer(argv[2], xLoc, yLoc, width, height, argv[8], theDomain, interp);
        }

    }  


    else if (strcmp(argv[1],"plot") == 0) {

	int xLoc, yLoc, width, height;
	if (argc < 9) {
	    cerr << "WARNING recorder display fileName? windowTitle? xLoc yLoc pixelsX pixelsY -columns colX1 colY1 -columns colX2 ...";
	    return TCL_ERROR;
	}    
	if (Tcl_GetInt(interp, argv[4], &xLoc) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[5], &yLoc) != TCL_OK)	
	    return TCL_ERROR;	      
	if (Tcl_GetInt(interp, argv[6], &width) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[7], &height) != TCL_OK)	
	    return TCL_ERROR;	      

	FilePlotter *thePlotter = new FilePlotter(argv[2], argv[3], xLoc, yLoc, width, height);
	(*theRecorder) = thePlotter;


	int endMarker = 8;

	int loc = 0;
	ID cols(0,16);
	while (endMarker < argc) {
	  if (strcmp(argv[endMarker],"-columns") == 0) {
	    if (argc < endMarker+2)
	      return TCL_ERROR;

	    int colX, colY;
	    endMarker++;
	    if (Tcl_GetInt(interp, argv[endMarker], &colX) != TCL_OK)	
	      return TCL_ERROR;	
	    endMarker++;
	    if (Tcl_GetInt(interp, argv[endMarker], &colY) != TCL_OK)	
	      return TCL_ERROR;	
	    endMarker++;
	    cols[loc] = colX;
	    cols[loc+1] = colY;
	    loc+=2;
	  } 
	  else
	    endMarker++;
	}
	thePlotter->setCol(cols);
    }


    else if (strcmp(argv[1],"increments") == 0) {

	int xLoc, yLoc, width, height;
	
	if (theAlgorithm == 0) {
	    cerr << "WARNING recorder increments - only allowed as algorithmRecorder";
	    return TCL_ERROR;
	}
	if (argc < 7) {
	    cerr << "WARNING recorder display windowTitle? xLoc yLoc pixelsX pixelsY ";
	    return TCL_ERROR;
	}    
	if (Tcl_GetInt(interp, argv[3], &xLoc) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[4], &yLoc) != TCL_OK)	
	    return TCL_ERROR;	      
	if (Tcl_GetInt(interp, argv[5], &width) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[6], &height) != TCL_OK)	
	    return TCL_ERROR;	      

	char *fileName = 0;
	bool displayRecord = false;
	int endMarker = 7;
	while (endMarker < argc) {
	  if ((strcmp(argv[endMarker],"-file") == 0) ||
	      (strcmp(argv[endMarker],"-file") == 0)) {
	    
	      if (argc < endMarker+1)
		return TCL_ERROR;
	      endMarker++;
	      fileName = argv[endMarker];
	      endMarker++;
	  } else if (strcmp(argv[endMarker],"-display") == 0) {
	    displayRecord = true;
	    endMarker++;				   
	  } 
	  else
	    endMarker++;
	}

	
	AlgorithmIncrements *thePlotter =  new AlgorithmIncrements(theAlgorithm, 
								   argv[2], xLoc, yLoc, width, height, 
								   displayRecord, fileName);
	(*theRecorder) = thePlotter;

    }




    // no recorder type specified yet exists
    else {
	cerr << "WARNING No recorder type exists ";
	cerr << "for recorder of type:" << argv[1];
    
	return TCL_ERROR;
    }    

    // check we instantiated a recorder .. if not ran out of memory
    if ((*theRecorder) == 0) {
	cerr << "WARNING ran out of memory - recorder " << argv[1]<< endl;
	return TCL_ERROR;
    } 
	/*
    // add the recorder to the domain
    if (theDomain.addRecorder(*theRecorder) < 0) {
	cerr << "WARNING could not add to domain - recorder " << argv[1]<< endl;
	delete theRecorder;
	return TCL_ERROR;
    } 
	*/
    // operation successfull
    return TCL_OK;
}


int 
TclAddRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	       char **argv, Domain &theDomain)
{
	Recorder *theRecorder;
	TclCreateRecorder(clientData, interp, argc, argv, theDomain, &theRecorder);

	if ((theRecorder == 0) || (theDomain.addRecorder(*theRecorder)) < 0) {
		cerr << "WARNING could not add to domain - recorder " << argv[1]<< endl;
		if (theRecorder == 0) 
			cerr << "could not create recorder\n";
		else
			delete theRecorder;
		return TCL_ERROR;
	} 
	return TCL_OK;
	
}


int 
TclAddAlgorithmRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	       char **argv, Domain &theDomain, EquiSolnAlgo *theAlgo)
{
	Recorder *theRecorder = 0;
	theAlgorithm = theAlgo;
	if (TclCreateRecorder(clientData, interp, argc, argv, theDomain,
			&theRecorder) == TCL_ERROR) {
		return TCL_ERROR;
	} else {
		// add the recorder to the domain, 
		// NOTE: will not be called with theALgo == 0
		// see ~/g3/SRC/tcl/commands.C file
		if (theRecorder == 0 || theAlgo->addRecorder(*theRecorder) < 0) {
			cerr << "WARNING could not add to algorithm - recorder " << argv[1]<< endl;
			delete theRecorder;
			return TCL_ERROR;
		} 
		return TCL_OK;
	}
}

