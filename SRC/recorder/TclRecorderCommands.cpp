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
                                                                        
// $Revision: 1.16 $
// $Date: 2003-02-20 08:00:23 $
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <EquiSolnAlgo.h>

// recorders
#include <NodeRecorder.h>
#include <EnvelopeNodeRecorder.h>
#include <EnvelopeElementRecorder.h>
#include <DriftRecorder.h>
#include <ElementRecorder.h>
#include <TclFeViewer.h>
#include <FilePlotter.h>
#include <AlgorithmIncrements.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <Node.h>
#include <Element.h>
#include <MeshRegion.h>
#include <GSA_Recorder.h>

#include <EquiSolnAlgo.h>

static EquiSolnAlgo *theAlgorithm =0;

int 
TclCreateRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
		  char **argv, Domain &theDomain, Recorder **theRecorder)
{
    // make sure at least one other argument to contain integrator
    if (argc < 2) {
	opserr << "WARNING need to specify a Recorder type\n"; 
	return TCL_ERROR;
    }    

    //
    // check argv[1] for type of Recorder, parse in rest of arguments
    // needed for the type of Recorder, create the object and add to Domain
    //
    (*theRecorder) = 0;

    // an Element Recorder or ElementEnvelope Recorder
    if ((strcmp(argv[1],"Element") == 0) || (strcmp(argv[1],"EnvelopeElement") == 0) 
	|| (strcmp(argv[1],"ElementEnvelope") == 0)) {

        /* KEEP - FOR LEGACY REASONS NEED TO KEEP THE FOLLOWING UGLY STUFF */
        int eleID;
        if (argc < 4) {
	    opserr << "WARNING recorder Element eleID1? eleID2? ...  <-time> "
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
	} else if (Tcl_GetInt(interp, argv[endEleIDs], &eleID) != TCL_OK) 
	  ;
	else
	  endEleIDs++;
	  
	numEle = endEleIDs-2;
	
	// create an ID to hold ele tags
        ID eleIDs(numEle, numEle+1); 

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
        /* ********************* END OF KEEP ****************************** */

	double dT = 0.0;
	bool echoTime = false;
	char *fileName = 0;
	int loc = endEleIDs;
	int flags = 0;
	int eleData = 0;

	while (flags == 0 && loc < argc) {

	  if ((strcmp(argv[loc],"-ele") == 0) ||
	      (strcmp(argv[loc],"-eles") == 0) ||
	      (strcmp(argv[loc],"-element") == 0)) {
      
	    // ensure no segmentation fault if user messes up
	    if (argc < loc+2) {
	      opserr << "WARNING recorder Element .. -ele tag1? .. - no ele tags specified\n";
	      return TCL_ERROR;
	    }
	    
	    //
	    // read in a list of ele until end of command or other flag
	    //
	    loc++;
	    int eleTag;
	    while (loc < argc && Tcl_GetInt(interp, argv[loc], &eleTag) == TCL_OK) {
	      eleIDs[numEle++] = eleTag;
	      loc++;

	    }

	    if (strcmp(argv[loc],"all") == 0) {
	      ElementIter &theEleIter = theDomain.getElements();
	      Element *theEle;
	      while ((theEle = theEleIter()) != 0)
		eleIDs[numEle++] = theEle->getTag();
	      loc++;
	    }

	  } else if (strcmp(argv[loc],"-eleRange") == 0) {
	    
	    // ensure no segmentation fault if user messes up
	    if (argc < loc+3) {
	      opserr << "WARNING recorder Element .. -eleRange start? end?  .. - no ele tags specified\n";
	      return TCL_ERROR;
	    }
	    
	    //
	    // read in start and end tags of two elements & add set [start,end]
	    //
	    
	    int start, end;
	    if (Tcl_GetInt(interp, argv[loc+1], &start) != TCL_OK) {
	      opserr << "WARNING recorder Element -eleRange start? end? - invalid start " << argv[loc+1] << endln;
	      return TCL_ERROR;
	    }      
	    if (Tcl_GetInt(interp, argv[loc+2], &end) != TCL_OK) {
	      opserr << "WARNING recorder Element -eleRange start? end? - invalid end " << argv[loc+2] << endln;
	      return TCL_ERROR;
	    }      
	    if (start > end) {
	      int swap = end;
	      end = start;
	      start = swap;
	    }

	    for (int i=start; i<=end; i++)
	      eleIDs[numEle++] = i;	    

	    loc += 3;
	  } 

	  else if (strcmp(argv[loc],"-region") == 0) {
	    // allow user to specif elements via a region

	    if (argc < loc+2) {
	      opserr << "WARNING recorder Element .. -region tag?  .. - no region specified\n";
	      return TCL_ERROR;
	    }
	    int tag;
	    if (Tcl_GetInt(interp, argv[loc+1], &tag) != TCL_OK) {
	      opserr << "WARNING recorder Element -region tag? - invalid tag " << argv[loc+1] << endln;
	      return TCL_ERROR;
	    }      
	    MeshRegion *theRegion = theDomain.getRegion(tag);
	    if (theRegion == 0) {
	      opserr << "WARNING recorder Element -region " << tag << " - region does not exist" << endln;
	      return TCL_OK;
	    }      
	    const ID &eleRegion = theRegion->getElements();
	    for (int i=0; i<eleRegion.Size(); i++)
	      eleIDs[numEle++] = eleRegion(i);

	    loc += 2;
	  } 

	  else if ((strcmp(argv[loc],"-time") == 0) || (strcmp(argv[loc],"-load") == 0)) { 
	    // allow user to specify const load
	    echoTime = true;
	    loc++;
	  } 

	  else if (strcmp(argv[loc],"-dT") == 0) {
	    // allow user to specify time step size for recording
	    loc++;
	    if (Tcl_GetDouble(interp, argv[loc], &dT) != TCL_OK)	
	      return TCL_ERROR;	
	    loc++;
	  } 
	  else if (strcmp(argv[loc],"-file") == 0) {
	    // allow user to specify load pattern other than current
	    loc++;
	    fileName = argv[loc];
	    loc++;
	  }
	  else {
	    // first unknown string then is assumed to start 
	    // element response request starts
	    eleData = loc;
	    flags = 1;
	  }
	}

	// if user has specified no element tags lets assume he wants them all
	if (numEle == 0) {
	  ElementIter &theEleIter = theDomain.getElements();
	  Element *theEle;
	  while ((theEle = theEleIter()) != 0)
	    eleIDs[numEle++] = theEle->getTag();
	}

	// now construct the recorder
	if (strcmp(argv[1],"Element") == 0) 
	  (*theRecorder) = new ElementRecorder(eleIDs, theDomain, &argv[eleData], 
					       argc-eleData, echoTime, dT, fileName);
	else

	  (*theRecorder) = new EnvelopeElementRecorder(eleIDs, theDomain, &argv[eleData], 
						       argc-eleData, dT, fileName);
    }
    
    // create a recorder to write nodal displacement quantities to a file
    else if ((strcmp(argv[1],"Node") == 0) || (strcmp(argv[1],"EnvelopeNode") == 0) 
	     || (strcmp(argv[1],"NodeEnvelope") == 0)) {	

      if (argc < 7) {
	    opserr << "WARNING recorder Node ";
	    opserr << "-node <list nodes> -dof <doflist> -file <fileName> -dT <dT> reponse";
	    return TCL_ERROR;
	}    

      char *responseID = 0;
      char *fileName = 0;

      int pos = 2;

      /* KEEP - FOR LEGACY REASONS NEED TO KEEP THE FOLLOWING UGLY STUFF */
      if ((strcmp(argv[pos],"-time") != 0) && (strcmp(argv[pos],"-load") != 0) &&
	  (strcmp(argv[pos],"-dT") !=  0) && (strcmp(argv[pos],"-node") != 0) &&
	  (strcmp(argv[pos],"-dof") != 0) && (strcmp(argv[pos],"-file") != 0)) {
	pos = 4;
	responseID = argv[3];
	fileName = argv[2];
      } 
      /* ********************** END OF KEEP ***************************  */

      int timeFlag = 0;
      int flags = 0;
      double dT = 0.0;
      int numNodes = 0;
	
      // create ID's to contain the node tags & the dofs
      ID theNodes(0,16);
      ID theDofs(0, 6);

      while (flags == 0 && pos < argc) {

	if (strcmp(argv[pos],"-time") == 0) {
	  timeFlag = 1;
	  pos++;
	}
	
	else if (strcmp(argv[pos],"-load") == 0) {
	  timeFlag = 2;      
	  pos++;
	}

	else if (strcmp(argv[pos],"-file") == 0) {
	  fileName = argv[pos+1];
	  pos += 2;
	}

	else if (strcmp(argv[pos],"-dT") == 0) {
	  pos ++;
	  if (Tcl_GetDouble(interp, argv[pos], &dT) != TCL_OK)	
	    return TCL_ERROR;		  
	  pos++;
	}

	else if ((strcmp(argv[pos],"-node") == 0) || 
		 (strcmp(argv[pos],"-nodes") == 0)) {
	  pos++;
	  
	  // read in the node tags or 'all' can be used
	  if (strcmp(argv[pos],"all") == 0) {
	    numNodes = theDomain.getNumNodes();
	    
	    NodeIter &theNodeIter = theDomain.getNodes();
	    Node *theNode;
	    int loc=0;
	    while ((theNode= theNodeIter()) != 0) {
	      int tag = theNode->getTag();
	      theNodes[loc++] = tag;
	    }
	    pos++;
	  } else {
	    int node;
	    for (int j=pos; j< argc; j++) 
	      if (Tcl_GetInt(interp, argv[pos], &node) != TCL_OK) 
		j = argc;
	      else {
		theNodes[numNodes] = node;
		numNodes++;
		pos++;
	      }
	  }
	} 

	else if (strcmp(argv[pos],"-nodeRange") == 0) {
	    
	  // ensure no segmentation fault if user messes up
	  if (argc < pos+3) {
	    opserr << "WARNING recorder Node .. -nodeRange start? end?  .. - no ele tags specified\n";
	    return TCL_ERROR;
	  }
	  
	  //
	  // read in start and end tags of two elements & add set [start,end]
	  //
	    
	  int start, end;
	  if (Tcl_GetInt(interp, argv[pos+1], &start) != TCL_OK) {
	    opserr << "WARNING recorder Node -nodeRange start? end? - invalid start " << argv[pos+1] << endln;
	    return TCL_ERROR;
	  }      
	  if (Tcl_GetInt(interp, argv[pos+2], &end) != TCL_OK) {
	    opserr << "WARNING recorder Node -nodeRange start? end? - invalid end " << argv[pos+2] << endln;
	    return TCL_ERROR;
	  }      
	  if (start > end) {
	    int swap = end;
	    end = start;
	    start = swap;
	  }
	  
	  for (int i=start; i<=end; i++)
	    theNodes[numNodes++] = i;	    
	  pos += 3;
	}

	else if (strcmp(argv[pos],"-region") == 0) {
	  // allow user to specif elements via a region
	  
	  if (argc < pos+2) {
	    opserr << "WARNING recorder Node .. -region tag?  .. - no region specified\n";
	    return TCL_ERROR;
	  }
	  int tag;
	  if (Tcl_GetInt(interp, argv[pos+1], &tag) != TCL_OK) {
	    opserr << "WARNING recorder Node -region tag? - invalid tag " << argv[pos+1] << endln;
	    return TCL_ERROR;
	  }      
	  MeshRegion *theRegion = theDomain.getRegion(tag);
	  if (theRegion == 0) {
	    opserr << "WARNING recorder Node -region " << tag << " - region does not exist" << endln;
	    return TCL_OK;
	  }      
	  const ID &nodeRegion = theRegion->getNodes();
	  for (int i=0; i<nodeRegion.Size(); i++)
	    theNodes[numNodes++] = nodeRegion(i);
	  
	  pos += 2;
	} 
	
	else if (strcmp(argv[pos],"-dof") == 0) {
	  pos++;
	  int numDOF = 0;
	  int dof;
	  for (int j=pos; j< argc; j++) 
	    if (Tcl_GetInt(interp, argv[pos], &dof) != TCL_OK) 
	      j = argc;
	    else {
	      theDofs[numDOF] = dof-1;  // -1 for c indexing of the dof's
	      numDOF++;
	      pos++;
	    }
	}
	
	else	 
	  flags = 1;
	}
     
      if (responseID == 0) {
	responseID  = argv[pos];
      }

      if (numNodes == 0) {
	NodeIter &theNodeIter = theDomain.getNodes();
	Node *theNode;
	while ((theNode= theNodeIter()) != 0) {
	  int tag = theNode->getTag();
	  theNodes[numNodes++] = tag;
	}
      }
	
      // construct the recorder
      if (strcmp(argv[1],"Node") == 0) 
	(*theRecorder) = new NodeRecorder(theDofs, theNodes, 
					  theDomain,
					  fileName, responseID, dT, timeFlag);
      else
	
	(*theRecorder) = new EnvelopeNodeRecorder(theDofs, theNodes, 
						  theDomain,
						  fileName, responseID, dT);
    } 

    // Create a recorder to write nodal drifts to a file
    else if (strcmp(argv[1],"Drift") == 0) {
      if (argc < 7) {
	opserr << "WARNING recorder Drift filename? <startFlag> ";
	opserr << "node1? node2? dof? perpDirn?";
	return TCL_ERROR;
      }    
      
      
      int flag = 0;
      if (strcmp(argv[3],"-time") == 0) 
	flag = 1;
      if (strcmp(argv[3],"-load") == 0)
	flag = 2;      
      
      int pos = 3;
      if (flag != 0) pos = 4;
      
      int node1, node2, dof, perpDirn;
      
      if (Tcl_GetInt(interp, argv[pos++], &node1) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetInt(interp, argv[pos++], &node2) != TCL_OK)	
	return TCL_ERROR;	      
      if (Tcl_GetInt(interp, argv[pos++], &dof) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetInt(interp, argv[pos++], &perpDirn) != TCL_OK)	
	return TCL_ERROR;

      // Subtract one from dof and perpDirn for C indexing
      (*theRecorder) = new DriftRecorder(node1, node2, dof-1, perpDirn-1,
					 theDomain,
					 argv[2], flag);
    } 


    // a recorder for the graphical display of the domain
    else if (strcmp(argv[1],"display") == 0) {

	int xLoc, yLoc, width, height;

	if (argc < 7) {
	    opserr << "WARNING recorder display title xLoc yLoc pixelsX pixelsY <-file fileName?>";
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


    else if (strcmp(argv[1],"plot") == 0) {

	int xLoc, yLoc, width, height;
	if (argc < 9) {
	    opserr << "WARNING recorder display fileName? windowTitle? xLoc yLoc pixelsX pixelsY -columns colX1 colY1 -columns colX2 ...";
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

	int loc = 8;

	double dT = 0.0;
	loc = 0;
	ID cols(0,16);
	int numCols = 0;
	while (loc < argc) {
	  if ((strcmp(argv[loc],"-columns") == 0) ||
	      (strcmp(argv[loc],"-cols") == 0) ||
	      (strcmp(argv[loc],"-col") == 0)) {
	    if (argc < loc+2)
	      return TCL_ERROR;

	    int colX, colY;
	    if (Tcl_GetInt(interp, argv[loc+1], &colX) != TCL_OK)	
	      return TCL_ERROR;	

	    if (Tcl_GetInt(interp, argv[loc+2], &colY) != TCL_OK)	
	      return TCL_ERROR;	

	    cols[numCols++] = colX;
	    cols[numCols++] = colY;
	    loc += 3;
	  } 
	  else if (strcmp(argv[loc],"-dT") == 0) {

	    if (Tcl_GetDouble(interp, argv[loc+1], &dT) != TCL_OK)	
	      return TCL_ERROR;	
	    loc += 2;	    
	  }
	  else
	    loc++;
	}

	FilePlotter *thePlotter = new FilePlotter(argv[2], argv[3], xLoc, yLoc, width, height, dT);
	(*theRecorder) = thePlotter;    
	thePlotter->setCol(cols);
    }


    else if (strcmp(argv[1],"increments") == 0) {

	int xLoc, yLoc, width, height;
	
	if (theAlgorithm == 0) {
	    opserr << "WARNING recorder increments - only allowed as algorithmRecorder";
	    return TCL_ERROR;
	}
	if (argc < 7) {
	    opserr << "WARNING recorder display windowTitle? xLoc yLoc pixelsX pixelsY ";
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
	int loc = 7;
	while (loc < argc) {
	  if ((strcmp(argv[loc],"-file") == 0) ||
	      (strcmp(argv[loc],"-file") == 0)) {
	    
	      if (argc < loc+1)
		return TCL_ERROR;
	      loc++;
	      fileName = argv[loc];
	      loc++;
	  } else if (strcmp(argv[loc],"-display") == 0) {
	    displayRecord = true;
	    loc++;				   
	  } 
	  else
	    loc++;
	}

	
	AlgorithmIncrements *thePlotter =  new AlgorithmIncrements(theAlgorithm, 
								   argv[2], xLoc, yLoc, width, height, 
								   displayRecord, fileName);
	(*theRecorder) = thePlotter;

    }

    else if (strcmp(argv[1],"GSA") == 0) {
	if (argc < 3) {
	  opserr << argc;
	  opserr << "WARNING recorder GSA -file fileName? -dT deltaT? - not enough arguments\n";
	  return TCL_ERROR;
	}    
	char *fileName = 0;
	double dT = 0.0;
	int loc = 2;

	while (loc < argc) {
	  if ((strcmp(argv[loc],"-file") == 0) ||
	      (strcmp(argv[loc],"-file") == 0)) {
	    fileName = argv[loc+1];
	    loc += 2;
	  }
	  else if (strcmp(argv[loc],"-dT") == 0) {
	    if (Tcl_GetDouble(interp, argv[loc+1], &dT) != TCL_OK)	
	      return TCL_ERROR;	      
	    loc += 2;
	  }
	  else
	    loc++;
	}

	GSA_Recorder *theR = new GSA_Recorder(theDomain, fileName, dT);
	(*theRecorder) = theR;
    }

    // no recorder type specified yet exists
    else {
	opserr << "WARNING No recorder type exists ";
	opserr << "for recorder of type:" << argv[1];
    
	return TCL_ERROR;
    }    

    // check we instantiated a recorder .. if not ran out of memory
    if ((*theRecorder) == 0) {
	opserr << "WARNING ran out of memory - recorder " << argv[1]<< endln;
	return TCL_ERROR;
    } 

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
		opserr << "WARNING could not add to domain - recorder " << argv[1]<< endln;
		if (theRecorder == 0) 
			opserr << "could not create recorder\n";
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
			opserr << "WARNING could not add to algorithm - recorder " << argv[1]<< endln;
			delete theRecorder;
			return TCL_ERROR;
		} 
		return TCL_OK;
	}
}

