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
                                                                        
// $Revision: 1.10 $
// $Date: 2001-08-30 22:13:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclModelBuilder.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/99
//
// Description: This file contains the class definition for TclModelBuilder.
// A TclModelBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3 
// framework. currently these elements include:
//	1) linear-elastic 2 and 3d beam-column elements
//	2) non-linear material truss
//	3) non-linear 2 and 3d fiber-beam-column elements

//
// What: "@(#) TclModelBuilder.cpp, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <NodalLoad.h>
#include <LoadPattern.h>

#include <SectionForceDeformation.h>
#include <SectionRepres.h>

#include <CrdTransf2d.h>
#include <CrdTransf3d.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <TclModelBuilder.h>
#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>
#include <MultiSupportPattern.h>

#include <Block2D.h>
#include <Block3D.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static Domain *theTclDomain =0;
static TclModelBuilder *theTclBuilder =0;
extern LoadPattern *theTclLoadPattern;
extern MultiSupportPattern *theTclMultiSupportPattern;
static int eleArgStart = 0;
static int nodeLoadTag = 0;
// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclModelBuilder_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
			char **argv);

int
TclModelBuilder_addElement(ClientData clientData, Tcl_Interp *interp,  int argc, 
			   char **argv);

int
TclModelBuilder_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
				    char **argv);

int
TclModelBuilder_addNDMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
			    char **argv);

int
TclModelBuilder_addSection(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv);
			    
int
TclModelBuilder_addPattern(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv);

int
TclModelBuilder_addSeries(ClientData clientData, Tcl_Interp *interp, int argc,   
			  char **argv);

int
TclModelBuilder_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, int argc,
				 char **argv);
int
TclModelBuilder_addHomogeneousBC_X(ClientData clientData, Tcl_Interp *interp, int argc,
				   char **argv);
int
TclModelBuilder_addHomogeneousBC_Y(ClientData clientData, Tcl_Interp *interp, int argc,
				   char **argv);
int
TclModelBuilder_addHomogeneousBC_Z(ClientData clientData, Tcl_Interp *interp, int argc,
				   char **argv);
int
TclModelBuilder_addEqualDOF_MP (ClientData clientData, Tcl_Interp *interp,
                                int argc, char **argv);

int
TclModelBuilder_addMP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      char **argv);

int
TclModelBuilder_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			     char **argv);
int
TclModelBuilder_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc,   
			     char **argv);
int
TclModelBuilder_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      char **argv);

int
TclModelBuilder_addImposedMotionSP(ClientData clientData, 
				   Tcl_Interp *interp, 
				   int argc,    
				   char **argv);	



int
TclModelBuilder_doBlock2D(ClientData clientData, Tcl_Interp *interp, int argc, 
			  char **argv);

int
TclModelBuilder_doBlock3D(ClientData clientData, Tcl_Interp *interp, int argc, 
			  char **argv);


int
TclModelBuilder_addRemoPatch(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,   
			     char **argv);  

int
TclModelBuilder_addRemoLayer(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,   
			     char **argv);   
			       
int
TclModelBuilder_addRemoFiber(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,    
			     char **argv);   

int
TclModelBuilder_addRemoGeomTransf(ClientData clientData, 
				  Tcl_Interp *interp, 
				  int argc,   
				  char **argv); 

				  

int
TclModelBuilder_addGroundMotion(ClientData clientData, 
				Tcl_Interp *interp, 
				int argc,    
				char **argv);

/// added by ZHY
int
TclModelBuilder_UpdateMaterialStage(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc,
				    char **argv);
			   
// REMO
extern int
TclModelBuilder_addPatch (ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv,
			  TclModelBuilder *theTclBuilder);

			  
extern int
TclModelBuilder_addFiber (ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv,
			  TclModelBuilder *theTclBuilder);
			  

extern int
TclModelBuilder_addReinfLayer (ClientData clientData, Tcl_Interp *interp,
			       int argc, char **argv,
			       TclModelBuilder *theTclBuilder);


extern int
TclModelBuilder_addGeomTransf(ClientData, Tcl_Interp *, int, char **,
			      Domain*, TclModelBuilder *);   
	
					 
					 
//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

// constructor: the constructor will add certain commands to the interpreter
TclModelBuilder::TclModelBuilder(Domain &theDomain, Tcl_Interp *interp, int NDM, int NDF)
  :ModelBuilder(theDomain), ndm(NDM), ndf(NDF), theInterp(interp)
{
  theUniaxialMaterials = new ArrayOfTaggedObjects(32);
  theNDMaterials = new ArrayOfTaggedObjects(32);
  theSections  = new ArrayOfTaggedObjects(32);
  theSectionRepresents = new ArrayOfTaggedObjects(32);  
  the2dGeomTransfs = new ArrayOfTaggedObjects(32);  
  the3dGeomTransfs = new ArrayOfTaggedObjects(32);  

  // call Tcl_CreateCommand for class specific commands
  Tcl_CreateCommand(interp, "node", TclModelBuilder_addNode,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "element", TclModelBuilder_addElement,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "uniaxialMaterial", TclModelBuilder_addUniaxialMaterial,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "nDMaterial", TclModelBuilder_addNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "section", TclModelBuilder_addSection,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "pattern", TclModelBuilder_addPattern,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "load", TclModelBuilder_addNodalLoad,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "mass", TclModelBuilder_addNodalMass,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fix", TclModelBuilder_addHomogeneousBC,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixX", TclModelBuilder_addHomogeneousBC_X,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixY", TclModelBuilder_addHomogeneousBC_Y,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixZ", TclModelBuilder_addHomogeneousBC_Z,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "sp", TclModelBuilder_addSP,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "imposedSupportMotion", 
		    TclModelBuilder_addImposedMotionSP,
		    (ClientData)NULL, NULL);  
  
  Tcl_CreateCommand(interp, "groundMotion", 
		    TclModelBuilder_addGroundMotion,
		    (ClientData)NULL, NULL);    

  Tcl_CreateCommand(interp, "equalDOF", TclModelBuilder_addEqualDOF_MP,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "mp", TclModelBuilder_addMP,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "block2D", TclModelBuilder_doBlock2D,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "block3D", TclModelBuilder_doBlock3D,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "patch", TclModelBuilder_addRemoPatch,
		    (ClientData)NULL, NULL);  

  Tcl_CreateCommand(interp, "layer", TclModelBuilder_addRemoLayer,
		    (ClientData)NULL, NULL);    
  
  Tcl_CreateCommand(interp, "fiber", TclModelBuilder_addRemoFiber,
		    (ClientData)NULL, NULL);    

  Tcl_CreateCommand(interp, "geomTransf", TclModelBuilder_addRemoGeomTransf,
		    (ClientData)NULL, NULL);    

  
  ///new command for elast2plast in Multi-yield plasticity, by ZHY
  Tcl_CreateCommand(interp, "updateMaterialStage", 
		    TclModelBuilder_UpdateMaterialStage,
		    (ClientData)NULL, NULL);
  
  // set the static pointers in this file
  theTclBuilder = this;
  theTclDomain = &theDomain;
  theTclLoadPattern = 0;
  theTclMultiSupportPattern = 0;  

  nodeLoadTag = 0;
  eleArgStart = 0;
}

TclModelBuilder::~TclModelBuilder()
{
  theUniaxialMaterials->clearAll();
  theNDMaterials->clearAll();
  theSections->clearAll(); 
  theSectionRepresents->clearAll();
  the2dGeomTransfs->clearAll();
  the3dGeomTransfs->clearAll();

  // free up memory allocated in the constructor
  delete theUniaxialMaterials;
  delete theNDMaterials;
  delete theSections;
  delete theSectionRepresents;
  delete the2dGeomTransfs;
  delete the3dGeomTransfs;

  // set the pointers to 0 
  theTclDomain =0;
  theTclBuilder =0;
  theTclLoadPattern =0;
  theTclMultiSupportPattern = 0;  
  
  // may possibly invoke Tcl_DeleteCommand() later
  Tcl_DeleteCommand(theInterp, "node");
  Tcl_DeleteCommand(theInterp, "element");
  Tcl_DeleteCommand(theInterp, "uniaxialMaterial");
  Tcl_DeleteCommand(theInterp, "nDMaterial");
  Tcl_DeleteCommand(theInterp, "section");
  Tcl_DeleteCommand(theInterp, "pattern");
  Tcl_DeleteCommand(theInterp, "load");
  Tcl_DeleteCommand(theInterp, "mass");
  Tcl_DeleteCommand(theInterp, "fix");
  Tcl_DeleteCommand(theInterp, "fixX");
  Tcl_DeleteCommand(theInterp, "fixY");
  Tcl_DeleteCommand(theInterp, "fixZ");
  Tcl_DeleteCommand(theInterp, "sp");
  Tcl_DeleteCommand(theInterp, "imposedSupportMotion");
  Tcl_DeleteCommand(theInterp, "groundMotion");
  Tcl_DeleteCommand(theInterp, "equalDOF");
  Tcl_DeleteCommand(theInterp, "mp");
  Tcl_DeleteCommand(theInterp, "block2D");
  Tcl_DeleteCommand(theInterp, "block3D");
  Tcl_DeleteCommand(theInterp, "patch");
  Tcl_DeleteCommand(theInterp, "layer");
  Tcl_DeleteCommand(theInterp, "fiber");
  Tcl_DeleteCommand(theInterp, "geomTransf");
  Tcl_DeleteCommand(theInterp, "updateMaterialStage");
}


//
// CLASS METHODS
//

int 
TclModelBuilder::buildFE_Model(void)
{
  // does nothing
  return 0;
}

int 
TclModelBuilder::getNDM(void) const
{
  return ndm;
}

int 
TclModelBuilder::getNDF(void) const
{
  return ndf;
}

int 
TclModelBuilder::addUniaxialMaterial(UniaxialMaterial &theMaterial)
{
  bool result = theUniaxialMaterials->addComponent(&theMaterial);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addUniaxialMaterial() - failed to add material: " << theMaterial;
    return -1;
  }
}


UniaxialMaterial *
TclModelBuilder::getUniaxialMaterial(int tag)
{
  TaggedObject *mc = theUniaxialMaterials->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  UniaxialMaterial *result = (UniaxialMaterial *)mc;
  return result;
}

int 
TclModelBuilder::addNDMaterial(NDMaterial &theMaterial)
{
  bool result = theNDMaterials->addComponent(&theMaterial);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addNDMaterial() - failed to add material: " << theMaterial;
    return -1;
  }
}


NDMaterial *
TclModelBuilder::getNDMaterial(int tag)
{
  TaggedObject *mc = theNDMaterials->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  NDMaterial *result = (NDMaterial *)mc;
  return result;
}

int 
TclModelBuilder::addSection(SectionForceDeformation &theSection)
{
  bool result = theSections->addComponent(&theSection);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addSection() - failed to add section: " << theSection;
    return -1;
  }
}



SectionForceDeformation *
TclModelBuilder::getSection(int tag)
{
  TaggedObject *mc = theSections->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  SectionForceDeformation *result = (SectionForceDeformation *)mc;
  return result;
}



int 
TclModelBuilder::addSectionRepres(SectionRepres &theSectionRepres)
{
  bool result = theSectionRepresents->addComponent(&theSectionRepres);

  if (result == true)
    return 0;
  else {
      cerr << "TclModelBuilder::addSectionRepres() - failed to add SectionRepres\n";
      return -1;
  }
}


SectionRepres *
TclModelBuilder::getSectionRepres(int tag)
{
  TaggedObject *mc = theSectionRepresents->getComponentPtr(tag);
  if (mc == 0) return 0;
  SectionRepres *result = (SectionRepres *)mc;
  return result;
}



int 
TclModelBuilder::addCrdTransf2d(CrdTransf2d &theCrdTransf)
{
  bool result = the2dGeomTransfs->addComponent(&theCrdTransf);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addCrdTransf() - failed to add crdTransf: " << theCrdTransf;
    return -1;
  }
}


int 
TclModelBuilder::addCrdTransf3d(CrdTransf3d &theCrdTransf)
{
  bool result = the3dGeomTransfs->addComponent(&theCrdTransf);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addCrdTransf() - failed to add crdTransf: " << theCrdTransf;
    return -1;
  }
}



CrdTransf2d *
TclModelBuilder::getCrdTransf2d(int tag)
{
  TaggedObject *mc = the2dGeomTransfs->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  CrdTransf2d *result = (CrdTransf2d *)mc;
  return result;
}


CrdTransf3d *
TclModelBuilder::getCrdTransf3d(int tag)
{
  TaggedObject *mc = the3dGeomTransfs->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  CrdTransf3d *result = (CrdTransf3d *)mc;
  return result;
}

//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

void printCommand(int argc, char **argv)
{
  cerr << "Input command: ";
  for (int i=0; i<argc; i++)
    cerr << argv[i] << " ";
  cerr << endl;
} 

int
TclModelBuilder_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
                        char **argv)
{

  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed" << endl;
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  // make sure corect number of arguments on command line
  if (argc < 2+ndm) {
    cerr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    cerr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }    

  Node *theNode = 0;

  // get the nodal id
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeTag\n";
    cerr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  // read in the coordinates and create the node
  double xLoc, yLoc, zLoc;
  if (ndm == 1) { 
    // create a node in 1d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      cerr << "WARNING invalid XCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc);
  } 

  else if (ndm == 2) { 
    // create a node in 2d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      cerr << "WARNING invalid XCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      cerr << "WARNING invalid YCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc,yLoc);
  } 

  else if (ndm == 3) { 
    // create a node in 3d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      cerr << "WARNING invalid XCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      cerr << "WARNING invalid YCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &zLoc) != TCL_OK) {
      cerr << "WARNING invalid ZCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc,yLoc,zLoc);
  } else {
      cerr << "WARNING invalid ndm\n";
      cerr << "node: " << nodeId << endl;;
      return TCL_ERROR;
  }

  if (theNode == 0) {
    cerr << "WARNING ran out of memory creating node\n";
    cerr << "node: " << nodeId << endl;
    return TCL_ERROR;
  }

  if (theTclDomain->addNode(theNode) == false) {
    cerr << "WARNING failed to add node to the domain\n";
    cerr << "node: " << nodeId << endl;
    delete theNode; // otherwise memory leak
    return TCL_ERROR;
  }

  // check for mass terms
  if (argc > 2+ndm) {
    if (strcmp(argv[2+ndm],"-mass") == 0) {
      if (argc < 3+ndm+ndf) {
        cerr << "WARNING incorrect number of nodal mass terms\n";
        cerr << "node: " << nodeId << endl;
        return TCL_ERROR;      
      }	
      Matrix mass(ndf,ndf);
      double theMass;
      for (int i=0; i<ndf; i++) {
	if (Tcl_GetDouble(interp, argv[i+3+ndm], &theMass) != TCL_OK) {
	  cerr << "WARNING invalid nodal mass term\n";
	  cerr << "node: " << nodeId << ", dof: " << i+1 << endl;
	  return TCL_ERROR;
	}
	mass(i,i) = theMass;
      }
      theNode->setMass(mass);      
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}





// the function for creating ne material objects and patterns is in a seperate file.
// this allows new material and patternobjects to be added without touching this file.
// does so at the expense of an extra procedure call.


extern int 
TclModelBuilderElementCommand(ClientData clientData, 
			      Tcl_Interp *interp, int argc,    
			      char **argv, 
			      Domain *theDomain, TclModelBuilder *theTclBuilder);
int
TclModelBuilder_addElement(ClientData clientData, Tcl_Interp *interp, 
			   int argc,    char **argv)
                          
{
  return TclModelBuilderElementCommand(clientData, interp, 
				       argc, argv, theTclDomain, theTclBuilder);
}


extern int
TclModelBuilderUniaxialMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, 
				    int argc, char **argv)
                          
{
  return TclModelBuilderUniaxialMaterialCommand(clientData, interp, 
						argc, argv, theTclBuilder);
}

extern int
TclModelBuilderNDMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				  char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addNDMaterial(ClientData clientData, Tcl_Interp *interp, 
			    int argc,    char **argv)
                          
{
  return TclModelBuilderNDMaterialCommand(clientData, interp, 
					  argc, argv, theTclBuilder);
}

extern int
TclModelBuilderSectionCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				  char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addSection(ClientData clientData, Tcl_Interp *interp, 
			    int argc,    char **argv)
                          
{
  return TclModelBuilderSectionCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}



extern int
TclPatternCommand(ClientData clientData, Tcl_Interp *interp, 
			   int argc, char **argv, Domain *theDomain);
			   
int
TclModelBuilder_addPattern(ClientData clientData, Tcl_Interp *interp, 
			   int argc, char **argv)
			  
{
  return TclPatternCommand(clientData, interp, argc, argv, theTclDomain);
}




extern int
TclGroundMotionCommand(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc,    
		       char **argv,
		       MultiSupportPattern *thePattern);

int
TclModelBuilder_addGroundMotion(ClientData clientData, Tcl_Interp *interp, 
			   int argc, char **argv)
			  
{
  return TclGroundMotionCommand(clientData, interp, argc, argv, 
				theTclMultiSupportPattern);
}


int
TclModelBuilder_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			 char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - load \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  NodalLoad *theLoad = 0;
  
  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: load nodeId " << ndf << " forces\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeId: " << argv[1];
    cerr << " - load nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // get the load vector
  Vector forces(ndf);
  for (int i=0; i<ndf; i++) {
    double theForce;
    if (Tcl_GetDouble(interp, argv[2+i], &theForce) != TCL_OK) {
      cerr << "WARNING invalid force " << i+1 << " - load " << nodeId;
      cerr << " " << ndf << " forces\n";
      return TCL_ERROR;
    } else
      forces(i) = theForce;
  }

  bool isLoadConst = false;
  int loadPatternTag = 123456789; // some pattern that will never be used!

  // allow some additional options at end of command
  int endMarker = 2+ndf;
  while (endMarker != argc) {
    if (strcmp(argv[endMarker],"-const") == 0) {
      // allow user to specify const load
      isLoadConst = true;
    } else if (strcmp(argv[endMarker],"-pattern") == 0) {
      // allow user to specify load pattern other than current
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

	cerr << "WARNING invalid patternTag - load " << nodeId << " ";
	cerr << ndf << " forces pattern patterntag\n";
	return TCL_ERROR;
      }
    }
    endMarker++;
  }

  // get the current pattern tag if no tag given in i/p
  if (loadPatternTag == 123456789)
    if (theTclLoadPattern == 0) {
	cerr << "WARNING no current load pattern - load " << nodeId;
	cerr << " " << ndf << " forces\n";
	return TCL_ERROR;
    } else 
	loadPatternTag = theTclLoadPattern->getTag();

  // create the load
  theLoad = new NodalLoad(nodeLoadTag, nodeId, forces, isLoadConst);
  if (theLoad == 0) {
    cerr << "WARNING ran out of memory for load  - load " << nodeId;
    cerr << " " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // add the load to the domain
  if (theTclDomain->addNodalLoad(theLoad, loadPatternTag) == false) {
    cerr << "WARNING TclModelBuilder - could not add load to domain ";
    printCommand(argc, argv);
    delete theLoad;
    return TCL_ERROR;
  }
  nodeLoadTag++;

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclModelBuilder_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc, 
                        char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - load \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: mass nodeId " << ndf << " mass values\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeId: " << argv[1];
    cerr << " - mass nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  Node *theNode = 0;
  theNode = theTclDomain->getNode(nodeId);

  if (theNode == 0)
  {
    cerr << "WARNING failed to get node pointer from the domain\n";
    cerr << "node: " << nodeId << endl;
    return TCL_ERROR;
  }

  // check for mass terms
  Matrix mass(ndf,ndf);
  double theMass;
  for (int i=0; i<ndf; i++) 
  {
     if (Tcl_GetDouble(interp, argv[i+2], &theMass) != TCL_OK) 
     {
	  cerr << "WARNING invalid nodal mass term\n";
	  cerr << "node: " << nodeId << ", dof: " << i+1 << endl;
	  return TCL_ERROR;
      }
      mass(i,i) = theMass;
  }
      
  theNode->setMass(mass);      

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}





int
TclModelBuilder_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, int argc,   
				 char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();
  int numSPs = theTclDomain->getNumSPs();

  // check number of arguments
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: fix nodeId " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      cerr << "WARNING invalid nodeId - fix nodeId " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // get the fixity condition and add the constraint if fixed
  for (int i=0; i<ndf; i++) {
    int theFixity;
    if (Tcl_GetInt(interp, argv[2+i], &theFixity) != TCL_OK) {
      cerr << "WARNING invalid fixity " << i+1 << " - load " << nodeId;
      cerr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } else {
      if (theFixity != 0) {
	// create a homogeneous constraint
	SP_Constraint *theSP = new SP_Constraint(numSPs, nodeId, i, 0.0);
	if (theSP == 0) {
	  cerr << "WARNING ran out of memory for SP_Constraint ";
	  cerr << "fix " << nodeId << " " << ndf << " [0,1] conditions\n";
	  return TCL_ERROR;
	}
	if (theTclDomain->addSP_Constraint(theSP) == false) {
	  cerr << "WARNING could not add SP_Constraint to domain - fix";
	  cerr << nodeId << " " << ndf << " [0,1] conditions\n";
	  delete theSP;
	  return TCL_ERROR;
	}
	numSPs++;      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}


int
TclModelBuilder_addHomogeneousBC_X(ClientData clientData, Tcl_Interp *interp, 
				   int argc, char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();
  int numSPs = theTclDomain->getNumSPs();

  // check number of arguments
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: fixX xLoc " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the xCrd of nodes to be constrained
  double xLoc;
  if (Tcl_GetDouble(interp, argv[1], &xLoc) != TCL_OK) {
      cerr << "WARNING invalid xCrd - fixX xLoc " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; i++) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      cerr << "WARNING invalid fixity " << i+1 << " - fixX " << xLoc;
      cerr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } 
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
    if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
      cerr << "WARNING invalid tol specified - fixX " << xLoc << endl;
      return TCL_ERROR;
    }       
  }

  NodeIter &theNodes = theTclDomain->getNodes();
  Node *theNode;

  // loop over all the nodes
  while ((theNode = theNodes()) != 0) {
    const Vector &theCrd = theNode->getCrds();
    double nodeX = theCrd(0);

    // add a single point constraint if Xcrd of node is within tol of xLoc
    if (fabs(nodeX - xLoc) < tol) {
      int nodeId = theNode->getTag();
      int theFixity = 0;

      // loop over all the ndf values valid for the node
      int numDOF = theNode->getNumberDOF();
      if (numDOF  < ndf) numDOF = ndf;

      for (int i=0; i<numDOF; i++) {
	theFixity = fixity(i);
	if (theFixity != 0) {
	  // create a homogeneous constraint
	  SP_Constraint *theSP = new SP_Constraint(numSPs, nodeId, i, 0.0);
	  if (theSP == 0) {
	    cerr << "WARNING ran out of memory for SP_Constraint at node " << nodeId;
	    cerr << " - fixX " << xLoc << " " << ndf << " [0,1] conditions\n";
	    return TCL_ERROR;
	  }
	  if (theTclDomain->addSP_Constraint(theSP) == false) {
	    cerr << "WARNING could not add SP_Constraint to domain for node " << nodeId;
	    cerr << " - fixX " << xLoc << " " << ndf << " [0,1] conditions\n";
	    delete theSP;
	    return TCL_ERROR;
	  }
	  numSPs++;      
	}
      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}




int
TclModelBuilder_addHomogeneousBC_Y(ClientData clientData, Tcl_Interp *interp, 
				   int argc, char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();
  int numSPs = theTclDomain->getNumSPs();

  // check number of arguments
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: fixY yLoc " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the yCrd of nodes to be constrained
  double yLoc;
  if (Tcl_GetDouble(interp, argv[1], &yLoc) != TCL_OK) {
      cerr << "WARNING invalid yCrd - fixY yLoc " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; i++) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      cerr << "WARNING invalid fixity " << i+1 << " - fixY " << yLoc;
      cerr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } 
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
    if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
      cerr << "WARNING invalid tol specified - fixY " << yLoc << endl;
      return TCL_ERROR;
    }       
  }

  NodeIter &theNodes = theTclDomain->getNodes();
  Node *theNode;

  // loop over all the nodes
  while ((theNode = theNodes()) != 0) {
    const Vector &theCrd = theNode->getCrds();
    if (theCrd.Size() > 1) {
      double nodeY = theCrd(1);

      // add a single point constraint if Xcrd of node is within tol of yLoc
      if (fabs(nodeY - yLoc) < tol) {

	int nodeId = theNode->getTag();
	int theFixity = 0;

	// loop over all the ndf values valid for the node
	int numDOF = theNode->getNumberDOF();
	if (numDOF  < ndf) numDOF = ndf;

	for (int i=0; i<numDOF; i++) {
	  theFixity = fixity(i);
	  if (theFixity != 0) {
	    // create a homogeneous constraint
	    SP_Constraint *theSP = new SP_Constraint(numSPs, nodeId, i, 0.0);
	    if (theSP == 0) {
	      cerr << "WARNING ran out of memory for SP_Constraint at node " << nodeId;
	      cerr << " - fixY " << yLoc << " " << ndf << " [0,1] conditions\n";
	      return TCL_ERROR;
	    }
	    if (theTclDomain->addSP_Constraint(theSP) == false) {
	      cerr << "WARNING could not add SP_Constraint to domain for node " << nodeId;
	      cerr << " - fixY " << yLoc << " " << ndf << " [0,1] conditions\n";
	      delete theSP;
	      return TCL_ERROR;
	    }
	    numSPs++;      
	  }
	}
      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclModelBuilder_addHomogeneousBC_Z(ClientData clientData, Tcl_Interp *interp, 
				   int argc, char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();
  int numSPs = theTclDomain->getNumSPs();

  // check number of arguments
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: fixZ zLoc " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the yCrd of nodes to be constrained
  double zLoc;
  if (Tcl_GetDouble(interp, argv[1], &zLoc) != TCL_OK) {
      cerr << "WARNING invalid zCrd - fixZ zLoc " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; i++) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      cerr << "WARNING invalid fixity " << i+1 << " - fixZ " << zLoc;
      cerr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } 
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
    if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
      cerr << "WARNING invalid tol specified - fixZ " << zLoc << endl;
      return TCL_ERROR;
    }       
  }

  NodeIter &theNodes = theTclDomain->getNodes();
  Node *theNode;

  // loop over all the nodes
  while ((theNode = theNodes()) != 0) {
    const Vector &theCrd = theNode->getCrds();
    if (theCrd.Size() > 2) {
      double nodeZ = theCrd(2);

      // add a single point constraint if Xcrd of node is within tol of zLoc
      if (fabs(nodeZ - zLoc) < tol) {

	int nodeId = theNode->getTag();
	int theFixity = 0;

	// loop over all the ndf values valid for the node
	int numDOF = theNode->getNumberDOF();
	if (numDOF  < ndf) numDOF = ndf;

	for (int i=0; i<numDOF; i++) {
	  theFixity = fixity(i);
	  if (theFixity != 0) {
	    // create a homogeneous constraint
	    SP_Constraint *theSP = new SP_Constraint(numSPs, nodeId, i, 0.0);
	    if (theSP == 0) {
	      cerr << "WARNING ran out of memory for SP_Constraint at node " << nodeId;
	      cerr << " - fixZ " << zLoc << " " << ndf << " [0,1] conditions\n";
	      return TCL_ERROR;
	    }
	    if (theTclDomain->addSP_Constraint(theSP) == false) {
	      cerr << "WARNING could not add SP_Constraint to domain for node " << nodeId;
	      cerr << " - fixZ " << zLoc << " " << ndf << " [0,1] conditions\n";
	      delete theSP;
	      return TCL_ERROR;
	    }
	    numSPs++;      
	  }
	}
      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}





int
TclModelBuilder_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - sp \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < 4) {
    cerr << "WARNING bad command - want: sp nodeId dofID value";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId;
  double value;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeId: " << argv[1] << " -  sp nodeId dofID value\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    cerr << "WARNING invalid dofId: " << argv[2] << " -  sp ";
    cerr << nodeId << " dofID value\n";
      return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    cerr << "WARNING invalid value: " << argv[3] << " -  sp ";
    cerr << nodeId << " dofID value\n";
      return TCL_ERROR;
  }

  bool isSpConst = false;
  int loadPatternTag = 123456789; // some pattern that will never be used!

  // allow some additional options at end of command
  theTclLoadPattern->getTag();
  int endMarker = 4;
  while (endMarker != argc) {
    if (strcmp(argv[endMarker],"-const") == 0) {
      // allow user to specify const load
      isSpConst = true;
    } else if (strcmp(argv[endMarker],"-pattern") == 0) {
      // allow user to specify load pattern other than current
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

	cerr << "WARNING invalid patternTag - load " << nodeId << " ";
	cerr << ndf << " forces pattern patterntag\n";
	return TCL_ERROR;
      }
    }  
    endMarker++;
  }

  // if load pattern tag has not changed - get the pattern tag from current one
  if (loadPatternTag == 123456789) {
    if (theTclLoadPattern == 0) {
      cerr << "WARNING no current pattern - sp " << nodeId << " dofID value\n";
      return TCL_ERROR;
    } else	
      loadPatternTag = theTclLoadPattern->getTag();
  }
  
  LoadPattern *thePattern = theTclDomain->getLoadPattern(loadPatternTag);
  SP_ConstraintIter &theSPs = thePattern->getSPs();
  int numSPs = 0;
  SP_Constraint *theSP2;
  while ((theSP2 = theSPs()) != 0)
      numSPs++;
  
  
  // create a homogeneous constraint
  SP_Constraint *theSP = new SP_Constraint(numSPs, nodeId, dofId, value, isSpConst);

  if (theSP == 0) {
    cerr << "WARNING ran out of memory for SP_Constraint ";
    cerr << " - sp " << nodeId << " dofID value\n";
    return TCL_ERROR;
  }
  if (theTclDomain->addSP_Constraint(theSP, loadPatternTag) == false) {
    cerr << "WARNING could not add SP_Constraint to domain ";
    printCommand(argc, argv);
    delete theSP;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclModelBuilder_addImposedMotionSP(ClientData clientData, 
				   Tcl_Interp *interp, 
				   int argc,   
				   char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - sp \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < 4) {
    cerr << "WARNING bad command - want: imposedSupportMotion nodeId dofID gMotionID\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId, gMotionID;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeId: " << argv[1];
    cerr << " - imposedSupportMotion nodeId dofID gMotionID\n";    
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    cerr << "WARNING invalid dofId: " << argv[2] << " -  imposedSupportMotion ";
    cerr << nodeId << " dofID gMotionID\n";    
      return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetInt(interp, argv[3], &gMotionID) != TCL_OK) {
    cerr << "WARNING invalid gMotionID: " << argv[3] << " -  imposedSupportMotion ";
    cerr << nodeId << " dofID gMotionID\n";
    return TCL_ERROR;
  }

  bool alt = false;
  if (argc == 5) {
    if (strcmp(argv[4],"-other") == 0) 
      alt = true;
  }

  MultiSupportPattern *thePattern = theTclMultiSupportPattern;
  int loadPatternTag = thePattern->getTag();
  
  GroundMotion *theGMotion = thePattern->getMotion(gMotionID);
  if (theGMotion == 0) {
    cerr << "WARNING no GroundMotion with tag: " << argv[3];
    cerr << " in current MultipleSupportPattern\n";
    return TCL_ERROR;      
  }
  
  SP_ConstraintIter &theSPs = thePattern->getSPs();
  int numSPs = 0;
  SP_Constraint *theSP2;
  while ((theSP2 = theSPs()) != 0)
      numSPs++;
  
  // create a new ImposedMotionSP
  SP_Constraint *theSP;
  if (alt == true) {
    theSP = new ImposedMotionSP1(numSPs, nodeId, dofId, 
				*theGMotion, false);
  }
  else {
    theSP = new ImposedMotionSP(numSPs, nodeId, dofId, 
				*theGMotion, false);
  }
  if (theSP == 0) {
    cerr << "WARNING ran out of memory for ImposedMotionSP ";
    cerr << " -  imposedSupportMotion ";
    cerr << nodeId << " " << dofId++ << " " << gMotionID << endl;
    return TCL_ERROR;
  }
  if (thePattern->addSP_Constraint(theSP) == false) {
    cerr << "WARNING could not add SP_Constraint to pattern ";
    printCommand(argc, argv);
    delete theSP;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}






int
TclModelBuilder_addEqualDOF_MP (ClientData clientData, Tcl_Interp *interp,
                                int argc, char **argv)
{
        // Ensure the destructor has not been called
        if (theTclBuilder == 0) {
	  cerr << "WARNING builder has been destroyed - equalDOF \n";
	  return TCL_ERROR;
        }

        // Check number of arguments
        if (argc < 4) {
	  cerr << "WARNING bad command - want: equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	  printCommand (argc, argv);
	  return TCL_ERROR;
        }

        // Read in the node IDs and the DOF
        int RnodeID, CnodeID, dofID;

        if (Tcl_GetInt (interp, argv[1], &RnodeID) != TCL_OK) {
	  cerr << "WARNING invalid RnodeID: " << argv[1]
	       << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	  return TCL_ERROR;
        }
        if (Tcl_GetInt (interp, argv[2], &CnodeID) != TCL_OK) {
	  cerr << "WARNING invalid CnodeID: " << argv[2]
	       << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	  return TCL_ERROR;
        }

        // The number of DOF to be coupled
        int numDOF = argc - 3;

        // The constraint matrix ... U_c = C_cr * U_r
        Matrix Ccr (numDOF, numDOF);
        Ccr.Zero();

        // The vector containing the retained and constrained DOFs
        ID rcDOF (numDOF);

        int i, j;
        // Read the degrees of freedom which are to be coupled
        for (i = 3, j = 0; i < argc; i++, j++) {
	  if (Tcl_GetInt (interp, argv[i], &dofID) != TCL_OK) {
	    cerr << "WARNING invalid dofID: " << argv[3]
		 << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	    return TCL_ERROR;
	  }

	  rcDOF (j) = --dofID;    // Decrement for C++ indexing
	  Ccr (j,j) = 1.0;
        }

        // Use this for the MP tag
        int numMPs = theTclDomain->getNumMPs();

        // Create the multi-point constraint
        MP_Constraint *theMP = new MP_Constraint (numMPs, RnodeID, CnodeID, Ccr, rcDOF, rcDOF);
        if (theMP == 0) {
	  cerr << "WARNING ran out of memory for equalDOF MP_Constraint ";
	  printCommand (argc, argv);
	  return TCL_ERROR;
        }

        // Add the multi-point constraint to the domain
        if (theTclDomain->addMP_Constraint (theMP) == false) {
	  cerr << "WARNING could not add equalDOF MP_Constraint to domain ";
	  printCommand(argc, argv);
	  delete theMP;
	  return TCL_ERROR;
        }

        return TCL_OK;
}



int
TclModelBuilder_addMP(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  cerr << "WARNING - TclModelBuilder_addMP() not yet implemented\n";
  return TCL_OK;
}

int
TclModelBuilder_doBlock2D(ClientData clientData, Tcl_Interp *interp, int argc,   
			  char **argv)
{

  int ndm = theTclBuilder->getNDM();
  if (ndm < 2) {
    cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    cerr << " : model dimension (ndm) must be at leat 2 " << endl;
    return TCL_ERROR;
  }

  if (argc < 8) {
    cerr << "WARNING incorrect numer of args :block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    return TCL_ERROR;
  }
  int numX, numY, startNodeNum, startEleNum;
  if (Tcl_GetInt (interp, argv[1], &numX) != TCL_OK) {
    cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid numX: " << argv[1] << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[2], &numY) != TCL_OK) {
    cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid numY: " << argv[2] << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[3], &startNodeNum) != TCL_OK) {
    cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid startNode: " << argv[3] << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[4], &startEleNum) != TCL_OK) {
    cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid startEle: " << argv[4] << endl;
    return TCL_ERROR;
  }


  static Matrix Coordinates(9,3);
  static ID     haveNode(9);
  Coordinates.Zero();
  for (int k=0; k<9; k++) haveNode(k) = -1;

  int numNodes = 4;
  if (argc == 10) {
    if (strcmp(argv[7],"-numEleNodes") == 0) 
      if (Tcl_GetInt (interp, argv[8], &numNodes) != TCL_OK) {
	cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
	cerr << " -numEleNodes numNodes?: invalid numNodes: " << argv[8] << endl;
	return TCL_ERROR;
      }
    if (numNodes != 4 && numNodes != 9) {
      cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
      cerr << "-numEleNodes numNodes?: invalid numNodes: " << argv[8] << " 4 or 9 only\n";
      return TCL_ERROR;
    }

    if (numNodes == 9) {
      if (((numX % 2) != 0) || ((numY % 2) != 0)) {
	cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
	cerr << "-numEleNodes 9: numX and numY MUST BOTH BE EVEN\n";
	return TCL_ERROR;
      }
    }
  }


  char *nodalInfo;
  if (numNodes == 4)
    nodalInfo = argv[7];
  else
    nodalInfo = argv[9];

  char **argvNodes;
  int  argcNodes;
  
  Tcl_SplitList(interp, nodalInfo, &argcNodes, &argvNodes);

  int ndf = theTclBuilder->getNDF();
  
  int count = 0;
  while (count < argcNodes) {
    if ((count + ndm + 1) >  argcNodes) {
      cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      cerr << " : invalid number of node args: " << argv[7] << endl;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR; 
    }
    int nodeTag;
    double value;
    if (Tcl_GetInt (interp, argvNodes[count], &nodeTag) != TCL_OK) {
      cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      cerr << " : invalid node tag: " << argvNodes[count] << endl;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR; 
    }
    if (nodeTag < 1 || nodeTag > 9) {
      cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      cerr << " : invalid node tag out of bounds [1,9]: " << argvNodes[count] << endl;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    for (int i=0; i<ndm; i++) {
      if (Tcl_GetDouble(interp, argvNodes[count+1+i], &value) != TCL_OK) {
	cerr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
	cerr << " : invalid node coordinate for node: " << argvNodes[count] << endl;
	Tcl_Free((char *)argvNodes);
	return TCL_ERROR;
      }
      Coordinates(nodeTag-1,i) = value;
      haveNode(nodeTag-1) = nodeTag;
    }      
    count += 1 + ndm;
  }

  Tcl_Free((char *)argvNodes);

  Block2D  theBlock(numX, numY, haveNode, Coordinates, numNodes);

  // create the nodes: (numX+1)*(numY+1) nodes to be created
  int nodeID = startNodeNum;
  int jj;
  for (jj=0; jj<=numY; jj++) {
    for (int ii=0; ii<=numX; ii++) {
      const Vector &nodeCoords = theBlock.getNodalCoords(ii,jj);
      double xLoc = nodeCoords(0);
      double yLoc = nodeCoords(1);
      double zLoc = nodeCoords(2);
      Node *theNode = 0;
      if (ndm == 2) {
	theNode = new Node(nodeID,ndf,xLoc, yLoc);
      } else if (ndm == 3) {
	theNode = new Node(nodeID,ndf,xLoc, yLoc, zLoc);
      } 

      if (theNode == 0) {
	cerr << "WARNING ran out of memory creating node\n";
	cerr << "node: " << nodeID << endl;
	return TCL_ERROR;
      }

      if (theTclDomain->addNode(theNode) == false) {
	cerr << "WARNING failed to add node to the domain\n";
	cerr << "node: " << nodeID << endl;
	delete theNode; // otherwise memory leak
	return TCL_ERROR;
      }

      nodeID++;
    }
  }
    
  // create the elements: numX*numY elements to be created if 4 node elements
  //                      numX/2 * numY /2 nodes to be v=created if 9 node elements
  char *eleType = argv[5];
  char *additionalEleArgs = argv[6];
  //  const ID &nodeTags = theBlock.getElementNodes(0,0);  
  //  int numNodes = nodeTags.Size();

  // assumes 15 is largest string for individual nodeTags
  count = 10 + strlen(eleType) + strlen(additionalEleArgs) + 15 * (numNodes+1);
  char *eleCommand = new char[count];
  int initialCount = 8 + strlen(eleType);

  int  eleID = startEleNum; 
  if (numNodes == 9) {
    numX /= 2;
    numY /= 2;
  }
    

  for (jj=0; jj<numY; jj++) {
    for (int ii=0; ii<numX; ii++) {
      count = initialCount;

      const ID &nodeTags = theBlock.getElementNodes(ii,jj);
      
      // create the string to be evaluated
      strcpy(eleCommand, "element ");
      strcpy(&eleCommand[8], eleType);
      count += sprintf(&eleCommand[count], " %d ", eleID);
      for (int i=0; i<numNodes; i++) {
	int nodeTag = nodeTags(i)+startNodeNum;
	count += sprintf(&eleCommand[count], " %d ", nodeTag);
      }
      strcat(eleCommand, additionalEleArgs);      

      // now to create the element we get the string eveluated
      if (Tcl_Eval(interp, eleCommand) != TCL_OK) {
          delete [] eleCommand;
	return TCL_ERROR;
      }
      eleID++;
    }
  }

  delete [] eleCommand;
  return TCL_OK;
}


int
TclModelBuilder_doBlock3D(ClientData clientData, Tcl_Interp *interp, int argc,   
			  char **argv)
{

  int ndm = theTclBuilder->getNDM();
  if (ndm < 3) {
    cerr << "WARNING block3D numX? numY? startNode? startEle? eleType? eleArgs?";
    cerr << " : model dimension (ndm) must be at leat 2 " << endl;
    return TCL_ERROR;
  }

  int numX, numY, numZ, startNodeNum, startEleNum;
  if (Tcl_GetInt (interp, argv[1], &numX) != TCL_OK) {
    cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid numX: " << argv[1] << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[2], &numY) != TCL_OK) {
    cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid numY: " << argv[2] << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[3], &numZ) != TCL_OK) {
    cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid numZ: " << argv[3] << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[4], &startNodeNum) != TCL_OK) {
    cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid startNode: " << argv[4] << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[5], &startEleNum) != TCL_OK) {
    cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    cerr << " : invalid startEle: " << argv[5] << endl;
    return TCL_ERROR;
  }

  static Matrix Coordinates(27,3);
  static ID     haveNode(27);
  Coordinates.Zero();
  for (int k=0; k<27; k++) haveNode(k) = -1;

  char *nodalInfo = argv[8];
  char **argvNodes;
  int  argcNodes;
  
  Tcl_SplitList(interp, nodalInfo, &argcNodes, &argvNodes);

  int ndf = theTclBuilder->getNDF();
  
  int count = 0;
  while (count < argcNodes) {
    if ((count + ndm + 1) > argcNodes) {
      cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      cerr << " : invalid number of node args: " << argv[8] << endl;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR; 
    }
    int nodeTag;
    double value;
    if (Tcl_GetInt (interp, argvNodes[count], &nodeTag) != TCL_OK) {
      cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      cerr << " : invalid node id in node args: " << argvNodes[count] << endl;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR; 
    }
    if (nodeTag < 1 || nodeTag > 27) {
      cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      cerr << " : node tag out of bounds [1, 27]: " << argvNodes[count] << endl;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    for (int i=0; i<ndm; i++) {
      if (Tcl_GetDouble(interp, argvNodes[count+1+i], &value) != TCL_OK) {
	cerr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
	cerr << " : invalid coordinate in node args: " << argvNodes[count] << endl;
	Tcl_Free((char *)argvNodes);
	return TCL_ERROR;
      }
      Coordinates(nodeTag-1,i) = value;
      haveNode(nodeTag-1) = nodeTag;
    }      
    count += 1 + ndm;
  }

  Tcl_Free((char *)argvNodes);

  Block3D  theBlock(numX, numY, numZ, haveNode, Coordinates);

  // create the nodes: (numX+1)*(numY+1) nodes to be created
  int nodeID = startNodeNum;
  int kk;
  for (kk=0; kk<=numZ; kk++) {
    for (int jj=0; jj<=numY; jj++) {
      for (int ii=0; ii<=numX; ii++) {
	const Vector &nodeCoords = theBlock.getNodalCoords(ii,jj,kk);
	double xLoc = nodeCoords(0);
	double yLoc = nodeCoords(1);
	double zLoc = nodeCoords(2);
	Node *theNode = 0;
	theNode = new Node(nodeID,ndf,xLoc, yLoc, zLoc);
	
	if (theNode == 0) {
	  cerr << "WARNING ran out of memory creating node\n";
	  cerr << "node: " << nodeID << endl;
	  return TCL_ERROR;
	}

	if (theTclDomain->addNode(theNode) == false) {
	  cerr << "WARNING failed to add node to the domain\n";
	  cerr << "node: " << nodeID << endl;
	  delete theNode; // otherwise memory leak
	  return TCL_ERROR;
	}
	
	nodeID++;
      }
    }
  }
    
  // create the elements: numX*numY elements to be created
  char *eleType = argv[6];
  char *additionalEleArgs = argv[7];
  const ID &nodeTags = theBlock.getElementNodes(0,0,0);  
  int numNodes = nodeTags.Size();

  // assumes 15 is largest string for individual nodeTags
  count = 10 + strlen(eleType) + strlen(additionalEleArgs) + 15 * (numNodes+1);
  char *eleCommand = new char[count];
  int initialCount = 8 + strlen(eleType);

  int  eleID = startEleNum;  
  for (kk=0; kk<numZ; kk++) {
    for (int jj=0; jj<numY; jj++) {
      for (int ii=0; ii<numX; ii++) {
	count = initialCount;

	const ID &nodeTags = theBlock.getElementNodes(ii,jj,kk);
      
	// create the string to be evaluated
	strcpy(eleCommand, "element ");
	strcpy(&eleCommand[8], eleType);
	count += sprintf(&eleCommand[count], " %d ", eleID);
	for (int i=0; i<numNodes; i++) {
	  int nodeTag = nodeTags(i)+startNodeNum;
	  count += sprintf(&eleCommand[count], " %d ", nodeTag);
	}
	strcat(eleCommand, additionalEleArgs);      
	
	// now to create the element we get the string eveluated
	if (Tcl_Eval(interp, eleCommand) != TCL_OK) {
        delete [] eleCommand;
	  return TCL_ERROR;
	}
	eleID++;
      }
    }
  }

  delete [] eleCommand;
  return TCL_OK;
}




int
TclModelBuilder_addRemoPatch(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  return TclModelBuilder_addPatch(clientData, interp, argc,argv,
				    theTclBuilder);
}

int
TclModelBuilder_addRemoFiber(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  return TclModelBuilder_addFiber(clientData, interp, argc,argv,
				  theTclBuilder);
}

int
TclModelBuilder_addRemoLayer(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  return TclModelBuilder_addReinfLayer(clientData, interp, argc,argv,
				       theTclBuilder);
}


					 
int
TclModelBuilder_addRemoGeomTransf(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  return TclModelBuilder_addGeomTransf(clientData, interp, argc,argv,
				       theTclDomain,
				       theTclBuilder);
}


/// added by ZHY
extern int 
TclModelBuilderUpdateMaterialStageCommand(ClientData clientData, 
					  Tcl_Interp *interp, 
					  int argc, 
					  char **argv, 
					  TclModelBuilder *theTclBuilder);
int
TclModelBuilder_UpdateMaterialStage(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc, 
				    char **argv)
{
  return TclModelBuilderUpdateMaterialStageCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}
