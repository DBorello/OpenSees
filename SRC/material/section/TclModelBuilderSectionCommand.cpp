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
                                                                        
// $Revision: 1.9 $
// $Date: 2001-09-06 00:23:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/TclModelBuilderSectionCommand.cpp,v $
                                                                        
                                                                        
// File: ~/material/section/TclModelBuilderSectionCommand.C
// 
// Written: rms, MHS 
// Created: 07/99
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the section command in the interpreter.
//
// What: "@(#) TclModelBuilderMaterialCommands.C, revA"

#include <TclModelBuilder.h>

#include <ElasticSection2d.h>
#include <ElasticSection3d.h>
#include <GenericSection1d.h>
#include <GenericSectionNd.h>
#include <SectionAggregator.h>
#include <FiberSection.h>
#include <FiberSection2d.h>
#include <FiberSection3d.h>
#include <FiberSectionRepr.h>

#include <ElasticPlateSection.h>
#include <ElasticMembranePlateSection.h>
#include <MembranePlateFiberSection.h>

#include <QuadPatch.h>
#include <CircPatch.h>
#include <QuadCell.h>
#include <StraightReinfLayer.h>
#include <CircReinfLayer.h>
#include <ReinfBar.h>

#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>

#include <string.h>
#include <fstream.h>

static void printCommand(int argc, char **argv)
{
    cerr << "Input command: ";
    for (int i=0; i<argc; i++)
	cerr << argv[i] << " ";
    cerr << endl;
} 

int
TclModelBuilder_addFiberSection (ClientData clientData, Tcl_Interp *interp, int argc,
				 char **argv, TclModelBuilder *theBuilder);

int
TclModelBuilder_addUCFiberSection (ClientData clientData, Tcl_Interp *interp, int argc,
				   char **argv, TclModelBuilder *theBuilder);


int
TclModelBuilderSectionCommand (ClientData clientData, Tcl_Interp *interp, int argc,
			       char **argv, TclModelBuilder *theTclBuilder)
{
    // Pointer to a section that will be added to the model builder
    SectionForceDeformation *theSection = 0;

    // Check argv[1] for section type
    if (strcmp(argv[1],"Elastic") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    //cerr << "Want: section Elastic tag? EA? EIz? <EIy? GJ?>" << endl;
		cerr << "Want: section Elastic tag? E? A? Iz? <Iy? G? J?>" << endl;
	    return TCL_ERROR;
	}
	
	int tag;
	double E, A, Iz, Iy, G, J;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid section Elastic tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
	    cerr << "WARNING invalid E" << endl;
	    cerr << "Elastic section: " << tag << endl;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &A) != TCL_OK) {
	    cerr << "WARNING invalid A" << endl;
	    cerr << "Elastic section: " << tag << endl;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[5], &Iz) != TCL_OK) {
	    cerr << "WARNING invalid Iz" << endl;
	    cerr << "Elastic section: " << tag << endl;	    	    
	    return TCL_ERROR;
	}	
	
	if (argc > 8) {
	    if (Tcl_GetDouble (interp, argv[6], &Iy) != TCL_OK) {
		cerr << "WARNING invalid EIy" << endl;
		cerr << "Elastic section: " << tag << endl;	    		
		return TCL_ERROR;
	    }
		
	    if (Tcl_GetDouble (interp, argv[7], &G) != TCL_OK) {
		cerr << "WARNING invalid G" << endl;
		cerr << "Elastic section: " << tag << endl;	    
		return TCL_ERROR;
	    }

		if (Tcl_GetDouble (interp, argv[8], &J) != TCL_OK) {
		cerr << "WARNING invalid J" << endl;
		cerr << "Elastic section: " << tag << endl;	    
		return TCL_ERROR;
	    }
		
	    // Parsing was successful, allocate the section
	    theSection = new ElasticSection3d (tag, E, A, Iz, Iy, G, J);
	}
	else
	    // Parsing was successful, allocate the section
	    theSection = new ElasticSection2d (tag, E, A, Iz);
    }	
	
    else if (strcmp(argv[1],"Generic1D") == 0 || strcmp(argv[1],"Generic1d") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: section Generic1D tag? 1DTag? code?" << endl;
	    return TCL_ERROR;
	}

	int tag, uniTag, code;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid section Generic1D tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[3], &uniTag) != TCL_OK) {
	    cerr << "WARNING invalid 1DTag" << endl;
	    cerr << "Generic1D section: " << tag << endl;	    
	    return TCL_ERROR;		
	}

	if (strcmp(argv[4],"Mz") == 0)
	    code = SECTION_RESPONSE_MZ;
	else if (strcmp(argv[4],"P") == 0)
	    code = SECTION_RESPONSE_P;
	else if (strcmp(argv[4],"Vy") == 0)
	    code = SECTION_RESPONSE_VY;
	else if (strcmp(argv[4],"My") == 0)
	    code = SECTION_RESPONSE_MY;
	else if (strcmp(argv[4],"Vz") == 0)
	    code = SECTION_RESPONSE_VZ;
	else if (strcmp(argv[4],"T") == 0)
	    code = SECTION_RESPONSE_T;
	else {
	    cerr << "WARNING invalid code" << endl;
	    cerr << "Generic1D section: " << tag << endl;
	    return TCL_ERROR;		
	}
		
	// Retrieve the uniaxial material from the model builder
	UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(uniTag);
	
	if (theMat == 0) {
	    cerr << "WARNING uniaxial material does not exist\n";
	    cerr << "uniaxial material: " << uniTag; 
	    cerr << "\nGeneric1D section: " << tag << endl;
	    return TCL_ERROR;
	}
	
	// Parsing was successful, allocate the section
	theSection = new GenericSection1d (tag, *theMat, code);
    }

    else if (strcmp(argv[1],"GenericND") == 0 || strcmp(argv[1],"GenericNd") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: section GenericNd tag? NDTag? code?" << endl;
	    return TCL_ERROR;
	}
	
	int tag, NDTag;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid section GenericNd tag" << endl;
	    return TCL_ERROR;		
	}
	
	if (Tcl_GetInt(interp, argv[3], &NDTag) != TCL_OK) {
	    cerr << "WARNING invalid NDTag" << endl;
	    cerr << "GenericNd section: " << tag << endl;	    
	    return TCL_ERROR;		
	}
	
	ID code(argc-4);
	
	int i,j;
	
	// Read in the code
	for (i = 4, j = 0; i < argc; i++, j++) {
	    if (strcmp(argv[i],"Mz") == 0)
		code(j) = SECTION_RESPONSE_MZ;
	    else if (strcmp(argv[i],"P") == 0)
		code(j) = SECTION_RESPONSE_P;
	    else if (strcmp(argv[i],"Vy") == 0)
		code(j) = SECTION_RESPONSE_VY;
	    else if (strcmp(argv[i],"My") == 0)
		code(j) = SECTION_RESPONSE_MY;
	    else if (strcmp(argv[i],"Vz") == 0)
		code(j) = SECTION_RESPONSE_VZ;
	    else if (strcmp(argv[i],"T") == 0)
		code(j) = SECTION_RESPONSE_T;
	    else {
		cerr << "WARNING invalid GenericND code" << endl;
		cerr << "\nGenericND section: " << tag << endl;
		return TCL_ERROR;		
	    }
	}
	
	// Retrieve the uniaxial material from the model builder
	NDMaterial *theMat = theTclBuilder->getNDMaterial(NDTag);
	
	if (theMat == 0) {
	    cerr << "WARNING nD material does not exist\n";
	    cerr << "nD material: " << NDTag; 
	    cerr << "\nGenericNd section: " << tag << endl;
	    return TCL_ERROR;
	}
	
	// Parsing was successful, allocate the section
	theSection = new GenericSectionNd (tag, *theMat, code);
    }	
    
    else if (strcmp(argv[1],"AddDeformation") == 0 || strcmp(argv[1],"Aggregator") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: section Aggregator tag? uniTag1? code1? ... <-section secTag?>" << endl;
	    return TCL_ERROR;
	}
	    
	int tag;
	int secTag;
	SectionForceDeformation *theSec = 0;
	    
	int nArgs = argc-3;
	
	for (int ii = 5; ii < argc; ii++) {
	    if (strcmp(argv[ii],"-section") == 0 && ++ii < argc) {
		if (Tcl_GetInt(interp, argv[ii], &secTag) != TCL_OK) {
		    cerr << "WARNING invalid Aggregator tag" << endl;
		    return TCL_ERROR;		
		}
		
		theSec = theTclBuilder->getSection(secTag);
		
		if (theSec == 0) {
		    cerr << "WARNING section does not exist\n";
		    cerr << "section: " << secTag; 
		    cerr << "\nsection Aggregator: " << tag << endl;
		    return TCL_ERROR;
		}
		
		nArgs -= 2;
	    }
	}
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid Aggregator tag" << endl;
	    return TCL_ERROR;		
	}

	int nMats = nArgs / 2;
	
	if (nArgs%2 != 0) {
	    cerr << "WARNING improper number of arguments for Aggregator" << endl;
	    return TCL_ERROR;
	}
	
	UniaxialMaterial **theMats = 0;
	ID codes(nMats);
	
	theMats = new UniaxialMaterial *[nMats];
	
	if (theMats == 0) {
	    cerr << "TclModelBuilderSection (Aggregator) -- unable to create uniaxial array" << endl;
	    return TCL_ERROR;
	}	
	
	int tagI;
	int i, j;
	
	for (i = 3, j = 0; j < nMats; i++, j++) {
	    if (Tcl_GetInt(interp, argv[i], &tagI) != TCL_OK) {
		cerr << "WARNING invalid Aggregator matTag" << endl;
		return TCL_ERROR;		
	    }
	    
	    theMats[j] = theTclBuilder->getUniaxialMaterial(tagI);
	    
	    if (theMats[j] == 0) {
		cerr << "WARNING uniaxial material does not exist\n";
		cerr << "uniaxial material: " << tagI; 
		cerr << "\nsection Aggregator: " << tag << endl;
		return TCL_ERROR;
	    }
	    
	    i++;
	    
	    if (strcmp(argv[i],"Mz") == 0)
		codes(j) = SECTION_RESPONSE_MZ;
	    else if (strcmp(argv[i],"P") == 0)
		codes(j) = SECTION_RESPONSE_P;
	    else if (strcmp(argv[i],"Vy") == 0)
		codes(j) = SECTION_RESPONSE_VY;
	    else if (strcmp(argv[i],"My") == 0)
		codes(j) = SECTION_RESPONSE_MY;
	    else if (strcmp(argv[i],"Vz") == 0)
		codes(j) = SECTION_RESPONSE_VZ;
	    else if (strcmp(argv[i],"T") == 0)
		codes(j) = SECTION_RESPONSE_T;
	    else {
		cerr << "WARNING invalid code" << endl;
		cerr << "\nsection Aggregator: " << tag << endl;
		return TCL_ERROR;		
	    }
	}
	
	if (theSec)
	    theSection = new SectionAggregator (tag, *theSec, nMats, theMats, codes);
	else
	    theSection = new SectionAggregator (tag, nMats, theMats, codes);
	
	delete [] theMats;
    }		
    
    else if (strcmp(argv[1],"Fiber") == 0 || strcmp(argv[1],"fiberSec") == 0)
	return TclModelBuilder_addFiberSection (clientData, interp, argc, argv,
						theTclBuilder);

    else if (strcmp(argv[1],"UCFiber") == 0)
	return TclModelBuilder_addUCFiberSection (clientData, interp, argc, argv,
						  theTclBuilder);

    else if (strcmp(argv[1],"ElasticPlateSection") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: section ElasticPlateSection tag? E? nu? h? " << endl;
	    return TCL_ERROR;
	}
	
	int tag;
	double E, nu, h;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid section ElasticPlateSection tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
	    cerr << "WARNING invalid E" << endl;
	    cerr << "ElasticPlateSection section: " << tag << endl;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &nu) != TCL_OK) {
	    cerr << "WARNING invalid nu" << endl;
	    cerr << "ElasticPlateSection section: " << tag << endl;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[5], &h) != TCL_OK) {
	    cerr << "WARNING invalid h" << endl;
	    cerr << "ElasticPlateSection section: " << tag << endl;	    	    
	    return TCL_ERROR;
	}	

	theSection = new ElasticPlateSection (tag, E, nu, h);
    }	

    else if (strcmp(argv[1],"ElasticMembranePlateSection") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: section ElasticMembranePlateSection tag? E? nu? h? <rho?>" << endl;
	    return TCL_ERROR;
	}
	
	int tag;
	double E, nu, h;
	double rho = 0.0;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid section ElasticMembranePlateSection tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
	    cerr << "WARNING invalid E" << endl;
	    cerr << "ElasticMembranePlateSection section: " << tag << endl;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &nu) != TCL_OK) {
	    cerr << "WARNING invalid nu" << endl;
	    cerr << "ElasticMembranePlateSection section: " << tag << endl;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[5], &h) != TCL_OK) {
	    cerr << "WARNING invalid h" << endl;
	    cerr << "ElasticMembranePlateSection section: " << tag << endl;	    	    
	    return TCL_ERROR;
	}	

	if (argc > 6 && Tcl_GetDouble (interp, argv[6], &rho) != TCL_OK) {
	    cerr << "WARNING invalid rho" << endl;
	    cerr << "ElasticMembranePlateSection section: " << tag << endl;	    	    
	    return TCL_ERROR;
	}

	theSection = new ElasticMembranePlateSection (tag, E, nu, h, rho);
    }	

    else if (strcmp(argv[1],"PlateFiber") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: section PlateFiber tag? matTag? h? " << endl;
	    return TCL_ERROR;
	}
	
	int tag, matTag;
	double  h;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid section PlateFiber tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
	    cerr << "WARNING invalid matTag" << endl;
	    cerr << "PlateFiber section: " << matTag << endl;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &h) != TCL_OK) {
	    cerr << "WARNING invalid h" << endl;
	    cerr << "PlateFiber section: " << tag << endl;	    	    
	    return TCL_ERROR;
	}	

	NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matTag);
	if (theMaterial == 0) {
	    cerr << "WARNING nD material does not exist\n";
	    cerr << "nD material: " << matTag; 
	    cerr << "\nPlateFiber section: " << tag << endl;
	    return TCL_ERROR;
	}

	theSection = new MembranePlateFiberSection( tag, h, *theMaterial );
    }	
    
    else {
	cerr << "WARNING unknown type of section: " << argv[2];
	cerr << "Valid types: Elastic, Generic1d, GenericNd, Aggregator, Fiber\n";
	return TCL_ERROR;
    }
    
    // Ensure we have created the Material, out of memory if got here and no section
    if (theSection == 0) {
	cerr << "WARNING ran out of memory creating section\n";
	cerr << argv[2] << endl;
	return TCL_ERROR;
    }
    
    // Now add the material to the modelBuilder
    if (theTclBuilder->addSection(*theSection) < 0) {
	cerr << "WARNING could not add section to the domain\n";
	cerr << *theSection << endl;
	delete theSection; // invoke the material objects destructor, otherwise mem leak
	return TCL_ERROR;
    }
    
    return TCL_OK;
}

static int currentSectionTag = 0;
    
int
buildSection (Tcl_Interp *interp, TclModelBuilder *theTclModelBuilder, int secTag);

int
TclModelBuilder_addFiberSection (ClientData clientData, Tcl_Interp *interp, int argc,
				 char **argv, TclModelBuilder *theTclModelBuilder)
{
    int secTag;
    int maxNumPatches = 30; 
    int maxNumReinfLayers = 30;
    
    if (argc < 4) 
	return TCL_ERROR;
    
    if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
	interp->result = "WARNING bad command - want: \nsection fiberSec secTag { \n\tpatch <patch arguments> \n\tlayer <layer arguments> \n}";
	return TCL_ERROR;
    }
    
    currentSectionTag = secTag;
      
    // create the fiber section representation (with the geometric information) 
      
    SectionRepres *fiberSectionRepr =
	new FiberSectionRepr(secTag, maxNumPatches, maxNumReinfLayers);  

    if (fiberSectionRepr == 0) {
	interp->result = "WARNING - ran out of memory to create section representation";
	return TCL_ERROR;
    }

    if (theTclModelBuilder->addSectionRepres(*fiberSectionRepr) < 0) {
	interp->result = "WARNING - cannot add section representation";
	return TCL_ERROR;
    }	

    int brace = 3; // Start of recursive parse

    // parse the information inside the braces (patches and reinforcing layers)
    if (Tcl_Eval(interp, argv[brace]) != TCL_OK) {
	cerr << "WARNING - error reading information in { } ";
	return TCL_ERROR;
    }

    // build the fiber section (for analysis)
    if (buildSection(interp, theTclModelBuilder, secTag) != TCL_OK) {
	cerr << "WARNING - error constructing the section ";
	return TCL_ERROR;
    }
    
//    currentSectionTag = 0;

    return TCL_OK;
}




// add patch to fiber section
int
TclModelBuilder_addPatch(ClientData clientData, Tcl_Interp *interp, int argc, 
			     char **argv, TclModelBuilder *theTclModelBuilder)
{
    // check if a section is being processed
    if (currentSectionTag == 0) {
	interp->result = "WARNING subcommand 'patch' is only valid inside a 'section' command";
	return TCL_ERROR;
    }	   
    
    // make sure at least one other argument to contain patch type
    if (argc < 2) {
	interp->result = "WARNING need to specify a patch type ";
	return TCL_ERROR;
    }    

    // check argv[1] for type of patch  and create the object
    if (strcmp(argv[1], "quad") == 0 || strcmp(argv[1], "quadr") == 0) {
	int numSubdivIJ, numSubdivJK, matTag, secTag;
	double vertexCoordY, vertexCoordZ;
	static Matrix vertexCoords(4,2);
	int j, argi;

	if (argc < 13) {
	    interp->result = "WARNING invalid number of parameters: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
	    return TCL_ERROR;
	}
  
	argi = 2;
      
      if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK)
      {
         interp->result = "WARNING invalid matTag: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
         return TCL_ERROR;
      }
      //cerr << "\n\tmatTag: " << matTag;

      if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK)
      {
         interp->result = "WARNING invalid numSubdivIJ: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
         return TCL_ERROR;
      }
      //cerr << "\n\tnumSubdivIJ: " << numSubdivIJ;
 
      if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK)
      {
         interp->result = "WARNING invalid numSubdivJK: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
         return TCL_ERROR;
      }
      //cerr << "\n\tnumSubdivJK: " << numSubdivJK;

      for (j=0; j < 4; j++)
      {
         //cerr << "\n\tVertexCoord: " << j;
         if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK)
         {
            interp->result = "WARNING invalid Coordinate y: ...yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
            return TCL_ERROR;
         }
         //cerr << "\n\t\tvertexCoordY: " << vertexCoordY; 

         if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK)
         {
            interp->result = "WARNING invalid Coordinate z: ...yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
            return TCL_ERROR;
         }
         //cerr << "\n\t\tvertexCoordZ: " << vertexCoordZ; 

         vertexCoords(j,0) = vertexCoordY;
         vertexCoords(j,1) = vertexCoordZ;
      }
       
      // get section representation
      secTag = currentSectionTag;
      
      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) 
      {
         interp->result = "WARNING cannot retrieve section";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection)
      {
         interp->result = "WARNING section invalid: patch can only be added to fiber sections";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      // create patch

      QuadPatch *patch = new QuadPatch(matTag, numSubdivIJ, numSubdivJK, vertexCoords);
      if (!patch)
      {
         interp->result = "WARNING cannot alocate patch";
         return TCL_ERROR;
      }

      //cerr << "\n\tpatch: " << *patch;
      
      // add patch to section representation

      int error = fiberSectionRepr->addPatch(*patch);
      delete patch;
      
      if (error)
      {
         interp->result = "WARNING cannot add patch to section";
         return TCL_ERROR;
      }  
  }
    
    
    // check argv[1] for type of patch  and create the object
    else if (strcmp(argv[1], "rect") == 0 || 
	     strcmp(argv[1], "rectangular") == 0) {
	
	int numSubdivIJ, numSubdivJK, matTag, secTag;
	double vertexCoordY, vertexCoordZ;
	static Matrix vertexCoords(4,2);
	int j, argi;

	if (argc < 9) {
	    interp->result = "WARNING invalid number of parameters: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertK zVertK";
	    return TCL_ERROR;
	}
  
	argi = 2;
      
	if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
	    interp->result = "WARNING invalid matTag: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK) {
	    interp->result = "WARNING invalid numSubdivIJ: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
	    return TCL_ERROR;
	}
 
	if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK) {
	    interp->result = "WARNING invalid numSubdivJK: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
	    return TCL_ERROR;
	}

	for (j=0; j < 2; j++) {
	    if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK) {
		interp->result = "WARNING invalid Coordinate y: ...yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
		return TCL_ERROR;
	    }

	    if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK) {
		interp->result = "WARNING invalid Coordinate z: ...yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL";
		return TCL_ERROR;
	    }

	    vertexCoords(j*2,0) = vertexCoordY;
	    vertexCoords(j*2,1) = vertexCoordZ;
	}

	vertexCoords(1,0) = vertexCoords(2,0);
	vertexCoords(1,1) = vertexCoords(0,1);	
	vertexCoords(3,0) = vertexCoords(0,0);
	vertexCoords(3,1) = vertexCoords(2,1);		
	    
      // get section representation
      secTag = currentSectionTag;
      
      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) {
         interp->result = "WARNING cannot retrieve section";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection) {
         interp->result = "WARNING section invalid: patch can only be added to fiber sections";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      // create patch

      QuadPatch *patch = new QuadPatch(matTag, numSubdivIJ, numSubdivJK, vertexCoords);
      if (!patch)
      {
         interp->result = "WARNING cannot alocate patch";
         return TCL_ERROR;
      }

      //cerr << "\n\tpatch: " << *patch;
      
      // add patch to section representation

      int error = fiberSectionRepr->addPatch(*patch);
      delete patch;
      
      if (error)
      {
         interp->result = "WARNING cannot add patch to section";
         return TCL_ERROR;
      }  
  }    
    
    
    
         
    else if (strcmp(argv[1], "circ") == 0) {
	int numSubdivRad, numSubdivCirc, matTag, secTag;
	double yCenter, zCenter;
	static Vector centerPosition(2);
	double intRad, extRad;
	double startAng, endAng;

	int argi;

      argi = 2;
      if (argc < 11)
      {
         interp->result = "WARNING invalid number of parameters: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
  
      if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK)
      {
         interp->result = "WARNING invalid matTag: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tmatTag: " << matTag;

      if (Tcl_GetInt(interp, argv[argi++], &numSubdivCirc) != TCL_OK)
      {
         interp->result = "WARNING invalid numSubdivCirc: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tnumSubdivCirc: " << numSubdivCirc;

      if (Tcl_GetInt(interp, argv[argi++], &numSubdivRad) != TCL_OK)
      {
         interp->result = "WARNING invalid numSubdivRad: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tnumSubdivRad: " << numSubdivRad;

      if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK)
      {
         interp->result = "WARNING invalid yCenter: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tyCenter: " << yCenter;

      if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK)
      {
         interp->result = "WARNING invalid zCenter: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tzCenter: " << zCenter;

      if (Tcl_GetDouble(interp, argv[argi++], &intRad) != TCL_OK)
      {
         interp->result = "WARNING invalid intRad: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tintRad: " << intRad;

      if (Tcl_GetDouble(interp, argv[argi++], &extRad) != TCL_OK)
      {
         interp->result = "WARNING invalid extRad: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\textRad: " << extRad;

      if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK)
      {
         interp->result = "WARNING invalid startAng: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tstartAngle: " << startAng;

      if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK)
      {
         interp->result = "WARNING invalid endAng: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tendAng: " << endAng;


      // get section 
      secTag = currentSectionTag;
      
      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) 
      {
         interp->result = "WARNING cannot retrieve section";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection)
      {
         interp->result = "WARNING section invalid: patch can only be added to fiber sections";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      centerPosition(0) = yCenter; 
      centerPosition(1) = zCenter; 
   
      // create patch

      CircPatch *patch = new CircPatch(matTag, numSubdivCirc, numSubdivRad,
                                       centerPosition, intRad, extRad, 
                                       startAng, endAng);
      if (!patch)
      {
         interp->result = "WARNING cannot alocate patch";
         return TCL_ERROR;
      }

      //cerr << "\n\tpatch: " << *patch;
      
      // add patch to section

      int error = fiberSectionRepr->addPatch(*patch);
      delete patch;
      
      if (error)
      {
         interp->result = "WARNING cannot add patch to section";
         return TCL_ERROR;
      }
   }

   else
   {
      interp->result = "WARNING patch type is not available";
      return TCL_ERROR;
   }
  
   return TCL_OK;
}



// add patch to fiber section
int
TclModelBuilder_addFiber(ClientData clientData, Tcl_Interp *interp, int argc, 
			 char **argv, TclModelBuilder *theTclModelBuilder)
{
    // check if a section is being processed
    if (currentSectionTag == 0) {
	interp->result = "WARNING subcommand 'fiber' is only valid inside a 'section' command";
	return TCL_ERROR;
    }	   
    
    // make sure at least one other argument to contain patch type
    if (argc < 5) {
	interp->result = "WARNING invalid num args: fiber yLoc zLoc area matTag";
	return TCL_ERROR;
    }    

    SectionRepres *sectionRepres = 
	theTclModelBuilder->getSectionRepres(currentSectionTag);
    
    if (sectionRepres == 0) {
	interp->result = "WARNING cannot retrieve section";
	return TCL_ERROR;
    }    
	
    if (sectionRepres->getType() != SEC_TAG_FiberSection) {
	interp->result = "WARNING section invalid: patch can only be added to fiber sections";
	return TCL_ERROR;
    }

    FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;
    int numFibers = fiberSectionRepr->getNumFibers();    
    
    int NDM = theTclModelBuilder->getNDM();  
    
    Fiber *theFiber =0;
      
    int matTag;
    double yLoc, zLoc, area;

    
    if (Tcl_GetDouble(interp, argv[1], &yLoc) != TCL_OK) {
         interp->result = "WARNING invalid yLoc: fiber yLoc zLoc area matTag";
         return TCL_ERROR;
     }    
    if (Tcl_GetDouble(interp, argv[2], &zLoc) != TCL_OK) {
         interp->result = "WARNING invalid zLoc: fiber yLoc zLoc area matTag";
         return TCL_ERROR;
     }        
    if (Tcl_GetDouble(interp, argv[3], &area) != TCL_OK) {
         interp->result = "WARNING invalid area: fiber yLoc zLoc area matTag";
         return TCL_ERROR;
     }            
    
    if (Tcl_GetInt(interp, argv[4], &matTag) != TCL_OK) {
         interp->result = "WARNING invalid matTag: fiber yLoc zLoc area matTag";
         return TCL_ERROR;
     }                
    
    UniaxialMaterial *material = theTclModelBuilder->getUniaxialMaterial(matTag);
    
    // creates 2d section      
    if (NDM == 2) {

	if (material == 0) {
	    interp->result = "WARNING invalid material ID for patch";
	    return TCL_ERROR;
	}   

	theFiber = new UniaxialFiber2d(numFibers, *material, area, zLoc);
	if (theFiber == 0) {
	    interp->result = "WARNING unable to allocate fiber ";
	    return TCL_ERROR;
	}    
    }

    else if (NDM == 3) {

      static Vector fiberPosition(2);
	fiberPosition(0) = yLoc;
	fiberPosition(1) = zLoc;
	    
	theFiber = new UniaxialFiber3d(numFibers, *material, area, fiberPosition);
	if (theFiber == 0) {
	    interp->result = "WARNING unable to allocate fiber ";
	    return TCL_ERROR;
	}    
    }

    else {
	interp->result = "WARNING fiber command for FiberSection only fo 2 or 3d ";
	return TCL_ERROR;
    }    
	
    // add patch to section representation
    int error = fiberSectionRepr->addFiber(*theFiber);

    if (error) {
	interp->result = "WARNING cannot add patch to section";
	return TCL_ERROR;
    }  

    return TCL_OK;
}




// add layers of reinforcing bars to fiber section
          
int
TclModelBuilder_addReinfLayer(ClientData clientData, Tcl_Interp *interp, int argc, 
				  char **argv, TclModelBuilder *theTclModelBuilder)
{
   //cerr << "\nreading layer:";

   // check if a section is being processed
   if (currentSectionTag == 0)
   {
      interp->result = "WARNING subcommand 'patch' is only valid inside a 'section' command";
      return TCL_ERROR;
   }	   

   // make sure at least one other argument to contain layer type
   if (argc < 2) 
   {
      interp->result = "WARNING need to specify a layer type ";
      return TCL_ERROR;
   }    

   // check argv[1] for type of layer and create the object
   if (strcmp(argv[1], "straight") == 0) 
   {
      if (argc < 9)
      {
         interp->result = "WARNING invalid number of parameters: layer straight matTag numReinfBars reinfBarArea yStartPt zStartPt yEndPt zEndPt";
         return TCL_ERROR;
      }

      int secTag, matTag, numReinfBars;
      double reinfBarArea;
      double yStartPt, zStartPt, yEndPt, zEndPt;
     
      int argi;

      argi = 2;
      
      if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK)
      {
         interp->result = "WARNING invalid matTag: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt";
         return TCL_ERROR;
      }
      //cerr << "\n\tmatTag: " << matTag;

      if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK)
      {
         interp->result = "WARNING invalid numReinfBars: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt";
         return TCL_ERROR;
      }
      //cerr << "\n\tnumReinfBars: " << numReinfBars;

      if (Tcl_GetDouble(interp, argv[argi++], &reinfBarArea) != TCL_OK)
      {
         interp->result = "WARNING invalid reinfBarArea: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt";
         return TCL_ERROR;
      }
      //cerr << "\n\treinfBarArea: " << reinfBarArea;

      if (Tcl_GetDouble(interp, argv[argi++], &yStartPt) != TCL_OK)
      {
         interp->result = "WARNING invalid yStartPt: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt";
         return TCL_ERROR;
      }
      //cerr << "\n\tyStartPt: " << yStartPt;
    
      if (Tcl_GetDouble(interp, argv[argi++], &zStartPt) != TCL_OK)
      {
         interp->result = "WARNING invalid zStartPt: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt";
         return TCL_ERROR;
      }
      //cerr << "\n\tzStartPt: " << zStartPt;
       
      if (Tcl_GetDouble(interp, argv[argi++], &yEndPt) != TCL_OK)
      {
         interp->result = "WARNING invalid yEndPt: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt";
         return TCL_ERROR;
      }
      //cerr << "\n\tyEndPt: " << yEndPt;
    
      if (Tcl_GetDouble(interp, argv[argi++], &zEndPt) != TCL_OK)
      {
         interp->result = "WARNING invalid zEndPt: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt";
         return TCL_ERROR;
      }
      
      //cerr << "\n\tzEndPt: " << zEndPt;
      
      // get section 
      secTag = currentSectionTag;

      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) 
      {
         interp->result = "WARNING cannot retrieve section";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection)
      {
         interp->result = "WARNING section invalid: patch can only be added to fiber sections";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;
 
      // create the reinforcing layer

      static Vector startPt(2);
      static Vector endPt(2);

      startPt(0)  = yStartPt;
      startPt(1)  = zStartPt;
      endPt(0) = yEndPt;
      endPt(1) = zEndPt;

      StraightReinfLayer *reinfLayer = new StraightReinfLayer (matTag,
                                                   numReinfBars, reinfBarArea,
                                                   startPt, endPt);
      if (!reinfLayer)
      {
         interp->result = "WARNING cannot alocate reinfLayer";
         return TCL_ERROR;
      }
      //cerr << "\nStraigthReinfLayer: " << *reinfLayer;

      // add reinfLayer to section
      int error = fiberSectionRepr->addReinfLayer(*reinfLayer);
      delete reinfLayer;
      
      if (error)
      {
         interp->result = "WARNING cannot add reinforcing layer to section";
         return TCL_ERROR;
      }
      
   }
   else if (strcmp(argv[1], "circ") == 0) 
   {
      if (argc < 8)
      {
         interp->result = "WARNING invalid number of parameters: layer circ matTag numReinfBars reinfBarArea yCenter zCenter arcRadius <startAng endAng>";
         return TCL_ERROR;
      }

      int secTag, matTag, numReinfBars;
      double reinfBarArea;
      double yCenter, zCenter, radius, startAng, endAng;
     
      int argi;

      argi = 2;
      
      if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK)
      {
         interp->result = "WARNING invalid matTag: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tmatTag: " << matTag;

      if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK)
      {
         interp->result = "WARNING invalid numReinfBars: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tnumReinfBars: " << numReinfBars;

      if (Tcl_GetDouble(interp, argv[argi++], &reinfBarArea) != TCL_OK)
      {
         interp->result = "WARNING invalid reinfBarArea: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\treinfBarArea: " << reinfBarArea;

      if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK)
      {
         interp->result = "WARNING invalid yCenter: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tyCenter: " << yCenter;
    
      if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK)
      {
         interp->result = "WARNING invalid zCenter: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tzCenter: " << zCenter;
       
      if (Tcl_GetDouble(interp, argv[argi++], &radius) != TCL_OK)
      {
         interp->result = "WARNING invalid radius: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng";
         return TCL_ERROR;
      }
      //cerr << "\n\tradius: " << radius;
    
	  bool anglesSpecified = false;

      if (argc > 9) {
		  if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK)
		  {
			 interp->result = "WARNING invalid startAng: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng";
			 return TCL_ERROR;
		  }
		  //cerr << "\n\tstartAng: " << startAng;

		  if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK)
		  {
			 interp->result = "WARNING invalid endAng: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng";
			 return TCL_ERROR;
		  }
		  //cerr << "\n\tendAng: " << endAng;

		  anglesSpecified = true;
	  }

      // get section 
      secTag = currentSectionTag;
      
      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) 
      {
         interp->result = "WARNING cannot retrieve section";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection)
      {
         interp->result = "WARNING section invalid: patch can only be added to fiber sections";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;
 
      // create the reinforcing layer

      static Vector center(2);

      center(0) = yCenter; 
      center(1) = zCenter; 

      CircReinfLayer *reinfLayer = 0;
	  if (anglesSpecified)
		  // Construct arc
		  reinfLayer = new CircReinfLayer (matTag, numReinfBars, reinfBarArea,
			center, radius, startAng, endAng);
	  else
		  // Construct circle
		  reinfLayer = new CircReinfLayer (matTag, numReinfBars, reinfBarArea,
			center, radius);

      if (!reinfLayer)
      {
         interp->result = "WARNING cannot alocate reinfLayer";
         return TCL_ERROR;
      }
      //cerr << "\nCircReinfLayer: " << *reinfLayer;

      // add reinfLayer to section
      int error = fiberSectionRepr->addReinfLayer(*reinfLayer);
      delete reinfLayer;
      
      if (error)
      {
         interp->result = "WARNING cannot add reinforcing layer to section";
         return TCL_ERROR;
      }
      
   }
   else
   {
      interp->result = "WARNING reinforcing layer type is not available";
      return TCL_ERROR;
   }
  
   return TCL_OK;
}    



// build the section
int 
buildSection (Tcl_Interp *interp, TclModelBuilder *theTclModelBuilder, int secTag)
{
   SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
   if (sectionRepres == 0) 
   {
      interp->result = "WARNING cannot retrieve section";
      return TCL_ERROR;
   }    
     
   if (sectionRepres->getType() == SEC_TAG_FiberSection)
   {
      // build the section
  
      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      int i, j, k;
      int numFibers;
      
      int numPatches;
      Patch **patch;

      int  numReinfLayers;
      ReinfLayer **reinfLayer;

      numPatches     = fiberSectionRepr->getNumPatches();
      patch          = fiberSectionRepr->getPatches();
      numReinfLayers = fiberSectionRepr->getNumReinfLayers();
      reinfLayer     = fiberSectionRepr->getReinfLayers(); 

      int numSectionRepresFibers = fiberSectionRepr->getNumFibers();
      Fiber **sectionRepresFibers = fiberSectionRepr->getFibers();
      
      numFibers = numSectionRepresFibers;
      for (i = 0; i < numPatches; i++)
         numFibers += patch[i]->getNumCells();
      
      for (i = 0; i < numReinfLayers; i++)
         numFibers += reinfLayer[i]->getNumReinfBars();
      
      //cerr << "\nnumFibers: " << numFibers;
      
      static Vector fiberPosition(2);
      int    matTag;
      
      ID     fibersMaterial(numFibers-numSectionRepresFibers);
      Matrix fibersPosition(2,numFibers-numSectionRepresFibers);
      Vector fibersArea(numFibers-numSectionRepresFibers);

      int  numCells;
      Cell **cell;
    
      k = 0;
      for (i = 0; i < numPatches; i++)
      {
         //cerr << "\nPatch :" << i;
      
         numCells   = patch[i]->getNumCells();
         matTag = patch[i]->getMaterialID();

         //cerr << "\nmatTag: " << matTag(k);

         cell = patch[i]->getCells();

         if (cell == 0)
         {
            interp->result = "WARNING out of run to create fibers";
            return TCL_ERROR;
         }    
         
         //cerr << "\n\tnumCells :" << numCells;
      
         for (j = 0; j < numCells; j++)
         {
	    fibersMaterial(k) = matTag;
            fibersArea(k)     = cell[j]->getArea();
            fiberPosition     = cell[j]->getCentroidPosition();

            fibersPosition(0,k) = fiberPosition(0);
	    fibersPosition(1,k) = fiberPosition(1);
	      
            k++;
         }
  
         for (j = 0; j < numCells; j++)
           delete cell[j];
  
         delete [] cell;
      }
         
      ReinfBar *reinfBar;
      int numReinfBars;

      for (i = 0; i < numReinfLayers; i++)
      {
         numReinfBars = reinfLayer[i]->getNumReinfBars();
         reinfBar     = reinfLayer[i]->getReinfBars();
         matTag  = reinfLayer[i]->getMaterialID();
   
         for (j = 0; j < numReinfBars; j++)
         {
	    fibersMaterial(k) = matTag; 
            fibersArea(k)     = reinfBar[j].getArea();
            fiberPosition     = reinfBar[j].getPosition();
     
	    fibersPosition(0,k) = fiberPosition(0);
	    fibersPosition(1,k) = fiberPosition(1);
	
            k++;
         }
         delete [] reinfBar;
      }

      UniaxialMaterial *material;
      
      int NDM = theTclModelBuilder->getNDM();   // dimension of the structure (1d, 2d, or 3d)


      Fiber **fiber = new Fiber *[numFibers];
      if (fiber == 0) {
	  interp->result = "WARNING unable to allocate fibers ";
	  return TCL_ERROR;
      }          
      
      // copy the section repres fibers
      for (i=0; i<numSectionRepresFibers; i++)
	  fiber[i] = sectionRepresFibers[i];

      // creates 2d section      


      if (NDM == 2)     
      {
	 k = 0;
	 for (i = numSectionRepresFibers; i < numFibers; i++)
	 {    
            material = theTclModelBuilder->getUniaxialMaterial(fibersMaterial(k));
            if (material == 0)
            {
               interp->result = "WARNING invalid material ID for patch";
               return TCL_ERROR;
            }   
	    
	    fiber[i] = new UniaxialFiber2d(k, *material, fibersArea(k), fibersPosition(0,k));
            if (!fiber[i]) 
            {
               interp->result = "WARNING unable to allocate fiber ";
               return TCL_ERROR;
            }    
   
            //cerr << *fiber[k];
	    k++;
	 }
	
	 SectionForceDeformation *section = new FiberSection2d(secTag, numFibers, fiber);
	 //SectionForceDeformation *section = new FiberSection(secTag, numFibers, fiber);
   
	 // Delete fibers
	 for (i = 0; i < numFibers; i++)
	   delete fiber[i];

         if (section == 0)
         {
            interp->result = "WARNING - cannot construct section";
            return TCL_ERROR;
         }
       
         if (theTclModelBuilder->addSection (*section) < 0)
         {
            interp->result = "WARNING - cannot add section";
            return TCL_ERROR;
         }

        //cerr << "section: " << *section;
     
      }
      else if (NDM == 3)     
      {

	 static Vector fiberPosition(2);
	 k = 0;
	 for (i = numSectionRepresFibers; i < numFibers; i++)
	 {    
            material = theTclModelBuilder->getUniaxialMaterial(fibersMaterial(k));
            if (material == 0)
            {
               interp->result = "WARNING invalid material ID for patch";
               return TCL_ERROR;
            }   
	    
	    fiberPosition(0) = fibersPosition(0,k);
	    fiberPosition(1) = fibersPosition(1,k);
	    
	    fiber[i] = new UniaxialFiber3d(k, *material, fibersArea(k), fiberPosition);
            if (!fiber[k]) 
            {
               interp->result = "WARNING unable to allocate fiber ";
               return TCL_ERROR;
            }    
	    k++;
            //cerr << *fiber[k];
	 }
	
	 //SectionForceDeformation *section = new FiberSection(secTag, numFibers, fiber);
	 SectionForceDeformation *section = new FiberSection3d(secTag, numFibers, fiber);
   
	 // Delete fibers
	 for (i = 0; i < numFibers; i++)
	   delete fiber[i];

         if (section == 0)
         {
            interp->result = "WARNING - cannot construct section";
            return TCL_ERROR;
         }
       
         if (theTclModelBuilder->addSection (*section) < 0)
         {
            interp->result = "WARNING - cannot add section";
            return TCL_ERROR;
         }

        //cerr << "section: " << *section;
       
      }
      else
      {
         cerr << "WARNING NDM = " << NDM << " is imcompatible with available frame elements";
         return TCL_ERROR;
      }

      // Delete fiber array
      delete [] fiber;

   }
   else 
   {
      interp->result = "WARNING section invalid: can only build fiber sections";
      return TCL_ERROR;
   }    

   return TCL_OK;
}




int
TclModelBuilder_addUCFiberSection (ClientData clientData, Tcl_Interp *interp, int argc,
				 char **argv, TclModelBuilder *theTclModelBuilder)
{
    int secTag;
    
    if (argc < 4) 
	return TCL_ERROR;
    
    if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
      interp->result = "could not read section tag\n";
      return TCL_ERROR;
    }

    currentSectionTag = secTag;

    // first create an empty FiberSection
    int NDM = theTclModelBuilder->getNDM();   // dimension of the structure (1d, 2d, or 3d)

    SectionForceDeformation *section = 0;
    FiberSection2d *section2d =0;
    FiberSection3d *section3d =0;

    if (NDM == 2) {
      section2d = new FiberSection2d(secTag, 0, 0);
      section = section2d;
      //SectionForceDeformation *section = new FiberSection(secTag, 0, 0);
    } else if (NDM == 3) {
      section3d = new FiberSection3d(secTag, 0, 0);
      section = section3d;
    } 

    if (section == 0) {
      return TCL_ERROR;
    }

    //
    // now parse the ouput file containing the fiber data, 
    // create fibers and add them to the section
    //

    // open the file
    char *fileName = argv[3];
    ifstream theFile;
    theFile.open(fileName, ios::in);
    if (!theFile) {
      cerr << "section UCFiber - could not open file named " << *fileName;
      return TCL_ERROR;
    } else {
      int foundStart = 0;
      static char garbage[100];

      // parse through until find start of fiber data
      while (foundStart == 0 && theFile >> garbage) 
	if (strcmp(garbage, "#FIBERS") == 0) 
	  foundStart = 1;

      if (foundStart == 0) {
	theFile.close();
	return TCL_ERROR;
      }

      // parse the fiber data until eof, creating a fiber and adding to section as go
      double ycoord, zcoord, area, prestrain;
      int matTag;
      int fiberCount = 0;
      
      while (theFile >> ycoord >> zcoord >> area >> prestrain >> garbage >> matTag) {

	UniaxialMaterial *theMaterial = theTclModelBuilder->getUniaxialMaterial(matTag);
	if (theMaterial == 0) {
	  cerr << "section UCFiber - no material exists with tag << " << matTag << endl;
	  return TCL_ERROR;
	}
	
	Fiber *theFiber = 0;
	if (NDM == 2) {
	  theFiber = new UniaxialFiber2d(fiberCount++, *theMaterial, area, zcoord);
	  if (theFiber != 0) {
	    section2d->addFiber(*theFiber);
	    delete theFiber;
	  }
	} else {
	  static Vector pos(2);
	  pos(0) = ycoord; pos(1) = zcoord;
	  theFiber = new UniaxialFiber3d(fiberCount++, *theMaterial, area, pos);
	  if (theFiber != 0) {
	    section3d->addFiber(*theFiber);
	    delete theFiber;
	  }
	}   
      }

      // close the file
      theFile.close();
    }

    // finally add the section to our modelbuilder
    if (theTclModelBuilder->addSection (*section) < 0) {
      interp->result = "WARNING - cannot add section";
      return TCL_ERROR;
    }


    return TCL_OK;
}
