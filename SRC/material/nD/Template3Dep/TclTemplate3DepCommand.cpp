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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-14 23:01:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/Template3Dep/TclTemplate3DepCommand.cpp,v $
                                                                        
// Written: fmk 
// Created: 04/01
//
// Description: This file contains the implementation of the TclModelBuilder_addTemplate3Dep() 
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>

#include <Domain.h>

#include <ErrorHandler.h>
#include <TclModelBuilder.h>

#include <Template3Dep.h>
#include <straint.h>
#include <stresst.h>

#include <YS.h>
#include <DP_YS.h>
#include <VM_YS.h>
#include <CAM_YS.h>

#include <PS.h>
#include <DP_PS.h>
#include <VM_PS.h>
#include <CAM_PS.h>

#include <EPState.h>

#include <EL_S.h>
#include <EL_LEeq.h>
#include <EL_NLEp.h>

#include <EL_T.h>
#include <EL_LEij.h>





// the functions to create the component objects (defined at eof)
YieldSurface     *EvaluateYieldSurfaceArgs(ClientData, Tcl_Interp *, char *tclString);
PotentialSurface *EvaluatePotentialSurfaceArgs(ClientData, Tcl_Interp *, char *tclString);
EPState          *EvaluateEPStateArgs(ClientData, Tcl_Interp *, char *tclString);
EvolutionLaw_S   *EvaluateEvolutionLawSArgs(ClientData, Tcl_Interp *, char *tclString);
EvolutionLaw_T   *EvaluateEvolutionLawTArgs(ClientData, Tcl_Interp *, char *tclString);


// little function to free memory after invoke Tcl_SplitList
//   note Tcl_Split list stores the array of pointers and the strings in 
//   one array, which is why Tcl_Free needs only be called on the array.
static void cleanup(char **argv) {
	  Tcl_Free((char *) argv);
}

Template3Dep *
TclModelBuilder_addTemplate3Dep(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  char **argv, TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // create some empty pointers which we fill in as parse the command line
  int tag =0;
  YieldSurface     *YS =0;
  PotentialSurface *PS =0;
  EPState          *EPS =0;
  EvolutionLaw_S   *ELS1 =0; 
  EvolutionLaw_S   *ELS2 =0; 
  EvolutionLaw_S   *ELS3 =0; 
  EvolutionLaw_S   *ELS4 =0; 
  EvolutionLaw_T   *ELT1 =0;
  EvolutionLaw_T   *ELT2 =0;
  EvolutionLaw_T   *ELT3 =0;
  EvolutionLaw_T   *ELT4 =0;

  int loc = eleArgStart;
  if (Tcl_GetInt(interp, argv[loc++], &tag) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep - invalid tag " << argv[loc] << endln;
    return 0;		
  }

  // parse rest of command, switching on the flags -YS, -PS, -EPS, -ELS, -ELT, ...
  while (loc < argc) {
    if ((strcmp(argv[loc],"-YS") == 0) || (strcmp(argv[loc],"-ys") == 0)) {
      YS = EvaluateYieldSurfaceArgs(clientData, interp, argv[loc+1]);
      if (YS == 0) {
	opserr << "nDMaterial Templated3Dep - could not create a YS from" << argv[loc+1] << endln;
	return 0;
      }
    }

    else if ((strcmp(argv[loc],"-PS") == 0) || (strcmp(argv[loc],"-ps") == 0)) {
      PS = EvaluatePotentialSurfaceArgs(clientData, interp, argv[loc+1]);
      if (PS == 0) {
	opserr << "nDMaterial Templated3Dep - could not create a PS from" << argv[loc+1] << endln;
	return 0;	
      }
    }    

   else if ((strcmp(argv[loc],"-EPS") == 0) || (strcmp(argv[loc],"-eps") == 0)) {
      EPS = EvaluateEPStateArgs(clientData, interp, argv[loc+1]);
      if (EPS == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an EPS from" << argv[loc+1] << endln;
	return 0;		
      }
    }    
    
    else if ((strcmp(argv[loc],"-ELS1") == 0) || (strcmp(argv[loc],"-els1") == 0)) {
      ELS1 = EvaluateEvolutionLawSArgs(clientData, interp, argv[loc+1]);
      if (ELS1 == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an ELS1 from" << argv[loc+1] << endln;
	return 0;		
      }
    }    

    else if ((strcmp(argv[loc],"-ELS2") == 0) || (strcmp(argv[loc],"-els2") == 0)) {
      ELS2 = EvaluateEvolutionLawSArgs(clientData, interp, argv[loc+1]);
      if (ELS2 == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an ELS2 from" << argv[loc+1] << endln;
	return 0;		
      }
    }    

    else if ((strcmp(argv[loc],"-ELS3") == 0) || (strcmp(argv[loc],"-els3") == 0)) {
      ELS3 = EvaluateEvolutionLawSArgs(clientData, interp, argv[loc+1]);
      if (ELS3 == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an ELS3 from" << argv[loc+1] << endln;
	return 0;		
      }
    }    

    else if ((strcmp(argv[loc],"-ELS4") == 0) || (strcmp(argv[loc],"-els4") == 0)) {
      ELS4 = EvaluateEvolutionLawSArgs(clientData, interp, argv[loc+1]);
      if (ELS4 == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an ELS4 from" << argv[loc+1] << endln;
	return 0;		
      }
    }    

    else if ((strcmp(argv[loc],"-ELT1") == 0) || (strcmp(argv[loc],"-elt1") == 0)) {
      ELT1 = EvaluateEvolutionLawTArgs(clientData, interp, argv[loc+1]);
      if (ELT1 == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an ELT1 from" << argv[loc+1] << endln;
	return 0;		
      }
    }    

    else if ((strcmp(argv[loc],"-ELT2") == 0) || (strcmp(argv[loc],"-elt2") == 0)) {
      ELT2 = EvaluateEvolutionLawTArgs(clientData, interp, argv[loc+1]);
      if (ELT2 == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an ELT2 from" << argv[loc+1] << endln;
	return 0;		
      }
    }    

    else if ((strcmp(argv[loc],"-ELT3") == 0) || (strcmp(argv[loc],"-elt3") == 0)) {
      ELT3 = EvaluateEvolutionLawTArgs(clientData, interp, argv[loc+1]);
      if (ELT3 == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an ELT3 from" << argv[loc+1] << endln;
	return 0;		
      }
    }    

    else if ((strcmp(argv[loc],"-ELT4") == 0) || (strcmp(argv[loc],"-elt4") == 0)) {
      ELT4 = EvaluateEvolutionLawTArgs(clientData, interp, argv[loc+1]);
      if (ELT4 == 0) {
	opserr << "nDMaterial Templated3Dep - could not create an ELT4 from" << argv[loc+1] << endln;
	return 0;		
      }
    }    

    else {
      opserr << "nDMaterial Templated3Dep - don't understand %s\n";
	return 0;	
    }

    // increment locator by 2 and do next one
    loc += 2;
  }

  // now depending on the objects types that are not NULL we use the appropriate 
  // constructor to construct or Template3Dep material object
  Template3Dep *theMaterial = 0;
  if ( (YS != 0) && (PS != 0) && (EPS != 0) && 
       (ELS1 != 0) && (ELS2 != 0) && (ELS3 != 0) && (ELS4 != 0) &&
       (ELT1 != 0) && (ELT2 != 0) && (ELT3 != 0) && (ELT4 != 0) )
    theMaterial = new Template3Dep(tag, YS, PS, EPS, ELS1, ELS2, ELS3, ELS4,
				   ELT1, ELT2, ELT3, ELT4);

  // constructor 0
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) && 
       (ELS1 == 0) && (ELS2 == 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 == 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, YS, PS, EPS);

  // constructor 1
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) && 
       (ELS1 != 0) && (ELS2 == 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 == 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, YS, PS, EPS, ELS1);

  // constructor 2
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) && 
       (ELS1 == 0) && (ELS2 == 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 != 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, YS, PS, EPS, ELT1);

  // constructor 3
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) && 
       (ELS1 != 0) && (ELS2 == 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 != 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, YS, PS, EPS, ELS1, ELT1);

  // constructor 4
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) && 
       (ELS1 != 0) && (ELS2 != 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 != 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, YS, PS, EPS, ELS1, ELS2, ELT1);

  // constructor 5
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) && 
       (ELS1 != 0) && (ELS2 != 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 != 0) && (ELT2 != 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, YS, PS, EPS, ELS1, ELS2, ELT1, ELT2);

  else
    opserr << "invalid number of args used to create a Template3Dep material\n";

    return theMaterial;
}


// Function - to create a YieldSurface
YieldSurface *EvaluateYieldSurfaceArgs(ClientData clientData, Tcl_Interp *interp, char *tclString)
{
  int argc;
  char **argv;
  
  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }

  if (argc == 0) 
    return 0;

  // now parse the list & construct the required object
  YieldSurface *YS = 0;

  // 1. Drucker-Prager Yield Surface
  if ((strcmp(argv[0],"-DP") == 0) || (strcmp(argv[0],"-dp") == 0)) {
    YS = new DPYieldSurface();
  }

  // 2. von Mises Yield Surface
  if ( (strcmp(argv[0],"-VM") == 0) || (strcmp(argv[0],"-vM") == 0) || (strcmp(argv[0],"-vm") == 0)) {
    YS = new VMYieldSurface();
  }
  
  // 3. Cam Clay yield surface
  if ((strcmp(argv[0],"-CC") == 0) || (strcmp(argv[0],"-cc") == 0)) {

    double mp = 1.2;
    if (argc > 1) {
      if (Tcl_GetDouble(interp, argv[1], &mp) != TCL_OK) {
	opserr << "invalid M: argv[1] for -PS CamClay M\n";
	return 0;		
      }
    }
    
    // create the object & return it
    YS = new CAMYieldSurface(mp);
  }

  cleanup(argv);
  return YS;
}

// Function - to create a PotentialSurface object
PotentialSurface *EvaluatePotentialSurfaceArgs(ClientData clientData, Tcl_Interp *interp, char *tclString)
{
  int argc;
  char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }

  // now parse the list & construct the required object
  PotentialSurface *PS = 0;
  //1. Drucker-Prager Potential Surface
  if ((strcmp(argv[0],"-DP") == 0) || (strcmp(argv[0],"-dp") == 0)) {
    double alpha = 0.0;
    if (argc > 1)
      if (Tcl_GetDouble(interp, argv[1], &alpha) != TCL_OK) {
	opserr << "nDMaterial Templated3Dep - invalid alpha " << argv[1] << endln;
	cleanup(argv);
	return 0;		
      }
    PS = new DPPotentialSurface(alpha);
  }

  //2. von-Mises Potential Surface
  if ((strcmp(argv[0],"-VM") == 0) || (strcmp(argv[0],"-vm") == 0)) {
    PS = new VMPotentialSurface();
  }

  // 3. CamClay potential surface
  else if ((strcmp(argv[0],"-CC") == 0) || (strcmp(argv[0],"-cc") == 0)) {

    double mp = 1.2;
    if (argc == 2) {
      if (Tcl_GetDouble(interp, argv[1], &mp) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep - invalid M " << argv[1] << endln;
	return 0;		
      }
    }
    
    // create the object & return it
    PS = new CAMPotentialSurface(mp);
  }


  cleanup(argv);
  return PS;
}

// Function - to read in a stress tensor
int EvaluateStressTensor(ClientData clientData, Tcl_Interp *interp, char *tclString,
			 stresstensor &stress)
{
  int argc;
  char **argv;
    
  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }
    
  //opserr << " argc " << argc << "  argv: " << *argv << endln;
  
  //Found a bug here!! Joey May 15, 2001
  //argc = 1 and only one value is passed in. 
  //But we need 9 values to initialze a tensor
  //double *values = new double[argc];
  double *values = new double[9];
  if (values == 0) {
    opserr << "nDMaterial Template3DEp -ESP stresstenor - out of memory for size: " << argc << endln;
    cleanup(argv);
    return -1;
  }

  int i;
  for (i=0; i<argc; i++) {
    if (Tcl_GetDouble(interp, argv[i], &values[i]) != TCL_OK) {
      opserr << "nDMaterial Templated3Dep -ESP stresstensor - invalid value: " << argv[i] << endln;
      cleanup(argv);
      //delete [] values;
      return -1;
    }


  }

  for (i=1; i<9; i++) {
     values[i] = 0;
     if ((i == 4)|| (i==8))
        values[i] = -1.0 * values[0];
  }
  values[0] = -1.0 * values[0];
  stresstensor newStress(values);
  
  //newStress.printshort("tcl:");

  stress = newStress;
  cleanup(argv);
  delete [] values;
  return 0;
}

// Function - to read in a strain tensor
int EvaluateStrainTensor(ClientData clientData, Tcl_Interp *interp, char *tclString,
			 straintensor &strain)
{
  int argc;
  char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }
  
  double *values = new double[argc];
  if (values == 0) {
    opserr << "nDMaterial Template3DEp -ESP straintensor - out of memory for size: " << argc << endln;
    cleanup(argv);
    return -1;
  }

  for (int i=0; i<argc; i++) {
    if (Tcl_GetDouble(interp, argv[i], &values[i]) != TCL_OK) {
      opserr << "nDMaterial Templated3Dep -ESP straintensor - invalid value: " << argv[i] << endln;
      cleanup(argv);
      //delete [] values;
      return -1;
    }
  }

  straintensor newStrain(values);
  strain = newStrain;
  cleanup(argv);
  delete [] values;
  return 0;
}

// Function - to create an EPState object
EPState *EvaluateEPStateArgs(ClientData clientData, Tcl_Interp *interp, char *tclString)
{
  int argc;
  char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }
  EPState *EPS = 0;

  double Eod, Ed, nu, rho;
  if (argc < 5) {

    cleanup(argv);
    return 0;
  }

  if (Tcl_GetDouble(interp, argv[0], &Eod) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep -EPS - invalid Eod " << argv[1] << endln;
    cleanup(argv);
    return 0;		
  }
  if (Tcl_GetDouble(interp, argv[1], &Ed) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep -EPS - invalid Ed " << argv[1] << endln;
    cleanup(argv);
    return 0;		
  }
  if (Tcl_GetDouble(interp, argv[2], &nu) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep -EPS - invalid nu " << argv[1] << endln;
    cleanup(argv);
    return 0;		
  }
  if (Tcl_GetDouble(interp, argv[3], &rho) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep -EPS - invalid rho " << argv[1] << endln;
  }
  int loc = 4;
  stresstensor stressp;
  straintensor strainp;
  straintensor Estrainp;
  straintensor Pstrainp;

  int NoS = 0;
  int NoD = 0;
  double *scalars = 0;
  stresstensor *tensors = 0;

  // switch on remaining args;
  while (loc < argc) {
    
    if ((strcmp(argv[loc],"-stressP") == 0) || (strcmp(argv[loc],"-stressp") == 0)) {
      if (EvaluateStressTensor(clientData, interp, argv[loc+1], stressp) < 0) {
	opserr << "nDMaterial Templated3Dep -EPS - invalid stressp " << argv[loc+1] << endln;
	cleanup(argv);
	return 0;
      }
      loc+=2;
    }
    else if ((strcmp(argv[loc],"-strainP") == 0) || (strcmp(argv[loc],"-strainp") == 0)) {
      if (EvaluateStrainTensor(clientData, interp, argv[loc+1], strainp) < 0) {
	opserr << "nDMaterial Templated3Dep -EPS - invalid strainp " << argv[loc+1] << endln;
	cleanup(argv);
	return 0;
      }
    }    
    else if ((strcmp(argv[loc],"-EstrainP") == 0) || (strcmp(argv[loc],"-Estrainp") == 0)) {
      if (EvaluateStrainTensor(clientData, interp, argv[loc+1], Estrainp) < 0) 
	opserr << "nDMaterial Templated3Dep -EPS - invalid Estrainp " << argv[loc+1] << endln;
	cleanup(argv);
	return 0;
      loc+=2;
    }
    else if ((strcmp(argv[loc],"-PstrainP") == 0) || (strcmp(argv[loc],"-Estrainp") == 0)) {
      if (EvaluateStrainTensor(clientData, interp, argv[loc+1], Pstrainp) < 0) {
	opserr << "nDMaterial Templated3Dep -EPS - invalid Pstrainp " << argv[loc+1] << endln;
	cleanup(argv);
	return 0;
      }
      loc+=2;
    }    
    else if ((strcmp(argv[loc],"-NOS") == 0) || (strcmp(argv[loc],"-nos") == 0)) {
      if (Tcl_GetInt(interp, argv[loc+1], &NoS) != TCL_OK) {
	opserr << "nDMaterial Templated3Dep -EPS - invalid NOS " << argv[loc+1] << endln;
	return 0;	
      }
      loc+=2;

      if (NoS > 0) {
	scalars = new double[NoS];
	for (int i=0; i<NoS; i++) {
	  if (Tcl_GetDouble(interp, argv[loc++], &scalars[i]) != TCL_OK) {
	    opserr << "nDMaterial Templated3Dep -EPS - invalid scalar: " << argv[loc-1] << endln;
	    return 0;		  
	  }
	}
      } else 
	NoS = 0;
    }
    else if ((strcmp(argv[loc],"-NOD") == 0) || (strcmp(argv[loc],"-nod") == 0)) {
      if (Tcl_GetInt(interp, argv[loc+1], &NoD) != TCL_OK) {
	opserr << "nDMaterial Templated3Dep -EPS - invalid NOD: " << argv[loc+1] << endln;
	return 0;
      }
      loc+=2;
      
      if (NoD > 0) {
	tensors = new stresstensor[NoD];
	if (tensors == 0) {
	  opserr << "nDMaterial Templated3Dep -EPS - invalid NOD " << argv[loc+1] << endln;
	  return 0;
	}
      } else 
	NoD = 0;
    } else
      loc++;
  }

  EPS = new EPState(Eod, Ed, nu, rho, stressp, strainp, Estrainp, Pstrainp, 
		    NoS, scalars, NoD, tensors);
  
  if (EPS == 0) {
    opserr << "nDMaterial Templated3Dep -EPS - out of memory\n";
  }

  cleanup(argv);

  //if (NoD > 0) 
  //   for (int i=0; i<=NoD; i++) {
  //     delete tensors;
  //   }
  if (NoS > 0) 
     delete [] scalars;  

  return EPS;
}



// Function - to create an EvolutionLaw_S object
EvolutionLaw_S *EvaluateEvolutionLawSArgs(ClientData clientData, Tcl_Interp *interp, char *tclString)
{
  int argc;
  char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }

  //1. Linear scalar evolution law: eq
  EvolutionLaw_S *ELS = 0;
  if ((strcmp(argv[0],"-Leq") == 0) || (strcmp(argv[0],"-leq") == 0)) {
    double alpha = 0.0;
    if (argc > 1)
      if (Tcl_GetDouble(interp, argv[1], &alpha) != TCL_OK) {
	opserr << "nDMaterial Templated3Dep - invalid alpha " << argv[1] << endln;
	cleanup(argv);
	return 0;		
      }
    ELS = new EvolutionLaw_L_Eeq(alpha);
  }

  //2. Nonlinear scalar evolution law: p for Cam clay
  //EvolutionLaw_S *ELS = 0;
  if ((strcmp(argv[0],"-NLp") == 0) || (strcmp(argv[0],"-nlp") == 0)) {
    double eod = 0.65, lambdad=0.19, kappad=0.06;
    if (argc ==3 ) {
      if (Tcl_GetDouble(interp, argv[1], &eod) != TCL_OK) {
	opserr << "nDMaterial Templated3Dep - invalid eo " << argv[1] << endln;
				
	//cleanup(argv);
	return 0;		
      }
      if (Tcl_GetDouble(interp, argv[2], &lambdad) != TCL_OK) {
	opserr << "nDMaterial Templated3Dep - invalid lambda " << argv[2] << endln;
	  
	//cleanup(argv);
	return 0;		
      }    
      if (Tcl_GetDouble(interp, argv[3], &kappad) != TCL_OK) {
	opserr << "nDMaterial Templated3Dep - invalid kappa " << argv[3] << endln;
	  
	cleanup(argv);
	return 0;		
      }

    }
    ELS = new EvolutionLaw_NL_Ep(eod, lambdad, kappad);
  }

  cleanup(argv);
  return ELS;
}

// Function - to create an EvolutionLaw_T object
EvolutionLaw_T *EvaluateEvolutionLawTArgs(ClientData clientData, Tcl_Interp *interp, char *tclString)
{
  int argc;
  char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }

  EvolutionLaw_T *ELT = 0;
  if ((strcmp(argv[0],"-Linear") == 0) || (strcmp(argv[0],"-linear") == 0)) {
    double alpha = 0.0;
    if (argc > 1)
      if (Tcl_GetDouble(interp, argv[1], &alpha) != TCL_OK) {
	opserr << "nDMaterial Templated3Dep - invalid alpha " << argv[1] << endln;
				
	cleanup(argv);
	return 0;		
      }
    ELT = new EvolutionLaw_L_Eij(alpha);
  }

  cleanup(argv);
  return ELT;
}


