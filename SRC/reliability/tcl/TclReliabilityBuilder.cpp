/* ****************************************************************** **

**    OpenSees - Open System for Earthquake Engineering Simulation    **

**          Pacific Earthquake Engineering Research Center            **

**                                                                    **

**                                                                    **

** (C) Copyright 2001, The Regents of the University of California    **

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

** Reliability module developed by:                                   **

**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **

**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **

**                                                                    **

** ****************************************************************** */

                                                                        

// $Revision: 1.2 $

// $Date: 2001-06-14 07:19:33 $

// $Source: /usr/local/cvs/OpenSees/SRC/reliability/tcl/TclReliabilityBuilder.cpp,v $





//

// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000

// Revised: haukaas 06/00 (core code)

//			haukaas 06/01 (made part of official OpenSees)

//



#include <stdlib.h>

#include <string.h>

#include <iostream.h>

#include <fstream.h>

#include <iomanip.h>



#include <Matrix.h>

#include <Vector.h>

#include <ID.h>

#include <ArrayOfTaggedObjects.h>



#include <Domain.h>



#include <ReliabilityDomain.h>

#include <RandomVariable.h>

#include <CorrelationCoefficient.h>

#include <LimitStateFunction.h>

#include <RandomVariablePositioner.h>

#include <NormalRV.h>

#include <LognormalRV.h>

#include <GammaRV.h>

#include <ShiftedExponentialRV.h>

#include <ShiftedRayleighRV.h>

#include <ExponentialRV.h>

#include <RayleighRV.h>

#include <UniformRV.h>

#include <BetaRV.h>

#include <Type1LargestValueRV.h>

#include <Type1SmallestValueRV.h>

#include <Type2LargestValueRV.h>

#include <Type3SmallestValueRV.h>

#include <ChiSquareRV.h>

#include <GumbelRV.h>

#include <WeibullRV.h>

#include <LaplaceRV.h>

#include <ParetoRV.h>

#include <GFunEvaluator.h>

#include <SensitivityEvaluator.h>

#include <StepSizeRule.h>

#include <SearchDirection.h>

#include <XuTransformation.h>

#include <NatafXuTransformation.h>

#include <FindDesignPoint.h>

#include <ReliabilityAnalysis.h>

#include <HLRFSearchDirection.h>

#include <ArmijoRule.h>

#include <FixedStepSizeRule.h>

#include <OpenSeesGFunEvaluator.h>

#include <OpenSeesSensitivityEvaluator.h>

#include <BasicGFunEvaluator.h>

#include <SensitivityByFiniteDifference.h>

#include <SensitivityByFiniteDifference.h>

#include <SearchWithStepSizeAndStepDirection.h>

#include <FORMAnalysis.h>

#include <SimulationAnalysis.h>

#include <RandomNumberGenerator.h>

#include <CStdLibRandGenerator.h>

#include <FindCurvatures.h>

#include <FirstPrincipalCurvature.h>

#include <CurvaturesBySearchAlgorithm.h>

#include <SORMAnalysis.h>

#include <SystemAnalysis.h>





#include <TclReliabilityBuilder.h>



//

// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER

//

static ReliabilityDomain *theReliabilityDomain = 0;

static Domain *theStructuralDomain = 0;



static GFunEvaluator *theGFunEvaluator = 0;

static SensitivityEvaluator *theSensitivityEvaluator = 0;

static StepSizeRule *theStepSizeRule = 0;

static SearchDirection *theSearchDirection = 0;

static XuTransformation *theXuTransformation = 0;

static RandomNumberGenerator *theRandomNumberGenerator = 0;

static FindDesignPoint *theFindDesignPoint = 0;

static FindCurvatures *theFindCurvatures = 0;

static ReliabilityAnalysis *theFORMAnalysis = 0;

static ReliabilityAnalysis *theSORMAnalysis = 0;

static ReliabilityAnalysis *theSimulationAnalysis = 0;

static SystemAnalysis *theSystemAnalysis = 0;



// 

// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER

//

int TclReliabilityModelBuilder_addRandomVariable(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addCorrelate(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addLimitState(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addRandomVariablePositioner(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addXuTransformation(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addRandomNumberGenerator(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addSearchDirection(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addStepSizeRule(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addgFunEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addSensitivityEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addFindDesignPoint(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addFindCurvatures(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addFORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addSORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addSystemAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_addSimulationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int TclReliabilityModelBuilder_analyzeReliability(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);







//

// CLASS CONSTRUCTOR & DESTRUCTOR

//



// constructor: the constructor will add certain commands to the interpreter

TclReliabilityBuilder::TclReliabilityBuilder(Domain &passedDomain, Tcl_Interp *interp)

{



  // call Tcl_CreateCommand for class specific commands

  Tcl_CreateCommand(interp, "randomVariable",	TclReliabilityModelBuilder_addRandomVariable,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "correlate", TclReliabilityModelBuilder_addCorrelate,

		    (ClientData)NULL, NULL); 

  Tcl_CreateCommand(interp, "limitState", TclReliabilityModelBuilder_addLimitState,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "randomVariablePositioner", 

		    TclReliabilityModelBuilder_addRandomVariablePositioner,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "findDesignPoint",	TclReliabilityModelBuilder_addFindDesignPoint,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "gFunEvaluator",	TclReliabilityModelBuilder_addgFunEvaluator,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "sensitivityEvaluator", 

		    TclReliabilityModelBuilder_addSensitivityEvaluator,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "stepSizeRule", 

		    TclReliabilityModelBuilder_addStepSizeRule,	

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "searchDirection",	TclReliabilityModelBuilder_addSearchDirection,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "xuTransformation",	TclReliabilityModelBuilder_addXuTransformation,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "findCurvatures",	TclReliabilityModelBuilder_addFindCurvatures,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "randomNumberGenerator",  

		    TclReliabilityModelBuilder_addRandomNumberGenerator,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "addFORMAnalysis",TclReliabilityModelBuilder_addFORMAnalysis,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "addSORMAnalysis",TclReliabilityModelBuilder_addSORMAnalysis,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "addSystemAnalysis",TclReliabilityModelBuilder_addSystemAnalysis,

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "addSimulationAnalysis",

		    TclReliabilityModelBuilder_addSimulationAnalysis,	

		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "analyzeReliability",  

		    TclReliabilityModelBuilder_analyzeReliability, 

		    (ClientData)NULL, NULL);





  // set the static pointers in this file

  theStructuralDomain	= &passedDomain;

  theReliabilityDomain	= new ReliabilityDomain();



//  ResultsContainer *theResultsContainer = new ResultsContainer();

//  theReliabilityDomain->addResultsContainer(theResultsContainer);



}



TclReliabilityBuilder::~TclReliabilityBuilder()

{

  cerr << "Terje needs to write code: TclReliabilityBuilder::~TclReliabilityBuilder()\n";

}





//

// CLASS METHODS

//



/*

int 

TclReliabilityBuilder::buildFE_Model(void)

{

  // does nothing

  return 0;

}

*/



ReliabilityDomain *

TclReliabilityBuilder::getReliabilityDomain()

{

	return theReliabilityDomain;

}



//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addRandomVariable(ClientData clientData,Tcl_Interp *interp,int argc,char **argv)

{

  RandomVariable *theRandomVariable = 0;

  int tag;

  double mean;

  double stdv;

  double startPt;

  double parameter1;

  double parameter2;

  double parameter3;

  double parameter4;

  int numberOfArguments = argc;



  if (numberOfArguments==5)  {   // (Use mean/stdv WITHOUT startPt)



	  // GET INPUT PARAMETER (integer)

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {

		cerr << "WARNING invalid input: tag \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[3], &mean) != TCL_OK) {

		cerr << "WARNING invalid input: mean \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[4], &stdv) != TCL_OK) {

		cerr << "WARNING invalid input: stdv \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	  if (strcmp(argv[1],"normal") == 0) {

		  theRandomVariable = new NormalRV(tag, mean, stdv);

	  }

	  else if (strcmp(argv[1],"lognormal") == 0) {

		  if (mean < 0.0  ||  stdv < 0.0) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			  theRandomVariable = new LognormalRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"gamma") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			  theRandomVariable = new GammaRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"shiftedExponential") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ShiftedExponentialRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"shiftedRayleigh") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ShiftedRayleighRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"exponential") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ExponentialRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"rayleigh") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new RayleighRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"uniform") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new UniformRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"beta") == 0) {

		  cerr << "WARNING: 'Beta' type random variable: use parameters to create!\n";

	  }

	  else if (strcmp(argv[1],"type1LargestValue") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new Type1LargestValueRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"type1SmallestValue") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new Type1SmallestValueRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"type2LargestValue") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new Type2LargestValueRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"type3SmallestValue") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new Type3SmallestValueRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"chiSquare") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ChiSquareRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"gumbel") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new GumbelRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"weibull") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new WeibullRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"laplace") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new LaplaceRV(tag, mean, stdv);

		  }

	  }

	  else if (strcmp(argv[1],"pareto") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ParetoRV(tag, mean, stdv);

		  }

	  }

	  else {

		cerr << "WARNING unrecognized type of random variable \n";

		cerr << "random variable: " << tag << endl;

	  }



	  if (theRandomVariable == 0) {

		cerr << "WARNING could not create random variable" << tag << endl;

		return TCL_ERROR;

	  }

  }



  if (numberOfArguments==6)  {   // (Use mean/stdv AND startPt)



	  // GET INPUT PARAMETER (integer)

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {

		cerr << "WARNING invalid input: tag \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[3], &mean) != TCL_OK) {

		cerr << "WARNING invalid input: mean \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[4], &stdv) != TCL_OK) {

		cerr << "WARNING invalid input: stdv \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[5], &startPt) != TCL_OK) {

		cerr << "WARNING invalid input: startPt \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	  if (strcmp(argv[1],"normal") == 0) {

		  theRandomVariable = new NormalRV(tag, mean, stdv, startPt);

	  }

	  else if (strcmp(argv[1],"lognormal") == 0) {

		  if (mean < 0.0  ||  stdv < 0.0) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			  theRandomVariable = new LognormalRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"gamma") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			  theRandomVariable = new GammaRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"shiftedExponential") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ShiftedExponentialRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"shiftedRayleigh") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ShiftedRayleighRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"exponential") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ExponentialRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"rayleigh") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new RayleighRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"uniform") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new UniformRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"beta") == 0) {

		  cerr << "WARNING: 'Beta' type random variable: use parameters to create!\n";

	  }

	  else if (strcmp(argv[1],"type1LargestValue") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new Type1LargestValueRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"type1SmallestValue") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new Type1SmallestValueRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"type2LargestValue") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new Type2LargestValueRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"type3SmallestValue") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new Type3SmallestValueRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"chiSquare") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ChiSquareRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"gumbel") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new GumbelRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"weibull") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new WeibullRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"laplace") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new LaplaceRV(tag, mean, stdv, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"pareto") == 0) {

		  if ( stdv < 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else {

			theRandomVariable = new ParetoRV(tag, mean, stdv, startPt);

		  }

	  }

	  else {

		cerr << "WARNING unrecognized type of random variable \n";

		cerr << "random variable: " << tag << endl;

	  }



	  if (theRandomVariable == 0) {

		cerr << "WARNING could not create random variable" << tag << endl;

		return TCL_ERROR;

	  }

  }



  if (numberOfArguments==7)  {  // (Use parameters WITHOUT startPt)

		// GET INPUT PARAMETER (integer)

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {

		cerr << "WARNING invalid input: tag \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[3], &parameter1) != TCL_OK) {

		cerr << "WARNING invalid input: parameter1 \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[4], &parameter2) != TCL_OK) {

		cerr << "WARNING invalid input: parameter2 \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[5], &parameter3) != TCL_OK) {

		cerr << "WARNING invalid input: parameter3 \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[6], &parameter4) != TCL_OK) {

		cerr << "WARNING invalid input: parameter4 \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	  if (strcmp(argv[1],"normal") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new NormalRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"lognormal") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new LognormalRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"gamma") == 0) {

		  if ( parameter1 <= 0.0  || parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new GammaRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"shiftedExponential") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ShiftedExponentialRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"shiftedRayleigh") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ShiftedRayleighRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"exponential") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ExponentialRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"rayleigh") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new RayleighRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"uniform") == 0) {

		  if ( parameter1 >= parameter2 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new UniformRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"beta") == 0) {

		  if ( parameter1 >= parameter2  ||  parameter3 <= 0.0  || parameter4 <= 0.0  ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new BetaRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"type1LargestValue") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new Type1LargestValueRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"type1SmallestValue") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new Type1SmallestValueRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"type2LargestValue") == 0) {

		  if ( parameter1 <= 0.0  || parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new Type2LargestValueRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"type3SmallestValue") == 0) {

		  if ( parameter2 <= 0.0  ||  parameter3 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new Type3SmallestValueRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"chiSquare") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ChiSquareRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"gumbel") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new GumbelRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"weibull") == 0) {

		  if ( parameter1 <= 0.0  ||  parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new WeibullRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"laplace") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new LaplaceRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else if (strcmp(argv[1],"pareto") == 0) {

		  if ( parameter1 <= 0.0  ||  parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ParetoRV(tag, parameter1, parameter2, parameter3, parameter4);

		  }

	  }

	  else {

		cerr << "WARNING unrecognized type of random variable \n";

		cerr << "random variable: " << tag << endl;

	  }



	  if (theRandomVariable == 0) {

		cerr << "WARNING could not create random variable" << tag << endl;

		return TCL_ERROR;

	  }

  }





  if (numberOfArguments==8)  {  // (Use parameters AND startPt)

		// GET INPUT PARAMETER (integer)

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {

		cerr << "WARNING invalid input: tag \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[3], &parameter1) != TCL_OK) {

		cerr << "WARNING invalid input: parameter1 \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[4], &parameter2) != TCL_OK) {

		cerr << "WARNING invalid input: parameter2 \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[5], &parameter3) != TCL_OK) {

		cerr << "WARNING invalid input: parameter3 \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[6], &parameter4) != TCL_OK) {

		cerr << "WARNING invalid input: parameter4 \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (double)

	  if (Tcl_GetDouble(interp, argv[7], &startPt) != TCL_OK) {

		cerr << "WARNING invalid input: startPt \n";

		return TCL_ERROR;

	  }



	  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	  if (strcmp(argv[1],"normal") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new NormalRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"lognormal") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new LognormalRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"gamma") == 0) {

		  if ( parameter1 <= 0.0  || parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new GammaRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"shiftedExponential") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ShiftedExponentialRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"shiftedRayleigh") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ShiftedRayleighRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"exponential") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ExponentialRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"rayleigh") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new RayleighRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"uniform") == 0) {

		  if ( parameter1 >= parameter2 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new UniformRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"beta") == 0) {

		  if ( parameter1 >= parameter2  ||  parameter3 <= 0.0  || parameter4 <= 0.0  ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new BetaRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"type1LargestValue") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new Type1LargestValueRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"type1SmallestValue") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new Type1SmallestValueRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"type2LargestValue") == 0) {

		  if ( parameter1 <= 0.0  || parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new Type2LargestValueRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"type3SmallestValue") == 0) {

		  if ( parameter2 <= 0.0  ||  parameter3 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new Type3SmallestValueRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"chiSquare") == 0) {

		  if ( parameter1 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ChiSquareRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"gumbel") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new GumbelRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"weibull") == 0) {

		  if ( parameter1 <= 0.0  ||  parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new WeibullRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"laplace") == 0) {

		  if ( parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new LaplaceRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else if (strcmp(argv[1],"pareto") == 0) {

		  if ( parameter1 <= 0.0  ||  parameter2 <= 0.0 ) {

			  cerr << "Invalid input parameters to random variable" << tag << endl;

		  }

		  else  {

			theRandomVariable = new ParetoRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);

		  }

	  }

	  else {

		cerr << "WARNING unrecognized type of random variable \n";

		cerr << "random variable: " << tag << endl;

	  }



	  if (theRandomVariable == 0) {

		cerr << "WARNING could not create random variable" << tag << endl;

		return TCL_ERROR;

	  }

  }





  // ADD THE OBJECT TO THE DOMAIN

  if (theReliabilityDomain->addRandomVariable(theRandomVariable) == false) {

	cerr << "WARNING failed to add random variable to the domain (wrong number of arguments?)\n";

	cerr << "random variable: " << tag << endl;

	delete theRandomVariable; // otherwise memory leak

	return TCL_ERROR;

  }



  return TCL_OK;

}

					   





//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addCorrelate(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

  CorrelationCoefficient *theCorrelationCoefficient = 0;

  int tag;

  int rv1;

  int rv2;

  double correlationValue;



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {

	cerr << "WARNING invalid input: tag \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[2], &rv1) != TCL_OK) {

	cerr << "WARNING invalid input: rv1 \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[3], &rv2) != TCL_OK) {

	cerr << "WARNING invalid input: rv2 \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (double)

  if (Tcl_GetDouble(interp, argv[4], &correlationValue) != TCL_OK) {

	cerr << "WARNING invalid input: correlationValue \n";

	return TCL_ERROR;

  }



  // CREATE THE OBJECT

  theCorrelationCoefficient = new CorrelationCoefficient(tag, rv1, rv2, correlationValue);



  if (theCorrelationCoefficient == 0) {

	cerr << "WARNING ran out of memory creating correlation coefficient \n";

	cerr << "correlation coefficient: " << tag << endl;

	return TCL_ERROR;

  }



  // ADD THE OBJECT TO THE DOMAIN

  if (theReliabilityDomain->addCorrelationCoefficient(theCorrelationCoefficient) == false) {

	cerr << "WARNING failed to add correlation coefficient to the domain\n";

	cerr << "correlation coefficient: " << tag << endl;

	delete theCorrelationCoefficient; // otherwise memory leak

	return TCL_ERROR;

  }



  return TCL_OK;

}







//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addLimitState(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

  LimitStateFunction *theLimitStateFunction = 0;

  int tag;



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {

	cerr << "WARNING invalid input: tag \n";

	return TCL_ERROR;

  }

  

  // CREATE THE OBJECT (passing on argv[2])

  theLimitStateFunction = new LimitStateFunction(tag, argv[2]);

  if (theLimitStateFunction == 0) {

	cerr << "WARNING ran out of memory creating limit-state function \n";

	cerr << "limit-state function: " << tag << endl;

	return TCL_ERROR;

  }



  // ADD THE OBJECT TO THE DOMAIN

  if (theReliabilityDomain->addLimitStateFunction(theLimitStateFunction) == false) {

	cerr << "WARNING failed to add limit-state function to the domain\n";

	cerr << "limit-state function: " << tag << endl;

	delete theLimitStateFunction; // otherwise memory leak

	return TCL_ERROR;

  }



  return TCL_OK;









/*



  LimitStateFunction *theLimitStateFunction = 0;

  int tag;

  int node;

  int dof;

  double displacementLimit;



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {

	cerr << "WARNING invalid input: tag \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[3], &node) != TCL_OK) {

	cerr << "WARNING invalid input: node \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[4], &dof) != TCL_OK) {

	cerr << "WARNING invalid input: dof \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (double)

  if (Tcl_GetDouble(interp, argv[5], &displacementLimit) != TCL_OK) {

	cerr << "WARNING invalid input: displacementLimit \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT

  if (strcmp(argv[1],"disp") == 0) {

	  theLimitStateFunction = new LimitStateFunction(tag, node, dof, displacementLimit);

  }

  else {

	cerr << "WARNING unrecognized type of limit-state function \n";

	cerr << "limit-state function: " << tag << endl;

  }



  if (theLimitStateFunction == 0) {

	cerr << "WARNING ran out of memory creating limit-state function \n";

	cerr << "limit-state function: " << tag << endl;

	return TCL_ERROR;

  }



  // ADD THE OBJECT TO THE DOMAIN

  if (theReliabilityDomain->addLimitStateFunction(theLimitStateFunction) == false) {

	cerr << "WARNING failed to add limit-state function to the domain\n";

	cerr << "limit-state function: " << tag << endl;

	delete theLimitStateFunction; // otherwise memory leak

	return TCL_ERROR;

  }



  return TCL_OK;

*/

}







//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addRandomVariablePositioner(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	RandomVariablePositioner *theRandomVariablePositioner = 0;

	int tag;

	int rvNumber;

	int tagOfObject;

	DomainComponent *theObject;





	// READ THE TAG NUMBER

	if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {

		cerr << "WARNING invalid input: tag \n";

		return TCL_ERROR;

	}



	// READ THE RANDOM VARIABLE NUMBER

	if (Tcl_GetInt(interp, argv[2], &rvNumber) != TCL_OK) {

		cerr << "WARNING invalid input: rvNumber \n";

		return TCL_ERROR;

	}



	// IF UNCERTAIN *ELEMENT* PROPERTY

	if (strcmp(argv[3],"element") == 0) {



		if (Tcl_GetInt(interp, argv[4], &tagOfObject) != TCL_OK) {

			cerr << "WARNING invalid input: tagOfObject \n";

			return TCL_ERROR;

		}



		theObject = (DomainComponent *)theStructuralDomain->getElement(tagOfObject);


		theRandomVariablePositioner = new RandomVariablePositioner(tag,rvNumber,theObject,&argv[5],argc-5);

	}





	// IF UNCERTAIN *LOAD*

	else if (strcmp(argv[3],"loadPattern") == 0) {



		if (Tcl_GetInt(interp, argv[4], &tagOfObject) != TCL_OK) {

			cerr << "WARNING invalid input: tagOfObject \n";

			return TCL_ERROR;

		}

		theObject = (DomainComponent *)theStructuralDomain->getLoadPattern(tagOfObject);


		theRandomVariablePositioner = new RandomVariablePositioner(tag,rvNumber,theObject,&argv[5],argc-5);

	}





	// ADD THE RANDOMVARIABLEPOSITIONER TO THE DOMAIN

	if (theReliabilityDomain->addRandomVariablePositioner(theRandomVariablePositioner) == false) {

		cerr << "WARNING failed to add random variable identificator to the domain\n";

		cerr << "randomvariableID: " << tag << endl;

		delete theRandomVariablePositioner; // otherwise memory leak

		return TCL_ERROR;

	}



	return TCL_OK;















/*

// THE OLD VERSION OF THIS FUNCTION:



  RandomVariablePositioner *theRandomVariablePositioner = 0;

  int tag;

  int rvNumber;

  int typeOfObject;

  int tagOfObject;

  int typeOfParameterInObject;



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {

	cerr << "WARNING invalid input: tag \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[2], &rvNumber) != TCL_OK) {

	cerr << "WARNING invalid input: rvNumber \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[3], &typeOfObject) != TCL_OK) {

	cerr << "WARNING invalid input: typeOfObject \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[4], &tagOfObject) != TCL_OK) {

	cerr << "WARNING invalid input: tagOfObject \n";

	return TCL_ERROR;

  }



  // GET INPUT PARAMETER (integer)

  if (Tcl_GetInt(interp, argv[5], &typeOfParameterInObject) != TCL_OK) {

	cerr << "WARNING invalid input: typeOfParameterInObject \n";

	return TCL_ERROR;

  }





  // CREATE THE OBJECT

  theRandomVariablePositioner = new RandomVariablePositioner(tag, rvNumber, typeOfObject, tagOfObject, typeOfParameterInObject);



  if (theRandomVariablePositioner == 0) {

	cerr << "WARNING ran out of memory creating random variable identificator \n";

	cerr << "randomVariableID: " << tag << endl;

	return TCL_ERROR;

  }



  // ADD THE OBJECT TO THE DOMAIN

  if (theReliabilityDomain->addRandomVariablePositioner(theRandomVariablePositioner) == false) {

	cerr << "WARNING failed to add random variable identificator to the domain\n";

	cerr << "randomvariableID: " << tag << endl;

	delete theRandomVariablePositioner; // otherwise memory leak

	return TCL_ERROR;

  }



  return TCL_OK;

*/

}









//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addRandomNumberGenerator(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT

  if (strcmp(argv[1],"CStdLibRandGenerator") == 0) {

	  theRandomNumberGenerator = new CStdLibRandGenerator();

  }

  else {

	cerr << "WARNING unrecognized type of RandomNumberGenerator \n";

	return TCL_ERROR;

  }



  if (theRandomNumberGenerator == 0) {

	cerr << "WARNING could not create theRandomNumberGenerator \n";

	return TCL_ERROR;

  }

  return TCL_OK;

}







//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addXuTransformation(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT

  if (strcmp(argv[1],"NatafXuTransformation") == 0) {

	  theXuTransformation = new NatafXuTransformation(theReliabilityDomain);

  }

  else {

	cerr << "WARNING unrecognized type of XuTransformation \n";

	return TCL_ERROR;

  }



  if (theXuTransformation == 0) {

	cerr << "WARNING could not create theXuTransformation \n";

	return TCL_ERROR;

  }

  return TCL_OK;

}









//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addSearchDirection(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	if (strcmp(argv[1],"HLRFSearchDirection") == 0) {

		theSearchDirection = new HLRFSearchDirection();

	}

	else {

		cerr << "WARNING unrecognized type of SearchDirection \n";

		return TCL_ERROR;

	}



	if (theSearchDirection == 0) {

		cerr << "WARNING could not create theSearchDirection \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}







//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addStepSizeRule(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	if (strcmp(argv[1],"ArmijoRule") == 0) {

		// Check that the necessary ingredients are present

		if (theGFunEvaluator == 0 ) {

			cerr << "Need theGFunEvaluator before a ArmijoRule can be created" << endl;

			return TCL_ERROR;

		}

		if (theGFunEvaluator == 0 ) {

			cerr << "Need theXuTransformation before a ArmijoRule can be created" << endl;

			return TCL_ERROR;

		}



		theStepSizeRule = new ArmijoRule(theGFunEvaluator,theXuTransformation);

	}

	else if (strcmp(argv[1],"FixedStepSizeRule") == 0) {



		double stepSize;



		// GET INPUT PARAMETER (double)

		if (Tcl_GetDouble(interp, argv[2], &stepSize) != TCL_OK) {

			cerr << "WARNING invalid input: stepSize \n";

			return TCL_ERROR;

		}



		theStepSizeRule = new FixedStepSizeRule(stepSize);

	}

	else {

		cerr << "WARNING unrecognized type of StepSizeRule \n";

		return TCL_ERROR;

	}



	if (theStepSizeRule == 0) {

		cerr << "WARNING could not create theStepSizeRule \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}







//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addgFunEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	if (strcmp(argv[1],"OpenSeesGFunEvaluator") == 0) {



	  int numIter;

	  // GET INPUT PARAMETER (integer)

	  if (Tcl_GetInt(interp, argv[2], &numIter) != TCL_OK) {

		cerr << "WARNING invalid input: numIter \n";

		return TCL_ERROR;

	  }



		theGFunEvaluator = new OpenSeesGFunEvaluator(numIter, interp, theReliabilityDomain);

	}

	else if (strcmp(argv[1],"BasicGFunEvaluator") == 0) {

		theGFunEvaluator = new BasicGFunEvaluator(interp, theReliabilityDomain);

	}

	else {

		cerr << "WARNING unrecognized type of GFunEvaluator \n";

		return TCL_ERROR;

	}

	

	if (theGFunEvaluator == 0) {

		cerr << "WARNING could not create the theGFunEvaluator \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}











//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addSensitivityEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	if (strcmp(argv[1],"SensitivityByFiniteDifference") == 0) {



		// Check that the necessary ingredients are present

		if (theGFunEvaluator == 0 ) {

			cerr << "Need theGFunEvaluator before a SensitivityByFiniteDifference can be created" << endl;

			return TCL_ERROR;

		}



		theSensitivityEvaluator = new SensitivityByFiniteDifference(theGFunEvaluator, theReliabilityDomain);

	}

	else if (strcmp(argv[1],"OpenSeesSensitivityEvaluator") == 0) {



		theSensitivityEvaluator = new OpenSeesSensitivityEvaluator(interp, theReliabilityDomain);

	}

	else {

		cerr << "WARNING unrecognized type of SensitivityEvaluator \n";

		return TCL_ERROR;

	}



	if (theSensitivityEvaluator == 0) {

		cerr << "WARNING could not create theSensitivityEvaluator \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}





//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addFindCurvatures(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	if (strcmp(argv[1],"onlyFirstPrincipal") == 0) {



		theFindCurvatures = new FirstPrincipalCurvature();

	}

	else if (strcmp(argv[1],"bySearchAlgorithm") == 0) {



		// Check that the necessary ingredients are present

		if (theFindDesignPoint == 0 ) {

			cerr << "Need theFindDesignPoint before a CurvaturesBySearchAlgorithm can be created" << endl;

			return TCL_ERROR;

		}



		int numberOfCurvatures;



		// GET INPUT PARAMETER (integer)

		if (Tcl_GetInt(interp, argv[2], &numberOfCurvatures) != TCL_OK) {

			cerr << "WARNING invalid input: numberOfCurvatures \n";

			return TCL_ERROR;

		}



		theFindCurvatures = new CurvaturesBySearchAlgorithm(numberOfCurvatures,theFindDesignPoint);

	}

	else {

		cerr << "WARNING unrecognized type of FindCurvatures \n";

		return TCL_ERROR;

	}



	if (theSensitivityEvaluator == 0) {

		cerr << "WARNING could not create theFindCurvatures \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}







//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addFindDesignPoint(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{







	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	if (strcmp(argv[1],"SearchWithStepSizeAndStepDirection") == 0) {



	int maxNumIter;

	double e1;

	double e2;



	// GET INPUT PARAMETER (integer)

	if (Tcl_GetInt(interp, argv[2], &maxNumIter) != TCL_OK) {

		cerr << "WARNING invalid input: maxNumIter \n";

		return TCL_ERROR;

	}



	// GET INPUT PARAMETER (double)

	if (Tcl_GetDouble(interp, argv[3], &e1) != TCL_OK) {

		cerr << "WARNING invalid input: e1 \n";

		return TCL_ERROR;

	}



	// GET INPUT PARAMETER (double)

	if (Tcl_GetDouble(interp, argv[4], &e2) != TCL_OK) {

		cerr << "WARNING invalid input: e2 \n";

		return TCL_ERROR;

	}



	// Check that the necessary ingredients are present

	if (theGFunEvaluator == 0 ) {

		cerr << "Need theGFunEvaluator before a FindDesignPoint can be created" << endl;

		return TCL_ERROR;

	}

	if (theSensitivityEvaluator == 0 ) {

		cerr << "Need theSensitivityEvaluator before a FindDesignPoint can be created" << endl;

		return TCL_ERROR;

	}

	if (theStepSizeRule == 0 ) {

		cerr << "Need theStepSizeRule before a FindDesignPoint can be created" << endl;

		return TCL_ERROR;

	}

	if (theSearchDirection == 0 ) {

		cerr << "Need theSearchDirection before a FindDesignPoint can be created" << endl;

		return TCL_ERROR;

	}

	if (theXuTransformation == 0 ) {

		cerr << "Need theXuTransformation before a FindDesignPoint can be created" << endl;

		return TCL_ERROR;

	}





	theFindDesignPoint = new SearchWithStepSizeAndStepDirection(

					maxNumIter, 

					e1, 

					e2,

					theGFunEvaluator,

					theSensitivityEvaluator,

					theStepSizeRule,

					theSearchDirection,

					theXuTransformation);

	}

	else {

		cerr << "WARNING unrecognized type of FindDesignPoint Algorithm \n";

		return TCL_ERROR;

	}



	if (theFindDesignPoint == 0) {

		cerr << "WARNING could not create theFindDesignPoint \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}









//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addFORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	if (theFindDesignPoint == 0 ) {

		cerr << "Need theFindDesignPoint before a FORMAnalysis can be created" << endl;

		return TCL_ERROR;

	}



	theFORMAnalysis 

		= new FORMAnalysis(theReliabilityDomain, theFindDesignPoint );



	if (theFORMAnalysis == 0) {

		cerr << "WARNING could not create theFORMAnalysis \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}



//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addSORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	if (theFindCurvatures == 0 ) {

		cerr << "Need theFindCurvatures before a SORMAnalysis can be created" << endl;

		return TCL_ERROR;

	}



	theSORMAnalysis 

		= new SORMAnalysis(theReliabilityDomain, theFindCurvatures );



	if (theSORMAnalysis == 0) {

		cerr << "WARNING could not create theSORMAnalysis \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}



//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addSystemAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{



	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT

	if (strcmp(argv[1],"allInSeries") == 0) {



		theSystemAnalysis = new SystemAnalysis(theReliabilityDomain);



		if (theSystemAnalysis == 0) {

			cerr << "WARNING could not create theSystemAnalysis \n";

			return TCL_ERROR;

		}

	}

	return TCL_OK;

}



//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_addSimulationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	int numberOfSimulations;

	double targetCOV;



	if (theXuTransformation == 0 ) {

		cerr << "Need theXuTransformation before a SimulationAnalyis can be created" << endl;

		return TCL_ERROR;

	}

	if (theGFunEvaluator == 0 ) {

		cerr << "Need theGFunEvaluator before a SimulationAnalyis can be created" << endl;

		return TCL_ERROR;

	}

	if (theRandomNumberGenerator == 0 ) {

		cerr << "Need theRandomNumberGenerator before a SimulationAnalyis can be created" << endl;

		return TCL_ERROR;

	}









		// GET INPUT PARAMETER (integer)

      if (Tcl_GetInt(interp, argv[2], &numberOfSimulations) != TCL_OK) {

		cerr << "WARNING invalid input: numberOfSimulations \n";

		return TCL_ERROR;

	  }



	// GET INPUT PARAMETER (double)

	if (Tcl_GetDouble(interp, argv[3], &targetCOV) != TCL_OK) {

		cerr << "WARNING invalid input: targetCOV \n";

		return TCL_ERROR;

	}





	theSimulationAnalysis 

		= new SimulationAnalysis(theReliabilityDomain, theXuTransformation, 

								theGFunEvaluator, theRandomNumberGenerator, 

								argv[1], numberOfSimulations, 

								targetCOV);



	if (theSimulationAnalysis == 0) {

		cerr << "WARNING could not create theSimulationAnalysis \n";

		return TCL_ERROR;

	}

	return TCL_OK;

}







//////////////////////////////////////////////////////////////////

int 

TclReliabilityModelBuilder_analyzeReliability(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)

{

	bool analysisIsPerformed = false;

	bool systemAnalysisIsPerformed = false;



	if (theFORMAnalysis != 0) {

		if (theFORMAnalysis->analyze() < 0) {

			cerr << "FORMAnalysis was NOT completed successfully." << endl;

		}

		else {

			cerr << "FORMAnalysis completed successfully." << endl;

			analysisIsPerformed = true;

		}

	}

	if (theSORMAnalysis != 0) {

		if (theSORMAnalysis->analyze() < 0) {

			cerr << "SORMAnalysis was NOT completed successfully." << endl;

		}

		else {

			cerr << "SORMAnalysis completed successfully." << endl;

			analysisIsPerformed = true;

		}

	}

	if (theSimulationAnalysis != 0) {

		if (theSimulationAnalysis->analyze() < 0) {

			cerr << "SimulationAnalysis was NOT completed successfully." << endl;

		}

		else {

			cerr << "SimulationAnalysis completed successfully." << endl;

			analysisIsPerformed = true;

		}

	}

	if (theSystemAnalysis != 0) {

		if (theSystemAnalysis->analyze() < 0) {

			cerr << "SystemAnalysis was NOT completed successfully." << endl;

		}

		else {

			cerr << "SystemAnalysis completed successfully." << endl;

			analysisIsPerformed = true;

			systemAnalysisIsPerformed = true;

		}

	}





	// Open output file

	ofstream outputFile( "reliabilityResults.txt", ios::out );

	



	if (analysisIsPerformed) {



		// Number of limit-state functions

		int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();



		// PRINT SUMMARY OF RESULTS

		outputFile << "#######################################################################" << endl;

		outputFile << "#  SUMMARY OF RESULTS                                                 #" << endl;

		outputFile << "#                                   Reliability         Probability   #" << endl;

		outputFile << "#                                   Index Beta          of Failure    #" << endl;

		outputFile << "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #" << endl;



		// Loop over number of limit-state functions

		int lsf;

		for (lsf=1; lsf<=numLsf; lsf++ ) {



			// Get limit-state function pointer

			LimitStateFunction *theLimitStateFunction = 

				theReliabilityDomain->getLimitStateFunctionPtr(lsf);



			// Print results for this limit-state function

			theLimitStateFunction->printSummaryOfResults(outputFile);

		}

		outputFile << "#   Explanations:                                                     #" << endl;

		outputFile << "#   SORM(1): Curvatures found by search algorithm.                    #" << endl;

		outputFile << "#   SORM(2): Curvatures found by point-fitting.                       #" << endl;

		outputFile << "#   SORM(3): Curvatures found by curvature-fitting.                   #" << endl;

		outputFile << "#                                                                     #" << endl;

		outputFile << "#######################################################################" << endl << endl << endl;



		// PRINT ECHO OF INPUT

		outputFile << "#######################################################################" << endl;

		outputFile << "#  SUMMARY OF INPUT                                                   #" << endl;

		outputFile << "#  (void in present version...)                                       #" << endl;

		outputFile << "#                                                                     #" << endl;

		outputFile << "#######################################################################" << endl << endl << endl;



		// PRINT COMPLETE RESULTS

		// Loop over number of limit-state functions

		for (lsf=1; lsf<=numLsf; lsf++ ) {



			// Get limit-state function pointer

			LimitStateFunction *theLimitStateFunction = 

				theReliabilityDomain->getLimitStateFunctionPtr(lsf);



			// Print results for this limit-state function

			theLimitStateFunction->printResults(outputFile);

		}



		// PRINT SYSTEM RESULTS

		if (systemAnalysisIsPerformed) {



			double lowerBound = theSystemAnalysis->getLowerBound();

			double upperBound = theSystemAnalysis->getUpperBound();

			

			outputFile << "#######################################################################" << endl;

			outputFile << "#                                                                     #" << endl;

			outputFile << "#  SYSTEM ANALYSIS RESULTS                                            #" << endl;

			outputFile << "#                                                                     #" << endl;

			outputFile << "#  Lower probability bound: ........................... " 

				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lowerBound

				<< "  #" << endl;

			outputFile << "#  Upper probability bound: ........................... " 

				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<upperBound

				<< "  #" << endl;

			outputFile << "#                                                                     #" << endl;

			outputFile << "#######################################################################" << endl << endl << endl;



















		}



		return TCL_OK;

	}

	else {

		outputFile << "#######################################################################" << endl;

		outputFile << "#                                                                     #" << endl;

		outputFile << "#  NO ANALYSIS WAS PERFORMED                                          #" << endl;

		outputFile << "#                                                                     #" << endl;

		outputFile << "#######################################################################" << endl << endl << endl;

		return TCL_ERROR;

	}

}









