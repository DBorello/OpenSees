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
// $Date: 2001-06-14 08:06:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/OpenSeesGFunEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <fstream.h>

#include <OpenSeesGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>

#include <tcl.h>
#include <string.h>


OpenSeesGFunEvaluator::OpenSeesGFunEvaluator(int passedNumIter, Tcl_Interp *passedTclInterp,
					ReliabilityDomain *passedReliabilityDomain)
:GFunEvaluator()
{
	numIter					= passedNumIter;
	theTclInterp			= passedTclInterp;
	theReliabilityDomain	= passedReliabilityDomain;
}

OpenSeesGFunEvaluator::~OpenSeesGFunEvaluator()
{
}


int
OpenSeesGFunEvaluator::evaluate_g(Vector passed_x)
{

	// Zero out the response in the structural domain to make ready for next analysis
	char theRevertToStartCommand[10] = "reset";
	Tcl_Eval( theTclInterp, theRevertToStartCommand );

	// Put random variables into the structural domain according to the RandomVariablePositioners
	int numberOfRandomVariablePositioners = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	RandomVariablePositioner *theRandomVariablePositioner;
	int rvNumber;
	int i;
	for ( i=1 ; i<=numberOfRandomVariablePositioners ; i++ )  {
		theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		rvNumber				= theRandomVariablePositioner->getRvNumber();
		theRandomVariablePositioner->update(passed_x(rvNumber-1));
	}

	// Run the structural analysis
//	char theAnalyzeCommand[12] = "analyze 100";
	char theAnalyzeCommand[12];
	sprintf(theAnalyzeCommand,"analyze %d",numIter);
	Tcl_Eval( theTclInterp, theAnalyzeCommand );
	
	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsf);

	// Tokenize and parse the limit-state function(s loop)
	char *theExpression = theLimitStateFunction->getExpression();
	char *lsf_forTokenizing = new char[500];
	strcpy(lsf_forTokenizing,theExpression);
	char lsf_expression[500] = "";
	char *dollarSign = "$";
	char *underscore = "_";
	char *tokenPtr = strtok( lsf_forTokenizing, " _");
	while ( tokenPtr != NULL ) {
		// If a basic random variable is detected
		if ( strcmp(tokenPtr, "x") == 0) {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);
			tokenPtr = strtok( NULL, " _");
			strcat(lsf_expression, tokenPtr);
			int basicRandomVariableNumber = atoi( tokenPtr );
			char tclAssignment[500];
			sprintf(tclAssignment , "set x_%d  %15.5f", basicRandomVariableNumber, passed_x(basicRandomVariableNumber-1) );
			Tcl_Eval( theTclInterp, tclAssignment);
		}
		// If a nodal displacement is detected
		else if ( strcmp(tokenPtr, "nd") == 0) {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);
			tokenPtr = strtok( NULL, " _");
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);
			int nodeNumber = atoi( tokenPtr );
			tokenPtr = strtok( NULL, " _");
			strcat(lsf_expression, tokenPtr);
			int direction = atoi( tokenPtr );

			char tclAssignment[500];
			sprintf(tclAssignment,"set nd_%d_%d [nodeDisp %d %d ]",nodeNumber,direction,nodeNumber,direction);
			
			Tcl_Eval( theTclInterp, tclAssignment);
		}
		// If a bending moment is detected
		else if (tokenPtr == "bm") {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);

		}
		// If a shear force is detected
		else if (tokenPtr == "sf") {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);
		}
		// If an axial force is detected
		else if (tokenPtr == "af") {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);
		}
		// If a stress is detected
		else if (tokenPtr == "stress") {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);
		}
		// If a strain is detected
		else if (tokenPtr == "strain") {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);
		}
		else {
			strcat(lsf_expression, tokenPtr);
		}
		tokenPtr = strtok( NULL, " _");
	}
	// Evaluate the limit-state function
	double gvalue = 0.0;
	Tcl_ExprDouble( theTclInterp, lsf_expression, &gvalue );

	delete lsf_forTokenizing;

	g = gvalue;

	return 0;

}

double
OpenSeesGFunEvaluator::get_g()
{
	return g;
}
