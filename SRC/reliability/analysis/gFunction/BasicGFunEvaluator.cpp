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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-08-02 18:13:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/BasicGFunEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <fstream.h>

#include <BasicGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>

#include <tcl.h>
#include <string.h>


BasicGFunEvaluator::BasicGFunEvaluator(Tcl_Interp *passedTclInterp, 
										ReliabilityDomain *passedReliabilityDomain)

:GFunEvaluator()
{
	theTclInterp			= passedTclInterp;
	theReliabilityDomain = passedReliabilityDomain;
}

BasicGFunEvaluator::~BasicGFunEvaluator()
{
}


int
BasicGFunEvaluator::evaluate_g(Vector passed_x)
{
	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsf);

	// Tokenize and parse the limit-state function(s loop)
	char *theExpression = theLimitStateFunction->getExpression();
	char *lsf_forTokenizing = new char[100];
	strcpy(lsf_forTokenizing,theExpression);
	char lsf_expression[100] = "";
	char *dollarSign = "$";
	char *underscore = "_";
	char *tokenPtr = strtok( lsf_forTokenizing, " _");
	while ( tokenPtr != NULL ) {
		// If a nodal basic random variable is detected
		if ( strcmp(tokenPtr, "x") == 0) {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr);
			strcat(lsf_expression, underscore);
			tokenPtr = strtok( NULL, " _");
			strcat(lsf_expression, tokenPtr);
			int basicRandomVariableNumber = atoi( tokenPtr );
			char tclAssignment[100];
			sprintf(tclAssignment , "set x_%d  %15.5f", basicRandomVariableNumber, passed_x(basicRandomVariableNumber-1) );
			Tcl_Eval( theTclInterp, tclAssignment);
		}
		else {
			strcat(lsf_expression, tokenPtr);
		}
		tokenPtr = strtok( NULL, " _");
	}
	// Evaluate the limit-state function
	double gvalue = 0.0;
	Tcl_ExprDouble( theTclInterp, lsf_expression, &gvalue );

	delete [] lsf_forTokenizing;

	g = gvalue;

	return 0;
}


double
BasicGFunEvaluator::get_g()
{
	return g;
}
