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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/OpenSeesSensitivityEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <OpenSeesSensitivityEvaluator.h>
#include <Vector.h>
#include <Matrix.h>
#include <SensitivityEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <GFunEvaluator.h>
#include <RandomVariable.h>

#include <fstream.h>
#include <tcl.h>
#include <string.h>


OpenSeesSensitivityEvaluator::OpenSeesSensitivityEvaluator(
					Tcl_Interp *passedTclInterp,
					ReliabilityDomain *passedReliabilityDomain)
:SensitivityEvaluator()
{
	theTclInterp = passedTclInterp;
	theReliabilityDomain = passedReliabilityDomain;

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	grad_g = new Vector(nrv);
}

OpenSeesSensitivityEvaluator::~OpenSeesSensitivityEvaluator()
{
	delete grad_g;
}




Vector
OpenSeesSensitivityEvaluator::get_grad_g()
{
	return (*grad_g);
}


int
OpenSeesSensitivityEvaluator::evaluate_grad_g(double gFunValue, Vector passed_x)
{
	////////////////////////////////////
	// Notation:
	// g: the limit-state function
	// The gradient of the limit-state function:
	//   dgdx(m,1) = [ dgdu(1,n) * dudx(n,m) ]^T   +  dgdxufixed(m,1)
	// m: number of random variables x(i)
	// n: number of response quantities u(i)
	////////////////////////////////////

	// Some initial declaractions
	char tclAssignment[500];
	double g;
	double onedgdxufixed;
	double onedgdu;
	double onedudx;
	
	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsf);
	char *theExpression = theLimitStateFunction->getExpression();

	// Declare some much used tokens
	char *dollarSign = "$";
	char *underscore = "_";

	// Declare strings for the strings to be evaluated by Tcl
	char lsf_expression[500] = "";
	char grad_lsf_expression[500] = "";

	// Create copies of the limit-state function to be tokenized
	char *lsf_copy1 = new char[500];
	char *lsf_copy2 = new char[500];
	strcpy(lsf_copy1,theExpression);
	strcpy(lsf_copy2,theExpression);

	// Have a counter to see how many response quantities are in the lsf
	int count_u = 0;
	
	// Tokenize the limit-state function and assign values to the found quantities. 
	char *tokenPtr1 = strtok( lsf_copy1, " _"); // read first token
	while ( tokenPtr1 != NULL ) {

		// If a basic random variable is detected
		if ( strcmp(tokenPtr1, "x") == 0) {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr1);
			strcat(lsf_expression, underscore);
			tokenPtr1 = strtok( NULL, " _");	 // read a token
			strcat(lsf_expression, tokenPtr1);
			int basicRandomVariableNumber = atoi( tokenPtr1 );
			sprintf(tclAssignment , "set x_%d  %15.5f", basicRandomVariableNumber, passed_x(basicRandomVariableNumber-1) );
			Tcl_Eval( theTclInterp, tclAssignment);
		}

		// If a nodal displacement is detected
		else if ( strcmp(tokenPtr1, "nd") == 0) {
			count_u++;
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr1);
			strcat(lsf_expression, underscore);
			tokenPtr1 = strtok( NULL, " _");  // read a token
			strcat(lsf_expression, tokenPtr1);
			strcat(lsf_expression, underscore);
			int nodeNumber = atoi( tokenPtr1 );
			tokenPtr1 = strtok( NULL, " _");  // read a token
			strcat(lsf_expression, tokenPtr1);
			int direction = atoi( tokenPtr1 );
			sprintf(tclAssignment,"set nd_%d_%d [nodeDisp %d %d ]",nodeNumber,direction,nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
		}
		else {
			strcat(lsf_expression, tokenPtr1);
		}

		tokenPtr1 = strtok( NULL, " _");  // read a token
	} 

	// Assign value of of limit-state function to 'lsfvalue' in Tcl 
	Tcl_ExprDouble( theTclInterp, lsf_expression, &g );
	sprintf(tclAssignment,"set lsfvalue %35.20f", g);
	Tcl_Eval( theTclInterp, tclAssignment);

	// Declaration of quantities used in the gradient computations
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	Vector dgdx(nrv);
	Vector dgdu(count_u);
	Matrix dudx(count_u, nrv);
	Vector dgdxufixed(nrv);  // this must be initialized to zero: is it???
	RandomVariable *theRandomVariable;

	// Tokenize the limit-state function and COMPUTE GRADIENTS
	count_u = 0;
	char *tokenPtr2 = strtok( lsf_copy2, " _"); // read first token
	while ( tokenPtr2 != NULL ) {

		// If a basic random variable is detected
		if ( strcmp(tokenPtr2, "x") == 0) {
			tokenPtr2 = strtok( NULL, " _");	 // read a token
			int basicRandomVariableNumber = atoi( tokenPtr2 );
			// Get the standard deviation of the random variable
			theRandomVariable = theReliabilityDomain->getRandomVariablePtr(basicRandomVariableNumber);
			double stdv = theRandomVariable->getStdv();
			// Assign a perturbed value
			sprintf(tclAssignment , "set x_%d  %35.20f", basicRandomVariableNumber, 
				(passed_x(basicRandomVariableNumber-1)+0.005*stdv) );
			Tcl_Eval( theTclInterp, tclAssignment);
			// Evaluate limit-state function again
			Tcl_ExprDouble( theTclInterp, lsf_expression, &g );
			sprintf(tclAssignment,"set lsfperturbed %35.20f", g);
			Tcl_Eval( theTclInterp, tclAssignment);
			// Compute the gradient 'dgdxufixed' by finite difference
			sprintf(tclAssignment , "($lsfperturbed-$lsfvalu)/%35.20f", (0.005*stdv) );
			Tcl_ExprDouble( theTclInterp, tclAssignment, &onedgdxufixed );
			dgdxufixed(basicRandomVariableNumber-1)=onedgdxufixed;
			// Make assignment back to its original value
			sprintf(tclAssignment , "set x_%d  %35.20f", basicRandomVariableNumber, passed_x(basicRandomVariableNumber-1) );
			Tcl_Eval( theTclInterp, tclAssignment);
		}

		// If a nodal displacement is detected
		else if ( strcmp(tokenPtr2, "nd") == 0) {
			count_u++;
			tokenPtr2 = strtok( NULL, " _");  // read a token
			int nodeNumber = atoi( tokenPtr2 );
			tokenPtr2 = strtok( NULL, " _");  // read a token
			int direction = atoi( tokenPtr2 );
			// Assign a perturbed value
			sprintf(tclAssignment,"set nd_%d_%d [ expr [nodeDisp %d %d ]*1.001 ]",nodeNumber,direction,nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
			// Evaluate the limit-state function again
			Tcl_ExprDouble( theTclInterp, lsf_expression, &g );
			sprintf(tclAssignment,"set lsfperturbed %35.20f", g);
			Tcl_Eval( theTclInterp, tclAssignment);
			// Compute the gradient 'dgdu' by finite difference
			sprintf(tclAssignment , "($lsfperturbed-$lsfvalue)/([nodeDisp %d %d ]*0.001)",nodeNumber,direction);
			Tcl_ExprDouble( theTclInterp, tclAssignment, &onedgdu );
			dgdu(count_u-1)=onedgdu;
			// Make assignment back to its original value
			sprintf(tclAssignment,"set nd_%d_%d [nodeDisp %d %d ]",nodeNumber,direction,nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
			// Assign one row in the 'dudx' matrix
			for (int i=1; i<=nrv; i++) {
				sprintf(tclAssignment , "set sens [sensNodeDisp %d %d %d ]",nodeNumber,direction,i);
				Tcl_Eval( theTclInterp, tclAssignment);
				sprintf(tclAssignment , "$sens ");
				Tcl_ExprDouble( theTclInterp, tclAssignment, &onedudx );
				dudx( (count_u-1), (i-1) ) = onedudx;
			}

		}

		tokenPtr2 = strtok( NULL, " _");  // read next token and go up and check the while condition again

	} 

	// Evaluate the product that yields 'dgdx'
	Vector temp;
	temp = dudx ^ dgdu;
	dgdx = temp + dgdxufixed;

	delete lsf_copy1;
	delete lsf_copy2;

	(*grad_g) = dgdx;

	return 0;

}
