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
                                                                        
// $Revision: 1.1 $
// $Date: 2001-06-13 05:06:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/SensitivityByFiniteDifference.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <SensitivityByFiniteDifference.h>
#include <Vector.h>
#include <SensitivityEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <GFunEvaluator.h>
#include <RandomVariable.h>
#include <fstream.h>


SensitivityByFiniteDifference::SensitivityByFiniteDifference(
					GFunEvaluator *passedGFunEvaluator,
					ReliabilityDomain *passedReliabilityDomain)
:SensitivityEvaluator()
{
	theGFunEvaluator = passedGFunEvaluator;
	theReliabilityDomain = passedReliabilityDomain;

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	grad_g = new Vector(nrv);
}

SensitivityByFiniteDifference::~SensitivityByFiniteDifference()
{
	delete grad_g;
}



Vector
SensitivityByFiniteDifference::get_grad_g()
{
	return (*grad_g);
}



int
SensitivityByFiniteDifference::evaluate_grad_g(double gFunValue, Vector passed_x)
{
	
	// Initial declarations
	int numberOfRandomVariables = passed_x.Size();
	Vector perturbed_x(numberOfRandomVariables);
	RandomVariable *theRandomVariable;
	int i;
	double h;
	double gFunValueAStepAhead;
	double stdv;


	// Set perturbation factor (perturbation = stdv/perturbationFactor)
	double perturbationFactor = 200.0;


	// For each random variable: perturb and run analysis again
	for ( i=0 ; i<numberOfRandomVariables ; i++ )
	{
		// Get random variable from domain
		theRandomVariable = theReliabilityDomain->getRandomVariablePtr(i+1);


		// Get the standard deviation
		stdv = theRandomVariable->getStdv();


		// Compute perturbation
		h = stdv/perturbationFactor;


		// Compute perturbed vector of random variables realization
		perturbed_x = passed_x;
		perturbed_x(i) = perturbed_x(i) + h;


		// Evaluate limit-state function
		double result = theGFunEvaluator->evaluate_g(perturbed_x);
		if (result < 0) {
			cerr << "SensitivityByFiniteDifference::evaluate_grad_g() - " << endl
				<< " could not evaluate limit-state function. " << endl;
			return -1;
		}
		gFunValueAStepAhead = theGFunEvaluator->get_g();


		// Compute the derivative by finite difference
		(*grad_g)(i) = (gFunValueAStepAhead - gFunValue) / h;

	}

	return 0;

}
