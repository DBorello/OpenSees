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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-14 23:01:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/stepSize/ArmijoRule.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//          haukaas 06/15/01 (notify the user if the step size is reduced)
//			haukaas 08/19/01 (modifications for Release 1.2 of OpenSees)
//

#include <ArmijoRule.h>
#include <StepSizeRule.h>
#include <XuTransformation.h>
#include <math.h>
#include <Vector.h>


ArmijoRule::ArmijoRule(	GFunEvaluator *passedGFunEvaluator,
						XuTransformation *passedXuTransformation,
						int passedMaxNumReductions,
						double passedFirstTrialStep)
:StepSizeRule()
{
	theGFunEvaluator = passedGFunEvaluator;
	theXuTransformation = passedXuTransformation;
	firstTrialStep = passedFirstTrialStep;
	maxNumReductions = passedMaxNumReductions;
	gFunValue = 0;
}

ArmijoRule::~ArmijoRule()
{
}



double 
ArmijoRule::getStepSize()
{
	return stepSize;
}



double 
ArmijoRule::getGFunValue()
{
	return gFunValue;
}




int
ArmijoRule::computeStepSize(Vector u, 
									Vector grad_G, 
									double G, 
									Vector d)
{
	// Check that the norm of the gradient is not zero
	double gradNorm = grad_G.Norm();
	if (gradNorm == 0.0) {
		opserr << "ArmijoRule::computeStepSize() - the norm " << endln
			<< " of the gradient is zero. " << endln;
		return -1;
	}

	// Compute factor 'c'
	double c = ( u.Norm() / gradNorm ) * 2.0 + 10.0;


	// Compute merit function
	double merit = 0.5 * u.Norm() * u.Norm() + c * fabs(G);

	
	// Set the first trial step size
	double trial_step_size = firstTrialStep;
	

	// Take a trial step in standard normal space
	Vector trial_u = u + d * trial_step_size;


	// Transform the trial point into original space
	double result = theXuTransformation->set_u(trial_u);
	if (result < 0) {
		opserr << "ArmijoRule::computeStepSize() - could not set " << endln
			<< " vector u in the xu-transformation. " << endln;
		return -1;
	}
	result = theXuTransformation->transform_u_to_x();
	if (result < 0) {
		opserr << "ArmijoRule::computeStepSize() - could not  " << endln
			<< " transform u to x. " << endln;
		return -1;
	}
	Vector trial_x = theXuTransformation->get_x();
	

	// Evaluate the limit-state function
	result = theGFunEvaluator->evaluate_g(trial_x);
	if (result < 0) {
		opserr << "ArmijoRule::computeStepSize() - could not  " << endln
			<< " evaluate the limit-state function. " << endln;
		return -1;
	}
	double trial_G = theGFunEvaluator->get_g();


	// Compute merit function again, based on the new limit-state function value
	double merit_new = 0.5 * trial_u.Norm() * trial_u.Norm() + c * fabs(trial_G);


	// If the merit function did not decrease, we try to 
	// repeatedly reduce the step size until it does. 
	int i = 1;
	while ( (merit_new > merit) && (i<(maxNumReductions+1)) ) {

		
		// Notify user that step sizes are being reduced
		opserr << "... the step size is being reduced by the Armijo rule ..." << endln;

		// Update merit function value
		merit = merit_new;


		// Cut the step size in half and try that instead
		trial_step_size = trial_step_size * 0.5;
		trial_u = u + d * trial_step_size;


		// Transform the trial point into original space
		double result = theXuTransformation->set_u(trial_u);
		if (result < 0) {
			opserr << "ArmijoRule::computeStepSize() - could not set " << endln
				<< " vector u in the xu-transformation. " << endln;
			return -1;
		}
		result = theXuTransformation->transform_u_to_x();
		if (result < 0) {
			opserr << "ArmijoRule::computeStepSize() - could not  " << endln
				<< " transform u to x. " << endln;
			return -1;
		}
		Vector trial_x = theXuTransformation->get_x();


		// Evaluate the limit-state function
		result = theGFunEvaluator->evaluate_g(trial_x);
		if (result < 0) {
			opserr << "ArmijoRule::computeStepSize() - could not  " << endln
				<< " evaluate the limit-state function. " << endln;
			return -1;
		}
		double trial_G = theGFunEvaluator->get_g();


		// Compute merit function again, based on the new limit-state function value
		merit_new = 0.5 * trial_u.Norm() * trial_u.Norm() + c * fabs(trial_G);


		// Report to user if the step reduction was not successful
		if ( (i==5) && (merit_new > merit) ) {
			opserr << "... tried to half step size " << maxNumReductions 
				<< " times before continuing ..." << endln;
		}

		
		// Increment counter
		i++;
 	
	}

	stepSize = trial_step_size;
	gFunValue = trial_G;

	return 0;

}
