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
// $Date: 2003-02-14 23:01:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/designPoint/SearchWithStepSizeAndStepDirection.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//			haukaas 08/19/01 (modifications for Release 1.2 of OpenSees)
//

#include <SearchWithStepSizeAndStepDirection.h>
#include <FindDesignPoint.h>
#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <XuTransformation.h>
#include <NatafXuTransformation.h>
#include <GFunEvaluator.h>
#include <SensitivityEvaluator.h>
#include <RandomVariable.h>
#include <CorrelationCoefficient.h>
#include <MatrixOperations.h>
#include <Matrix.h>
#include <Vector.h>
#include <GammaRV.h>



SearchWithStepSizeAndStepDirection::SearchWithStepSizeAndStepDirection(
					int passedMaxNumberOfIterations, 
					double passedConvergenceCriterionE1, 
					double passedConvergenceCriterionE2,
					GFunEvaluator *passedGFunEvaluator,
					SensitivityEvaluator *passedSensitivityEvaluator,
					StepSizeRule *passedStepSizeRule,
					SearchDirection *passedSearchDirection,
					XuTransformation *passedXuTransformation)
:FindDesignPoint()
{
	maxNumberOfIterations			= passedMaxNumberOfIterations;
	convergenceCriterionE1			= passedConvergenceCriterionE1;
	convergenceCriterionE2			= passedConvergenceCriterionE2;
	theGFunEvaluator				= passedGFunEvaluator;
	theSensitivityEvaluator			= passedSensitivityEvaluator;
	theStepSizeRule					= passedStepSizeRule;
	theSearchDirection				= passedSearchDirection;
	theXuTransformation				= passedXuTransformation;
}



SearchWithStepSizeAndStepDirection::~SearchWithStepSizeAndStepDirection()
{
}


int
SearchWithStepSizeAndStepDirection::findDesignPoint(
											Vector *passedStartPoint, 
											ReliabilityDomain *passedReliabilityDomain)
{
	// Set the start point (a data member of this class)
	startPoint = passedStartPoint;

	// Set the reliability domain (a data member of this class)
	theReliabilityDomain = passedReliabilityDomain;

	// RUN the generic algorithm
	int result = doTheActualSearch(false);

	return result;
}



int
SearchWithStepSizeAndStepDirection::approachDesignPointThroughOrthogonalSubspace(
											Vector *passedDesignPoint_u,
											int passedNumberOfAxes,
											Vector *passedPrincipalAxes, 
											Vector *passedStartPoint, 
											ReliabilityDomain *passedReliabilityDomain)
{
	// Set number of axes in the long vector
	numberOfAxes = passedNumberOfAxes;

	// Set the start point (a data member of this class)
	startPoint = passedStartPoint;

	// Set the start point (a data member of this class)
	designPoint_uStar = passedDesignPoint_u;

	// Set vector of principal axes (a data member of this class)
	principalAxesPtr = passedPrincipalAxes;

	// Set the reliability domain (a data member of this class)
	theReliabilityDomain = passedReliabilityDomain;

	// RUN the generic algorithm
	int result = doTheActualSearch(true);

	return result;
}



int
SearchWithStepSizeAndStepDirection::doTheActualSearch(	bool doProjectionToOrthogonalSubspace)
{
	// Get number of random variables and correlation coefficients from reliability domain
	int numberOfRandomVariables = theReliabilityDomain->getNumberOfRandomVariables();
	
	// Declaration of data used in the algorithm
	int j;
	int zeroFlag;
	x(numberOfRandomVariables);
	u(numberOfRandomVariables);
	uSecondLast(numberOfRandomVariables);
	Vector uNew(numberOfRandomVariables);
	alpha (numberOfRandomVariables);
	gamma (numberOfRandomVariables);
	alphaSecondLast(numberOfRandomVariables);
	double gFunctionValue;
	Vector gradientOfgFunction(numberOfRandomVariables);
	Vector gradientInStandardNormalSpace(numberOfRandomVariables);
	double normOfGradient =0;
	double stepSize;
	Matrix jacobian_x_u(numberOfRandomVariables,numberOfRandomVariables);
	Vector u_minus_alpha_u_alpha;
	double alpha_times_u;
	double NormOf_u_minus_alpha_u_alpha;
	double possibleNewGFunValue;
	double criterion1, criterion2;


	// Get starting point
	x = *startPoint;


	// Transform starting point into standard normal space
	int result = theXuTransformation->set_x(x);
	if (result < 0) {
		opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not set x in the xu-transformation." << endln;
		return -1;
	}


	result = theXuTransformation->transform_x_to_u();
	if (result < 0) {
		opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not transform from x to u." << endln;
		return -1;
	}
	u = theXuTransformation->get_u();


	// Loop to find design point
	i = 1;
	while ( i <= maxNumberOfIterations )
	{

		// Transform from u to x space
		result = theXuTransformation->set_u(u);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not set u in the xu-transformation." << endln;
			return -1;
		}

		result = theXuTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not transform from u to x and compute Jacobian." << endln;
			return -1;
		}
		x = theXuTransformation->get_x();
		jacobian_x_u = theXuTransformation->getJacobian_x_u();


		// Evaluate limit-state function unless it has been done in 
		// a trial step by the "stepSizeAlgorithm"
		if (possibleNewGFunValue == gFunctionValue || i==1) {
			result = theGFunEvaluator->evaluate_g(x);
			if (result < 0) {
				opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
					<< " could not evaluate limit-state function. " << endln;
				return -1;
			}
			gFunctionValue = theGFunEvaluator->get_g();
		}
		else {
			gFunctionValue = possibleNewGFunValue;
		}


		// Provide user with information about the limit-state function value
		// (only in step one)
		opserr << " ITERATION #" << i;
		if (i==1) {
			opserr << ", limit-state function value: " << gFunctionValue << endln;
		}


		// Gradient in original space
		result = theSensitivityEvaluator->evaluate_grad_g(gFunctionValue,x);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute gradients of the limit-state function. " << endln;
			return -1;
		}
		gradientOfgFunction = theSensitivityEvaluator->get_grad_g();


		// Check if all components of the vector is zero
		zeroFlag = 0;
		for (j=0; j<gradientOfgFunction.Size(); j++) {
			if (gradientOfgFunction[j] != 0.0) {
				zeroFlag = 1;
			}
		}
		if (zeroFlag == 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " all components of the gradient vector is zero. " << endln;
			return -1;
		}


		// Gradient in standard normal space
		gradientInStandardNormalSpace = jacobian_x_u ^ gradientOfgFunction;


		// Set scale parameter
		if (i == 1)
		{
			Go = gFunctionValue;
		}


		// Compute the norm of the gradient in standard normal space
		normOfGradient = gradientInStandardNormalSpace.Norm();


		// Check that the norm is not zero
		if (normOfGradient == 0.0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " the norm of the gradient is zero. " << endln;
			return -1;
		}


		// Compute alpha-vector
		alpha = gradientInStandardNormalSpace *  ( (-1) / normOfGradient );


		// Check convergence
		alpha_times_u = alpha ^ u;
		u_minus_alpha_u_alpha = u - alpha * alpha_times_u;
		NormOf_u_minus_alpha_u_alpha = u_minus_alpha_u_alpha.Norm();
		criterion1 = fabs(gFunctionValue / Go);
		criterion2 = NormOf_u_minus_alpha_u_alpha;


		// Inform user about convergence status 
		if (i != 1) {
			opserr << ", convergence checks: (" << criterion1 << ") (" << criterion2 << ")" << endln;
		}


		// Postprocessing if the analysis converged
		if ( ( criterion1 < convergenceCriterionE1 ) & 
			( criterion2 < convergenceCriterionE2 ) ) 
		{
			// Compute the gamma vector
			MatrixOperations theMatrixOperations(jacobian_x_u);


			result = theMatrixOperations.computeTranspose();
			if (result < 0) {
				opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
					<< " could not compute transpose of jacobian matrix. " << endln;
				return -1;
			}
			Matrix transposeOfJacobian_x_u = theMatrixOperations.getTranspose();


			Matrix jacobianProduct = jacobian_x_u * transposeOfJacobian_x_u;


			theMatrixOperations.setMatrix(jacobianProduct);
			result = theMatrixOperations.computeSquareRoot();
			if (result < 0) {
				opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
					<< " could not compute square root of matrix elements. " << endln;
				return -1;
			}
			Matrix squareRootOfJacobianProduct = theMatrixOperations.getSquareRoot();
			

			Matrix D_prime(numberOfRandomVariables,numberOfRandomVariables);
			for (j=0; j<numberOfRandomVariables; j++) {
				D_prime(j,j) = squareRootOfJacobianProduct(j,j);
			}


			Matrix jacobian_u_x = theXuTransformation->getJacobian_u_x();

			Vector tempProduct = jacobian_u_x ^ alpha;

			gamma = D_prime ^ tempProduct;

			return 1;
		}


		// Store 'u' and 'alpha' at the second last iteration point
		uSecondLast = u;
		alphaSecondLast = alpha;


		// Determine search direction
		result = theSearchDirection->computeSearchDirection(
			u, gFunctionValue, gradientInStandardNormalSpace );
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute search direction. " << endln;
			return -1;
		}
		searchDirection = theSearchDirection->getSearchDirection();


		// Determine step size
		result = theStepSizeRule->computeStepSize(
			u, gradientInStandardNormalSpace, gFunctionValue, searchDirection);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute step size. " << endln;
			return -1;
		}
		stepSize = theStepSizeRule->getStepSize();


		// For efficiency: get g-function value from the step size algorithm.
		// If it is different from the original then a trial step was actually performed. 
		possibleNewGFunValue = theStepSizeRule->getGFunValue();


		// Determine new iteration point (take the step)
		u += searchDirection * stepSize;

		// Possibly do projection to subspace
		if (doProjectionToOrthogonalSubspace) {
			doProjection(u, uNew);
			u = uNew;
		}

		// Increment the loop parameter
		i++;

	}


	// Print a message if max number of iterations was reached
	opserr << "Maximum number of iterations was reached before convergence." << endln;

	return 0;
}




int
SearchWithStepSizeAndStepDirection::doProjection(Vector uOld, Vector uNew)
{

	// THIS METHOD IS CURRENTLY UNDER DEVELOPMENT

	// Initial declarations
	double dNorm;
	double temp;

	// Number of random variables
	int nrv = uOld.Size();

	// Declaration of the 'current' orthogonality vector 'd' (and the whole vector of d's)
	Vector d(nrv);
	Vector principalAxes = *principalAxesPtr;

	// Other auxiliary declarations
	Vector uDifference(nrv);
	Vector correction(nrv);

	for (int i=0; i<numberOfAxes; i++) {

		// Extract the current orthogonality vector
		for (int j=0; j<nrv; j++) {
			d(j) = principalAxes((i*nrv)+j);
		}

		// Norm of the orthogonality vector
		temp = d.Norm();
		dNorm = temp*temp;

		// Compute the formula
		uDifference = *designPoint_uStar - uOld;

		double dotProduct = d ^ uDifference;

		correction = (1.0/dNorm * dotProduct) * d;
		
		uNew = uOld + correction;
		uOld = uNew;
	}

	return 0;
}


Vector
SearchWithStepSizeAndStepDirection::get_x()
{
	return x;
}

Vector
SearchWithStepSizeAndStepDirection::get_u()
{
	return u;
}

Vector
SearchWithStepSizeAndStepDirection::get_alpha()
{
	return alpha;
}

Vector
SearchWithStepSizeAndStepDirection::get_gamma()
{
	return gamma;
}

int
SearchWithStepSizeAndStepDirection::getNumberOfIterations()
{
	return i;
}

Vector
SearchWithStepSizeAndStepDirection::getSecondLast_u()
{
	return uSecondLast;
}

Vector
SearchWithStepSizeAndStepDirection::getSecondLast_alpha()
{
	return alphaSecondLast;
}

Vector
SearchWithStepSizeAndStepDirection::getLastSearchDirection()
{
	return searchDirection;
}

double
SearchWithStepSizeAndStepDirection::getFirstGFunValue()
{
	return Go;
}






