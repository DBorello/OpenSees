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
// $Date: 2001-06-14 08:06:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SimulationAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <SimulationAnalysis.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <LimitStateFunction.h>
#include <XuTransformation.h>
#include <NatafXuTransformation.h>
#include <GFunEvaluator.h>
#include <BasicGFunEvaluator.h>
#include <RandomNumberGenerator.h>
#include <RandomVariable.h>
#include <NormalRV.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixOperations.h>
#include <NormalRV.h>
#include <math.h>
#include <stdlib.h>


SimulationAnalysis::SimulationAnalysis(	ReliabilityDomain *passedReliabilityDomain,
										XuTransformation *passedXuTransformation,
										GFunEvaluator *passedGFunEvaluator,
										RandomNumberGenerator *passedRandomNumberGenerator,
										char *passedPointToSampleAround,
										int passedNumberOfSimulations,
										double passedTargetCOV)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	theXuTransformation = passedXuTransformation;
	theGFunEvaluator = passedGFunEvaluator;
	theRandomNumberGenerator = passedRandomNumberGenerator;
	pointToSampleAround = new char[250];
	strcpy(pointToSampleAround,passedPointToSampleAround);
	numberOfSimulations = passedNumberOfSimulations;
	targetCOV = passedTargetCOV;

}


SimulationAnalysis::~SimulationAnalysis()
{
}



int 
SimulationAnalysis::analyze(void)
{

	// Alert the user that the simulation analysis has started
	cerr << "Simulation Analysis is running ... " << endl;


	// Declaration of some of the data used in the algorithm
	double gFunctionValue;
	int result;
	int I, i, j, k;
	double det_covariance;
	double phi;
	double h;
	double q_bar, pf, cov, beta, q;
	double variance_of_q_bar;
	double cov_of_q_bar;
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	Matrix covariance(numRV, numRV);
	Matrix chol_covariance(numRV, numRV);
	Matrix inv_covariance(numRV, numRV);
	Vector point(numRV);
	Vector startValues(numRV);
	Vector x(numRV);
	Vector z(numRV);
	Vector u(numRV);
	Vector randomArray(numRV);
	RandomVariable *aRandomVariable;
	LimitStateFunction *theLimitStateFunction = 0;
	NormalRV *aStdNormRV = 0;
	aStdNormRV = new NormalRV(1,0.0,1.0,0.0);

	
	// Check if computer ran out of memory
	if (aStdNormRV==0) {
		cerr << "SimulationAnalysis::analyze() - out of memory while instantiating internal objects." << endl;
		return -1;
	}

	
	// Establish covariance matrix
	for (i=0;  i<numRV;  i++) {
		for (j=0;  j<numRV;  j++) {
			if (i==j) {
				covariance(i,j) = 1.0;
			}
			else
				covariance(i,j) = 0.0;
		}
	}


	// Create object to do matrix operations on the covariance matrix
	MatrixOperations *theMatrixOperations = 0;
	theMatrixOperations = new MatrixOperations(covariance);
	if (theMatrixOperations == 0) {
		cerr << "SimulationAnalysis::analyze() - could not create" << endl
			<< " the object to perform matrix operations." << endl;
		return -1;
	}


	// Cholesky decomposition of covariance matrix
	result = theMatrixOperations->computeLowerCholesky();
	if (result < 0) {
		cerr << "SimulationAnalysis::analyze() - could not compute" << endl
			<< " the Cholesky decomposition of the covariance matrix." << endl;
		return -1;
	}
	chol_covariance = theMatrixOperations->getLowerCholesky();


	// Inverse of covariance matrix
	result = theMatrixOperations->computeInverse();
	if (result < 0) {
		cerr << "SimulationAnalysis::analyze() - could not compute" << endl
			<< " the inverse of the covariance matrix." << endl;
		return -1;
	}
	inv_covariance = theMatrixOperations->getInverse();


	// Compute the determinant, knowing that this is a diagonal matrix
	result = theMatrixOperations->computeTrace();
	if (result < 0) {
		cerr << "SimulationAnalysis::analyze() - could not compute" << endl
			<< " the trace of the covariance matrix." << endl;
		return -1;
	}
	det_covariance = theMatrixOperations->getTrace();
	

	// Initializations
	double sum_q;
	double sum_q_squared;


	// Pre-compute some factors to minimize computations inside simulation loop
	double pi = 3.14159265358979;
	double factor1 = 1.0 / ( pow((2.0*pi),((double)numRV/2.0)) );
	double factor2 = 1.0 / ( pow((2.0*pi),((double)numRV/2.0)) * sqrt(det_covariance) );


	// Number of limit-state functions
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();


	// Loop over number of limit-state functions
	for (int lsf=1; lsf<=numLsf; lsf++ ) {


		// Inform the user which limit-state function is being evaluated
		cerr << "Limit-state function number: " << lsf << endl;


		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);


		// Get the limit-state function pointer
		theLimitStateFunction = 0;
		lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
		theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
		if (theLimitStateFunction == 0) {
			cerr << "SimulationAnalysis::analyze() - could not find" << endl
				<< " limit-state function with tag #" << lsf << "." << endl;
			return -1;
		}


		// Determine the point to simulate around
		if ( strcmp(pointToSampleAround, "aroundStartPoint") == 0 ) {
			for ( j=1; j<=numRV; j++)
			{
				aRandomVariable = theReliabilityDomain->getRandomVariablePtr(j);
				startValues(j-1) = aRandomVariable->getStartValue();
			}
			result = theXuTransformation->set_x(startValues);
			if (result < 0) {
				cerr << "SimulationAnalysis::analyze() - could not " << endl
					<< " set the x-vector for xu-transformation. " << endl;
				return -1;
			}
			result = theXuTransformation->transform_x_to_u();
			if (result < 0) {
				cerr << "SimulationAnalysis::analyze() - could not " << endl
					<< " transform x to u. " << endl;
				return -1;
			}
			point = theXuTransformation->get_u();
		}
		else if ( strcmp(pointToSampleAround, "aroundOrigin") == 0 ) {
			for ( j=0; j<numRV; j++)
			{
				point(j) = 0.0;
			}
		}
		else if ( strcmp(pointToSampleAround, "aroundDesignPoint") == 0 ) {
			point = theLimitStateFunction->designPoint_u_inStdNormalSpace;
		}
		else {
			cerr << "SimulationAnalysis::analyze() - " << pointToSampleAround << endl
				<< " is not a valid point to sample around." << endl;
			return -1;
		}

		sum_q = 0.0;
		sum_q_squared = 0.0;
		k = 1;
		cov_of_q_bar = 999.0;

		while( k<=numberOfSimulations  && cov_of_q_bar>targetCOV ) {

			// Create array of standard normal random numbers
			result = theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numRV);
			if (result < 0) {
				cerr << "SimulationAnalysis::analyze() - could not generate" << endl
					<< " random numbers for simulation." << endl;
				return -1;
			}
			randomArray = theRandomNumberGenerator->getGeneratedNumbers();


			// Compute the point in standard normal space
			u = point + chol_covariance * randomArray;


			// Transform into original space
			result = theXuTransformation->set_u(u);
			if (result < 0) {
				cerr << "SimulationAnalysis::analyze() - could not " << endl
					<< " set the u-vector for xu-transformation. " << endl;
				return -1;
			}

			
			result = theXuTransformation->transform_u_to_x();
			if (result < 0) {
				cerr << "SimulationAnalysis::analyze() - could not " << endl
					<< " transform u to x. " << endl;
				return -1;
			}
			x = theXuTransformation->get_x();


			// Evaluate limit-state function
			result = theGFunEvaluator->evaluate_g(x);
			if (result < 0) {
				cerr << "SimulationAnalysis::analyze() - could not " << endl
					<< " evaluate limit-state function. " << endl;
				return -1;
			}
			gFunctionValue = theGFunEvaluator->get_g();


			// Collect result of sampling
			if (gFunctionValue < 0.0) {
				I = 1;
			}
			else {
				I = 0;
			}

			// Compute values of joint distributions at the u-point
			phi = factor1 * exp( -0.5 * (u ^ u) );
			Vector temp1 = inv_covariance ^ (u-point);
			double temp2 = temp1 ^ (u-point);
			h   = factor2 * exp( -0.5 * temp2 );


			// Update sums
			q = I * phi / h;
			sum_q = sum_q + q;
			sum_q_squared = sum_q_squared + q*q;

			if (sum_q > 0 && k>5 ) {
				// Compute coefficient of variation (of pf)
				q_bar = 1.0/(double)k * sum_q;

				variance_of_q_bar = 1.0/(double)k * 
					( 1.0/(double)k * sum_q_squared - (sum_q/(double)k)*(sum_q/(double)k));
				cov_of_q_bar = sqrt(variance_of_q_bar) / q_bar;
			}
			else {
				cov_of_q_bar = 999.0;
			}

			k++;
		}
		// Step 'k' back a step now that we went out
		k--;

		if (sum_q > 0) {
			// Compute probability of failure and reliability index
			pf = q_bar;
			cov = cov_of_q_bar;
			beta = -aStdNormRV->getInverseCDFvalue(pf);
		}
		else {
			pf = 0;
			cov = 0;
			beta = 0;
		}

		// Store results
		theLimitStateFunction->SimulationAnalysisPerformed=true;
		theLimitStateFunction->SimulationReliabilityIndexBeta=beta;
		theLimitStateFunction->SimulationProbabilityOfFailure_pfsim=pf;
		theLimitStateFunction->CoefficientOfVariationOfPfFromSimulation=cov;
		theLimitStateFunction->NumberOfSimulations=k;
	}

	delete aStdNormRV;

	return 0;
}

