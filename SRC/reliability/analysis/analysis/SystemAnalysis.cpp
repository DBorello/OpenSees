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
// $Date: 2001-06-14 08:06:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SystemAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <LimitStateFunction.h>
#include <MatrixOperations.h>
#include <NormalRV.h>
#include <math.h>
#include <fstream.h>

SystemAnalysis::SystemAnalysis(	ReliabilityDomain *passedReliabilityDomain)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
}


SystemAnalysis::~SystemAnalysis()
{
}


double 
SystemAnalysis::getLowerBound(void)
{
	return lowerBound;
}



double 
SystemAnalysis::getUpperBound(void)
{
	return upperBound;
}



int 
SystemAnalysis::analyze(void)
{
cerr << "System Reliability Analysis is running ... " << endl;

	// Initial declarations
	double beta;
	double pf1;
	Vector alpha;
	LimitStateFunction *theLimitStateFunction;
	NormalRV aStdNormalRV(1, 0.0, 1.0, 0.0);

	// Number of limit-state functions
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();

	// Number of random variables
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	// Allocate vectors to store ALL the betas and alphas
	Vector allBetas(numLsf);
	Vector allPf1s(numLsf);
	Vector allAlphas(numLsf*nrv);

	// Loop over number of limit-state functions and collect results
	int i;
	for (i=0; i<numLsf; i++ ) {

		// Get FORM results from the limit-state function
		theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(i+1);
		beta = theLimitStateFunction->FORMReliabilityIndexBeta;
		pf1 = theLimitStateFunction->FORMProbabilityOfFailure_pf1;
		alpha = theLimitStateFunction->normalizedNegativeGradientVectorAlpha;

		// Put FORM results into vector of all betas and alphas
		allBetas(i) = beta;
		allPf1s(i) = pf1;
		for (int j=0; j<nrv; j++ ) {
			allAlphas((i*nrv)+j) = alpha(j);
		}
	}

	// Compute vector of 'rhos', that is, dot products between the alphas
	// in order of the rows (without diagonal): r12, r13, r23 for the case of 3 r.v.s
	int rhosDim = (int)(0.5*(nrv*nrv-nrv));
	Vector rhos(rhosDim);
	int counter = 0;
	double dotProduct;
	for (i=0; i<numLsf; i++ ) {
		for (int j=i+1; j<numLsf; j++ ) {
			dotProduct = 1.0;
			for ( int k=0; k<nrv; k++ ) {
				dotProduct=dotProduct*allAlphas((i*nrv)+k)*allAlphas((j*nrv)+k);
			}
			rhos(counter) = dotProduct;
			counter++;
		}
	}

	// Compute the bi-variate normal distribution for all the pairs
	Vector Pmn(rhosDim);
	counter = 0;
	double beta1;
	double beta2;
	int m;
	for (m=0; m<numLsf; m++ ) {
		for (int n=m+1; n<numLsf; n++ ) {
			beta1 = allBetas(m);
			beta2 = allBetas(n);
			// Numerical integration with Simpsons rule
			double a = 0.0;				// Interval start
			double b = rhos(counter);	// Interval end
			double integral = 0.0;		// The result of integration
			int n_2 = 100;				// Half the number of intervals
			double h = b-a;				// Interval length
			double fa = functionToIntegrate(a,beta1,beta2);
			double fb = functionToIntegrate(b,beta1,beta2);
			double sum_fx2j = 0.0;
			double sum_fx2j_1 = 0.0;
			for (int j=1;  j<=n_2;  j++) {
				sum_fx2j = sum_fx2j + functionToIntegrate(((double)(j*2)*h),beta1,beta2);
				sum_fx2j_1 = sum_fx2j_1 + functionToIntegrate(((double)(j*2-1)*h),beta1,beta2);
			}
			sum_fx2j = sum_fx2j - functionToIntegrate(b,beta1,beta2);
			integral = h/3.0*(fa + 2.0*sum_fx2j + 4.0*sum_fx2j_1 + fb);
			Pmn(counter) = aStdNormalRV.getCDFvalue(-beta1)
				*aStdNormalRV.getCDFvalue(-beta2) + integral;
			counter++;
		}
	}

	// Compute lower probability bound
	counter = 0;
	double temp1;
	double rowSum;
	lowerBound = 0.0;
	for ( m=0; m<numLsf; m++ ) {
		rowSum = 0.0;
		for ( i=0; i<numLsf; i++ ) {
			for ( int j=i+1; j<numLsf; j++ ) {
				if (j==m) {
					rowSum += Pmn(counter);
				}
				counter++;
			}
		}
		if ((allPf1s(i)-rowSum)>0.0) {
			temp1 = (allPf1s(i)-rowSum);
		}
		else {
			temp1 = 0.0;
		}
		lowerBound += temp1;
	}

	// Compute upper probability bound;
	double largeNegativeNumber = -9999999.0;
	upperBound = 0.0;
	double theMaxValue = largeNegativeNumber;
	for ( m=0; m<numLsf; m++ ) {
		rowSum = 0.0;
		counter = 0;
		for ( i=0; i<numLsf; i++ ) {
			for ( int j=i+1; j<numLsf; j++ ) {
				if (j==m) {
					if (Pmn(counter)>theMaxValue) {
						theMaxValue = Pmn(counter);
					}
				}
				counter++;
			}
		}
		if (theMaxValue == largeNegativeNumber) {
			theMaxValue = 0.0;
		}
		upperBound += allPf1s(m) - theMaxValue;

	}

	return 0;
}



double
SystemAnalysis::functionToIntegrate(double rho, double beta1, double beta2)
{
	double pi = 3.14159265358979;
	return 1.0/(2.0*pi*sqrt(1.0-rho*rho)) 
		* exp(-(beta1*beta1+beta2*beta2-2*rho*beta1*beta2)
		/(2.0*(1.0-rho*rho)));
}
