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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/FORMAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <FORMAnalysis.h>
#include <FindDesignPoint.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixOperations.h>
#include <NormalRV.h>
#include <RandomVariable.h>
#include <math.h>

FORMAnalysis::FORMAnalysis(	ReliabilityDomain *passedReliabilityDomain,
							FindDesignPoint *passedFindDesignPoint)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	theFindDesignPoint = passedFindDesignPoint;
}


FORMAnalysis::~FORMAnalysis()
{
}



int 
FORMAnalysis::analyze(void)
{

	// Alert the user that the FORM analysis has started
	cerr << "FORM Analysis is running ... " << endl;


	// Declare variables used in this method
	Vector xStar;
	Vector uStar;
	Vector alpha;
	Vector gamma;
	int i;
	double Go;
	Vector uSecondLast;
	Vector alphaSecondLast;
	Vector lastSearchDirection;
	double beta;
	double pf1;
	int lsf;
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();
	Vector startValues(numRV);
	RandomVariable *aRandomVariable;
	LimitStateFunction *theLimitStateFunction;
	NormalRV *aStdNormRV=0;
	aStdNormRV = new NormalRV(1,0.0,1.0,0.0);


	// Check if computer ran out of memory
	if (aStdNormRV==0) {
		cerr << "FORMAnalysis::analyze() - out of memory while instantiating internal objects." << endl;
		return -1;
	}


	// Establish the user-given starting point in the original space
	for ( int j=1; j<=numRV; j++ )
	{
		aRandomVariable = 0;
		aRandomVariable = theReliabilityDomain->getRandomVariablePtr(j);
		if (aRandomVariable == 0) {
			cerr << "FORMAnalysis::analyze() - could not find" << endl
				<< " random variable with tag #" << j << "." << endl;
			return -1;
		}
		startValues(j-1) = aRandomVariable->getStartValue();
	}


	// Loop over number of limit-state functions and perform FORM analysis
	for (lsf=1; lsf<=numLsf; lsf++ ) {


		// Inform the user which limit-state function is being evaluated
		cerr << "Limit-state function number: " << lsf << endl;


		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);


		// Get the limit-state function pointer
		theLimitStateFunction = 0;
		lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
		theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
		if (theLimitStateFunction == 0) {
			cerr << "FORMAnalysis::analyze() - could not find" << endl
				<< " limit-state function with tag #" << lsf << "." << endl;
			return -1;
		}


		// Find the design point
		if (theFindDesignPoint->findDesignPoint(&startValues, theReliabilityDomain) < 0){
			cerr << "FORMAnalysis::analyze() - failed while finding the" << endl
				<< " design point for limit-state function number " << lsf << "." << endl;
			return -1;
		}

		
		// Get results from the "find desingn point algorithm"
		xStar = theFindDesignPoint->get_x();
		uStar = theFindDesignPoint->get_u();
		alpha = theFindDesignPoint->get_alpha();
		gamma = theFindDesignPoint->get_gamma();
		i        = theFindDesignPoint->getNumberOfIterations();
		Go    = theFindDesignPoint->getFirstGFunValue();
		uSecondLast		= theFindDesignPoint->getSecondLast_u();
		alphaSecondLast	= theFindDesignPoint->getSecondLast_alpha();
		lastSearchDirection	= theFindDesignPoint->getLastSearchDirection();


		// Postprocessing
		beta = alpha ^ uStar;
		pf1 = 1.0 - aStdNormRV->getCDFvalue(beta);


		// Store the results
		theLimitStateFunction->FORMAnalysisPerformed				= true;
		theLimitStateFunction->FORMReliabilityIndexBeta				= beta;
		theLimitStateFunction->FORMProbabilityOfFailure_pf1			= pf1;
		theLimitStateFunction->designPoint_x_inOriginalSpace		= xStar;
		theLimitStateFunction->designPoint_u_inStdNormalSpace		= uStar;
		theLimitStateFunction->normalizedNegativeGradientVectorAlpha= alpha;
		theLimitStateFunction->importanceVectorGamma				= gamma;
		theLimitStateFunction->numberOfIterationsToFindDesignPoint	= i;
		theLimitStateFunction->GFunValueAtStartPt					= Go;
		theLimitStateFunction->secondLast_u							= uSecondLast;
		theLimitStateFunction->secondLastAlpha						= alphaSecondLast;
		theLimitStateFunction->lastSearchDirection					= lastSearchDirection;

	}

	delete aStdNormRV;

	return 0;
}

