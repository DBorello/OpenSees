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
// $Date: 2001-06-14 08:06:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunction.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#ifndef LimitStateFunction_h
#define LimitStateFunction_h

#include <ReliabilityDomainComponent.h>
#include <Vector.h>
#include <fstream.h>

class LimitStateFunction : public ReliabilityDomainComponent
{

public:
	LimitStateFunction(	int tag, 
						char *expression);
	~LimitStateFunction();
	void Print(ostream &s, int flag =0);
	char *getExpression();
	void printResults(ofstream &outputFile);
	void printSummaryOfResults(ofstream &outputFile);

	
	// Flags to check if analyses are performed
	bool FORMAnalysisPerformed;
	bool SimulationAnalysisPerformed;
	bool EvaluateLimitStateFunctionAtStartPointPerformed;
	bool PointFittingSORMAnalysisPerformed;
	bool CurvatureFittingSORMAnalysisPerformed;
	bool CurvaturesFromSearchAlgorithmSORMAnalysisPerformed;

	// GFunValueAtStartPt analysis
	double GFunValueAtStartPt;

	// FORM analysis:
	double FORMReliabilityIndexBeta;
	double FORMProbabilityOfFailure_pf1;
	Vector designPoint_x_inOriginalSpace;
	Vector designPoint_u_inStdNormalSpace;
	Vector normalizedNegativeGradientVectorAlpha;
	Vector importanceVectorGamma;
	int numberOfIterationsToFindDesignPoint;
	
	// From Simulation analysis:
	double SimulationReliabilityIndexBeta;
	double SimulationProbabilityOfFailure_pfsim;
	double CoefficientOfVariationOfPfFromSimulation;
	int NumberOfSimulations;
	
	// From SORM analysis:
	double SORMCurvatureFittingBetaBreitung;
	double SORMCurvatureFittingPf2Breitung;
	double SORMPointFittingBetaBreitung;
	double SORMPointFittingPf2Breitung;
	double SORMUsingSearchBetaBreitung;
	double SORMUsingSearchPf2Breitung;
	Vector lastSearchDirection;
	int numberOfCurvatauresUsed;
	Vector secondLast_u;
	Vector secondLastAlpha;




protected:

private:
	int tag;
	char *expression;

};

#endif
