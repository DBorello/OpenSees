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
// $Date: 2001-06-13 05:06:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/designPoint/SearchWithStepSizeAndStepDirection.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#ifndef SearchWithStepSizeAndStepDirection_h
#define SearchWithStepSizeAndStepDirection_h

#include <FindDesignPoint.h>
#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <XuTransformation.h>
#include <GFunEvaluator.h>
#include <SensitivityEvaluator.h>
#include <Matrix.h>
#include <Vector.h>
#include <ReliabilityDomain.h>

class SearchWithStepSizeAndStepDirection : public FindDesignPoint
{

public:

	// Constructor and destructor
	SearchWithStepSizeAndStepDirection(
					int passedMaxNumberOfIterations, 
					double passedConvergenceCriterionE1, 
					double passedConvergenceCriterionE2,
					GFunEvaluator *passedGFunEvaluator,
					SensitivityEvaluator *passedSensitivityEvaluator,
					StepSizeRule *passedStepSizeRule,
					SearchDirection *passedSearchDirection,
					XuTransformation *passedXuTransformation);
	~SearchWithStepSizeAndStepDirection();
	
	int findDesignPoint(Vector *startPoint, 
						ReliabilityDomain *theReliabilityDomain);
	int approachDesignPointThroughOrthogonalSubspace(
						Vector *designPoint,
						int numberOfAxesInVector,
						Vector *principalAxes, 
						Vector *startPoint, 
						ReliabilityDomain *theReliabilityDomain);
	Vector get_x();
	Vector get_u();
	Vector get_alpha();
	Vector get_gamma();
	int getNumberOfIterations();
	Vector getSecondLast_u();
	Vector getSecondLast_alpha();
	Vector getLastSearchDirection();
	double getFirstGFunValue();

protected:

private:	

	// The reliability domain and tools for the analysis
	ReliabilityDomain *theReliabilityDomain;
	GFunEvaluator *theGFunEvaluator;
	SensitivityEvaluator *theSensitivityEvaluator;
	StepSizeRule *theStepSizeRule;
	SearchDirection *theSearchDirection;
	XuTransformation *theXuTransformation;

	// Private member functions to do the job
	int doTheActualSearch(bool doProjection);
	int doProjection(Vector uOld, Vector uNew);

	// Data members set when the object is created
	int maxNumberOfIterations;
	double convergenceCriterionE1;
	double convergenceCriterionE2;

	// Data members where the results are to be stored
	Vector x;
	Vector u;
	Vector alpha;
	Vector gamma;
	Vector uSecondLast;
	Vector alphaSecondLast;
	int i;
	Vector searchDirection;
	double Go;

	// Data members set through the call when a job is to be done
	Vector *startPoint;
	Vector *designPoint_uStar;
	int numberOfAxes;
	Vector *principalAxesPtr;

};

#endif
