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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/SensitivityByFiniteDifference.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#ifndef SensitivityByFiniteDifference_h
#define SensitivityByFiniteDifference_h

#include <SensitivityEvaluator.h>
#include <Vector.h>
#include <ReliabilityDomain.h>
#include <GFunEvaluator.h>


class SensitivityByFiniteDifference : public SensitivityEvaluator
{

public:
	SensitivityByFiniteDifference(GFunEvaluator *passedGFunEvaluator,
				ReliabilityDomain *passedReliabilityDomain);
	~SensitivityByFiniteDifference();

	int		evaluate_grad_g(double gFunValue, Vector passed_x);
	Vector	get_grad_g();

protected:

private:
	Vector *grad_g;
	GFunEvaluator *theGFunEvaluator;
	ReliabilityDomain *theReliabilityDomain;

};

#endif
