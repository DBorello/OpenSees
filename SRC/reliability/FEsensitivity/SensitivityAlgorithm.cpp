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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-09-05 17:54:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/SensitivityAlgorithm.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), July 2001
// Revised:	haukaas 08/19/01 (modifications for Release 1.2 of OpenSees)
//

#include <SensitivityAlgorithm.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <EquiSolnAlgo.h>
#include <ReliabilityDomain.h>
#include <RandomVariablePositioner.h>


SensitivityAlgorithm::SensitivityAlgorithm(ReliabilityDomain *passedReliabilityDomain,
										   EquiSolnAlgo *passedAlgorithm,
										   SensitivityIntegrator *passedSensitivityIntegrator,
										   bool passedTrueIfPathDependent)
{
	// The reliability domain is needed to get hold 
	// of the random variable positioners:
	theReliabilityDomain = passedReliabilityDomain;

	// The finite element equation solution algorithm is 
	// needed to get hold of the system of equations (SOE):
	theAlgorithm = passedAlgorithm;

	// The sensitivity integrator is needed to assemble the 
	// new right-hand side of the system of equations:
	theSensitivityIntegrator = passedSensitivityIntegrator;

	// Flag to say whether this algorithm is for path-
	// dependent problems or not
	trueIfPathDependent = passedTrueIfPathDependent;

}




SensitivityAlgorithm::~SensitivityAlgorithm()
{
}



int 
SensitivityAlgorithm::computeGradients(void)
{
	
	// Get pointer to the system of equations (SOE)
	LinearSOE *theSOE = theAlgorithm->getLinearSOEptr();

	// Get pointer to incremental integrator
	IncrementalIntegrator *theIncInt = theAlgorithm->getIncrementalIntegratorPtr();

	// Form current tangent at converged state
	if (theIncInt->formTangent(CURRENT_TANGENT) < 0){
		cerr << "WARNING SensitivityAlgorithm::computeGradients() -";
		cerr << "the Integrator failed in formTangent()\n";
		return -1;
	}

	// Get number of random variables and random variable positioners
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	int npos = theReliabilityDomain->getNumberOfRandomVariablePositioners();

	// Initial declarations
	RandomVariablePositioner *theRandomVariablePositioner;

	// Loop where one gradient vector is computed per iteration
	for ( int gradNumber=1; gradNumber<=nrv; gradNumber++ )  {

		// Clear sensitivity flags
		for ( int i=1; i<=npos; i++ ) {
			theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
			theRandomVariablePositioner->setSensitivityFlag(0);
		}

		// Loop over r.v. positioners to check which are "active" this time
		// and set flag for phase 1
		for ( int rvposNumber=1; rvposNumber<=npos; rvposNumber++ ) {

			// Get the random variable positioner and its rv#
			theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(rvposNumber);
			int rvNumber = theRandomVariablePositioner->getRvNumber();

			// Check if rv# of this rv positioner is the same as the gradient#
			if ( rvNumber==gradNumber ) {

				// Set sensitivity flag so that this one contributes to the RHS
				theRandomVariablePositioner->setSensitivityFlag(1);

			} // End if rv# == gradient#

		} // End loop over r.v. positioners to check if it contributes

		// Form new right-hand side
		theSensitivityIntegrator->formRightHandSide();

		// Solve the system of equation with the new right-hand side
		theSOE->solve();

		// Save 'v' to the nodes for a "sensNodeDisp node? dof?" command
		theSensitivityIntegrator->saveGradient( theSOE->getX(), gradNumber, nrv );

		// Commit unconditional history variables for path-dependent problems
		if (trueIfPathDependent) {
			theSensitivityIntegrator->commitGradient(gradNumber);
		}

	} // End loop for each gradient

    return 0;
}

bool 
SensitivityAlgorithm::isPathDependent(void)
{
	return (trueIfPathDependent);
}