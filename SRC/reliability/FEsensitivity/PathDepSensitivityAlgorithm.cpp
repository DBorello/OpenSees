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
// $Date: 2001-07-31 22:11:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/PathDepSensitivityAlgorithm.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), July 2001
// Revised: 
//

#include <PathDepSensitivityAlgorithm.h>
#include <SensitivityAlgorithm.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <EquiSolnAlgo.h>
#include <ReliabilityDomain.h>
#include <RandomVariablePositioner.h>


PathDepSensitivityAlgorithm::PathDepSensitivityAlgorithm(ReliabilityDomain *passedReliabilityDomain,
										   EquiSolnAlgo *passedAlgorithm,
										   SensitivityIntegrator *passedSensitivityIntegrator)
:SensitivityAlgorithm()
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
}




PathDepSensitivityAlgorithm::~PathDepSensitivityAlgorithm()
{
}



int 
PathDepSensitivityAlgorithm::computeGradients(void)
{
	
	// Get pointer to the system of equations (SOE)
	LinearSOE *theSOE = theAlgorithm->getLinearSOEptr();

	// Get number of random variables and random variable positioners
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	int npos = theReliabilityDomain->getNumberOfRandomVariablePositioners();

	// Initial declarations
	RandomVariablePositioner *theRandomVariablePositioner;
	Vector vectorOfActivePositioners(nrv);

	// Loop where one gradient vector is computed per iteration
	for ( int gradNumber=1; gradNumber<=nrv; gradNumber++ )  {

		// Clear sensitivity flags
		for ( int i=1; i<=npos; i++ ) {
			theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
			theRandomVariablePositioner->setSensitivityFlag(0);
		}

		// Loop over r.v. positioners to check which are "active" this time
		// and set flag for phase 1
		int numberOfActivePositioners = 0;
		for ( int rvposNumber=1; rvposNumber<=npos; rvposNumber++ ) {

			// Get the random variable positioner and its rv#
			theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(rvposNumber);
			int rvNumber = theRandomVariablePositioner->getRvNumber();

			// Check if rv# of this rv positioner is the same as the gradient#
			if ( rvNumber==gradNumber ) {

				// Pick up the number of this random variable positioner
				vectorOfActivePositioners(numberOfActivePositioners) = rvposNumber;
				numberOfActivePositioners++;

				// Get the random variable positioner
				theRandomVariablePositioner = theReliabilityDomain->
					getRandomVariablePositionerPtr(rvposNumber);

				// Set sensitivity flag so that this one contributes to the RHS
				theRandomVariablePositioner->setSensitivityFlag(1);

			} // End if rv# == gradient#

		} // End loop over r.v. positioners to check if it contributes

		// Form new right-hand side
		theSensitivityIntegrator->formRightHandSide();

		// Solve the system of equation with the new right-hand side
		theSOE->solve();

		// Commit 'v' to the nodes for a "sensNodeDisp node? dof?" command
		theSensitivityIntegrator->saveGradient( theSOE->getX(), gradNumber, nrv );

		// Commit unconditional history variables for path-dependent problems
		theSensitivityIntegrator->commitGradient(gradNumber);

	} // End loop for each gradient

    return 0;
}

