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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/StaticSensitivityIntegrator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), July 2001
// Revised: 
//

#include <SensitivityIntegrator.h>
#include <StaticSensitivityIntegrator.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>

StaticSensitivityIntegrator::StaticSensitivityIntegrator(AnalysisModel *theModel, LinearSOE *theLinSOE)
:SensitivityIntegrator()
{
	// The analysis model is needed to get hold of dof's, elements, etc:
    theAnalysisModel = theModel;

	// The system of equations (SOE) is needed to zero out and 
	// assemble the right-hand side:
    theSOE = theLinSOE;
}


StaticSensitivityIntegrator::~StaticSensitivityIntegrator()
{
}




int
StaticSensitivityIntegrator::formRightHandSide(void)
{
	// In a non-linear static analysis the response is found 
	// by solving the equation:
	//             Pint = Pext
	// Then, the sensitivity vector 'v' is found by solving the equation:
	//       Kt * v = dPext/dh - dPint/dh 
	// where Kt is the tangent of the converged FE analysis. 
	// The main task here is to assemble this right-hand side. 


	// Zero out the old right-hand side
    theSOE->zeroB();


	// Loop through the FE_Elements and add dPint/dh contributions
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
		theSOE->addB( (-1)*elePtr->gradient(0),  elePtr->getID());
	}

	// For now, the dPext/dh term is not included. Coming soon. 

	return 0;

}



int
StaticSensitivityIntegrator::saveGradient(const Vector &v, int gradNum, int numGrads)
{
    DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();
    DOF_Group 	*dofPtr;

    while ( (dofPtr = theDOFGrps() ) != 0)  {
		dofPtr->setGradient(v,gradNum,numGrads);
	}
    

    return 0;
}





int 
StaticSensitivityIntegrator::commitGradient(int gradNumber)
{
    
	// Loop through the FE_Elements and set unconditional sensitivities
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
		elePtr->gradient(gradNumber);
	}

    return 0;
}













