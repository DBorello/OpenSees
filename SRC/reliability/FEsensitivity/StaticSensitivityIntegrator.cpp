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
// $Date: 2001-08-20 00:37:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/StaticSensitivityIntegrator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), July 2001
// Revised: haukaas 08/19/01 (modifications for Release 1.2 of OpenSees)
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
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>


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


	// Loop through the loadPatterns and add the dPext/dh contributions
	Vector oneDimVectorWithOne(1);
	oneDimVectorWithOne(0) = 1.0;
	ID oneDimID(1);
	Node *aNode;
	DOF_Group *aDofGroup;
	int nodeNumber, dofNumber, relevantID, i, sizeRandomLoads, numRandomLoads;
	/////////////////////////////////////////////////
	LoadPattern *loadPatternPtr;
	Domain *theDomain = theAnalysisModel->getDomainPtr();
	LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    while((loadPatternPtr = thePatterns()) != 0) {

		const Vector &randomLoads = loadPatternPtr->gradient(true,0);
		sizeRandomLoads = randomLoads.Size();
		if (sizeRandomLoads == 1) {
			// No random loads in this load pattern
		}
		else {
			// Random loads: add contributions to the 'B' vector
			numRandomLoads = (int)(sizeRandomLoads/2);
			for (i=0; i<numRandomLoads*2; i=i+2) {
				nodeNumber = (int)randomLoads(i);
				dofNumber = (int)randomLoads(i+1);
				aNode = theDomain->getNode(nodeNumber);
				aDofGroup = aNode->getDOF_GroupPtr();
				const ID &anID = aDofGroup->getID();
				relevantID = anID(dofNumber-1);
				oneDimID(0) = relevantID;
				theSOE->addB(oneDimVectorWithOne, oneDimID);
			}
		}
	}

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













