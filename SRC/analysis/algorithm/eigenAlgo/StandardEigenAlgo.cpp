/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2001-11-19 22:41:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/eigenAlgo/StandardEigenAlgo.cpp,v $
                                                                        
// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition of StandardEigenAlgo.
// StandardEigenAlgo is a class which performs a eigen solution algorithm
// to solve standard eigenvalue equations. It is not expected that 
// this class will have subclasses.

#include <StandardEigenAlgo.h>
#include <AnalysisModel.h>
#include <EigenAnalysis.h>
#include <EigenIntegrator.h>
#include <EigenSOE.h>
#include <iostream.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>

StandardEigenAlgo::StandardEigenAlgo()
  :EigenAlgorithm(EigenALGORITHM_TAGS_Standard)
{
  // do nothing here.
}

StandardEigenAlgo::~StandardEigenAlgo()
{
  // do nothing here.
}

int 
StandardEigenAlgo::solveCurrentStep(int numModes)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  EigenSOE *theSOE = this->getEigenSOEptr();
  EigenIntegrator *theIntegrator = this->getEigenIntegratorPtr();
  
  if ((theModel == 0) || (theIntegrator == 0) || (theSOE == 0)) {

    g3ErrorHandler->warning("%s -- setLinks() has not been called",
			    "StandardEigenAlgo::solverCurrentStep()");
    return -1;
  }
  
  if (theIntegrator->formK() < 0) {
    g3ErrorHandler->warning("%s -- the Integrator failed in formK()",
			    "StandardEigenAlgo::solverCurrentStep()");
    return -2;
  }
  
  if (theSOE->solve(numModes) < 0) {
    g3ErrorHandler->warning("%s -- the EigenSOE failed in solve()",
			    "StandardEigenAlgo::solverCurrentStep()");
    return -4;
  }
  
  // now set the eigenvalues and eigenvectors in the model
  theModel->setNumEigenvectors(numModes);
  Vector theEigenvalues(numModes);
  for (int i = 1; i <= numModes; i++) {
    theEigenvalues[i-1] = theSOE->getEigenvalue(i);
    theModel->setEigenvector(i, theSOE->getEigenvector(i));
  }    
  theModel->setEigenvalues(theEigenvalues);
  
  return 0;
}

int 
StandardEigenAlgo::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
StandardEigenAlgo::recvSelf(int cTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
StandardEigenAlgo::Print(ostream &s, int flag)
{
  s << "\tStandardEigenAlgo\n";
}
