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
// $Date: 2001-07-18 16:24:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/KrylovNewton.cpp,v $

// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// KrylovNewton.  KrylovNewton is a class which uses a Krylov
// subspace accelerator on the modified Newton method.
// The accelerator is described by Carlson and Miller in
// "Design and Application of a 1D GWMFE Code"
// from SIAM Journal of Scientific Computing (Vol. 19, No. 3,
// pp. 728-765, May 1998)

#include <KrylovNewton.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <fstream.h>

// Constructor
KrylovNewton::KrylovNewton(int theTangentToUse)
:EquiSolnAlgo(EquiALGORITHM_TAGS_KrylovNewton),
 theTest(0), tangent(theTangentToUse),
 r(0), v(0), Av(0),
 numTests(0), numEqns(0)
{

}

KrylovNewton::KrylovNewton(ConvergenceTest &theT, int theTangentToUse)
:EquiSolnAlgo(EquiALGORITHM_TAGS_KrylovNewton),
 theTest(&theT), tangent(theTangentToUse),
 r(0), v(0), Av(0),
 numTests(0), numEqns(0)
{

}

// Destructor
KrylovNewton::~KrylovNewton()
{
  if (r != 0) {
    for (int i = 0; i < numTests; i++)
      delete r[i];
    delete [] r;
  }

  if (v != 0) {
    for (int i = 0; i < numTests; i++)
      delete v[i];
    delete [] v;
  }

  if (Av != 0) {
    for (int i = 0; i < numTests; i++)
      delete Av[i];
    delete [] Av;
  }
}

void 
KrylovNewton::setTest(ConvergenceTest &newTest)
{
  theTest = &newTest;
}

int 
KrylovNewton::solveCurrentStep(void)
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass
  AnalysisModel   *theAnaModel = this->getAnalysisModelPtr();
  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
  LinearSOE  *theSOE = this->getLinearSOEptr();
  
  if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
      || (theTest == 0)){
    cerr << "WARNING KrylovNewton::solveCurrentStep() - setLinks() has";
    cerr << " not been called - or no ConvergenceTest has been set\n";
    return -5;
  }	
  
  numTests = 2 + theTest->getMaxNumTests();
  numEqns  = theSOE->getNumEqn();
  
  if (r == 0) {
    r = new Vector*[numTests];
    for (int i = 0; i < numTests; i++)
      r[i] = new Vector(numEqns);
  }

  if (v == 0) {
    v = new Vector*[numTests];
    for (int i = 0; i < numTests; i++)
      v[i] = new Vector(numEqns);
  }

  if (Av == 0) {
    Av = new Vector*[numTests];
    for (int i = 0; i < numTests; i++)
      Av[i] = new Vector(numEqns);
  }

  // set itself as the ConvergenceTest objects EquiSolnAlgo
  theTest->setEquiSolnAlgo(*this);
  if (theTest->start() < 0) {
    cerr << "KrylovNewton::solveCurrentStep() -";
    cerr << "the ConvergenceTest object failed in start()\n";
    return -3;
    }
  
  // Evaluate system residual R(y_0)
  if (theIntegrator->formUnbalance() < 0) {
    cerr << "WARNING KrylovNewton::solveCurrentStep() -";
    cerr << "the Integrator failed in formUnbalance()\n";	
    return -2;
  }
  
  // Evaluate system Jacobian J = R'(y)|y_0
  if (theIntegrator->formTangent(tangent) < 0){
    cerr << "WARNING KrylovNewton::solveCurrentStep() -";
    cerr << "the Integrator failed in formTangent()\n";
    return -1;
  }    
  
  // Solve for residual f(y_0) = J^{-1} R(y_0)
  if (theSOE->solve() < 0) {
    cerr << "WARNING KrylovNewton::solveCurrentStep() -";
    cerr << "the LinearSysOfEqn failed in solve()\n";	
    return -3;
  }

  // Update y ... y_1 = y_0 + v_1 --> (y_0 is zero) --> y_1 = v_1
  // --> (v_1 = f(y_0)) --> y_1 = f(y_0)
  *(r[0]) = theSOE->getX();

  // Store first update v_1 = r_0
  *(v[1]) = *(r[0]);

  int result = -1;
  int k = 1;

  do {

    // Update system with v_k
    if (theIntegrator->update(*(v[k])) < 0) {
      cerr << "WARNING KrylovNewton::solveCurrentStep() -";
      cerr << "the Integrator failed in update()\n";	
      return -4;
    }	

    // Evaluate system residual R(y_k)
    if (theIntegrator->formUnbalance() < 0) {
      cerr << "WARNING KrylovNewton::solveCurrentStep() -";
      cerr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }
    
    // Solve for residual f(y_k) = J^{-1} R(y_k)
    if (theSOE->solve() < 0) {
      cerr << "WARNING KrylovNewton::solveCurrentStep() -";
      cerr << "the LinearSysOfEqn failed in solve()\n";	
      return -3;
    }

    // Store residual r_k \equiv f(y_k)
    *(r[k]) = theSOE->getX();

    // Compute Av_k = f(y_{k-1}) - f(y_k) = r_{k-1} - r_k
    // NOTE: Not using Av[0] so that notation follows paper
    *(Av[k]) = *(r[k-1]);
    Av[k]->addVector(1.0, *(r[k]), -1.0);

    // Solve least squares A w_{k+1} = r_k
    if (this->leastSquares(k) < 0) {
      cerr << "WARNING KrylovNewton::solveCurrentStep() -";
      cerr << "the Integrator failed in leastSquares()\n";
      return -1;
    }		    

    result = theTest->test();
    this->record(k++);

  } while (result == -1);
  
  if (result == -2) {
    cerr << "KrylovNewton::solveCurrentStep() -";
    cerr << "the ConvergenceTest object failed in test()\n";
    return -3;
  }
  
  // note - if postive result we are returning what the convergence
  // test returned which should be the number of iterations
  return result;
}

ConvergenceTest *
KrylovNewton::getTest(void)
{
  return theTest;
}

int
KrylovNewton::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
KrylovNewton::recvSelf(int cTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
KrylovNewton::Print(ostream &s, int flag)
{
  if (flag == 0)
    s << "KrylovNewton\n";
}

int
KrylovNewton::leastSquares(int k)
{
  // NOTE: Not using Av[0] so that notation follows paper
  const int maxK = 30;

  if (k > maxK) {
    cerr << "WARNING KrylovNewton::leastSquares() -";
    cerr << "Number of iterations exceeds matrix limits";
    return -2;
  }

  static double workA[maxK*maxK];
  static double workb[maxK];
  static double workc[maxK];

  Matrix A(workA,k,k);
  Vector b(workb,k);
  Vector c(workc,k);

  int i,j;

  // Form normal equations (Av_i,Av_j) c_j = (Av_i,r_k)
  for (i = 1; i <= k; i++) {
    for (j = 1; j <= k; j++)
      A(i-1,j-1) = *(Av[i]) ^ *(Av[j]);
    b(i-1) = *(Av[i]) ^ *(r[k]);
  }

  // Solve normal equations
  if (A.Solve(b,c) < 0) {
    cerr << "WARNING KrylovNewton::leastSquares() -";
    cerr << "solution of normal equations failed in Matrix::Solve()\n";
    return -1;
  }

  // v_{k+1} = w_{k+1} + q_{k+1}
  *(v[k+1]) = *(r[k]);

  double cj;
  for (j = 1; j <= k; j++) {

    cj = c(j-1);

    // Compute w_{k+1} = c_1 v_1 + ... + c_k v_k
    v[k+1]->addVector(1.0, *(v[j]), cj);

    // Compute least squares residual q_{k+1} = r_k - c_1 Av_1 - ... - c_k Av_k
    v[k+1]->addVector(1.0, *(Av[j]), -cj);
  }

  return 0;
}
