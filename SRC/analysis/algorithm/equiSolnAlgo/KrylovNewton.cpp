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

// $Revision: 1.3 $
// $Date: 2001-08-27 22:14:01 $
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
 v(0), Av(0), AvData(0), rData(0), work(0),
 maxTests(0), numEqns(0), numVecs(0), lwork(0)
{

}

KrylovNewton::KrylovNewton(ConvergenceTest &theT, int theTangentToUse)
:EquiSolnAlgo(EquiALGORITHM_TAGS_KrylovNewton),
 theTest(&theT), tangent(theTangentToUse),
 v(0), Av(0), AvData(0), rData(0), work(0),
 maxTests(0), numEqns(0), numVecs(0), lwork(0)
{

}

// Destructor
KrylovNewton::~KrylovNewton()
{
  if (v != 0) {
    for (int i = 0; i < numVecs; i++)
      delete v[i];
    delete [] v;
  }

  if (Av != 0) {
    for (int i = 0; i < numVecs; i++)
      delete Av[i];
    delete [] Av;
  }

  if (AvData != 0)
    delete [] AvData;

  if (rData != 0)
    delete [] rData;

  if (work != 0)
    delete [] work;
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

  // Get size information from convergence test and SOE
  maxTests = theTest->getMaxNumTests();
  numVecs  = 2 + maxTests;
  numEqns  = theSOE->getNumEqn();

  if (v == 0) {
    v = new Vector*[numVecs];
    for (int i = 0; i < numVecs; i++)
      v[i] = new Vector(numEqns);
  }

  if (Av == 0) {
    Av = new Vector*[numVecs];
    for (int i = 0; i < numVecs; i++)
      Av[i] = new Vector(numEqns);
  }

  if (AvData == 0)
    AvData = new double [maxTests*numEqns];

  if (rData == 0)
    // The LAPACK least squares subroutine overwrites the RHS vector
    // with the solution vector ... these vectors are not the same
    // size, so we need to use the max size
    rData = new double [(numEqns > maxTests) ? numEqns : maxTests];

  // Length of work vector should be >= 2*min(numEqns,maxTests)
  // See dgels subroutine documentation
  lwork = 2 * ((numEqns < maxTests) ? numEqns : maxTests);
  
  if (work == 0)
    work = new double [lwork];

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

  // Not using Av[0], for the subspace vectors, so
  // use it for the residual vector
  Vector &r = *(Av[0]);

  // Update y ... y_1 = y_0 + v_1 --> (y_0 is zero) --> y_1 = v_1
  // --> (v_1 = f(y_0)) --> y_1 = f(y_0)
  r = theSOE->getX();

  // Store first update v_1 = r_0
  *(v[1]) = r;

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

    // Compute Av_k = f(y_{k-1}) - f(y_k) = r_{k-1} - r_k
    // NOTE: Not using Av[0] so that notation follows paper
    *(Av[k]) = r;
    r = theSOE->getX();
    Av[k]->addVector(1.0, r, -1.0);

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
KrylovNewton::recvSelf(int cTag, Channel &theChannel, 
					   FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
KrylovNewton::Print(ostream &s, int flag)
{
  s << "KrylovNewton\n";
}

#ifdef _WIN32

extern "C" int _stdcall DGELS(char *T, unsigned int *SZ, int *M, int *N, int *NRHS,
			      double *A, int *LDA, double *B, int *LDB,
			      double *WORK, int *LWORK, int *INFO);

#else

extern "C" int dgels_(char *T, int *M, int *N, int *NRHS,
		      double *A, int *LDA, double *B, int *LDB,
		      double *WORK, int *LWORK, int *INFO);

#endif

int
KrylovNewton::leastSquares(int k)
{
  int i,j;

  // Put subspace vectors into AvData
  Matrix A(AvData, numEqns, k);
  for (i = 0; i < k; i++) {
    Vector &Ai = *(Av[i+1]);
    for (j = 0; j < numEqns; j++)
      A(j,i) = Ai(j);
  }

  // Reference to the residual vector
  Vector &r = *(Av[0]);

  // Put residual vector into rData (need to save r for later!)
  Vector B(rData, numEqns);
  B = r;
  
  // No transpose
  char *trans = "N";

  // The number of right hand side vectors
  int nrhs = 1;

  // Leading dimension of the right hand side vector
  int ldb = (numEqns > k) ? numEqns : k;

  // Subroutine error flag
  int info = 0;

  // Call the LAPACK least squares subroutine
#ifdef _WIN32
  unsigned int sizeC = 1;
  DGELS(trans, &sizeC, &numEqns, &k, &nrhs, AvData, &numEqns, rData, &ldb, work, &lwork, &info);
#else
  dgels_(trans, &numEqns, &k, &nrhs, AvData, &numEqns, rData, &ldb, work, &lwork, &info);
#endif
  
  // Check for error returned by subroutine
  if (info < 0) {
    cerr << "WARNING KrylovNewton::leastSquares() - \n";
    cerr << "error code " << info << " returned by LAPACK dgels\n";
    return info;
  }
  
  // v_{k+1} = w_{k+1} + q_{k+1}
  *(v[k+1]) = r;
  
  // Compute the correction vector
  double cj;
  for (j = 1; j <= k; j++) {
    
    // Solution to least squares is written to rData
    cj = rData[j-1];
    
    // Compute w_{k+1} = c_1 v_1 + ... + c_k v_k
    v[k+1]->addVector(1.0, *(v[j]), cj);
    
    // Compute least squares residual q_{k+1} = r_k - c_1 Av_1 - ... - c_k Av_k
    v[k+1]->addVector(1.0, *(Av[j]), -cj);
  }
  
  return 0;
}

