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
                                                                        
// $Revision: 1.2 $
// $Date: 2001-11-20 02:45:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/SymBandEigenSolver.cpp,v $

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for
// SymBandEigenSolver, which computes the eigenvalues and eigenvectors
// of a symmetric banded matrix using the LAPACK subroutine dsbevx.

#include <SymBandEigenSolver.h>
#include <f2c.h>
#include <math.h>
#include <stdio.h>
#include <AnalysisModel.h>
#include <DOF_GrpIter.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <Integrator.h>

SymBandEigenSolver::SymBandEigenSolver()
:EigenSolver(EigenSOLVER_TAGS_SymBandEigenSolver),
 theSOE(0), numModes(0), eigenvalue(0), eigenvector(0), eigenV(0)
{

}

SymBandEigenSolver::~SymBandEigenSolver()
{
  if (eigenvalue != 0)
    delete [] eigenvalue;
  if (eigenvector != 0)
    delete [] eigenvector;
  if (eigenV != 0)
    delete eigenV;
}

#ifdef _WIN32

extern "C" int _stdcall DSBEVX(char *jobz, unsigned int *a, char *range, unsigned int *b, 
							   char *uplo, unsigned int *c, int *n, int *kd,
			       double *ab, int *ldab, double *q, int *ldq,
			       double *vl, double *vu, int *il, int *iu, double *abstol,
			       int *m, double *w, double *z, int *ldz,
			       double *work, int *iwork, int *ifail, int *info);

#else

extern "C" int dsbevx_(char *jobz, char *range, char *uplo, int *n, int *kd,
		       double *ab, int *ldab, double *q, int *ldq,
		       double *vl, double *vu, int *il, int *iu, double *abstol,
		       int *m, double *w, double *z, int *ldz,
		       double *work, int *iwork, int *ifail, int *info);

#endif

int
SymBandEigenSolver::solve(int nModes)
{
  if (theSOE == 0) {
    g3ErrorHandler->warning("%s -- no EigenSOE has been set yet",
			    "SymBandEigenSolver::solve()");
    return -1;
  }
  
  // Set number of modes
  numModes = nModes;

  // Number of equations
  int n = theSOE->size;

  // Check for quick return
  if (numModes < 1) {
    numModes = 0;
    return 0;
  }

  // Simple check
  if (numModes > n)
    numModes = n;

  // Allocate storage for eigenvalues
  if (eigenvalue != 0)
    delete [] eigenvalue;
  eigenvalue = new double [n];

  // Real work array (see LAPACK dsbevx subroutine documentation)
  double *work = new double [7*n];

  // Integer work array (see LAPACK dsbevx subroutine documentation)
  int *iwork = new int [5*n];

  // Leading dimension of eigenvectors
  int ldz = n;

  // Allocate storage for eigenvectors
  if (eigenvector != 0)
    delete [] eigenvector;
  eigenvector = new double [ldz*numModes];

  // Number of superdiagonals
  int kd = theSOE->numSuperD;

  // Matrix data
  double *ab = theSOE->A;

  // Leading dimension of the matrix
  int ldab = kd + 1;

  // Leading dimension of q
  int ldq = n;

  // Orthogonal matrix used in reduction to tridiagonal form
  // (see LAPACK dsbevx subroutine documentation)
  double *q = new double [ldq*n];

  // Index ranges [1,numModes] of eigenpairs to compute
  int il = 1;
  int iu = numModes;

  // Compute eigenvalues and eigenvectors
  char *jobz = "V";

  // Selected eigenpairs are based on index range [il,iu]
  char *range = "I";

  // Upper triagle of matrix is stored
  char *uplo = "U";
  
  // Return value
  int *ifail = new int [n];
  int info = 0;

  // Number of eigenvalues returned
  int m = 0;

  // Not used
  double vl = 0.0;
  double vu = 1.0;

  // Not used ... I think!
  double abstol = -1.0;

  // Call the LAPACK eigenvalue subroutine
#ifdef _WIN32
  unsigned int sizeC = 1;
  DSBEVX(jobz, &sizeC, range, &sizeC, uplo, &sizeC, &n, &kd, ab, &ldab,
	 q, &ldq, &vl, &vu, &il, &iu, &abstol, &m,
	 eigenvalue, eigenvector, &ldz, work, iwork, ifail, &info);
#else
  dsbevx_(jobz, range, uplo, &n, &kd, ab, &ldab,
	  q, &ldq, &vl, &vu, &il, &iu, &abstol, &m,
	  eigenvalue, eigenvector, &ldz, work, iwork, ifail, &info);
#endif

  delete [] q;
  delete [] work;
  delete [] iwork;
  delete [] ifail;

  if (info < 0) {
    g3ErrorHandler->warning("%s -- invalid argument number %d passed to LAPACK dsbevx",
			    "SymBandEigenSolver::solve()", -info);
    return info;
  }

  if (info > 0) {
    g3ErrorHandler->warning("%s -- LAPACK dsbevx returned error code %d",
			    "SymBandEigenSolver::solve()", info);
    return -info;
  }

  if (m < numModes) {
    g3ErrorHandler->warning("%s -- LAPACK dsbevx only computed %d eigenvalues, %d were requested",
			    "SymBandEigenSolver::solve()", m, numModes);
    numModes = m;
  }

  theSOE->factored = true;

  return 0;
}

int
SymBandEigenSolver::setEigenSOE(SymBandEigenSOE &theBandSOE)
{
  theSOE = &theBandSOE;

  return 0;
}

const Vector &
SymBandEigenSolver::getEigenvector(int mode)
{
  if (mode < 1 || mode > numModes) {
    g3ErrorHandler->warning("%s -- mode %d is out of range (1 - %d)",
			    "SymBandEigenSolver::getEigenvector()", mode, numModes);
    eigenV->Zero();
    return *eigenV;  
  }
  
  int size = theSOE->size;

  int index = (mode - 1) * size;
  
  Vector &vec = *eigenV;
  if (eigenvector != 0) {
    for (int i = 0; i < size; i++) {
      vec(i) = eigenvector[index++];
    }	
  }
  else {
    g3ErrorHandler->warning("%s -- eigenvectors not yet computed",
			    "SymBandEigenSolver::getEigenvector()");
    eigenV->Zero();
  }      
  
  return *eigenV;  
}

double
SymBandEigenSolver::getEigenvalue(int mode)
{
  if (mode < 1 || mode > numModes) {
    g3ErrorHandler->warning("%s -- mode %d is out of range (1 - %d)",
			    "SymBandEigenSolver::getEigenvalue()", mode, numModes);
    return 0.0;
  }
  
  if (eigenvalue != 0)
    return eigenvalue[mode-1];
  else {
    g3ErrorHandler->warning("%s -- eigenvalues not yet computed",
			    "SymBandEigenSolver::getEigenvalue()");
    return 0.0;
  }      
}

int
SymBandEigenSolver::setSize()
{
  int size = theSOE->size;    

  if (eigenV == 0 || eigenV->Size() != size) {
    if (eigenV != 0)
      delete eigenV;
    
    eigenV = new Vector(size);
    if (eigenV == 0 || eigenV->Size() != size) {
      g3ErrorHandler->warning("%s -- ran out of memory for eigenvector of size %d",
			      "SymBandEigenSolver::setSize()", size);
      return -2;	    
    }
  }
  
  return 0;
}

int    
SymBandEigenSolver::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
SymBandEigenSolver::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  return 0;
}
