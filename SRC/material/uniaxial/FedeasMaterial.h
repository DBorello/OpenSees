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
// $Date: 2001-07-29 22:59:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FedeasMaterial.h,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasMaterial. FedeasMaterial provides a FORTRAN interface
// for programming uniaxial material models, using the subroutine
// interface from the FEDEAS ML1D library.

#ifndef FedeasMaterial_h
#define FedeasMaterial_h

// Pre-defined values for the matType variable
#define FEDEAS_Bond1		1
#define FEDEAS_Bond2		2
#define FEDEAS_Concrete1	3
#define FEDEAS_Concrete2	4
#define FEDEAS_Concrete3	5
#define FEDEAS_Hardening	6
#define FEDEAS_Hysteretic1	7
#define FEDEAS_Hysteretic2	8
#define FEDEAS_Steel1		9
#define FEDEAS_Steel2		10
// Add more #defines as needed

#include <UniaxialMaterial.h>

class FedeasMaterial : public UniaxialMaterial
{
 public:
  FedeasMaterial(int tag, int classTag, int type, int numHV, int numData);
  FedeasMaterial(int classTag, int type, int numHV, int numData);
  virtual ~FedeasMaterial();
  
  virtual int setTrialStrain(double strain, double strainRate = 0.0);
  virtual int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0);
  virtual double getStrain(void);
  virtual double getStress(void);
  virtual double getTangent(void);
  
  virtual int commitState(void);
  virtual int revertToLastCommit(void);    
  virtual int revertToStart(void);        
  
  // WARNING -- if you wish to override any method in this base class, you must
  // also override the getCopy method to return a pointer to the derived class!!!
  virtual UniaxialMaterial *getCopy(void);
  
  virtual int sendSelf(int commitTag, Channel &theChannel);  
  virtual int recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker);    
  
  virtual void Print(ostream &s, int flag = 0);
  
 protected:
  // Invokes the FORTRAN subroutine
  virtual int invokeSubroutine(int ist);
  
  double *data;		// Material parameters array
  double *hstv;		// History array: first half is committed, second half is trial
  
  int numData;		// Number of material parameters
  int numHstv;		// Number of history variables
  
  double epsilonP;	// Committed strain
  double sigmaP;	// Committed stress
  
 private:
  int matType;		// Material type ... tells which subroutine to call
  
  double epsilon;	// Trial strain
  double sigma;		// Trial stress
  double tangent;	// Trial tangent
};

#endif
