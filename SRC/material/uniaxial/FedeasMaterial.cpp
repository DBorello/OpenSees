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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FedeasMaterial.cpp,v $
                                                                        
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasMaterial. FedeasMaterial provides a FORTRAN interface
// for programming uniaxial material models, using the subroutine
// interface from the FEDEAS ML1D library.

#include <G3Globals.h>
#include <FedeasMaterial.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <string.h>
#include <stdlib.h>

FedeasMaterial::FedeasMaterial(int tag, int classTag, int type,
			       int nhv, int ndata)
  :UniaxialMaterial(tag,classTag),
   data(0), hstv(0), numData(ndata), numHstv(nhv),
   epsilonP(0.0), sigmaP(0.0), matType(type),
   epsilon(0.0), sigma(0.0), tangent(0.0)
{
  if (numHstv < 0)
    numHstv = 0;
  
  if (numHstv > 0) {
    // Allocate history array
    hstv = new double[2*numHstv];
    if (hstv == 0)
      g3ErrorHandler->fatal("%s -- failed to allocate history array -- type %d",
			    "FedeasMaterial::FedeasMaterial", matType);
  }
  
  if (numData < 0)
    numData = 0;
  
  if (numData > 0) {
    // Allocate material parameter array
    data = new double[numData];
    if (data == 0)
      g3ErrorHandler->fatal("%s -- failed to allocate data array -- type %d",
			    "FedeasMaterial::FedeasMaterial", matType);
  }
  
  // Initialize history variables
  this->revertToStart();
}

FedeasMaterial::FedeasMaterial(int classTag, int type, int nhv, int ndata)
  :UniaxialMaterial(0,classTag),
   data(0), hstv(0), numData(ndata), numHstv(nhv),
   epsilonP(0.0), sigmaP(0.0), matType(type),
   epsilon(0.0), sigma(0.0), tangent(0.0)
{
  if (numHstv < 0)
    numHstv = 0;
  
  if (numHstv > 0) {
    // Allocate history array
    hstv = new double[2*numHstv];
    if (hstv == 0)
      g3ErrorHandler->fatal("%s -- failed to allocate history array -- type %d",
			    "FedeasMaterial::FedeasMaterial", matType);
  }
  
  if (numData < 0)
    numData = 0;
  
  if (numData > 0) {
    // Allocate material parameter array
    data = new double[numData];
    if (data == 0)
      g3ErrorHandler->fatal("%s -- failed to allocate data array -- type %d",
			    "FedeasMaterial::FedeasMaterial", matType);
  }
  
  // Initialize history variables
  this->revertToStart();
}

FedeasMaterial::~FedeasMaterial()
{
  if (hstv != 0)
    delete [] hstv;
  
  if (data != 0)
    delete [] data;
}

int
FedeasMaterial::setTrialStrain(double strain, double strainRate)
{
  // Store the strain
  epsilon = strain;
  
  // Tells subroutine to do normal operations for stress and tangent
  int ist = 1;
  
  // Call the subroutine
  return this->invokeSubroutine(ist);
}

int
FedeasMaterial::setTrial(double strain, double &stress, double &stiff, double strainRate)
{
  // Store the strain
  epsilon = strain;
  
  // Tells subroutine to do normal operations for stress and tangent
  int ist = 1;
  
  // Call the subroutine
  int res = this->invokeSubroutine(ist);
  
  stress = sigma;
  stiff = tangent;
  
  return res;
}

double
FedeasMaterial::getStrain(void)
{
  return epsilon;
}

double
FedeasMaterial::getStress(void)
{
  return sigma;
}

double
FedeasMaterial::getTangent(void)
{
  return tangent;
}

int
FedeasMaterial::commitState(void)
{
  // Set committed values equal to corresponding trial values
  for (int i = 0; i < numHstv; i++)
    hstv[i] = hstv[i+numHstv];
  
  epsilonP = epsilon;
  sigmaP = sigma;
  
  return 0;
}

int
FedeasMaterial::revertToLastCommit(void)
{
  // Set trial values equal to corresponding committed values
  for (int i = 0; i < numHstv; i++)
    hstv[i+numHstv] = hstv[i];
  
  epsilon = epsilonP;
  sigma = sigmaP;
  
  return 0;
}

int
FedeasMaterial::revertToStart(void)
{
  // Set all trial and committed values to zero
  for (int i = 0; i < 2*numHstv; i++)
    hstv[i] = 0.0;
  
  epsilonP = 0.0;
  sigmaP = 0.0;
  
  return 0;
}

// WARNING -- if you wish to override any method in this base class, you must
// also override the getCopy method to return a pointer to the derived class!!!
UniaxialMaterial*
FedeasMaterial::getCopy(void)
{
  FedeasMaterial *theCopy = 
    new FedeasMaterial(this->getTag(), this->getClassTag(), matType, numHstv, numData);
  
  // Copy history variables
  int i;
  for (i = 0; i < 2*numHstv; i++)
    theCopy->hstv[i] = hstv[i];
  
  for (i = 0; i < numData; i++)
    theCopy->data[i] = data[i];
  
  theCopy->epsilonP = epsilonP;
  theCopy->sigmaP = sigmaP;
  
  return theCopy;
}

int 
FedeasMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static ID idData(4);
  
  idData(0) = this->getTag();
  idData(1) = numHstv;
  idData(2) = numData;
  idData(3) = matType;
  
  res += theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) 
    cerr << "FedeasMaterial::sendSelf() - failed to send ID data\n";
  
  Vector vecData(numHstv+numData+2);
  
  int i, j;
  // Copy only the committed history variables into vector
  for (i = 0; i < numHstv; i++)
    vecData(i) = hstv[i];
  
  // Copy material properties into vector
  for (i = 0, j = numHstv; i < numData; i++, j++)
    vecData(j) = data[i];
  
  vecData(j++) = epsilonP;
  vecData(j++) = sigmaP;
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) 
    cerr << "FedeasMaterial::sendSelf() - failed to send Vector data\n";
  
  return res;
}

int
FedeasMaterial::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static ID idData(4);
  
  res += theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    cerr << "FedeasMaterial::recvSelf() - failed to receive ID data\n";
    return res;
  }
  
  this->setTag(idData(0));
  numHstv = idData(1);
  numData = idData(2);
  matType = idData(3);
  
  Vector vecData(numHstv+numData+2);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    cerr << "FedeasMaterial::recvSelf() - failed to receive Vector data\n";
    return res;
  }
  
  int i, j;
  // Copy committed history variables from vector
  for (i = 0; i < numHstv; i++)
    hstv[i] = vecData(i);
  
  // Copy material properties from vector
  for (i = 0, j = numHstv; i < numData; i++, j++)
    data[i] = vecData(j);
  
  epsilonP = vecData(j++);
  sigmaP   = vecData(j++);
  
  return res;
}

void
FedeasMaterial::Print(ostream &s, int flag)
{
  s << "FedeasMaterial, type: ";
	
  switch (matType) {
  case FEDEAS_Bond1:
    s << "Bond1" << endl;
    break;
  case FEDEAS_Bond2:
    s << "Bond2" << endl;
    break;
  case FEDEAS_Concrete1:
    s << "Concrete1" << endl;
    break;
  case FEDEAS_Concrete2:
    s << "Concrete2" << endl;
    break;
  case FEDEAS_Concrete3:
    s << "Concrete3" << endl;
    break;
  case FEDEAS_Hardening:
    s << "Hardening" << endl;
    break;
  case FEDEAS_Hysteretic1:
    s << "Hysteretic1" << endl;
    break;
  case FEDEAS_Hysteretic2:
    s << "Hysteretic2" << endl;
    break;
  case FEDEAS_Steel1:
    s << "Steel1" << endl;
    break;
  case FEDEAS_Steel2:
    s << "Steel2" << endl;
    break;
    // Add more cases as needed
    
  default:
    s << "Material identifier = " << matType << endl;
    break;
  }
}

#ifdef _WIN32

extern "C" int _stdcall BOND_1(double *matpar, double *hstvP, double *hstv,
			       double *strainP, double *stressP, double *dStrain,
			       double *tangent, double *stress, int *ist);

extern "C" int _stdcall BOND_2(double *matpar, double *hstvP, double *hstv,
			       double *strainP, double *stressP, double *dStrain,
			       double *tangent, double *stress, int *ist);

extern "C" int _stdcall CONCRETE_1(double *matpar, double *hstvP, double *hstv,
				   double *strainP, double *stressP, double *dStrain,
				   double *tangent, double *stress, int *ist);

extern "C" int _stdcall CONCRETE_2(double *matpar, double *hstvP, double *hstv,
				   double *strainP, double *stressP, double *dStrain,
				   double *tangent, double *stress, int *ist);

extern "C" int _stdcall CONCRETE_3(double *matpar, double *hstvP, double *hstv,
				   double *strainP, double *stressP, double *dStrain,
				   double *tangent, double *stress, int *ist);

extern "C" int _stdcall HARD_1(double *matpar, double *hstvP, double *hstv,
			       double *strainP, double *stressP, double *dStrain,
			       double *tangent, double *stress, int *ist);

extern "C" int _stdcall HYSTER_1(double *matpar, double *hstvP, double *hstv,
				 double *strainP, double *stressP, double *dStrain,
				 double *tangent, double *stress, int *ist);

extern "C" int _stdcall HYSTER_2(double *matpar, double *hstvP, double *hstv,
				 double *strainP, double *stressP, double *dStrain,
				 double *tangent, double *stress, int *ist);

extern "C" int _stdcall STEEL_1(double *matpar, double *hstvP, double *hstv,
				double *strainP, double *stressP, double *dStrain,
				double *tangent, double *stress, int *ist);

extern "C" int _stdcall STEEL_2(double *matpar, double *hstvP, double *hstv,
				double *strainP, double *stressP, double *dStrain,
				double *tangent, double *stress, int *ist);

// Add more declarations as needed

#define bond_1_			BOND_1
#define bond_2_			BOND_2
#define concrete_1_		CONCRETE_1
#define concrete_2_		CONCRETE_2
#define concrete_3_		CONCRETE_3
#define hard_1_			HARD_1
#define hyster_1_		HYSTER_1
#define hyster_2_		HYSTER_2
#define steel_1_		STEEL_1
#define steel_2_		STEEL_2

#else

extern "C" int bond_1_(double *matpar, double *hstvP, double *hstv,
		       double *strainP, double *stressP, double *dStrain,
		       double *tangent, double *stress, int *ist);

extern "C" int bond_2_(double *matpar, double *hstvP, double *hstv,
		       double *strainP, double *stressP, double *dStrain,
		       double *tangent, double *stress, int *ist);

extern "C" int concrete_1_(double *matpar, double *hstvP, double *hstv,
			   double *strainP, double *stressP, double *dStrain,
			   double *tangent, double *stress, int *ist);

extern "C" int concrete_2_(double *matpar, double *hstvP, double *hstv,
			   double *strainP, double *stressP, double *dStrain,
			   double *tangent, double *stress, int *ist);

extern "C" int concrete_3_(double *matpar, double *hstvP, double *hstv,
			   double *strainP, double *stressP, double *dStrain,
			   double *tangent, double *stress, int *ist);

extern "C" int hard_1_(double *matpar, double *hstvP, double *hstv,
		       double *strainP, double *stressP, double *dStrain,
		       double *tangent, double *stress, int *ist);

extern "C" int hyster_1_(double *matpar, double *hstvP, double *hstv,
			 double *strainP, double *stressP, double *dStrain,
			 double *tangent, double *stress, int *ist);

extern "C" int hyster_2_(double *matpar, double *hstvP, double *hstv,
			 double *strainP, double *stressP, double *dStrain,
			 double *tangent, double *stress, int *ist);

extern "C" int steel_1_(double *matpar, double *hstvP, double *hstv,
			double *strainP, double *stressP, double *dStrain,
			double *tangent, double *stress, int *ist);

extern "C" int steel_2_(double *matpar, double *hstvP, double *hstv,
			double *strainP, double *stressP, double *dStrain,
			double *tangent, double *stress, int *ist);

// Add more declarations as needed

#endif

int
FedeasMaterial::invokeSubroutine(int ist)
{
  // Compute strain increment
  double dEpsilon = epsilon-epsilonP;
  
  switch (matType) {
  case FEDEAS_Bond1:
    //bond_1_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon,
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Bond2:
    //bond_2_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon,
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Concrete1:
    //concrete_1_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Concrete2:
    //concrete_2_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Concrete3:
    //concrete_3_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Hardening:
    //hard_1_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Hysteretic1:
    //hyster_1_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Hysteretic2:
    //hyster_2_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Steel1:
    //steel_1_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
    //	&sigma, &tangent, &ist);
    break;
    
  case FEDEAS_Steel2:
    //steel_2_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
    //	&sigma, &tangent, &ist);
    break;
    
    // Add more cases as needed
  default:
    g3ErrorHandler->fatal("%s -- unknown material type",
			  "FedeasMaterial::invokeSubroutine");
    return -1;
  }
  
  return 0;
}

