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
                                                                        
// $Revision: 1.4 $
// $Date: 2002-03-05 16:48:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/CableMaterial.cpp,v $
                                                                        
// Written: Charles Chadwell 
// Created: 07/01
//
// Description: This file contains the class definition for 
// CableMaterial. CableMaterial provides the abstraction
// of an elastic uniaxial material,
//
// The input parameters are the Prestress, E, Effective Self Weight (gravity component of 
// Weight per volume transverse to the cable), and Length of the cable.
//
// The cable Force Displacement is converted to Stress Strain material for use 
// with the truss element.  The stress strain ranges from slack (large strain at zero 
// stress) to taught (linear with modulus E).  The material has no history and is not
// path dependent.
//
//
// What: "@(#) CableMaterial.cpp, revA"

#include <CableMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>


CableMaterial::CableMaterial(int tag, double PRESTRESS, double e, double UNIT_WEIGHT_EFF, double L_Element)
:UniaxialMaterial(tag,MAT_TAG_CableMaterial),
 Ps(PRESTRESS), E(e), Mue(UNIT_WEIGHT_EFF), 
 L(L_Element), trialStrain(0.0), commitStrain(0.0)
{

}

CableMaterial::CableMaterial()
:UniaxialMaterial(0,MAT_TAG_CableMaterial),
 Ps(0.0), E(0.0), Mue(0.0), 
 L(0.0), trialStrain(0.0), commitStrain(0.0)
{

}

CableMaterial::~CableMaterial()
{
  // does nothing
}

int 
CableMaterial::setTrialStrain(double strain, double strainRate)
{
    trialStrain     = strain;
    return 0;
}

double 
CableMaterial::getStress(void)
{
    
	// Check if out side the range of inportance 
	double testStress, dP, curstrain, e0;
	int i = 0;
    
	// Perameters for bisection
	double L_bound = 0, U_bound, middle = 0;
	
	
    if (trialStrain < 0) U_bound = Ps;
	else
	{
		U_bound = E * trialStrain + Ps;
        testStress = U_bound;
    }
    
	// Check if slack in cable has been taken out and it is a bar
    e0 = Mue*Mue*L*L/(24*Ps*Ps) - Ps/E;
	if (trialStrain > 0 && abs(trialStrain - evalStress((trialStrain - e0)*E)) < 10e-9) return (trialStrain - e0)*E;

	// Check if all slack
	if (trialStrain < - Ps/E*10.0) return 0.0; 

	// if stress is in between then do itterations -- Bisection
    dP = U_bound - L_bound;
    
	while (abs(dP)/U_bound > 0.00000001 && i < 100)
	{
		middle = .5 * (U_bound + L_bound);
        curstrain = evalStress(middle);
        
		if (curstrain <= trialStrain) 
		{
			L_bound = middle;
		}
		else 
		{
			U_bound = middle;
		}
        dP = U_bound - L_bound;
	    i++;
	}
	
	
	// if it did not converge - return near zero stress
    if (i == 100) 
	{
		return 0.0;
    }
	else return middle;
}

double 
CableMaterial::evalStress(double stress)
{
    double strainG, strainE;
    
	// Should never be zero or less than zero
	if (stress <= 0) {return -10;}
		
    // Elastic Part
	strainE = 1 / E * (stress - Ps) * (1 + Mue * Mue * L * L / (24 * stress));
    // Geometric Part
	strainG = 1 / 24 * Mue * Mue * L * L * (1 / (Ps * Ps) - 1 / (stress * stress));
    return strainE + strainG;
}


double 
CableMaterial::getTangent(void) 
{
	
	double derivE, derivG, stress;

	stress = getStress();

	if (stress <= 0.0) {return 0;}
	
	// Elastic Part
	derivE = 1 / E * (1. - Mue * Mue * L * L / (24. * stress * stress) * (1. - 2. * Ps / stress));
    // Geometric Part
	derivG = 1 / 12. * Mue * Mue * L * L / (stress * stress * stress);
    
	if (derivE + derivG != 0.0)
		return 1.0 / (derivE + derivG);
	else 
		return 1e-8;
};
 
int 
CableMaterial::commitState(void)
{
    commitStrain     = trialStrain;
    return 0;
}

int 
CableMaterial::revertToLastCommit(void)
{
    trialStrain     = commitStrain;
    return 0;
}

int 
CableMaterial::revertToStart(void)
{
    commitStrain = 0.0;
	trialStrain = 0.0;
    return 0;
}

UniaxialMaterial *
CableMaterial::getCopy(void)
{
    CableMaterial *theCopy = new CableMaterial(this->getTag(), Ps, E, Mue, L);
    theCopy->trialStrain = trialStrain;
    return theCopy;
}

int 
CableMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = Ps;
  data(2) = E;
  data(3) = Mue;
  data(4) = L;
  data(5) = commitStrain;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    cerr << "CableMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
CableMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      cerr << "CableMaterial::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag(data(0));
    Ps   = data(1);
    E = data(2);
    Mue     = data(3);
    L = data(4);
    commitStrain = data(5);
	trialStrain      = commitStrain;
;
  }
    
  return res;
}

void 
CableMaterial::Print(ostream &s, int flag)
{
   s << "CableMaterial tag: " << this->getTag() << endl;
    s << "  E: " << E << " Prestress: " << Ps << endl;
}

//int
//CableMaterial::setParameter(char **argv, int argc, Information &info)
//{
//	if (strcmp(argv[0],"E") == 0) {
//		info.theType = DoubleType;
///		return 1;
//	}
//	else if (strcmp(argv[0],"eta") == 0) {
//		info.theType = DoubleType;
//		return 2;
//	}
//	else
//		return -1;
//}
//
//int 
//CableMaterial::updateParameter(int parameterID, Information &info)
//{
//	switch(parameterID) {
//	case -1:
//		return -1;
//	case 1:
//		E = info.theDouble;
//		return 0;
//	case 2:
//		eta = info.theDouble;
//		return 0;
//	default:
//		return -1;
//	}
//}

 
double 
CableMaterial::abs(double value)
{
	if (value < 0) return -value;
	else return value;
}
