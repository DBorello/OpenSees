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
// $Date: 2003-08-12 23:01:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FractureMaterial.cpp,v $
                                                      
// Written: Patxi

// Description: This file contains the class definition for 
// FractureMaterial.  FractureMaterial wraps a UniaxialMaterial
// and imposes fatigue limits.

#include <stdlib.h>

#include <FractureMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

FractureMaterial::FractureMaterial(int tag, UniaxialMaterial &material,
				   double dmax, double nf, double e0, double fe)
  :UniaxialMaterial(tag,MAT_TAG_Fracture), theMaterial(0), strain(0.0),
   Tfailed(false), Cfailed(false)
{
  D   = 0; // Damage index
  X   = 0; // Range in consideration
  Y   = 0; // Previous Adjacent Range
  A   = 0; // Strain at first  cycle peak/valley
  B   = 0; // Strain at second cycle peak/valley
  R1F = 0; // Flag for first  cycle count
  R2F = 0; // Flag for second cycle count
  CS  = 0; // Current Slope
  PS  = 0; // Previous slope
  EP  = 0; // Previous Strain
  FF  = 0; // Failure Flag
  SF  = 0; // Start Flag - for initializing the very first strain
  PF  = 0; // Peak Flag --> Did we reach a peak/valley at current strain?

  Dmax = dmax;
  Nf  = nf;
  E0  = e0;
  FE  = fe;
  b   =  log10(E0) - log10(Nf)*FE; //Theoretical Intercept

  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "FractureMaterial::FractureMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

FractureMaterial::FractureMaterial()
  :UniaxialMaterial(0,MAT_TAG_Fracture), theMaterial(0), strain(0.0),
   Tfailed(false), Cfailed(false)
{
  D   = 0; // Damage index
  X   = 0; // Range in consideration
  Y   = 0; // Previous Adjacent Range
  A   = 0; // Strain at first  cycle peak/valley
  B   = 0; // Strain at second cycle peak/valley
  R1F = 0; // Flag for first  cycle count
  R2F = 0; // Flag for second cycle count
  CS  = 0; // Current Slope
  PS  = 0; // Previous slope
  EP  = 0; // Previous Strain
  FF  = 0; // Failure Flag
  SF  = 0; // Start Flag - for initializing the very first strain
  PF  = 0; // Peak Flag --> Did we reach a peak/valley at current strain?

  Dmax = 0;
  Nf   = 0;
  E0   = 0;
  FE   = 0;
  b    =  0;
}

FractureMaterial::~FractureMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

static int sign(double a) {
  if (a < 0)
    return -1;
  else
    return 1;
}

int 
FractureMaterial::setTrialStrain(double eps, double strainRate)
{
  if (Cfailed)
    return 0;

  strain = eps;
  return theMaterial->setTrialStrain(strain, strainRate);
}

double 
FractureMaterial::getStress(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getStress();
}

double 
FractureMaterial::getTangent(void)
{
  if (Tfailed)
    //return 0.0;
    return 1.0e-8*theMaterial->getInitialTangent();
  else
    return theMaterial->getTangent();
}

double 
FractureMaterial::getDampTangent(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}



double 
FractureMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
FractureMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
FractureMaterial::commitState(void)
{	
  if (Cfailed == false) {
 
    //Now check to see if we have reached failure although  
    // not at a peak or valley, assume
    // a 1/2 cycle to the current point
    if (SF == 1) {
      double Ytemp = fabs(strain - B);
      double Dtemp = D +  0.5 / pow(10.0, ((log10(Ytemp)-b)/FE ) );
      if (Dtemp > Dmax) {
	D = Dtemp;
	Tfailed = true;
      }
    }
    
    if (Tfailed == false) {
      if (SF == 0) {
	PS = strain;
	A  = strain;
	SF = 1;
      }
      
      CS = strain - EP;         // Determine Current Slope
      
      // Identify Peak or Valley (slope Change)
      if ((sign(PS) != sign(CS)) && (FF == 0)) {
	
	X =  fabs(strain - B);   // Range of current Cycle
	
	if (R1F == 0) {            // After Reaching First Peak
	  Y = fabs(strain - A);
	  B = strain;
	  R1F = 1;
	}
	
	if ((R2F == 0) && (R1F == 1)) { // Mark the first cycle  
	  A = B;
	  B = strain;
	  R2F = 1;
	  D = 0.5 / pow(10.0, ((log10(Y)-b)/FE ) );
	  PF = 1; // We reached a peak
	  Y = X;
	} else if (X < Y) {
	  // Do Nothing; keep counting, not a peak or 
	  // valley by definition
	  ;
	}
	else {
	  A = B;
	  B = strain;
	  D = D + 0.5 / pow(10.0, ((log10(Y)-b)/FE ));
	  PF = 1; // We reached a peak
	  //plot(r(i), B, 'ko','MarkerSize',20)
	  Y = X;
	}
      }
      
      PS = CS;       // Previous Slope
      EP = strain;   // Keep track of previous strain
      
      if (D >= Dmax) {
	Tfailed = true;
      }
    }
  } else {
    theMaterial->commitState();
  }

  Cfailed = Tfailed;
    
  return  0;
}

int 
FractureMaterial::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int 
FractureMaterial::revertToStart(void)
{
  Cfailed = false;
  Tfailed = false;
  
  return theMaterial->revertToStart();
}

UniaxialMaterial *
FractureMaterial::getCopy(void)
{
  FractureMaterial *theCopy = 
    new FractureMaterial(this->getTag(), *theMaterial, Dmax, Nf, E0, FE);
  
  theCopy->Cfailed = Cfailed;
  theCopy->Tfailed = Tfailed;
  
  return theCopy;
}

int 
FractureMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
FractureMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
FractureMaterial::Print(OPS_Stream &s, int flag)
{
  s << "FractureMaterial tag: " << this->getTag() << endln;
  s << "\tMaterial: " << theMaterial->getTag() << endln;
  s << "\tD: " << D << " Dmax: " << Dmax << endln;
  s << "Nf: " << Nf <<  " E0: " << E0 << " FE: " << FE << " b: " << b << endln;
}
