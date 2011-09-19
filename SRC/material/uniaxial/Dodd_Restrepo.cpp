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
                                                                        
// $Revision: 1.20 $
// $Date: 2006/01/03 23:52:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Dodd_Restrepo.cpp,v $
                                                                        
// Written: fmk
// Created: Sep 2011
//
// Description: This file contains the class definition for
// Dodd_Restrepo. Dodd_Restrepo class calls the fortran routine


#include <elementAPI.h>

#include <OPS_Globals.h>
#include "Dodd_Restrepo.h"
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

 
int static numDoddRestrepo = 0;

void *
OPS_Dodd_Restrepo(void)
{
  if (numDoddRestrepo == 0) {
    numDoddRestrepo++;
    opserr << "Dodd_Restrepo unaxial material - Written by L.L. Dodd & J. Restepo\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[9];
  int numData;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 8 || numArgs > 10) {
    opserr << "WARNING wrong # args: uniaxialMaterial $tag $Fy $Fsu $ESH $ESU $Youngs $ESHI $FSHI <$OmegaFac $Conv>" << endln;
    return 0;
  }
    
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticPP tag" << endln;
    return 0;
  }

  numData = numArgs-1;
  dData[7] = 1.0; // OmegFac = 1.0
  dData[8] = 1.0; // Conv = 1.0

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E & ep\n";
    return 0;	
  }

  theMaterial = new Dodd_Restrepo(iData[0], dData[0], dData[1], dData[2],
				  dData[3], dData[4], dData[5], dData[6],
				  dData[7], dData[8]);       
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ElasticPPCpp\n";
    return 0;
  }
  
  return theMaterial;
}


#ifdef _WIN32

//#define steel STEEL

extern "C" int STEEL(double *Es, double *EpsLast, double *FpsLast, double *YpTanLast, 
		       double *EpsOld, double *Fy, double *Epy, double * EpSH, double *Epsu, 
		       double *Fpsu, double *Youngs, double *SHPower, double *Epr, double *Fpr, 
		       double *Epa, double *Fpa, double *Epo, double *EpoMax, double *EpsuSh, 
		       double *YoungsUn, double *Power, int *BFlag, int *LMR, double *EprM, 
		       double *FprM, double *EpaM, double *FpaM, double *YpTanM, double *PowerM, 
		       double * Eps, double *Fps, double *Fs, double *YpTan, double *YTan, double *OmegFac);

// Add more declarations as needed

#else

#define STEEL steel_

extern "C" int STEEL(double *Es, double *EpsLast, double *FpsLast, double *YpTanLast, 
			double *EpsOld, double *Fy, double *Epy, double * EpSH, double *Epsu, 
			double *Fpsu, double *Youngs, double *SHPower, double *Epr, double *Fpr, 
			double *Epa, double *Fpa, double *Epo, double *EpoMax, double *EpsuSh, 
			double *YoungsUn, double *Power, int *BFlag, int *LMR, double *EprM, 
			double *FprM, double *EpaM, double *FpaM, double *YpTanM, double *PowerM, 
			double * Eps, double *Fps, double *Fs, double *YpTan, double *YTan, double *OmegFac);

#endif



Dodd_Restrepo::Dodd_Restrepo(int tag, 
			     double fy, 
			     double fsu, 
			     double esh, 
			     double eSU,
			     double youngs, 
			     double eshi, 
			     double fSHI, 
			     double omegaFac,
			     double conv)

  :UniaxialMaterial(tag,0),
   Fy(fy), Fsu(fsu), ESH(esh), ESU(eSU), Youngs(youngs), 
   ESHI(eshi), FSHI(fSHI), Conv(conv), OmegaFac(omegaFac)
{
  if (OmegaFac < 0.6 || OmegaFac > 1.3) 
    OmegaFac = 1.0;
  
  Epy       = Fy/Youngs;
  EpSH      = log(1+ESH/Conv);
  Epsu      = log(1+ESU/Conv);
  Fpsu      = Fsu*(1+ESU/Conv);
  EpsuSh[0] =  Epsu;
  EpsuSh[1] = -Epsu;
  YoungsUn  = Youngs;

  LMR       = 0;
  BFlag[0]  = 0;
  BFlag[1]  = 0;

  Epa[0]    = 0.0;
  Epa[1]    = 0.0;
  EpaM[0]   = 0.0;
  EpaM[1]   = 0.0;
  Epo[0]   = 0.0;
  Epo[1]    = 0.0;
  EpoMax    = 0.0;
  Epr[0]    = 0.0;
  Epr[1]    = 0.0;
  EprM[0]   = 0.0;
  EprM[1]   = 0.0;

  //
  // Calculate the power term for the strain hardening branch
  //

  double C1       ; // Temporary constant
  double EpSHI    ; // Intermediate strain hardening curve natural strain
  double FpSH     ; // True stress at initiation of strain hardening curve
  double FpSHI    ; // Intermediate strain hardening curve true stress

  EpSHI   = log(1+ESHI/Conv);
  FpSH    = Fy  *(1+ESH/Conv);
  FpSHI   = FSHI*(1+ESHI/Conv);
  C1      = FpSH-Fpsu+Fpsu*(Epsu-EpSH);
  SHPower = log((FpSHI+Fpsu*(Epsu-EpSHI)-Fpsu)/C1) / log((Epsu-EpSHI)/(Epsu-EpSH));

  tStrain = 0.0;
  tTangent = Youngs;
  tStress = 0.0;

  STEEL(&tStrain, &EpsLast, &FpsLast, &YpTanLast, 
	  &EpsOld, &Fy, &Epy,  &EpSH, &Epsu, 
	  &Fpsu, &Youngs, &SHPower, Epr, Fpr, 
	  Epa, Fpa, Epo, &EpoMax, EpsuSh,&YoungsUn, Power, BFlag, &LMR, EprM,
	  FprM, EpaM, FpaM, YpTanM, PowerM,
	  &Eps, &Fps, &Fs, &YpTan, &YTan, &OmegaFac);
  
  this->commitState();
}

UniaxialMaterial *
Dodd_Restrepo::getCopy(void)
{
  return new Dodd_Restrepo(this->getTag(), Fy, Fsu, ESH, ESU, Youngs, ESHI, FSHI, OmegaFac, Conv);
}

int
Dodd_Restrepo::setTrialStrain(double strain, double strainRate)
{

  if (fabs(strain-tStrain) > DBL_EPSILON) {

    tStrain = strain;

    STEEL(&tStrain, &EpsLast, &FpsLast, &YpTanLast, 
	    &EpsOld, &Fy, &Epy,  &EpSH, &Epsu, 
	    &Fpsu, &Youngs, &SHPower, Epr, Fpr, 
	    Epa, Fpa, Epo, &EpoMax, EpsuSh,&YoungsUn, Power, BFlag, &LMR, EprM,
	    FprM, EpaM, FpaM, YpTanM, PowerM,
	    &Eps, &Fps, &Fs, &YpTan, &YTan, &OmegaFac);

    tStress = Fs;
    tTangent = YTan;
  }

  return 0;
}

int
Dodd_Restrepo::setTrial(double strain, double &stress, double &stiff, double strainRate)
{
  opserr << "Dodd_Restrepo::setTrialStrain()2: strain: " << strain << endln;
  
  if (fabs(strain-tStrain) > DBL_EPSILON) {
    // Store the strain
    tStrain = strain;

    return STEEL(&tStrain, &EpsLast, &FpsLast, &YpTanLast, 
		      &EpsOld, &Fy, &Epy,  &EpSH, &Epsu, 
		      &Fpsu, &Youngs, &SHPower, Epr, Fpr, 
		      Epa, Fpa, Epo, &EpoMax, EpsuSh, 
		      &YoungsUn, Power, BFlag, &LMR, EprM, 
		      FprM, EpaM, FpaM, YpTanM, PowerM, 
		      &Eps, &Fps, &Fs, &YpTan, &YTan, &OmegaFac);


  }
  return 0;
  tStress = Fs;
  tTangent = YTan;
    
  stress = tStress;;
  stiff = tTangent;
  
  return 0;
}

double
Dodd_Restrepo::getStrain(void)
{
  return tStrain;
}

double
Dodd_Restrepo::getStress(void)
{
  return tStress;
}

double
Dodd_Restrepo::getTangent(void)
{
  return tTangent;
}

double
Dodd_Restrepo::getInitialTangent(void)
{
  return Youngs;
}

int
Dodd_Restrepo::commitState(void)
{
  if (cStrain != tStrain) {
    EpsOld    = EpsLast;
    EpsLast   = Eps;
    FpsLast   = Fps;
    YpTanLast = YpTan;
  }
  cStrain = tStrain;
  cStress = tStress;
  cTangent = tTangent;
  
  return 0;
}

int
Dodd_Restrepo::revertToLastCommit(void)
{
  tStrain = cStrain;
  tStress = cStress;
  tTangent = cTangent;

  return 0;
}

int
Dodd_Restrepo::revertToStart(void)
{
  double C1       ; // Temporary constant
  double EpSHI    ; // Intermediate strain hardening curve natural strain
  double ESH      ; // Engineering coordinate strain hardening strain
  double ESHI     ; // Intermediate strain hardening curve engineering strain
  double FpSH     ; // True stress at initiation of strain hardening curve
  double FpSHI    ; // Intermediate strain hardening curve true stress

  Epy       = Fy/Youngs;
  EpSH      = log(1+ESH/Conv);
  Epsu      = log(1+ESU/Conv);
  Fpsu      = Fsu*(1+ESU/Conv);
  EpsuSh[0] =  Epsu;
  EpsuSh[1] = -Epsu;
  YoungsUn  = Youngs;

  LMR       = 0;
  BFlag[0]  = 0;
  BFlag[1]  = 0;

  Epa[0]    = 0.0;
  Epa[1]    = 0.0;
  EpaM[0]   = 0.0;
  EpaM[1]   = 0.0;
  Epo[0]   = 0.0;
  Epo[1]    = 0.0;
  EpoMax    = 0.0;
  Epr[0]    = 0.0;
  Epr[1]    = 0.0;
  EprM[0]   = 0.0;
  EprM[1]   = 0.0;

  //
  // Calculate the power term for the strain hardening branch
  //

  EpSHI   = log(1+ESHI/Conv);
  FpSH    = Fy  *(1+ESH/Conv);
  FpSHI   = FSHI*(1+ESHI/Conv);
  C1      = FpSH-Fpsu+Fpsu*(Epsu-EpSH);
  SHPower = log((FpSHI+Fpsu*(Epsu-EpSHI)-Fpsu)/C1) / log((Epsu-EpSHI)/(Epsu-EpSH));

  tStrain = 0.0;
  tTangent = Youngs;
  tStress = 0.0;
    
  this->commitState();

  return 0;
}

int 
Dodd_Restrepo::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  return res;
}

int
Dodd_Restrepo::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  return res;
}

void
Dodd_Restrepo::Print(OPS_Stream &s, int flag)
{
  s << "Dodd_Restrepo: " << this->getTag() << endln;
  s << "tStrain: " << tStrain << "  tStress: " << tStress << " tTangent: " << tTangent << endln;
}






