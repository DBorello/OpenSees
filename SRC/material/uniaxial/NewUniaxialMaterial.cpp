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
// $Date: 2001-08-20 03:45:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NewUniaxialMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class implementation for 
// NewUniaxialMaterial.

#include <NewUniaxialMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

NewUniaxialMaterial::NewUniaxialMaterial(int tag)
  :UniaxialMaterial(tag,MAT_TAG_NewUniaxial),
   Tstrain(0.0), Tstress(0.0), Ttangent(0.0)
{

}

NewUniaxialMaterial::NewUniaxialMaterial()
  :UniaxialMaterial(0,MAT_TAG_NewUniaxial),
   Tstrain(0.0), Tstress(0.0), Ttangent(0.0)
{

}

NewUniaxialMaterial::~NewUniaxialMaterial()
{

}

int 
NewUniaxialMaterial::setTrialStrain(double strain, double strainRate)
{
  return 0;
}

double 
NewUniaxialMaterial::getStress(void)
{
  return Tstress;
}

double 
NewUniaxialMaterial::getTangent(void)
{
  return Ttangent;
}

double 
NewUniaxialMaterial::getStrain(void)
{
  return Tstrain;
}

int 
NewUniaxialMaterial::commitState(void)
{
  return 0;
}

int 
NewUniaxialMaterial::revertToLastCommit(void)
{
  return 0;
}

int 
NewUniaxialMaterial::revertToStart(void)
{
  return 0;
}

UniaxialMaterial *
NewUniaxialMaterial::getCopy(void)
{
  NewUniaxialMaterial *theCopy = new NewUniaxialMaterial(this->getTag());

  return theCopy;
}

int 
NewUniaxialMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
NewUniaxialMaterial::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
NewUniaxialMaterial::Print(ostream &s, int flag)
{
  return;
}


