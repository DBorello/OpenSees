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
// $Date: 2002-01-19 16:47:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/MinMaxMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// MinMaxMaterial.  MinMaxMaterial wraps a UniaxialMaterial
// and imposes min and max strain limits.

#include <MinMaxMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <G3Globals.h>

MinMaxMaterial::MinMaxMaterial(int tag, UniaxialMaterial &material,
							   double min, double max)
:UniaxialMaterial(tag,MAT_TAG_MinMax), theMaterial(0),
minStrain(min), maxStrain(max), Tfailed(false), Cfailed(false)
{
	theMaterial = material.getCopy();

	if (theMaterial == 0)
		g3ErrorHandler->fatal("%s -- failed to get copy of material",
			"MinMaxMaterial::MinMaxMaterial");
}

MinMaxMaterial::MinMaxMaterial()
:UniaxialMaterial(0,MAT_TAG_MinMax), theMaterial(0),
minStrain(0.0), maxStrain(0.0), Tfailed(false), Cfailed(false)
{

}

MinMaxMaterial::~MinMaxMaterial()
{
	if (theMaterial)
		delete theMaterial;
}

int 
MinMaxMaterial::setTrialStrain(double strain, double strainRate)
{
	if (Cfailed)
		return 0;

	if (strain >= maxStrain || strain <= minStrain) {
		Tfailed = true;
		return 0;
	}
	else {
		Tfailed = false;
		return theMaterial->setTrialStrain(strain, strainRate);
	}
}

double 
MinMaxMaterial::getStress(void)
{
	if (Tfailed)
		return 0.0;
	else
		return theMaterial->getStress();
}

double 
MinMaxMaterial::getTangent(void)
{
	if (Tfailed)
		return 0.0;
	else
		return theMaterial->getTangent();
}

double 
MinMaxMaterial::getDampTangent(void)
{
   	if (Tfailed)
		return 0.0;
	else
		return theMaterial->getDampTangent();
}

double 
MinMaxMaterial::getSecant(void)
{
	if (Tfailed)
		return 0.0;
	else
		return theMaterial->getSecant();
}

double 
MinMaxMaterial::getStrain(void)
{
    return theMaterial->getStrain();
}

double 
MinMaxMaterial::getStrainRate(void)
{
    return theMaterial->getStrainRate();
}

int 
MinMaxMaterial::commitState(void)
{	
	Cfailed = Tfailed;

	// Check if failed at current step
	if (Tfailed)
		return 0;
	else
		return theMaterial->commitState();
}

int 
MinMaxMaterial::revertToLastCommit(void)
{
	// Check if failed at last step
	if (Cfailed)
		return 0;
	else
		return theMaterial->revertToLastCommit();
}

int 
MinMaxMaterial::revertToStart(void)
{
	Cfailed = false;
	Tfailed = false;

	return theMaterial->revertToStart();
}

UniaxialMaterial *
MinMaxMaterial::getCopy(void)
{
    MinMaxMaterial *theCopy = 
		new MinMaxMaterial(this->getTag(), *theMaterial, minStrain, maxStrain);
        
	theCopy->Cfailed = Cfailed;
	theCopy->Tfailed = Tfailed;

	return theCopy;
}

int 
MinMaxMaterial::sendSelf(int cTag, Channel &theChannel)
{
	return -1;
}

int 
MinMaxMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
    return -1;
}

void 
MinMaxMaterial::Print(ostream &s, int flag)
{
    s << "MinMaxMaterial tag: " << this->getTag() << endl;
    s << "\tMaterial: " << theMaterial->getTag() << endl;
	s << "\tMin strain: " << minStrain << endl;
	s << "\tMax strain: " << maxStrain << endl;
}
