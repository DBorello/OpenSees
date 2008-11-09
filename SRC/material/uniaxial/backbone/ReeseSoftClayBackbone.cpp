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
// $Date: 2008-11-09 06:07:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/ReeseSoftClayBackbone.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// ReeseSoftClayBackbone

#include <ReeseSoftClayBackbone.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>

ReeseSoftClayBackbone::ReeseSoftClayBackbone(int tag, double p, double y):
  HystereticBackbone(tag,BACKBONE_TAG_ReeseSoftClay),
  pu(p), y50(y)
{

}

ReeseSoftClayBackbone::ReeseSoftClayBackbone():
  HystereticBackbone(0,BACKBONE_TAG_ReeseSoftClay),
  pu(0.0), y50(0.0)
{
  
}

ReeseSoftClayBackbone::~ReeseSoftClayBackbone()
{
  
}

double
ReeseSoftClayBackbone::getTangent (double strain)
{
  int signStrain = (strain > 0.0) ? 1 : -1;
  strain = signStrain*strain;

  double minStrain = 0.001*y50;
  double tangent;
  if (strain > 8*y50)
    tangent = 0.001*pu/y50;
  else if (strain > minStrain)
    tangent = pu/(6*y50)*pow(y50/strain,0.66667);
  else
    tangent = pu*0.5*pow(0.001,0.33333)/minStrain;

  return tangent;
}

double
ReeseSoftClayBackbone::getStress (double strain)
{
  int signStrain = (strain > 0.0) ? 1 : -1;
  strain = signStrain*strain;

  double minStrain = 0.001*y50;
  double stress;
  if (strain > 8*y50)
    stress = pu;
  else if (strain > minStrain)
    stress = pu*0.5*pow(strain/y50,0.33333);
  else
    stress = pu*0.5*pow(0.001,0.33333)/minStrain*strain;

  return signStrain*stress;
}

double
ReeseSoftClayBackbone::getEnergy (double strain)
{
  return 0.0;
}

double
ReeseSoftClayBackbone::getYieldStrain(void)
{
  return 0.0;
}

HystereticBackbone*
ReeseSoftClayBackbone::getCopy(void)
{
  ReeseSoftClayBackbone *theCopy =
    new ReeseSoftClayBackbone (this->getTag(), pu, y50);
  
  return theCopy;
}

void
ReeseSoftClayBackbone::Print(OPS_Stream &s, int flag)
{
  s << "ReeseSoftClayBackbone, tag: " << this->getTag() << endln;
  s << "\tpu: " << pu << endln;
  s << "\ty50: " << y50 << endln;
}

int
ReeseSoftClayBackbone::setVariable (char *argv)
{
  return -1;
}

int
ReeseSoftClayBackbone::getVariable (int varID, double &theValue)
{
  return -1;
}

int
ReeseSoftClayBackbone::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(3);
  
  data(0) = this->getTag();
  data(1) = pu;
  data(2) = y50;
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ReeseSoftClayBackbone::sendSelf -- could not send Vector" << endln;

    return res;
  }
  
  return res;
}

int
ReeseSoftClayBackbone::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(3);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ReeseSoftClayBackbone::recvSelf -- could not receive Vector" << endln;

    return res;
  }
  
  this->setTag(int(data(0)));
  pu = data(1);
  y50 = data(2);
  
  return res;
}
