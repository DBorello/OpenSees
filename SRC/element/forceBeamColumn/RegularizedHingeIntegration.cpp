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
// $Date: 2008-04-29 17:48:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/RegularizedHingeIntegration.cpp,v $

// Theory Reference
// ----------------
// Scott, M.H. and Hamutcuoglu, O.M. "Numerically consistent regularization
// of force-based frame elements." International Journal for Numerical
// Methods in Engineering. In press, Approved April 2008.

#include <RegularizedHingeIntegration.h>
#include <ElementalLoad.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

RegularizedHingeIntegration::RegularizedHingeIntegration(BeamIntegration &bi,
							 double lpi, double lpj,
							 double epsi, double epsj):
  BeamIntegration(BEAM_INTEGRATION_TAG_RegularizedHinge),
  beamInt(0), lpI(lpi), lpJ(lpj), epsI(epsi), epsJ(epsj)
{
  beamInt = bi.getCopy();
  if (beamInt == 0) {
    opserr << "RegularizedHingeIntegration::RegularizedHingeIntegration -- failed to get copy of BeamIntegration" << endln;
  }
}

RegularizedHingeIntegration::RegularizedHingeIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_RegularizedHinge),
  beamInt(0), lpI(0.0), lpJ(0.0), epsI(0.0), epsJ(0.0)
{

}

RegularizedHingeIntegration::~RegularizedHingeIntegration()
{
  if (beamInt != 0)
    delete beamInt;
}

void
RegularizedHingeIntegration::getSectionLocations(int numSections, double L,
						 double *xi)
{
  beamInt->getSectionLocations(numSections-2, L, xi);

  double oneOverL = 1.0/L;

  xi[numSections-2] = epsI*oneOverL;
  xi[numSections-1] = 1-epsJ*oneOverL;
}

void
RegularizedHingeIntegration::getSectionWeights(int numSections, double L,
					       double *wt)
{
  beamInt->getSectionWeights(numSections-2, L, wt);

  double oneOverL = 1.0/L;

  double betaI = lpI*oneOverL;
  wt[numSections-2] = wt[0]-betaI;
  wt[0] = betaI;

  double betaJ = lpJ*oneOverL;
  wt[numSections-1] = wt[numSections-3]-betaJ;
  wt[numSections-3] = betaJ;
}

BeamIntegration*
RegularizedHingeIntegration::getCopy(void)
{
  return new RegularizedHingeIntegration(*beamInt, lpI, lpJ, epsI, epsJ);
}

int
RegularizedHingeIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(4);

  data(0) = lpI;
  data(1) = lpJ;
  data(2) = epsI;
  data(3) = epsJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "RegularizedHingeIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
RegularizedHingeIntegration::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  static Vector data(4);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "RegularizedHingeIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  lpI = data(0);
  lpJ = data(1);
  epsI = data(2);
  epsJ = data(3);

  return 0;
}

void
RegularizedHingeIntegration::Print(OPS_Stream &s, int flag)
{
  s << "RegularizedHinge" << endln;
  s << " lpI = " << lpI;
  s << " lpJ = " << lpJ << endln;
  s << " epsI = " << epsI;
  s << " epsJ = " << epsJ << endln;

  beamInt->Print(s, flag);

  return;
}
