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

// $Revision: 1.3 $
// $Date: 2003-02-21 22:27:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/HingeRadauBeamIntegration2d.cpp,v $

#include <HingeRadauBeamIntegration2d.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

HingeRadauBeamIntegration2d::HingeRadauBeamIntegration2d(double e,
							 double a,
							 double i,
							 double lpi,
							 double lpj):
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadau2d),
  E(e), A(a), I(i), lpI(lpi), lpJ(lpj)
{
  // Nothing to do
}

HingeRadauBeamIntegration2d::HingeRadauBeamIntegration2d():
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadau2d),
  E(0.0), A(0.0), I(0.0), lpI(0.0), lpJ(0.0)
{

}

HingeRadauBeamIntegration2d::~HingeRadauBeamIntegration2d()
{
  // Nothing to do
}

void
HingeRadauBeamIntegration2d::getSectionLocations(int numSections, double L,
						 double *xi)
{
  xi[0] = 0.0;
  xi[1] = 1.0;
  for (int i = 2; i < numSections; i++)
    xi[i] = 0.0;
}

void
HingeRadauBeamIntegration2d::getSectionWeights(int numSections, double L,
					       double *wt)
{
  double oneOverL = 1.0/L;

  wt[0] = lpI*oneOverL;
  wt[1] = lpJ*oneOverL;
  for (int i = 2; i < numSections; i++)
    wt[i] = 1.0;
}

int
HingeRadauBeamIntegration2d::addElasticFlexibility(double L, Matrix &fElastic)
{
  double oneOverL = 1.0/L;

  double Le = L-lpI-lpJ;
  fElastic(0,0) += Le/(E*A);

  double x[4];
  double w[4];
  
  static const double eight3 = 8.0/3.0;
  static const double oneOverRoot3 = 1.0/sqrt(3.0);

  double oneOverEI = 1.0/(E*I);
  
  x[0] = eight3*lpI;
  w[0] = 3.0*lpI;

  x[1] = L-eight3*lpJ;
  w[1] = 3.0*lpJ;

  Le = L-4.0*(lpI+lpJ);

  x[2] = 4.0*lpI + 0.5*Le*(1.0-oneOverRoot3);
  w[2] = 0.5*Le;

  x[3] = 4.0*lpI + 0.5*Le*(1.0+oneOverRoot3);
  w[3] = w[2];

  double tmp = 0.0;
  double xL, xL1, wt;
  for (int i = 0; i < 4; i++) {
    xL  = x[i]*oneOverL;
    xL1 = xL-1.0;
    wt = w[i]*oneOverEI;
    fElastic(1,1) += xL1*xL1*wt;
    fElastic(2,2) += xL*xL*wt;
    tmp           += xL*xL1*wt;
  }
  fElastic(1,2) += tmp;
  fElastic(2,1) += tmp;

  return -1;
}

void
HingeRadauBeamIntegration2d::addElasticDeformations(ElementalLoad *theLoad,
						    double loadFactor,
						    double L, double *v0)
{
  return;
}

BeamIntegration*
HingeRadauBeamIntegration2d::getCopy(void)
{
  return new HingeRadauBeamIntegration2d(E, A, I, lpI, lpJ);
}

int
HingeRadauBeamIntegration2d::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(5);

  data(0) = E;
  data(1) = A;
  data(2) = I;
  data(3) = lpI;
  data(4) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "HingeRadauBeamIntegration2d::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
HingeRadauBeamIntegration2d::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  static Vector data(5);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "HingeRadauBeamIntegration2d::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  E   = data(0);
  A   = data(1);
  I   = data(2);
  lpI = data(3);
  lpJ = data(4);

  return 0;
}
  
