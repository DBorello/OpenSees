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
// $Date: 2003-02-14 03:44:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/HingeRadauTwoBeamIntegration2d.cpp,v $

#include <HingeRadauTwoBeamIntegration2d.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

HingeRadauTwoBeamIntegration2d::HingeRadauTwoBeamIntegration2d(double e,
							       double a,
							       double i,
							       double lpi,
							       double lpj):
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadauTwo2d),
  E(e), A(a), I(i), lpI(lpi), lpJ(lpj)
{
  // Nothing to do
}

HingeRadauTwoBeamIntegration2d::HingeRadauTwoBeamIntegration2d():
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadauTwo2d),
  E(0.0), A(0.0), I(0.0), lpI(0.0), lpJ(0.0)
{

}

HingeRadauTwoBeamIntegration2d::~HingeRadauTwoBeamIntegration2d()
{
  // Nothing to do
}

void
HingeRadauTwoBeamIntegration2d::getSectionLocations(int numSections, double L,
						    double *xi)
{
  double two3oneOverL = (2.0/3.0)/L;

  xi[0] = 0.0;
  xi[1] = lpI*two3oneOverL;
  xi[2] = 1.0-lpJ*two3oneOverL;
  xi[3] = 1.0;
  for (int i = 4; i < numSections; i++)
    xi[i] = 0.0;
}

void
HingeRadauTwoBeamIntegration2d::getSectionWeights(int numSections, double L,
						  double *wt)
{
  double oneOverL = 1.0/L;

  wt[0] = 0.25*lpI*oneOverL;
  wt[1] = 3.0*wt[0];
  wt[3] = 0.25*lpJ*oneOverL;
  wt[2] = 3.0*wt[3];
  for (int i = 4; i < numSections; i++)
    wt[i] = 1.0;
}

int
HingeRadauTwoBeamIntegration2d::addElasticFlexibility(double L, Matrix &fElastic)
{
  double oneOverL = 1.0/L;

  // Length of elastic interior
  double Le = L-lpI-lpJ;
  double LoverEA  = Le/(E*A);
  double Lover3EI = Le/(3*E*I);
  double Lover6EI = 0.5*Lover3EI;
  
  // Elastic flexibility of element interior
  static Matrix fe(2,2);
  fe(0,0) = fe(1,1) =  Lover3EI;
  fe(0,1) = fe(1,0) = -Lover6EI;
  
  // Equilibrium transformation matrix
  static Matrix B(2,2);
  double betaI = lpI*oneOverL;
  double betaJ = lpJ*oneOverL;
  B(0,0) = 1.0 - betaI;
  B(1,1) = 1.0 - betaJ;
  B(0,1) = -betaI;
  B(1,0) = -betaJ;
  
  // Transform the elastic flexibility of the element
  // interior to the basic system
  static Matrix ftmp(2,2);
  ftmp.addMatrixTripleProduct(0.0, B, fe, 1.0);

  fElastic(0,0) += LoverEA;
  fElastic(1,1) += ftmp(0,0);
  fElastic(1,2) += ftmp(0,1);
  fElastic(2,1) += ftmp(1,0);
  fElastic(2,2) += ftmp(1,1);

  return -1;
}

void
HingeRadauTwoBeamIntegration2d::addElasticDeformations(ElementalLoad *theLoad,
						       double loadFactor,
						       double L, double *v0)
{
  return;
}

BeamIntegration*
HingeRadauTwoBeamIntegration2d::getCopy(void)
{
  return new HingeRadauTwoBeamIntegration2d(E, A, I, lpI, lpJ);
}

int
HingeRadauTwoBeamIntegration2d::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(5);

  data(0) = E;
  data(1) = A;
  data(2) = I;
  data(3) = lpI;
  data(4) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    g3ErrorHandler->warning("HingeRadauTwoBeamIntegration2d::sendSelf() - %s\n",
			    "failed to send Vector data");
    return -1;
  }    

  return 0;
}

int
HingeRadauTwoBeamIntegration2d::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  static Vector data(5);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    g3ErrorHandler->warning("HingeRadauTwoBeamIntegration2d::recvSelf() - %s\n",
			    "failed to receive Vector data");
    return -1;
  }
  
  E   = data(0);
  A   = data(1);
  I   = data(2);
  lpI = data(3);
  lpJ = data(4);

  return 0;
}
