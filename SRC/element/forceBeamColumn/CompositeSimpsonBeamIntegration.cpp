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

// $Revision$
// $Date$
// $Source$

#include <CompositeSimpsonBeamIntegration.h>

CompositeSimpsonBeamIntegration::CompositeSimpsonBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_CompositeSimpson)
{
  // Nothing to do
}

CompositeSimpsonBeamIntegration::~CompositeSimpsonBeamIntegration()
{
  // Nothing to do
}

BeamIntegration*
CompositeSimpsonBeamIntegration::getCopy(void)
{
  return new CompositeSimpsonBeamIntegration();
}

void
CompositeSimpsonBeamIntegration::getSectionLocations(int numSections, double L,
						double *xi)
{
  // Check that num sections is odd
  if (numSections % 2 == 1) {
    int numIntervals = (numSections+1)/2; // Num intervals is even
    double h = 1.0/numIntervals;
    xi[0] = 0.0;
    xi[numSections-1] = 1.0;
    for (int i = 1; i < numSections-1; i++)
      xi[i] = h*i;    
  }
  else {
    opserr << "CompositeSimpson, numSections must be odd (" << numSections << " was input)" << endln;
  }
}

void
CompositeSimpsonBeamIntegration::getSectionWeights(int numSections, double L,
					      double *wt)
{
  // Check that num sections is odd
  if (numSections % 2 == 1) {
    int numIntervals = (numSections+1)/2; // Num intervals is even
    double h = 1.0/numIntervals;
    wt[0] = h/3.0;
    wt[numSections-1] = h/3.0;
    for (int i = 1; i < numSections; i += 2)
      wt[i] = 4*h/3.0;  
    for (int i = 2; i < numSections-1; i += 2)
      wt[i] = 2*h/3.0;  
  }
  else {
    opserr << "CompositeSimpson, numSections must be odd (" << numSections << " was input)" << endln;
  }
}

void
CompositeSimpsonBeamIntegration::Print(OPS_Stream &s, int flag)
{
  s << "CompositeSimpson" << endln;
}
