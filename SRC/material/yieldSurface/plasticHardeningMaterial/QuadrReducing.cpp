// QuadrReducing.cpp: implementation of the QuadrReducing class.
//
//////////////////////////////////////////////////////////////////////

#include "QuadrReducing.h"
#include <stdlib.h>

#define MAT_TAG_EXPON -1
#define DEBG 0
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

QuadrReducing::QuadrReducing(int tag, double kp0, double kp_half)
:PlasticHardeningMaterial(tag,MAT_TAG_EXPON),
  Kp0(kp0), Kp_half(kp_half)
{
  a = 2*(Kp0 - 2*Kp_half);
  b = -1*a - Kp0;
  c = Kp0;
}

QuadrReducing::~QuadrReducing()
{

}


double QuadrReducing::getTrialPlasticStiffness()
{
  double K = a*val_trial*val_trial + b*val_trial + c;
  return K;
}


void QuadrReducing::Print(ostream &s, int flag)
{
	s << "QuadrReducing, Tag = " << getTag() << endl;
	s << "Kp0 = " << Kp0 << endl;
	s << "Kp_half = " <<  Kp_half << endl;
}

PlasticHardeningMaterial *QuadrReducing::getCopy(void)
{
 	PlasticHardeningMaterial *theMat = new QuadrReducing(getTag(), Kp0, Kp_half);
    return theMat;
}

