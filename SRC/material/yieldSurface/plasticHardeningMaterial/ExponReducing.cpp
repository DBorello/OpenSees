// ExponReducing.cpp: implementation of the ExponReducing class.
//
//////////////////////////////////////////////////////////////////////

#include "ExponReducing.h"
#include <stdlib.h>

#define MAT_TAG_EXPON -1
#define DEBG 0
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ExponReducing::ExponReducing(int tag, double kp0, double alfa)
:PlasticHardeningMaterial(tag,MAT_TAG_EXPON),
  Kp0(kp0), alpha(alfa), resFactor(0.0)
{
}

ExponReducing::ExponReducing(int tag, double kp0, double alfa, double min_fact)
:PlasticHardeningMaterial(tag,MAT_TAG_EXPON),
  Kp0(kp0), alpha(alfa), resFactor(min_fact)
{
//	cout << "ResFact = " <<  res_fact << endl; cin.get();
}


ExponReducing::~ExponReducing()
{

}


double ExponReducing::getTrialPlasticStiffness()
{
	double K ;//= Kp0*exp(-1*val_trial*alpha);

	// if x0 and Kp0 is a const:
	// K = Kp0(1.0  - exp(-alpha*x0 + alpha*val_trial));	
	// K = Kp0*(1.0 - exp(-alpha + alpha*val_trial));
	
	// for pinching type stuff
	K = Kp0*(1 - exp(-1*alpha*val_trial));

	if(sFactor != 1.0)
		K = Kp0*sFactor;
	
	if(K < (Kp0*resFactor))
		K = Kp0*resFactor;

//	cout << "K = " << K << ", sFactor = " << sFactor << endl;
	
	if(K <0.0)
	{
		cerr << "Ri = " << val_trial << ", Factor = " << K/Kp0 << ", res_fact = " << resFactor << endl;
		cin.get();
	}
	
	return K;
}


void ExponReducing::Print(ostream &s, int flag)
{
	s << "MultiLinear, Tag = " << getTag() << endl;
	s << "Kp0 = " << Kp0 << endl;
	s << "Alpha = " <<  alpha << endl;
}

PlasticHardeningMaterial *ExponReducing::getCopy(void)
{
 	PlasticHardeningMaterial *theMat = new ExponReducing(getTag(), Kp0, alpha, resFactor);
    return theMat;
}

