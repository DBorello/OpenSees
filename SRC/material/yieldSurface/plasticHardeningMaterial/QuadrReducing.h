// QuadrReducing.h: interface for the QuadrReducing class.
//
//////////////////////////////////////////////////////////////////////


#ifndef QuadrReducing_h
#define QuadrReducing_h

#include "PlasticHardeningMaterial.h"
#include <math.h>

class QuadrReducing : public PlasticHardeningMaterial
{
public:
	QuadrReducing(int tag, double kp0, double kp_half);
	virtual ~QuadrReducing();
	
	double getTrialPlasticStiffness();
    PlasticHardeningMaterial *getCopy(void);
    void Print(OPS_Stream &s, int flag =0);

  private:
  double Kp0;
  double Kp_half;
  double a, b, c;
};

#endif
