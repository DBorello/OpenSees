//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              19AUg2003
//# UPDATE HISTORY:    Sept 2003
//#
//#
//===============================================================================

#ifndef LogWEnergy_CPP
#define LogWEnergy_CPP

#include <LogWEnergy.h>

//================================================================================
// Normal constructor
//================================================================================
LogWEnergy::LogWEnergy(double E_in, double nu_in )
{
   E = E_in;
   nu = nu_in;
   G = E / 2.0 / (1.0 + nu);
   K = E / 3.0 / (1.0 - 2.0 * nu);
}

//LogWEnergy::LogWEnergy( )
//{
//
//}

//================================================================================
// Normal destructor
//================================================================================
LogWEnergy::~LogWEnergy( )
{

}

//================================================================================
//create a clone of itself
//================================================================================
WEnergy * LogWEnergy::newObj( )
  {
    LogWEnergy  *new_WEnergy = new LogWEnergy( E,  nu);
    return new_WEnergy;
  }

const double LogWEnergy::getE()
  {
    return E;
  }

const double LogWEnergy::getnu()
  {
    return nu;
  }

//================================================================================
// w
//================================================================================
const double LogWEnergy::wE(const double &J_in, const Vector &lambda_wave_in)
  {
    double w_iso = G * ( log(lambda_wave_in(0)) * log(lambda_wave_in(0))
                                  + log(lambda_wave_in(1)) * log(lambda_wave_in(1))
                                  + log(lambda_wave_in(2)) * log(lambda_wave_in(2)) );
    double w_vol = 0.5 * K * (log(J_in) * log(J_in));
    double w_total = w_iso + w_vol;
    return w_total;
  }

//================================================================================
// d(iso)w / d(lambda)
//================================================================================
const Vector LogWEnergy::disowOdlambda(const Vector &lambda_wave_in)
  {
    Vector disowOverdlambda(3);
    disowOverdlambda(0) = 2.0 * G * log(lambda_wave_in(0)) / lambda_wave_in(0);
    disowOverdlambda(1) = 2.0 * G * log(lambda_wave_in(1)) / lambda_wave_in(1);
    disowOverdlambda(2) = 2.0 * G * log(lambda_wave_in(2)) / lambda_wave_in(2);
    return disowOverdlambda;
  }

//================================================================================
// d2(iso)w / d(lambda)2
//================================================================================
const Vector LogWEnergy::d2isowOdlambda2(const Vector &lambda_wave_in)
  {
    Vector d2isowOverdlambda2(3);
    d2isowOverdlambda2(0) = -2.0 * G * (1.0-log(lambda_wave_in(0))) / lambda_wave_in(0) / lambda_wave_in(0);
    d2isowOverdlambda2(1) = -2.0 * G * (1.0-log(lambda_wave_in(1))) / lambda_wave_in(1) / lambda_wave_in(1);
    d2isowOverdlambda2(2) = -2.0 * G * (1.0-log(lambda_wave_in(2))) / lambda_wave_in(2) / lambda_wave_in(2);
    return d2isowOverdlambda2;
  }

//================================================================================
// d(vol)w / dJ
//================================================================================
const double LogWEnergy::dvolwOdJ(const double &J_in)
{
//   printf("J=%lf\n",J_in);
   double dcolwOverdJ = K * log(J_in) / J_in;
//   printf("dW/dJ=%lf\n",dcolwOverdJ);
   return  dcolwOverdJ;
}

//================================================================================
// d2(vol)w / dJ2
//================================================================================
const double LogWEnergy::d2volwOdJ2(const double &J_in)
{
//   printf("J=%lf\n",J_in);
   double d2colwOverdJ2 = K * (1.0 - log(J_in) ) / J_in / J_in ;
//   printf("d2W/dJ2=%lf\n",d2colwOverdJ2);
   return  d2colwOverdJ2;
}

//================================================================================
// friend ostream functions for output
//================================================================================
//OPS_Stream& operator<< (OPS_Stream& os, const LogWEnergy & WEnergy)
//{
//    os << "Logarithmic Strain Energy Function: " << endln;
//    return os;
//}

#endif

