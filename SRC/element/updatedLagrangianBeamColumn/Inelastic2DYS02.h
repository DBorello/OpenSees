#ifndef Inelastic2DYS02_H
#define Inelastic2DYS02_H

#include "InelasticYS2DGNL.h"

class Inelastic2DYS02 : public InelasticYS2DGNL
{
public:
  Inelastic2DYS02(int tag, double A, double E, double Iz,
		  int Nd1, int Nd2,
		  YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
		  double pow1, double pow2, int rf_algo=-1, 
		  bool islinear=false, double rho=0.0);
  
  ~Inelastic2DYS02();
  int commitState(void);
  
 protected:
  void getLocalStiff(Matrix &K);
  
 private:
  void getLocalStiff(Matrix &K, double iFactor);
  
 private:
  double A;
  double E;
  double Iz;
  double powKmin;
  double powKunload;
  double fpeak;
  double peakRot;
  double powKmin_orig;
  double rotInElas;
  double forceInElas;
  double sumMi_2;
  double sumMi_1;
  double iFactor, iFactor_hist;
  double IminFactor, IminFactor_hist;
  double sumNatrDisp, sumNatrDisp_hist;
  bool   loading, loading_hist;
  double delPmax1, delPmax2;

  bool   didPlastify_hist; // what is this variable for .. unused in code
};

#endif

