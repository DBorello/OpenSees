//================================================================================
//# COPY LEFT and RIGHT:                                                         #
//# Commercial    use    of    this  program without express permission of the   #
//# University  of  California, is strictly encouraged. Copyright and Copyleft   #
//# are covered by the following clause:                                         #
//#                                                                              #
//# Woody's license:                                                             #
//# ``This    source    code is Copyrighted in U.S., by the The Regents of the   #
//# University  of  California,  for  an indefinite period, and anybody caught   #
//# using  it  without  our  permission,  will be mighty good friends of ourn,   #
//# cause  we  don't give a darn. Hack it. Compile it. Debug it. Run it. Yodel   #
//# it. Enjoy it. We wrote it, that's all we wanted to do.'' bj                  #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Rounded Mohr Coulomb Potential Surface                    #
//# CLASS:                                                                       #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++                                                       #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic jeremic@ucdavis.edu                         #
//#                    Zhao Cheng,                                               #
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic                                 #
//#                                                                              #
//#                                                                              #
//# DATE:              12 Feb. 2003                                              #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION: Functions for rounded Mohr-Coulomb potential function     #
//#                                                                              #
//================================================================================


#ifndef RMC01_PS_CPP
#define RMC01_PS_CPP

#include "RMC01_PS.h"
    
//================================================================================
//create a colne of itself
//================================================================================

PotentialSurface * RMC01PotentialSurface::newObj() 
{  
     PotentialSurface  *new_YS = new RMC01PotentialSurface();
     return new_YS;
}



//================================================================================
// tensor dQ/dsigma_ij                        
//================================================================================

tensor RMC01PotentialSurface::dQods(const EPState *EPS) const 
{
  tensor dQoverds( 2, def_dim_2, 0.0);
  double p = EPS->getStress().p_hydrostatic();
  double q = EPS->getStress().q_deviatoric();
  double theta = EPS->CurrentStress.theta(); 
  double temp_phi = EPS->getScalarVar(1)*3.14159265358979/180.0;
  double temp_cohesive = EPS->getScalarVar(2);
  tensor DpoDs = EPS->getStress().dpoverds();
  tensor DqoDs = EPS->getStress().dqoverds();
  tensor DthetaoDs = EPS->getStress().dthetaoverds();
  double a1 = sin(temp_phi);
  double a2 = temp_cohesive*cos(temp_phi);
  double e = (3.0-a1)/(3.0+a1);
  double Frou = g_0(temp_phi, e);
  double a3 = sqrt(6.0)*a2/3.0;
  double dQoverdp = -a1;
  double dQoverdq = a3*Frou;
  double dQoverdtheta = a3*q*g_prime(theta, e);

  dQoverds = DpoDs  * dQoverdp +
             DqoDs  * dQoverdq +
			 DthetaoDs * dQoverdtheta;

  return dQoverds;

}


//================================================================================
// tensor d2Q/dsigma_ij_2  
//================================================================================

tensor RMC01PotentialSurface::d2Qods2( const EPState *EPS ) const 
{
  tensor d2Qoverds2( 2, def_dim_2, 0.0);

  double p = EPS->getStress().p_hydrostatic();
  double q = EPS->getStress().q_deviatoric();
  double theta = EPS->CurrentStress.theta(); 
  double temp_phi = EPS->getScalarVar(1)*3.14159265358979/180.0;
  double temp_cohesive = EPS->getScalarVar(2);
  tensor DpoDs = EPS->getStress().dpoverds();
  tensor DqoDs = EPS->getStress().dqoverds();
  tensor DthetaoDs = EPS->getStress().dthetaoverds();
  double a1 = sin(temp_phi);
  double a2 = temp_cohesive*cos(temp_phi);
  double e = (3.0-a1)/(3.0+a1);
  double Frou = g_0(temp_phi, e);
  double a3 = sqrt(6.0)*a2/3.0;
  double dQoverdp = -a1;
  double dQoverdq = a3*Frou;
  double dQoverdtheta = a3*q*g_prime(theta, e);

  d2Qoverds2 = DpoDs  * dQoverdp +
             DqoDs  * dQoverdq +
			 DthetaoDs * dQoverdtheta;

  return d2Qoverds2;

}


//================================================================================
// friend ostream functions for output
//================================================================================

ostream& operator<< (ostream& os, const RMC01PotentialSurface & PS)
{
    os << "ROUNDED MC Potential Surface Parameters: " << endln;
    return os;
}


#endif

