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
//# CLASS:             RMC01YieldSurface                                         #
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
//# SHORT EXPLANATION: Functions for rounded Mohr-Coulomb yield function         #
//#                                                                              #
//================================================================================

#ifndef RMC01_YS_CPP
#define RMC01_YS_CPP

#include "RMC01_YS.h"

//================================================================================
//create a clone of itself
//================================================================================

YieldSurface * RMC01YieldSurface::newObj() 
{  
     YieldSurface  *new_YS = new RMC01YieldSurface();
     return new_YS;
}

//================================================================================
//  Yield criterion evaluation function F(EPState)
//================================================================================

double RMC01YieldSurface::f(const EPState *EPS) const {
    double p = EPS->getStress().p_hydrostatic();
    double q = EPS->CurrentStress.q_deviatoric();
	double theta = EPS->CurrentStress.theta(); 
    double temp_phi = EPS->getScalarVar(1)*3.14159265358979/180.0;
	double temp_cohesive = EPS->getScalarVar(2);
    double a1 = sin(temp_phi);
	double a2 = temp_cohesive*cos(temp_phi);
	double e = (3.0-a1)/(3.0+a1);
	double Frou = g_0(theta, e);
	double a3 = sqrt(6.0)*a2/3.0;
    double f = -a1*p+a3*q*Frou-a2;
    return f;
}

//================================================================================
// tensor dF/dsigma_ij  ( eq. 5.206 in Chen )
//================================================================================

tensor RMC01YieldSurface::dFods(const EPState *EPS) const {
  
  tensor dFoverds( 2, def_dim_2, 0.0);

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
  double dFoverdp = -a1;
  double dFoverdq = a3*Frou;
  double dFoverdtheta = a3*q*g_prime(theta, e);

  dFoverds = DpoDs  * dFoverdp +
             DqoDs  * dFoverdq +
			 DthetaoDs * dFoverdtheta;

  return dFoverds;

}



//================================================================================
// double xi_s1 = dF/dS1 = dF/dalfa1 = I1  Derivative in terms of first scalar var 
//================================================================================

double RMC01YieldSurface::xi_s1( const EPState *EPS ) const {
    double p = EPS->getStress().p_hydrostatic();
    double q = EPS->CurrentStress.q_deviatoric();
	double theta = EPS->CurrentStress.theta(); 
    double temp_phi = EPS->getScalarVar(1)*3.14159265358979/180.0;
	double temp_cohesive = EPS->getScalarVar(2);
    double a1 = sin(temp_phi);
	double a2 = temp_cohesive*cos(temp_phi);
	double e = (3.0-a1)/(3.0+a1);
	double Frou = g_0(theta, e);
	double a3 = sqrt(6.0)*a2/3.0;
    double temp1 = -p*cos(temp_phi);
	double temp2 = -a3*temp_cohesive*a1*q*Frou;
	double temp3 = temp_cohesive*a1;
	double temp = temp1+temp2+temp3;
	return temp;
}

//================================================================================
// double xi_s2 = dF/dS2 = dF/k = -1.0  Derivative in terms of second scalar var 
//================================================================================

double RMC01YieldSurface::xi_s2( const EPState *EPS ) const {
//    double p = EPS->getStress().p_hydrostatic();
    double q = EPS->CurrentStress.q_deviatoric();
	double theta = EPS->CurrentStress.theta(); 
    double temp_phi = EPS->getScalarVar(1)*3.14159265358979/180.0;
	double temp_cohesive = EPS->getScalarVar(2);
    double a1 = sin(temp_phi);
	double a2 = temp_cohesive*cos(temp_phi);
	double e = (3.0-a1)/(3.0+a1);
	double Frou = g_0(theta, e);
	double a3 = sqrt(6.0)*a2/3.0;
	double a4 = cos(temp_phi);
	double temp1 = a3*a4*q*Frou;
	double temp = temp1-a4;
	return temp;
}


//================================================================================
// double xi_t1 = dF/dt1  Derivative in terms of first tensorial var 
//================================================================================

//tensor RMC01YieldSurface::xi_t1( ) const 
//{


//}

//================================================================================
OPS_Stream& operator<< (OPS_Stream& os, const RMC01YieldSurface & YS)
{
   os << "Rounded Mohr Coulomb Surface Parameters: " << endln;
   return os;

}

#endif

