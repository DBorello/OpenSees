///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Drucker - Prager  yield criterion                         #
//# CLASS:             DPPotentialSurface                                        #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 03 '93                                             #
//# UPDATE HISTORY:    August 08 '00                                             #
//# UPDATE HISTORY:    20Aug2004 ZC added kinematic hardening part               #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

#ifndef DP_PS_CPP
#define DP_PS_CPP

#include "DP_PS.h"



//================================================================================
// Normal constructor
//================================================================================

//DPPotentialSurface::DPPotentialSurface( double a2d = 0.3) : alfa2(a2d) 
//{
//
//} 


//================================================================================
// Copy constrstructor  === not necessary for no pointer member
//================================================================================

DPPotentialSurface::DPPotentialSurface(const DPPotentialSurface &DPPS ) {

    alfa2 =  DPPS.getalfa2();

}


//================================================================================
//create a colne of itself
//================================================================================

PotentialSurface * DPPotentialSurface::newObj() {  

     PotentialSurface  *new_PS = new DPPotentialSurface( this->getalfa2() );
     return new_PS;

}



//================================================================================
// tensor dQ/dsigma_ij = a*delta_ij + 0.5*sqrt(2.0)*(s_ij-alpha_ij)/|s_ij-alpha_ij|
//================================================================================

tensor DPPotentialSurface::dQods(const EPState *EPS) const {
    stresstensor dQods;
    
	tensor KroneckerI("I", 2, def_dim_2);
    //double temp1 =  EPS->getScalarVar(1);    
    double temp1 = getalfa2();
    
    stresstensor alpha;
    stresstensor Tnsr0;
    stresstensor sigma = EPS->getStress();
    double halfRt2 = 0.5 * sqrt(2.0);
    int nod = EPS->getNTensorVar();
    if ( nod >=1 )  { //May not have kinematic hardening
      alpha = EPS->getTensorVar(1); 
    }
    stresstensor sigma_bar = sigma - alpha;   
    stresstensor s_bar = sigma_bar.deviator();
    stresstensor temp3 = s_bar("ij") * s_bar("ij");
    temp3.null_indices();
    double temp4 = temp3.trace();
    temp4 = sqrt(temp4);
    double eps = pow( d_macheps(), 0.5 );
    if ( fabs(temp4) > eps )  {
      Tnsr0 = s_bar * (1.0/temp4);
    }
    	    
    dQods = KroneckerI*temp1 + Tnsr0*halfRt2; 
        
    return dQods;
}


//================================================================================
// tensor d2Q/dsigma_ij_2
//================================================================================

tensor DPPotentialSurface::d2Qods2(const EPState *EPS) const {
    tensor d2Qods2(4, def_dim_4, 0.0);
    
	tensor KroneckerI("I", 2, def_dim_2);
	tensor T1 = KroneckerI("ij")*KroneckerI("mn");
	T1.null_indices();
	tensor T2 = (T1.transpose0110()+T1.transpose0111())*0.5 - T1*(1.0/3.0);
	
    //double temp1 =  EPS->getScalarVar(1);    
    double temp1 = getalfa2();
    
    stresstensor alpha;
    stresstensor sigma = EPS->getStress();
    double halfRt2 = 0.5 * sqrt(2.0);
    int nod = EPS->getNTensorVar();
    if ( nod >=1 )  { //May not have kinematic hardening
      alpha = EPS->getTensorVar(1); 
    }
    stresstensor sigma_bar = sigma - alpha;   
    stresstensor s_bar = sigma_bar.deviator();
    stresstensor temp3 = s_bar("ij") * s_bar("ij");
    temp3.null_indices();
    double temp4 = temp3.trace();
    temp4 = sqrt(temp4);
    tensor temp5 = s_bar("ij")*s_bar("mn");
    double eps = pow( d_macheps(), 0.5 );
    if ( fabs(temp4) > eps )  {
      d2Qods2 =T1 * (1.0/temp4) - temp5*(1.0/(temp4*temp4*temp4));
    }

    return d2Qods2;
}

// For Consistent Algorithm, Z Cheng, Jan 2004
tensor DPPotentialSurface::d2Qodsds1(const EPState *EPS) const 
{  
  tensor I("I", 2, def_dim_2);
  tensor d2Qoverdsds1 = I;
  return d2Qoverdsds1;
}

//================================================================================

double DPPotentialSurface::getalfa2() const {

    return alfa2; 

}

//================================================================================
OPS_Stream& operator<< (OPS_Stream& os, const DPPotentialSurface &PS)
{
   os << "Drucker-Prager Potential Surface Parameters: " << endln;
   os << "alfa2 = " << PS.getalfa2() << endln;
   return os;
}


#endif
















 
