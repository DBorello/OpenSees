
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Mazari-Dafalias  yield criterion                          #
//#                      (Ref. Geotechnique                    			 #
//#                      V.47 No.2 255-272, 1997)                                #
//# CLASS:             MDYieldSurface01                                          #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 03 '00                                             #
//# UPDATE HISTORY:    December 13, '00                                             #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================


#ifndef MD_YS_CPP
#define MD_YS_CPP

#include "MD_YS.h"
#include <basics.h>

#include <G3Globals.h>

//================================================================================
// Normal constructor
//================================================================================

MDYieldSurface::MDYieldSurface( )  { } 


//================================================================================
//create a colne of itself
//================================================================================

YieldSurface * MDYieldSurface::newObj() {  

     YieldSurface  *new_YS = new MDYieldSurface();
     return new_YS;

}

//================================================================================
// Copy constrstructor
//================================================================================
//
//MDYieldSurface::MDYieldSurface(MDYieldSurface &MDYS ) { }
//


//================================================================================
//  Yield criterion evaluation function F(EPState)
//================================================================================

double MDYieldSurface::f(const EPState *EPS) const
{
  stresstensor S = EPS->getStress().deviator();
  //cout << "s " << S;

  double p =EPS->getStress().p_hydrostatic();
  //cout << "p " << p;

  stresstensor alpha = EPS->getTensorVar( 1 );
  //cout << "alpha " << alpha;
  
  double m = EPS->getScalarVar(1);

  stresstensor temp1 = S - alpha*p;
  temp1.null_indices();
  stresstensor temp2 = temp1("ij") * temp1("ij");  

  double temp3 = sqrt( temp2.trace() );

  temp3 = temp3 - sqrt(2.0/3.0) * m * p;        
  //printf("\n========Inside f  temp3 = %.4f x = %.4f\n ", temp3, x);

  return temp3;
}


//================================================================================
// tensor dF/dsigma_ij
//================================================================================

tensor MDYieldSurface::dFods(const EPState *EPS) const {
  
  tensor dFoverds( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();
  //S.reportshort("S");

  double p = EPS->getStress().p_hydrostatic();
  //printf("Here we go!  p %f\n", p);
  	    
  stresstensor alpha = EPS->getTensorVar( 1 ); // getting alpha_ij from EPState
  //alpha.reportshort("alpha");
  
  //stresstensor n = EPS->getTensorVar( 3 );     // getting n_ij from EPState
  //n.reportshort("n");
  
  //-------------------------------------------------
  // might be moved to Evolution Law
    stresstensor r = S * (1.0 / p);
    //r.reportshort("r");
    stresstensor r_bar = r - alpha;
    //r_bar.reportshort("r_bar"); 
    stresstensor norm2 = r_bar("ij") * r_bar("ij");
    double norm = sqrt( norm2.trace() );
    
    //cout << "d_macheps " << d_macheps() << endln;

    stresstensor n;
    if ( norm >= d_macheps() ){ 
      n = r_bar *(1.0 / norm );
    }
    else {
      g3ErrorHandler->fatal("MDYieldSurface::dFods  |n_ij| = 0, divide by zero! Program exits.");
      exit(-1);
    }
    //EPS->setTensorVar( 3, n); //update n_ij//
  //-------------------------------------------------


  double m = EPS->getScalarVar( 1 );

  
  //tensorial multiplication produces 1st-order tensor
  //tensor temp = n("ij") * n("ij");
  //double temp1 = temp.trace();
  //printf("==== n_ij*n_ij %e\n", temp1);

  //!!Very important:  N = n_pq * alpha_pq +sqrt(2/3)*m (always) = n_pq*r_pq(not hold when not on the yield surface)
  tensor temp = n("ij") * alpha("ij");
  double N = temp.trace() + sqrt(2.0/3.0)*m; 
  //printf("    N =  %e\n", N);

  dFoverds =  n - I2 * N*(1.0/3.0);

  return dFoverds;

}


//================================================================================
// double xi_s1 = dF/dm = -(2/3)^0.5 p  Derivative in terms of first scalar var 
//================================================================================

double MDYieldSurface::xi_s1( const EPState *EPS ) const {

    //double p = EPS->getStress().Iinvariant1()/3.0;
    double p = EPS->getStress().p_hydrostatic();
    return -1.0*sqrt(2.0/3.0) * p;

}

//================================================================================
// double xi_t1 = dF/dt1 = dF/dalpha= -p*n  Derivative in terms of first tensorial var 
//================================================================================

tensor MDYieldSurface::xi_t1( const EPState *EPS) const {
  tensor dFoverds( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();

  double p = EPS->getStress().p_hydrostatic();
  	    
  stresstensor alpha = EPS->getTensorVar( 1 ); // getting alpha_ij from EPState
  
  stresstensor r = S * (1.0 / p);
  stresstensor r_bar = r - alpha;
  stresstensor norm2 = r_bar("ij") * r_bar("ij");
  double norm = sqrt( norm2.trace() );

  stresstensor n;
  if ( norm >= d_macheps() ){ 
    n = r_bar *(1.0 / norm );
  }
  else {
    g3ErrorHandler->fatal("MDYieldSurface::dFods  |n_ij| = 0, divide by zero! Program exits.");
    exit(-1);
  }
    
  return  n  * (-1.0)* p;
}


ostream& operator<< (ostream& os, const MDYieldSurface & YS)
{
   os << "Manzari-Dafalias Yield Surface Parameters: " << endln;
   return os;
}

#endif

