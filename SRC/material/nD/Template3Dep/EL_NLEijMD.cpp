/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
#                                                                                #
# CLASS:             EvolutionLaw_NL_EijMD (on plastic strain)                     #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              09-13-2000                                                  #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: This is a nonlinear evolution law for the evoltion of a     #
#                    tensorial variable alpha which depends on plastic strain    #
#                    i.e. dalpha_ij = 2/3*ha*dE_ij -Cr*de_eq*alpha_ij( Amstrong- #
//                   Frederick Model)                                            #
//================================================================================
*/

#ifndef EL_NLEIJMD_CPP
#define EL_NLEIJMD_CPP

#include "EL_NLEijMD.h"
#include <basics.h>
    

//================================================================================
// Default constructor
//================================================================================
EvolutionLaw_NL_EijMD::EvolutionLaw_NL_EijMD( double eod,
                    double ad, 
    		    double Mcd,    
    		    double Med,    
    		    double Lambdad,
    		    double ec_refd, 
    		    double p_refd, 
    		    double kc_bd, 
    		    double kc_dd, 
    		    double ke_bd, 
    		    double ke_dd, 
    		    double hod,	 
    		    double Cmd,
    		    double ed,    
    		    double Aod,    
    		    double Fmaxd,   
    	            double Cfd):    
  eo(eod), a(ad), Mc(Mcd), Me(Med), Lambda(Lambdad), ec_ref(ec_refd), p_ref(p_refd),
  kc_b(kc_bd), kc_d(kc_dd), ke_b(ke_bd), ke_d(ke_dd), ho(hod), Cm(Cmd),
  e(eod), D(0.0), Ao(Aod), Fmax(Fmaxd), Cf(Cfd)
  { } 

//================================================================================
// Copy constructor
//================================================================================

EvolutionLaw_NL_EijMD::EvolutionLaw_NL_EijMD(const EvolutionLaw_NL_EijMD &MDE ) {

    this->a      = MDE.geta();
    this->Mc     = MDE.getMc();
    this->Me     = MDE.getMe();
    this->Lambda = MDE.getLambda();
    this->ec_ref = MDE.getec_ref();
    this->p_ref=  MDE.getp_ref();
    this->kc_b =  MDE.getkc_b();
    this->kc_d =  MDE.getkc_d();
    this->ke_b =  MDE.getke_b();
    this->ke_d =  MDE.getke_d();
    this->ho   =  MDE.getho();
    this->eo   =  MDE.geteo();
    this->e   =  MDE.gete();
    this->Cm   =  MDE.getCm();
    //D is also needed in PS, so it is copied to 2nd cell of scalar var array
    this->D    =  MDE.getD();  
    this->Ao   =  MDE.getAo();
    this->Fmax =  MDE.getFmax(); 
    this->Cf   =  MDE.getCf();
    this->F    =  MDE.getF(); 
}


//================================================================================
//  Create a clone of itself 
//================================================================================
EvolutionLaw_NL_EijMD * EvolutionLaw_NL_EijMD::newObj() {
    
    //EvolutionLaw_T *newEL = new EvolutionLaw_NL_EijMD( *this );
    EvolutionLaw_NL_EijMD *newEL = new EvolutionLaw_NL_EijMD( *this );
    
    return newEL;

}
    
//================================================================================
// Evaluating h_s = h b_ij = h(alpha_ij_theta_b - alpha_ij) (For the evaluation of Kp)
//================================================================================

tensor EvolutionLaw_NL_EijMD::h_t( EPState *EPS, PotentialSurface *PS)
{
    //=========================================================================
    //calculate  n_ij
    stresstensor S = EPS->getStress().deviator();
    double p = EPS->getStress().p_hydrostatic();
    stresstensor alpha = EPS->getTensorVar( 1 );  // alpha_ij
   
    stresstensor r = S * (1.0 / p);
    //r.reportshort("r");
    stresstensor r_bar = r - alpha;
    stresstensor norm2 = r_bar("ij") * r_bar("ij");
    double norm = sqrt( norm2.trace() );
    
    stresstensor n;
    if ( norm >= d_macheps()){ 
      n = r_bar *(1.0 / norm );
    }
    else {
      g3ErrorHandler->fatal("EvolutionLaw_NL_EijMD::h_t  |n_ij| = 0, divide by zero! Program exits.");
      exit(-1);
    }
    //cout << "nij = " << n;
    
    
    //=========================================================================
    //calculating b_ij for Kp and d_ij for updating D

    //Calculate the state parameters xi 
    double e = gete(); // Stored in MD ELT
    
    if (p < 0.0)
    {
       g3ErrorHandler->fatal("EvolutionLaw_NL_EijMD::h_t   p < 0, Program exits.");
       exit(-1);
    }
    double ec = getec_ref() - getLambda() * log( p/getp_ref() );

    double xi = e - ec;
    cout << "e = " << e << " ec = " << ec << " xi = " << xi << endln;

    //Calculating the lode angle theta
    double J2_bar = r_bar.Jinvariant2();
    double J3_bar = r_bar.Jinvariant3();
    double tempd = 3.0*pow(3.0, 0.5)/2.0*J3_bar/ pow( J2_bar, 1.5);
    //cout << "theta -1 =" << tempd << endln;
    
    if (tempd > 1.0 ) tempd = 1.0; //bug. if tempd = 1.00000000003, acos gives NaN
    if (tempd < -1.0 ) tempd = -1.0;
    
    double theta = acos( tempd ) / 3.0;
    //cout << "theta = " << theta << endln;
    
    //calculate the alpha_theta_b and alpha_theta_d
    double m = EPS->getScalarVar(1);
    double c = getMe() / getMc();

    double cd = getke_d() / getkc_d();
    //stresstensor alpha_theta_d = n("ij") * (g_WW(theta, c) * Mc + g_WW(theta, cd) * kc_d * xi - m) * pow(2.0/3.0, 0.5);
    stresstensor alpha_theta_d = n("ij") * (g_A(theta, c) * Mc + g_A(theta, cd) * kc_d * xi - m) * pow(2.0/3.0, 0.5);
    //cout << "alpha_theta_d " << alpha_theta_d<<" g_WW(theta, c) "<< g_WW(theta, c) << endln;

    double cb = getke_b() / getkc_b();
    if ( xi > 0.0 ) xi = 0.0;  // < -xi > for alpha_theta_b
    //stresstensor alpha_theta_b = n("ij") * (g_WW(theta, c) * Mc - g_WW(theta, cb) * kc_b * xi - m) * pow(2.0/3.0, 0.5);
    stresstensor alpha_theta_b = n("ij") * (g_A(theta, c) * Mc - g_A(theta, cb) * kc_b * xi - m) * pow(2.0/3.0, 0.5);
    alpha_theta_b.null_indices();

    //=========================================================================
    // calculating h
    stresstensor b;
    b =  alpha_theta_b - alpha;
    b.null_indices();
    stresstensor d;
    d =  alpha_theta_d - alpha;
    d.null_indices();

    //double alpha_c_b = g_WW(0.0, c) * Mc + g_WW(0.0, cb) * kc_b * (-xi) - m;
    double alpha_c_b = g_A(0.0, c) * Mc + g_A(0.0, cb) * kc_b * (-xi) - m;
    double b_ref = 2.0 * pow(2.0/3.0, 0.5) * alpha_c_b;
    tensor temp1 = b("ij") * n("ij");
    double bn = temp1.trace();

    double h = getho() * fabs(bn) / ( b_ref - fabs(bn) ); 
    //cout << "ho =" << getho()  << "   h =" << h << endln;
    //cerr << "\nbn =" << bn << endln; 

    //Updating D and F---need to fine tune
    //temp1 = d("ij") * n("ij");
    //double dn = temp1.trace();
    //cout << "bn =" << bn << "  dn =" << dn << endln; 
    
    //// Calculating A
    //stresstensor F = EPS->getTensorVar( 2 );   // getting  F_ij from EPState
    //temp1 = F("ij") * n("ij");
    //double temp = temp1.trace();
    //if (temp < 0)   temp = 0;   
    //double A = getAo()*(1.0 + temp);

        
    tensor ht = h * b;

    return ht;

}


//================================================================================
//  Print vars defined in Linear Evolution Law
//================================================================================
void EvolutionLaw_NL_EijMD::print()
{
    cout << (*this);
}
    
//================================================================================
// prints Manzari-Dafalia EvolutionLaw's contents 
//================================================================================
ostream& operator<< (ostream& os, const EvolutionLaw_NL_EijMD & MDEL)
{
    os.unsetf( ios::scientific );
    os.precision(5);

    //os.width(10);       
    os << endln << "Manzari-Dafalias Evolution Law's parameters:" << endln;
    os << "a = " << MDEL.geta() << "; ";
    os << "Mc = " << MDEL.getMc() << "; ";
    //os.width(10);       
    os << "Me = "  << MDEL.getMe() << "; ";
    //os.width(10);       
    os << "Lambda = " << MDEL.getLambda() << "; ";
    //os.width(10);       
    os << "ec_ref = " << MDEL.getec_ref() << "; ";
    //os.width(10);       
    os << "p_ref = " << MDEL.getp_ref() << "kPa"  << "; " << endln;

    //os.width(10);       
    os << "kc_b = " << MDEL.getkc_b() << "; ";
    //os.width(10);       
    os << "kc_d = " << MDEL.getkc_d() << "; ";
    //os.width(10);       
    os << "ke_b = " << MDEL.getke_b() << "; ";
    //os.width(10);       
    os << "ke_d = " << MDEL.getke_d() << "; " << endln;

    //os.width(10);       
    //os << "h = " << MDEL.h << "; ";
    //os.width(10);       
    os << "ho = " << MDEL.getho() << "; ";
    //os.width(10);       
    os << "Cm = " << MDEL.getCm() << "; " << endln;

    //os.width(10);       
    os << "D = " << MDEL.getD() << "; ";
    //os.width(10);       
    os << "Ao = " << MDEL.getAo() << "; ";
    //os.width(10);       
    os << "Fmax = " << MDEL.getFmax() << "; ";
    //os.width(10);       
    os << "Cf = " << MDEL.getCf() << "; " << endln; 
    os << "F = " << MDEL.getF() << endln; 
           
    return os;
}  


//================================================================================
//  Initialize some  vars in EPState				        
//================================================================================
//void EvolutionLaw_L_EijMD::InitVars(EPState  *EPS) 

int EvolutionLaw_NL_EijMD::updateEeDm(EPState *EPS, double st_vol, double dLamda)
{
    int err = 0;
    // Updating E_Young corresponding to current stress state
    double p_atm = 100.0; //Kpa atmospheric pressure
    double p = EPS->getStress().p_hydrostatic();//Pc???
    
    double E = EPS->getEo();//* pow( (p/p_atm), geta());
    EPS->setE( E );

    // Updating e
    //double st_vol = st.p_hydrostatic();
    double e = gete();
    double de = -(1.0 + geteo())*st_vol;
    e = e + de;
    this->e = e;
    EPS->setScalarVar(3, e); // e also stored in scalar array's 3nd cell for PS

    //Updating D and m if dLamda != 0.0
    //if (dLamda != 0.0)
    //{

       //=========================================================================
       //calculate  n_ij
       stresstensor S = EPS->getStress().deviator();
       stresstensor alpha = EPS->getTensorVar( 1 );  // alpha_ij

       stresstensor r = S * (1.0 / p);
       //r.reportshort("r");
       stresstensor r_bar = r - alpha;
       stresstensor norm2 = r_bar("ij") * r_bar("ij");
       double norm = sqrt( norm2.trace() );
       
       stresstensor n;
       if ( norm >= d_macheps()){ 
         n = r_bar *(1.0 / norm );
       }
       else {
         g3ErrorHandler->fatal("EvolutionLaw_L_EijMD::dFods  |n_ij| = 0, divide by zero! Program exits.");
         err += 1;
	 exit(-1);
       }
       	         
       cout << " alpha = " << alpha;
       //cout << "\n **** n =     " << n;


       //Calculate the state parameters xi 

       if (p < 0.0)
       {
          g3ErrorHandler->fatal("EvolutionLaw_NL_EijMD::updateEeDm   p < 0, Program exits.");
          exit(-1);
       }
       
       double ec = getec_ref() - getLambda() * log( p/getp_ref() );       
       double xi = e - ec;
       
       //=========================================================================
       // Update m
       double m = EPS->getScalarVar(1);
       double Cm = getCm(); 
       double eo = geteo(); 
       double dm = dLamda * Cm * ( 1.0 + eo ) * D;
       m = m + dm;
       EPS->setScalarVar(1, m);
       //cout  << endln << "dm = " << dm << endln;

       //=========================================================================
       // Update F
       stresstensor dF;
       stresstensor F = EPS->getTensorVar( 2 );   // getting  F_ij from EPState
       double D = getD();
       //if (D != EPS->getScalarVar(2)) {
       //  g3ErrorHandler->warning("EvolutionLaw_L_EijMD::updateEeDm  D values contradict:%f %f ", D, EPS->getScalarVar(2));
       //	 err += 1;
       //  //exit(-1);
       //}         

       if ( D > 0.0 ) D = 0.0;
       dF =  dLamda * getCf() * (-D) * ( getFmax() * n("ij") + F("ij") );
       //cout << "dF" << dF;       
       F = F ;//- dF;
       EPS->setTensorVar(2, F);
       this->F = F;
       cout << " \n F =          " << F << endln;       

       tensor temp_tensor = F("ij") * n("ij");
       double temp = temp_tensor.trace();
       if (temp < 0)   temp = 0;
       double A = getAo()*(1.0 + temp);
       
       //Calculating the lode angle theta
       double J2_bar = r_bar.Jinvariant2();
       double J3_bar = r_bar.Jinvariant3();
       double tempd = 3.0*pow(3.0, 0.5)/2.0*J3_bar/ pow( J2_bar, 1.5);
       
       if (tempd > 1.0 ) tempd = 1.0; //bug. if tempd = 1.00000000003, acos gives nan
       if (tempd < -1.0 ) tempd = -1.0;
       
       double theta = acos( tempd ) / 3.0;
       cout << "theta = " << theta << endln;
       
       //Calculate alpha_theta_d
       double c = getMe() / getMc();
       
       double cd = getke_d() / getkc_d();
       
       //Testing g_WW
       //cout << "g_A(0, c) " << g_A(0.0, 0.0167) << endln;
       //cout << "g_A(60, c) " << g_A(1.0472, 0.0167) << endln;
       double alpha_theta_dd = g_A(theta, c) * Mc + g_A(theta, cd) * kc_d * xi - m;
       stresstensor alpha_theta_d = n("ij") * alpha_theta_dd * pow(2.0/3.0, 0.5);
       cout << " cd "<< cd << " c " << c << "alpha_theta_d_norm " << alpha_theta_dd << endln;
       cout << "alpha_theta_d " << alpha_theta_d<<" g_A(theta, cd) "<< g_A(theta, cd) << " Mc " << Mc << endln;
       		        
       stresstensor d;
       d =  alpha_theta_d - alpha;
       d.null_indices();
       EPS->setTensorVar(3, d); //d is also stored in Tensor array for 2nd derivative eval in PS 

       cout << " d = "  << d << " n = "  << n << endln;    
       //=========================================================================
       // Updating D --does not depend on dLamda
       tensor temp1 = d("ij") * n("ij");
       temp1.null_indices();
       double dn =  temp1.trace();
       
       // Restrictions 
       if ( xi > 0.0 && dn < 0.0)
         dn = 0.0;
       
       double newD = dn * A;
       EPS->setScalarVar(2, newD); // D also stored in scalar array's 2nd cell for PS
       this->D = newD;
       cout << " dn = "  << dn << " A = " << A << endln << endln;
       cout << "   ***** D =  " << newD << endln;    
       //cout << "alpha_theta_d " << alpha_theta_d<<" g_WW(theta, c) "<< g_WW(theta, c);
    
    //}
    
    return err;
}   


//================================================================================
//  Set initial value of D once the current stress hit the yield surface     	
//  for MD Hardening model only    						
//================================================================================
//void EvolutionLaw_L_EijMD::setInitD(EPState  *EPS) 
//
//void EvolutionLaw_NL_EijMD::updateDF(EPState  *EPS) 
//{
//
//}

//================================================================================
// Private accessory functions

//================================================================================
double EvolutionLaw_NL_EijMD::geta() const 
{
    return a;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getMc() const 
{
    return Mc;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getMe() const 
{
    return Me;
}


//================================================================================
double EvolutionLaw_NL_EijMD::getLambda() const
{
    return Lambda;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getec_ref() const
{
    return ec_ref;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getp_ref() const 
{
    return p_ref;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getkc_b() const
{
    return kc_b;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getkc_d() const  
{
    return kc_d;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getke_b() const
{
    return ke_b;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getke_d() const
{
    return ke_d;
}

//================================================================================
//double EvolutionLaw_NL_EijMD::geth() const
//{       
//    return h;
//}

//================================================================================
double EvolutionLaw_NL_EijMD::getho() const
{      
    return ho;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getCm() const
{       
    return Cm;
}

//================================================================================
double EvolutionLaw_NL_EijMD::geteo() const
{       
    return eo;
}

//================================================================================
double EvolutionLaw_NL_EijMD::gete() const
{       
    return e;
}


//================================================================================
double EvolutionLaw_NL_EijMD::getAo() const
{       
    return Ao;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getD() const
{       
    return D;
}


//================================================================================
stresstensor EvolutionLaw_NL_EijMD::getF() const
{       
    return F;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getFmax() const
{       
    return Fmax;
}

//================================================================================
double EvolutionLaw_NL_EijMD::getCf() const
{       
    return Cf;
}


//================================================================================
//  Interpolation function No.1  -- Agyris: g_A(theta, e) 
//================================================================================

double EvolutionLaw_NL_EijMD::g_A(double theta, double e) 
{
    double temp = 2.0 * e;
    temp  = temp / ( (1.0 + e) - (1.0 - e) * cos(3.0*theta) ); 

    return temp;
}


//================================================================================
//  Interpolation function No.2  -- Willan-Warkne: g_WW(theta, e)
//================================================================================
//g_WW(60, 0.0167) =  +3.3818e-02?????

double EvolutionLaw_NL_EijMD::g_WW(double theta, double e) 
{
    //Rotate 60 degree to produce 1 at 0.0, c at 60 degrees
    theta = theta - 3.1415926/3;
    double g1 = 4.0*( 1.0 - e*e ) * cos(theta) * cos(theta) + pow(2.0*e-1.0, 2.0);
    double d1 = 2.0*( 1.0 - e*e ) * cos( theta );
    double d2 = ( 2.0*e-1.0 ) * pow( (4.0*(1.0-e*e)*cos(theta)*cos(theta) + 5.0*e*e - 4.0*e), 0.5);
    double temp =( d1 + d2 ) / g1; 

    return temp;
}

#endif

