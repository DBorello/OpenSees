// Inelastic2DYS02.cpp
//////////////////////////////////////////////////////////////////////

#include "Inelastic2DYS02.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Inelastic2DYS02::Inelastic2DYS02(int tag, double a, double e, double iz,
                                 int Nd1, int Nd2,
				 YieldSurface_BC *ysEnd1,
				 YieldSurface_BC *ysEnd2,
				 double pow_k_min, double pow_k_unload,
				 int rf_algo, bool islinear, double rho
				 )
  :InelasticYS2DGNL (tag, Nd1, Nd2, ysEnd1, ysEnd2,
		     rf_algo, islinear, rho),
  A(a), E(e), Iz(iz), powKmin(pow_k_min), powKunload(pow_k_unload),
  fpeak(0.0), peakRot(0.0),  powKmin_orig(pow_k_min),
  rotInElas(-1.0), forceInElas(-1.0), sumMi_2(0.0), sumMi_1(0.0),
  iFactor(1.0), iFactor_hist(1.0), IminFactor(1.0), IminFactor_hist(1.0),
  sumNatrDisp(0.0), sumNatrDisp_hist(0.0), loading(true), loading_hist(true),
  delPmax1(0.0), delPmax2(0.0)
{
  massDof = A*L*rho;
  massDof = massDof/2;
  
}

Inelastic2DYS02::~Inelastic2DYS02()
{
	//does nothing
}


int Inelastic2DYS02::commitState()
{
	sumMi_2 = fabs(eleForce_hist(2)) + fabs(eleForce_hist(5));
	this->InelasticYS2DGNL::commitState();
	sumMi_1 = fabs(eleForce_hist(2)) + fabs(eleForce_hist(5));

	double diff = sumMi_1 - sumMi_2;
	iFactor_hist = iFactor;
	 IminFactor_hist =  IminFactor;
	 sumNatrDisp_hist =  sumNatrDisp;
	 loading_hist = loading;

	double delp1 = ys1->hModel->getTrialPlasticStrains(0);
	double delp2 = ys2->hModel->getTrialPlasticStrains(0);
	 
	if(delPmax1 < delp1)
		delPmax1 = delp1;
		
	if(delPmax2 < delp2)
		delPmax2 = delp2;

	cerr   << delp1
	<< " " << delPmax1
	<< " " << delp2
	<< " " << delPmax2
	<< endl;
	
	        
//	cout << "diff = " << diff; (diff < 0) ? cout << ", unloading\n" : cout << ", loading\n";

	return 0;
}

// very simple element
// provides the elastic stiffness based on the average degraded
// Iz at each ends
void Inelastic2DYS02::getLocalStiff(Matrix &K)
{
    iFactor = 1.0;
	double mp = ys1->getCap(0);
	double normalzF = forceInElas;

Vector tmp(6);
	getTrialNaturalDisp(tmp);
	sumNatrDisp =  fabs(tmp(2)) + fabs(tmp(5));
//    double diff =  sumNatrDisp - sumNatrDisp_hist;
//    cout << "sumNd = " << sumNatrDisp << ", sumNd_hist = " << sumNatrDisp_hist << endl;
		
Vector df(6), testForce(6);
	  getLocalStiff(K, 1);
     addInternalGeomStiff(K);
     getIncrNaturalDisp(tmp);
     df = K*tmp;
     testForce = eleForce_hist + df;

    double sumtest = fabs(testForce(2)) + fabs(testForce(5));
    
// 	double diff = sumMi_1 - sumMi_2;
	double diff = sumtest - sumMi_1;
//	cout << "sumtest = " << sumtest << ", i-1 = " <<  sumMi_1 << endl;
//	cout << "test force = " << testForce;
//	cout << "Force hist = " << this->eleForce_hist;
/*
	if(end1Plastify_hist || end2Plastify_hist)
	{
		diff =  sumNatrDisp - sumNatrDisp_hist;
	}
*/

	int maxEnd = 2;
	if(fabs(testForce(5)) > fabs(testForce(2)))
		 maxEnd = 5;
	 
	if(sign(testForce(maxEnd)) !=  sign(eleForce_hist(maxEnd)))
		loading = loading_hist;
	else
	{
		(diff >= 0) ? loading = true : loading = false;
	}
//    if(loading)
//		cout << "Loading\n"; else cout << "unloading \n"; cin.get();
			
	// cout << "Base = " << rotInElas << endl;
	// calculated only once
	if( (end1Damage && rotInElas==-1.0) || (end2Damage && rotInElas==-1.0) )
	{
		cout << " Initialized \n";
		end1Damage ? cout << "e1 \n" : cout << "!e1 \n";
		end2Damage ? cout << "e2 \n" : cout << "!e2 \n";
		rotInElas = fabs(tmp(2)) + fabs(tmp(5));
		forceInElas = fabs(eleForce_hist(2)) + fabs(eleForce_hist(5));

		cout << "BaseRot = " << rotInElas << ", Base Force = " << forceInElas << endl; // cin.get();
		cout << "mp = " << mp << endl;

	}//initialize it now
 	
	if(rotInElas == -1)
	{
		if( fabs(eleForce_hist(2)*eleForce_hist(5)) < 1e-10) // single curv
			normalzF = mp;
		else
			normalzF = mp*2; //approx
	}

	double normalzR = rotInElas;

	if(normalzR ==-1.0)
		normalzR = 1.0;
	/*if(normalzF ==-1.0)
		normalzF = 1.0;*/


	double sumRot = fabs(tmp(2)) + fabs(tmp(5));
//	cout << "r_curr before = " << m_avg;
	sumRot = sumRot/normalzR;
//	cout << ", r_curr after = " << m_avg << endl;

	if(peakRot < sumRot)
		peakRot = sumRot;

//!!
double delP_ult = 0.03;
	peakRot = delPmax1/delP_ult;
				
//	cout << "Curr Rotation = " << m_avg << ", Peak = " << fpeak << ", Base = " << rotInElas << endl;

// 	if(end1Damage)
// 		cout << "End1 Damage \n";
// 	if(end2Damage)
// 		cout << "End2 Damage \n";
//
// 		if(end1Damage || end2Damage) cin.get();
	// assuming that the 2 surfaces are same and symmetric (for now)
	// lets get some pinching in the system

// peak-force based
	
//	cout << "normalzF = " << normalzF << ", Base Force = " << forceInElas << endl;

	// assuming double curvature
	double sumForce = fabs(eleForce_hist(2)) + fabs(eleForce_hist(5));
		// m_avg = m_avg/2; // fix this!!

		// non dimensionalize
		sumForce = sumForce/normalzF;
	if(fpeak < sumForce)
		fpeak = sumForce;

double offset = 0.0; //1.5;
        if(peakRot <= offset) // || end1Plastify_hist || end2Plastify_hist )
                iFactor = 1;
        else {
		// degrade Kun0 also
	  double I0factor = 1/pow(peakRot,powKunload);
	  I0factor = pow((1 - peakRot),powKunload);
	  /* double dy = 0.31;
	     double x = peakRot;
	     
	     double y = dy*(delP_ult - x)/delP_ult;
	     y = y + (1 - dy);
	     I0factor = */
	  double i0f = I0factor;
	  if(I0factor < 0.5)
	    I0factor = 0.5;
	  
	  //		if( peakRot < 2.0)
	  //			IminFactor = I0factor;
	  //		else
	  //			IminFactor = 1/pow(peakRot, powKmin);
	  
	  //		if(peakRot < 10.0)
	  //			IminFactor = I0factor;
	  //		else
	  //			IminFactor = pow(I0factor, powKmin);
	  
	  //		powKmin = exp(0.02*peakRot)*powKmin_orig;
	  //		IminFactor = pow(I0factor, pow(1/I0factor, powKmin));
	  //		IminFactor = 1/pow(peakRot, powKmin);
	  
	  IminFactor = 1/exp(powKmin*peakRot) + 0.05;
	  double ductility = 12;
	  //double x = peakRot - offset;
	  //		IminFactor = pow((1 - peakRot),powKmin)*i0f; //(1 - pow(x, powKmin)/pow(ductility,powKmin))*I0factor;
	  // IminFactor = 1/(exp(powKmin*x));
	  // IminFactor = 1/exp(powKmin*x)*I0factor;
	  //IminFactor = I0factor/pow(peakRot,powKmin);
	  
	  //		cerr << "powKmin = " << powKmin << ", peakRot = " << peakRot << ", Imin= " << IminFactor << endl;
	  
	  
	  if(IminFactor < 0.02)
	    IminFactor = 0.02;
	  
	  iFactor = (sumForce*(I0factor - IminFactor))/fpeak + IminFactor;
	  
	  if(loading)  // iFactor = IminFactor;
	    iFactor = (IminFactor + iFactor)/2;
	  
	  //iFactor = (sumForce*(1 - IminFactor))/fpeak + IminFactor;
	  // cout << "peakRot = " << peakRot << ", Imin = " << IminFactor << ", Ifactor = " << iFactor << endl;
	  // cout << "sumForce = " << sumForce << ", Peak-Force = " << fpeak << endl << endl;
	}
	
	//	cout << "Rot i = " << tmp(2) << ", Rot j = " << tmp(5) << endl;
	
	// assuming some degrading parameter independent of damage
	// just based on peak
	
	getLocalStiff(K, iFactor);
	
	/*	double iz = Iz*iFactor;
		
		double	EIbyL = E*iz/L;
		
		K(0, 1) = K(0, 2) = K(0, 4) = K(0, 5)=0;
		K(1, 0) = K(1, 3) =0;
		K(2, 0) = K(2, 3) =0;
		K(3, 1) = K(3, 2) = K(3, 4) = K(3, 5)=0;
		K(4, 0) = K(4, 3) =0;
		K(5, 0) = K(5, 3) =0;

		K(0,0) = K(3,3) = (A/iz)*(EIbyL);
		K(0,3) = K(3,0) = (-A/iz)*(EIbyL);
		K(1,1) = K(4,4) = (12/(L*L))*(EIbyL);
		K(1,4) = K(4,1) = (-12/(L*L))*(EIbyL);
		K(1,2) = K(2,1) = K(1,5) = K(5,1) = (6/L)*(EIbyL);
		K(2,4) = K(4,2) = K(4,5) = K(5,4) = (-6/L)*(EIbyL);
		K(2,2) = K(5,5) = 4*(EIbyL);
		K(2,5) = K(5,2) = 2*(EIbyL);
	*/
}//getLocalStiff


void Inelastic2DYS02::getLocalStiff(Matrix &K, double i_factor)
{
	double iz = Iz*i_factor;
    double	EIbyL = E*iz/L;

    K(0, 1) = K(0, 2) = K(0, 4) = K(0, 5)=0;
    K(1, 0) = K(1, 3) =0;
    K(2, 0) = K(2, 3) =0;
    K(3, 1) = K(3, 2) = K(3, 4) = K(3, 5)=0;
    K(4, 0) = K(4, 3) =0;
    K(5, 0) = K(5, 3) =0;

	K(0,0) = K(3,3) = (A/iz)*(EIbyL);
	K(0,3) = K(3,0) = (-A/iz)*(EIbyL);
	K(1,1) = K(4,4) = (12/(L*L))*(EIbyL);
	K(1,4) = K(4,1) = (-12/(L*L))*(EIbyL);
	K(1,2) = K(2,1) = K(1,5) = K(5,1) = (6/L)*(EIbyL);
	K(2,4) = K(4,2) = K(4,5) = K(5,4) = (-6/L)*(EIbyL);
	K(2,2) = K(5,5) = 4*(EIbyL);
	K(2,5) = K(5,2) = 2*(EIbyL);


}
