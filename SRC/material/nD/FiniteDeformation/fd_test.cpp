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

// fd_test.cpp

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <iostream.h>

#include <Vector.h>
#include <Tensor.h>
#include <BJvector.h>
#include <BJtensor.h>

#include "FiniteDeformationElastic3D.h"


// standard C++ includes
#include <stdlib.h>
#include <iostream.h>

#include <G3Globals.h>
#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>



// init the global variabled defined in G3Globals.h
ErrorHandler *g3ErrorHandler =0;
double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;



// init the global variabled defined in OPS_Globals.h
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream &opserr = sserr;




#include "FiniteDeformationElastic3D.h"
#include "W.h"
#include "LogWEnergy.h"
#include "NeoHookeanWEnergy.h"
#include "SimoPisterWEnergy.h"


int main()
{
 printf("\n\n\n*** Finite Deformations: T E S T ***\n\n");

double rho_in = 1000.0;
double E_in = 206.0e9;
double nu_in = 0.33;
printf("   rho = %.12e\n",rho_in);
printf("   Ey = %.12e\n",E_in);
printf("   nu = %.12e\n",nu_in);
//LogWEnergy  thisFDW( E_in, nu_in );
NeoHookeanWEnergy  thisFDW( E_in, nu_in );
//SimoPisterWEnergy  thisFDW( E_in, nu_in );

FiniteDeformationElastic3D  thisFDstate( 0, 0, &thisFDW, rho_in);

double gammastart = 0.0;
double gammaend = 1.001;
double deltagamma = 0.1;
printf("   gammastart  = %.12e\n",gammastart);
printf("   gammaend  = %.12e\n",gammaend);
printf("   deltagamma  = %.12e\n",deltagamma);

double gamma = 0.0;
Tensor aStress(2, def_dim_2, 0.0);

for ( gamma = gammastart ; gamma <= gammaend ; gamma = gamma + deltagamma )
{
	printf("\n   gamma  = %.12e\n",gamma);

// /***********************************************************************************/
// /*   Simple shear                                                                  */
// /***********************************************************************************/
	double F11 = 1.0;    double F12 = gamma;  double F13 = 0.0;
   	double F21 = 0.0;    double F22 = 1.0;    double F23 = 0.0;
	double F31 = 0.0;    double F32 = 0.0;    double F33 = 1.0;

//
//**************************************************************************************/
//*   Pure Extension                                                                   */
//**************************************************************************************/
// 	 double F11 = 1.0+gamma;  double F12 = 0.0;        double F13 = 0.0;
// 	 double F21 = 0.0;        double F22 = 1.00001;     double F23 = 0.0;
// 	 double F31 = 0.0;        double F32 = 0.0;        double F33 =1.00002;


//
///************************************************************************************/
///*   Pure compression                                                               */
///************************************************************************************/
//  double F11 = 1.0-gamma;     double F12 = 0.0;              double F13 = 0.0;
//  double F21 = 0.0;           double F22 = 0.999999;          double F23 = 0.0;
//  double F31 = 0.0;           double F32 = 0.0;              double F33 = 1.000001;


//
///************************************************************************************/
///*   Triaxial compression                                                           */
///************************************************************************************/
//  double dilatancy = 0.2;
//  double F11 = 1.0-gamma; double F12 = 0.0;                     double F13 = 0.0;
//  double F21 = 0.0;       double F22 = 0.99999+gamma*dilatancy; double F23 = 0.0;
//  double F31 = 0.0;       double F32 = 0.0;                     double F33 = 1.0+gamma*dilatancy;


  	double F_values[] = {  F11,  F12,  F13,
                               F21,  F22,  F23,
                               F31,  F32,  F33 };

 	      Tensor thisf(2, def_dim_2, F_values);
              thisf.print("F","\nDeformation Gradient:");
              thisFDstate.setTrialF(thisf);
              Tensor thisC = thisFDstate.getC();
              thisC.print("C","\nTensor of C");
              double thisJ = thisFDstate.getJ();
              printf("\nJ = %lf\n", thisJ);
//              Tensor thisStrain = thisFDstate.getStrainTensor();
//              Tensor thisStress = thisFDstate.getStressTensor();
//              thisStrain.print("Strain","\nGreen Strain");
//              thisStress.print("Stress","\nThe 2nd PK Stress");
//              Vector thislambda =  thisFDstate.getlambda();
//              printf("\nlambda 1 = %lf, lambda 2 = %lf, lambda 3 = %lf\n", thislambda(0), thislambda(1), thislambda(2));
//              Vector thislambda_wave =  thisFDstate.getlambda_wave();
//              printf("\nlambda_w 1 = %lf, lambda_w 2 = %lf, lambda_w 3 = %lf\n", thislambda_wave(0), thislambda_wave(1), thislambda_wave(2));
//              Vector thisWa = thisFDstate.wa();
//              printf("\nWa 1 = %lf, Wa 2 = %lf, Wa 3 = %lf\n", thisWa(0), thisWa(1), thisWa(2));
//              Tensor thisyab = thisFDstate.Yab();
//              thisyab.print("yab","Yab");
//              Tensor thisisoStiffTensor = thisFDstate.FDisoStiffness();
//              thisisoStiffTensor.print("Kiso","\nISO Tangent Tensor");

//              Tensor thisStiffTensor = thisFDstate.getTangentTensor();

//              thisStiffTensor.print("K","\nTangent Tensor");
//              Tensor thisvolStiffTensor = thisFDstate.FDvolStiffness();
//              thisvolStiffTensor.print("K","\nTangent Tensor");
//              Tensor thisStrain = thisFDstate.getStrainTensor();
//              thisStrain.print("E","\nGreen Strain Tensor:");
//              Tensor thisPK2Stress = thisFDstate.getPK2StressTensor();
//              thisPK2Stress.print("S","\n2nd PK Stress Tensor:");
//              Tensor thisFPKStress = thisFDstate.getPK1StressTensor();
//              thisPK1Stress.print("P","\n1st PK Stress Tensor:");
              stresstensor tStress = thisFDstate.getCauchyStressTensor();
              tStress.print("S2","\n");
}
        return 1;
}


