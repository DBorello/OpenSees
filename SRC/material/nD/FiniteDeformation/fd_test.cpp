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
#include "OgdenWEnergy.h"
#include "MooneyRivlinWEnergy.h"
#include "OgdenSimoWEnergy.h"
#include "MooneyRivlinSimoWEnergy.h"

int main()
{
 printf("\n\n\n*** Finite Deformations: T E S T ***\n\n");

double rho_in = 0.0;
double E_in = 11.83e5;
double nu_in = 0.4;

////********** for Odegen
//int Nogden = 3;
//double cr[3];
//double mur[3];
//cr[0] = 6.3e5;
//cr[1] = 0.012e5;
//cr[2] = -0.1e5;
//mur[0] = 1.3;
//mur[1] = 5.0;
//mur[2] = -2.0;

////********** for Mooney-Rivlin
//double c1_in = 1.8484375e5;
//double c2_in =  0.2640625e5;

printf("   rho = %.12e\n",rho_in);
printf("   Ey = %.12e\n",E_in);
printf("   nu = %.12e\n",nu_in);

//LogWEnergy  thisFDW( E_in, nu_in );
NeoHookeanWEnergy  thisFDW( E_in, nu_in );
//SimoPisterWEnergy  thisFDW( E_in, nu_in );
//OgdenWEnergy  thisFDW( E_in, nu_in, Nogden, cr, mur);
//MooneyRivlinWEnergy  thisFDW( E_in, nu_in, c1_in, c2_in );
//OgdenSimoWEnergy  thisFDW( E_in, nu_in, Nogden, cr, mur );
//MooneyRivlinSimoWEnergy  thisFDW( E_in, nu_in, c1_in, c2_in );


FiniteDeformationElastic3D  thisFDstate( 0, 0, &thisFDW, rho_in);

double gammastart = 0.0;
double gammaend = 0.0011;
double deltagamma = 0.001;
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
 double F11 = 1.0;    double F12 = 0.0;    double F13 = gamma;
 double F21 = 0.0;    double F22 = 1.0;    double F23 = 0.0;
 double F31 = 0.0;    double F32 = 0.0;    double F33 = 1.0;

//
//**************************************************************************************/
//*   Pure Extension                                                                   */
//**************************************************************************************/
//   double F11 = 1.0+gamma;  double F12 = 0.0;        double F13 = 0.0;
//   double F21 = 0.0;        double F22 = 1.0;        double F23 = 0.0;
//   double F31 = 0.0;        double F32 = 0.0;        double F33 = 1.0;

//
//**************************************************************************************/
//*   Simple Extension                                                                   */
//**************************************************************************************/
//   double F11 = 1.0+gamma;  double F12 = 0.0;        double F13 = 0.0;
//   double F21 = 0.0;        double F22 = 1.0/sqrt(1.0+gamma);        double F23 = 0.0;
//   double F31 = 0.0;        double F32 = 0.0;        double F33 = 1.0/sqrt(1.0+gamma);

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
//  double F21 = 0.0;       double F22 = 1.0+gamma*dilatancy; double F23 = 0.0;
//  double F31 = 0.0;       double F32 = 0.0;                     double F33 = 1.0+gamma*dilatancy;


///************************************************************************************/
///*   Pure Volumn Change                                                           */
///************************************************************************************/
//  double F11 = 1.0+gamma;  double F12 = 0.0;                   double F13 = 0.0;
//  double F21 = 0.0;        double F22 = 1.0+gamma;             double F23 = 0.0;
//  double F31 = 0.0;        double F32 = 0.0;                   double F33 = 1.0+gamma;


    double F_values[] = {  F11,  F12,  F13,
                               F21,  F22,  F23,
                               F31,  F32,  F33 };

         Tensor thisf(2, def_dim_2, F_values);
              thisf.print("F","\n");
              thisFDstate.setTrialF(thisf);
//              Tensor thisC = thisFDstate.getC();
//              thisC.print("C","\n");
//              double thisJ = thisFDstate.getJ();
//              printf("\nJ = %lf\n", thisJ);


              Tensor thisPK2Stress = thisFDstate.getStressTensor();
              thisPK2Stress.print("PK2","\n2nd PK Stress Tensor:");
//              Tensor thisFPKStress = thisFDstate.getPK1StressTensor();
//              thisPK1Stress.print("P","\n1st PK Stress Tensor:");
              stresstensor tStress = thisFDstate.getCauchyStressTensor();
              tStress.print("Sigma","Cauchy Stress tensor:");
}
        return 1;
}

