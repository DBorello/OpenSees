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
#include <math.h>

#include <Matrix.h>
#include <Vector.h>
#include <Tensor.h>
#include <BJtensor.h>
#include <BJvector.h>
#include <ID.h>

#include <Channel.h>
#include <OPS_Globals.h>

#include <W.h>

#include <FiniteDeformationElastic3D.h>


//-----------------------------------------------------------------------------------------------------------------------------------------------
FiniteDeformationElastic3D::FiniteDeformationElastic3D(int tag,
                                                       int classTag,
                                                       WEnergy *wEnergy_in,
                                                       double rho_in= 0.0)
:NDMaterial(tag, classTag), rho(rho_in)
{
    if ( wEnergy_in )
    {
       W = wEnergy_in->newObj();
       E = W->getE();
       nu = W->getnu();

       Tensor F0("I", 2, def_dim_2);
       F = F0;
    }
    else
    {
      opserr << "FiniteDeformationElastic3D:: FiniteDeformationElastic3D failed to construct the W Energy\n";
      exit(-1);
    }
}

FiniteDeformationElastic3D::FiniteDeformationElastic3D(int tag,
                                                       WEnergy *wEnergy_in,
                                                       double rho_in = 0.0)
:NDMaterial(tag, ND_TAG_FiniteDeformationElastic3D), rho(rho_in)
{
    if ( wEnergy_in)
    {
       W = wEnergy_in->newObj();
       E = W->getE();
       nu = W->getnu();

       Tensor F0("I", 2, def_dim_2);
       F = F0;
    }
    else
    {
      opserr << "FiniteDeformationElastic3D:: FiniteDeformationElastic3D failed to construct the W Energy\n";
      exit(-1);
    }
}

FiniteDeformationElastic3D::FiniteDeformationElastic3D(int tag,
                                                       WEnergy *wEnergy_in)
:NDMaterial(tag, ND_TAG_FiniteDeformationElastic3D), rho(0.0)
{
    if ( wEnergy_in )
    {
       W = wEnergy_in->newObj();
       E = W->getE();
       nu = W->getnu();

       Tensor F0("I", 2, def_dim_2);
       F = F0;

    }
    else
    {
     opserr << "FiniteDeformationElastic3D:: FiniteDeformationElastic3D failed to construct the W Energy\n";
      exit(-1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------------
FiniteDeformationElastic3D::FiniteDeformationElastic3D( )
:NDMaterial(0, 0), rho(0.0)
{
    W = 0;
    E = 0.0;
    nu = 0.0;

    Tensor F0("I", 2, def_dim_2);
    F = F0;
}
//------------------------------------------------------------------------------------------------------------------------------------------------
FiniteDeformationElastic3D::~FiniteDeformationElastic3D()
{
   if (W)
     delete W;
}
//---------------------------------------------------------------------------------------------------------------------------------------
double FiniteDeformationElastic3D::getE(void)
{
   return E;
}
//------------------------------------------------------------------------------------------------------------------------------------------------
double FiniteDeformationElastic3D::getnu(void)
{
   return nu;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------
double FiniteDeformationElastic3D::getRho(void)
{
   return rho;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------
WEnergy *FiniteDeformationElastic3D::getWEnergy(void)
{
   return W;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::setTrialF(const Tensor &f)
{
   F = f;
   int ret = 0;
   ret = this->ComputeTrials();
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::setTrialF(const Tensor &f, const Tensor &d)
{
   F = f;
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::setTrialFIncr(const Tensor &f)
{
   F += f;
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::setTrialFIncr(const Tensor &f, const Tensor &d)
{
   F += f;
   return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor FiniteDeformationElastic3D::getF(void)
{
  return F;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor FiniteDeformationElastic3D::getC(void)
{
   F = this->getF();
   tensor C = (F.transpose11())("ij")*F("jk");
   F.null_indices();
   (F.transpose11()).null_indices();
   C.null_indices();
   C.symmetrize11();

   return C;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
const double FiniteDeformationElastic3D::getJ(void)
{
   F = this->getF();
   double J0 = F.determinant( );
   return J0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
const Vector FiniteDeformationElastic3D::getlambda(void)
{
  tensor C = this->getC();
  tensor eigtensor = C.eigenvalues();
  Vector lambda(3);

  double lambda1 = sqrt(eigtensor.cval(1));
  double lambda2 = sqrt(eigtensor.cval(2));
  double lambda3 = sqrt(eigtensor.cval(3));

  double diff12 = fabs(lambda1-lambda2);
  double diff23 = fabs(lambda2-lambda3);
  double diff31 = fabs(lambda3-lambda1);

  double perturbation = pow( d_macheps(), (0.5) );

//  if (diff12 < perturbation ) lambda1 += lambda1*perturbation;
//  if (diff23 < perturbation ) lambda2 += lambda2*perturbation;
//  if (diff31 < perturbation ) lambda3 += lambda3*perturbation;

  if (diff12 < perturbation ) lambda1 += 0.0;
  if (diff23 < perturbation ) lambda2 += lambda2*perturbation;
  if (diff31 < perturbation ) lambda3 -= lambda3*perturbation;

  if (lambda1 == 0.0 ) lambda1 = perturbation;
  if (lambda2 == 0.0 ) lambda2 = perturbation;
  if (lambda3 == 0.0 ) lambda3 = perturbation;

  lambda(0) = lambda1;
  lambda(1) = lambda2;
  lambda(2) = lambda3;

  return lambda;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------
const Vector FiniteDeformationElastic3D::getlambda_wave(void)
{
  Vector lambda(3);
  lambda = this->getlambda();
  double J = this->getJ();
  double JJJ = pow(J, -0.33333333333333333333333333333);
  Vector lambda_wave(3);
  lambda_wave(0) = lambda(0) *JJJ;
  lambda_wave(1) = lambda(1) *JJJ;
  lambda_wave(2) = lambda(2) *JJJ;

  return lambda_wave;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
const Vector FiniteDeformationElastic3D::wa(void)
{
  Vector Wa(3);
  Vector lambda_wave(3);
  lambda_wave = this->getlambda_wave();
  Vector disowOverlambda = W->disowOdlambda(lambda_wave);
  double temp = disowOverlambda(0) * lambda_wave(0) +
                disowOverlambda(1) * lambda_wave(1) +
                disowOverlambda(2) * lambda_wave(2) ;
  temp = temp * (-0.3333333333333333333333333333);
  Wa(0) = temp + disowOverlambda(0) * lambda_wave(0);
  Wa(1) = temp + disowOverlambda(1) * lambda_wave(1);
  Wa(2) = temp + disowOverlambda(2) * lambda_wave(2);
  return Wa;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor FiniteDeformationElastic3D::Yab(void)
{
        Tensor Y(2, def_dim_2, 0.0);
        Tensor I_ij("I", 2, def_dim_2);
        Vector lambda_wave(3);
        lambda_wave = this->getlambda_wave();
        Tensor  d2 = W->d2isowOdlambda1dlambda2(lambda_wave);
        Vector  d1 = W->disowOdlambda(lambda_wave);
        Vector tempi(3);
        double tempd = d1(0)*lambda_wave(0) + d1(1)*lambda_wave(1) + d1(2)*lambda_wave(2) ;
        double tempcd = 0.0;
        for (int i=0; i<3; i++)
        {
          tempi(i) = 0.0;
          for (int j=0; j<3; j++)
          {
              tempi(i) += d2.cval(i+1,j+1) * lambda_wave(i) * lambda_wave(j);
              tempcd   += d2.cval(i+1,j+1) * lambda_wave(i) * lambda_wave(j);
          }
        }
        for(int a=1; a<=3; a++)
        {
          for(int b=1; b<=3; b++)
          {
              Y.val(a,b) = d1(a-1)*I_ij.cval(a,b)*lambda_wave(b-1) + d2.cval(a,b)*lambda_wave(a-1)*lambda_wave(b-1) -
                           (  tempi(a-1) + tempi(b-1) + d1(a-1)*lambda_wave(a-1) + d1(b-1)*lambda_wave(b-1) ) / 3.0 +
                           ( tempcd + tempd ) / 9.0;
          }
        }
        return Y;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor FiniteDeformationElastic3D::FDisoStiffness(void)
{
  Tensor I_ij("I", 2, def_dim_2);
  Tensor I_ijkl( 4, def_dim_4, 0.0 );
  I_ijkl = I_ij("ij") * I_ij("kl");
  I_ij.null_indices();
  Tensor I_ikjl( 4, def_dim_4, 0.0 );
  I_ikjl = I_ijkl.transpose0110();
  Tensor I_iljk( 4, def_dim_4, 0.0 );
  I_iljk = I_ijkl.transpose0111();
  Tensor I4s = (I_ikjl+I_iljk)*0.5;
  Tensor  tempI = I4s - I_ijkl;

  Vector lambda(3);
    lambda = this->getlambda();
  Vector lambda_wave(3);
    lambda_wave = this->getlambda_wave();
  tensor C = this->getC( );
  double J = this->getJ( );

  Tensor Cinv = C.inverse( );
  Cinv.symmetrize11();

  double lambda1 = lambda(0);
  double lambda2 = lambda(1);
  double lambda3 = lambda(2);

  double I1 = lambda1*lambda1 + lambda2*lambda2 + lambda3*lambda3;

  double d1 = 2.0*lambda1*lambda1*lambda1*lambda1
                      - I1*lambda1*lambda1 + J*J/(lambda1*lambda1);
  double d2 = 2.0*lambda2*lambda2*lambda2*lambda2
                      - I1*lambda2*lambda2 + J*J/(lambda2*lambda2);
  double d3 = 2.0*lambda3*lambda3*lambda3*lambda3
                      - I1*lambda3*lambda3 + J*J/(lambda3*lambda3);
//  printf("lambda1,2,3=%e,%e,%e\n",lambda1,lambda2,lambda3);
//  printf("d1,d2,d3=%e,%e,%e\n",d1,d2,d3);
//  double test100 =1.0/d1 + 1.0/d2 + 1.0/d3;
//  printf("1/d1+1/d2+1/d3=%lf\n",test100);

  Tensor M1 = ( C - I_ij*(I1-lambda1*lambda1) + Cinv*(J*J/(lambda1*lambda1)) ) * (1.0/d1);
  Tensor M2 = ( C - I_ij*(I1-lambda2*lambda2) + Cinv*(J*J/(lambda2*lambda2)) ) * (1.0/d2);
  Tensor M3 = ( C - I_ij*(I1-lambda3*lambda3) + Cinv*(J*J/(lambda3*lambda3)) ) * (1.0/d3);

  double d1p = 4.0 *lambda1*lambda1*lambda1*lambda1 - I1*lambda1*lambda1
                   - J*J /(lambda1*lambda1);
  double d2p = 4.0 *lambda2*lambda2*lambda2*lambda2 - I1*lambda2*lambda2
                   - J*J /(lambda2*lambda2);
  double d3p = 4.0 *lambda3*lambda3*lambda3*lambda3 - I1*lambda3*lambda3
                   - J*J /(lambda3*lambda3);

  Tensor CinvCinv = Cinv("ij") * Cinv("kl") ;
  Cinv.null_indices(); CinvCinv.null_indices();

  Tensor ICinv = ( CinvCinv.transpose0110() + CinvCinv.transpose0111() ) * (0.5);

  Tensor CinvCinv_ICinv = CinvCinv - ICinv;

  Tensor Cm1M1M1Cm1 = Cinv("ij")*M1("kl") + M1("ij")*Cinv("kl");
  Cinv.null_indices(); M1.null_indices(); Cm1M1M1Cm1.null_indices();

  Tensor Cm1M2M2Cm1 = Cinv("ij")*M2("kl") + M2("ij")*Cinv("kl");
  Cinv.null_indices(); M2.null_indices(); Cm1M2M2Cm1.null_indices();

  Tensor Cm1M3M3Cm1 = Cinv("ij")*M3("kl") + M3("ij")*Cinv("kl");
  Cinv.null_indices(); M3.null_indices(); Cm1M3M3Cm1.null_indices();


  Tensor dM1M1d = I_ij("ij")*M1("kl") + M1("ij")*I_ij("kl");
  I_ij.null_indices(); M1.null_indices(); dM1M1d.null_indices();

  Tensor dM2M2d = I_ij("ij")*M2("kl") + M2("ij")*I_ij("kl");
  I_ij.null_indices(); M2.null_indices(); dM2M2d.null_indices();

  Tensor dM3M3d = I_ij("ij")*M3("kl") + M3("ij")*I_ij("kl");
  I_ij.null_indices(); M3.null_indices(); dM3M3d.null_indices();

  Tensor M1M1 = M1("ij") * M1("kl");
  M1.null_indices(); M1M1.null_indices();
  Tensor M2M2 = M2("ij") * M2("kl");
  M2.null_indices(); M2M2.null_indices();
  Tensor M3M3 = M3("ij") * M3("kl");
  M3.null_indices(); M3M3.null_indices();

  Tensor calM1 = ( tempI + (CinvCinv_ICinv -Cm1M1M1Cm1)*(J*J/(lambda1*lambda1)) + dM1M1d*(lambda1*lambda1) - M1M1*d1p ) *(1.0/d1);
//  calM1.print("M1","\n");
  Tensor calM2 = ( tempI + (CinvCinv_ICinv -Cm1M2M2Cm1)*(J*J/(lambda2*lambda2)) + dM2M2d*(lambda2*lambda2) - M2M2*d2p ) *(1.0/d2);
//  calM2.print("M2","\n");
  Tensor calM3 = ( tempI + (CinvCinv_ICinv -Cm1M3M3Cm1)*(J*J/(lambda3*lambda3)) + dM3M3d*(lambda3*lambda3) - M3M3*d3p ) *(1.0/d3);
//  calM3.print("M3","\n");


//---------------------------------------------------------------------------
// T           E           S           T          I          N          G
//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//  printf("\n  T E S T I N G     (M1, M2, M3) \n");
//
//---------------------------------------------------------------------------
//  tensor temp2 = M1 + M2 + M3;
//  temp2.print("M1pM2pM3","\nM1pM2pM3");
//  Cinv.print("Cinv","\nCinv");
//  if ( temp2 == Cinv )
//  {
//    printf(" M1+M2+M3 == Cinv  TRUE (OK) ( up to sqrt(macheps())) \n");
//  }
//  else
//  {
//    printf("\a M1+M2+M3 == Cinv    NOTTRUE  ( up to sqrt(macheps())) \n");
//    printf("\a M1+M2+M3 == Cinv    NOTTRUE  ( up to sqrt(macheps())) \n");
//    temp2.print("temp2","\n temp2");
//    Cinv.print("Cinv","\n Cinv");
//  }
//---------------------------------------------------------------------------
//  tensor temp3 = M1*lambda1*lambda1 + M2*lambda2*lambda2 + M3*lambda3*lambda3;
//  temp3.print("temp3","\ntemp3");
//  I_ij.print("KroneckerI","\nKroneckerI");
//  if ( temp3 == I_ij )
//  {
//    printf(" temp3 == KroneckerI  TRUE (OK) ( up to sqrt(macheps())) \n");
//  }
//  else
//  {
//    printf("\a  temp3 == I_ij    NOTTRUE  ( up to sqrt(macheps())) \n");
//    printf("\a  temp3 == I_ij    NOTTRUE  ( up to sqrt(macheps())) \n");
//    temp3.print("temp3","\n temp3");
//    I_ij.print("I","\n KroneckerI");
//  }
//test//---------------------------------------------------------------------------
//  tensor temp4 =   M1*lambda1*lambda1*lambda1*lambda1
//                 + M2*lambda2*lambda2*lambda2*lambda2
//                 + M3*lambda3*lambda3*lambda3*lambda3;
//  temp4.print("temp4","\ntemp4");
//  C.print("C","\nC");
//  if ( temp4 == C )
//  {
//    printf(" temp4 == C  TRUE (OK) ( up to sqrt(macheps())) \n");
//  }
//  else
//  {
//    printf("\a  temp4 == C    NOTTRUE  ( up to sqrt(macheps())) \n");
//    printf("\a  temp4 == C    NOTTRUE  ( up to sqrt(macheps())) \n");
//    temp4.print("temp4","\n temp4");
//    C.print("C","\n C");
//  }
//
//
//---------------------------------------------------------------------------
//  printf("\n  T E S T I N G     (calM1, calM2, calM3) \n");
//
//
//---------------------------------------------------------------------------
//  tensor calMsum = calM1*-1. + calM2*-1. + calM3*-1.;
//  if ( calMsum == ICinv )
//  {
//    printf(" calMsum == ICinv  TRUE (OK) ( up to sqrt(macheps())) \n");
//  }
//  else
//  {
//    printf("\a calMsum == ICinv    NOTTRUE  ( up to sqrt(macheps())) \n");
//    printf("\a calMsum == ICinv    NOTTRUE  ( up to sqrt(macheps())) \n");
//    calMsum.print("calMsum","\n calMsum");
//    ICinv.print("ICinv","\n ICinv");
//  }
//
//
//---------------------------------------------------------------------------
//  tensor ZERO(4,def_dim_4,0.0);
//  ZERO.print("z","\ntensor ZERO (4-th order with value assignment)");
//  tensor temp5 =   calM1*(lambda1*lambda1)
//                 + calM2*(lambda2*lambda2)
//                 + calM3*(lambda3*lambda3)
//                 + M1("ij")*M1("kl")*(lambda1*lambda1)
//                 + M2("ij")*M2("kl")*(lambda2*lambda2)
//                 + M3("ij")*M3("kl")*(lambda3*lambda3);
//  if ( temp5 == ZERO )
//  {
//    printf(" temp5 == ZERO     TRUE (OK) ( up to sqrt(macheps())) \n");
//  }
//  else
//  {
//    printf("\a temp5 == ZERO    NOTTRUE  ( up to sqrt(macheps())) \n");
//    printf("\a temp5 == ZERO    NOTTRUE  ( up to sqrt(macheps())) \n");
//    temp5.print("temp5","\n temp5");
//  }
//    temp5.print("temp5","\ntensor temp5 (should be zero!)");
//
//---------------------------------------------------------------------------
//  tensor temp6 =   calM1*(lambda1*lambda1*lambda1*lambda1)
//                 + calM2*(lambda2*lambda2*lambda2*lambda2)
//                 + calM3*(lambda3*lambda3*lambda3*lambda3)
//                 + M1("ij")*M1("kl")*2.0*(lambda1*lambda1*lambda1*lambda1)
//                 + M2("ij")*M2("kl")*2.0*(lambda2*lambda2*lambda2*lambda2)
//                 + M3("ij")*M3("kl")*2.0*(lambda3*lambda3*lambda3*lambda3);
//  if ( temp6 == I4s )
//  {
//    printf(" temp6 == I4s     TRUE (OK) ( up to sqrt(macheps())) \n");
//  }
//  else
//  {
//    printf("\a temp6 == I4s     NOTTRUE  ( up to sqrt(macheps())) \n");
//    printf("\a temp6 == I4s     NOTTRUE  ( up to sqrt(macheps())) \n");
//    temp6.print("temp6","\n temp6");
//    I4s.print("I4s","\n I4s");
//  }
//

  Vector Wa = this->wa();
//  printf("wa1,wa2,wa3=%lf,%lf,%lf\n",Wa(0),Wa(1),Wa(2));
  Tensor yab = this->Yab();
//  yab.print("yab","\n");
  Tensor L_iso_1 = ( calM1*Wa(0) + calM2*Wa(1) + calM3*Wa(2) ) * 2.0;
//  L_iso_1.print("Li1","\n");
  Tensor L_iso_2 =  M1("ij") * M1("kl") * yab.cval(1,1)  + M1("ij") * M2("kl") * yab.cval(1,2)  + M1("ij") * M3("kl") * yab.cval(1,3)  +
                    M2("ij") * M1("kl") * yab.cval(2,1)  + M2("ij") * M2("kl") * yab.cval(2,2)  + M2("ij") * M3("kl") * yab.cval(2,3)  +
                    M3("ij") * M1("kl") * yab.cval(3,1)  + M3("ij") * M2("kl") * yab.cval(3,2)  + M3("ij") * M3("kl") * yab.cval(3,3);
//  L_iso_2.print("Li2","\n");
  Tensor L_iso = L_iso_1 + L_iso_2 ;

  return L_iso;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor FiniteDeformationElastic3D::FDvolStiffness(void)
{
   tensor C = this->getC();
   Tensor Cinv = C.inverse( );
   Cinv.symmetrize11();
   Tensor CinvCinv = Cinv("ij")*Cinv("kl") ;
   Cinv.null_indices(); CinvCinv.null_indices();
   Tensor ICinv = ( CinvCinv.transpose0110() + CinvCinv.transpose0111() ) * (0.5);
   double J = this->getJ();
   double dWdJ = W->dvolwOdJ(J);
   double d2WdJ2 = W->d2volwOdJ2(J);
   double wj = d2WdJ2*J*J + J*dWdJ;

   Tensor L_vol = CinvCinv*wj - ICinv *2.0*J*dWdJ  ;

   return L_vol;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor& FiniteDeformationElastic3D::getTangentTensor(void)
{
    return Stiffness;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor
&FiniteDeformationElastic3D::getInitialTangentTensor(void)
{
    tensor I2("I", 2, def_dim_2);
    tensor I_ijkl = I2("ij")*I2("kl");
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    static tensor L0;
    L0 = I_ijkl*( E*nu / ( (1.0+nu)*(1.0 - 2.0*nu) ) ) + I4s*( E / (1.0 + nu) );

    return L0;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor FiniteDeformationElastic3D::getStrainTensor(void)
{
   return thisGreenStrain;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor FiniteDeformationElastic3D::getStressTensor(void)
{
   return thisPK2Stress;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor FiniteDeformationElastic3D::getPK1StressTensor(void)
{
//   F = this->getF();
   stresstensor thisSPKStress = this->getStressTensor();
   stresstensor thisFPKStress = thisSPKStress("ij") * (F.transpose11())("jk") ;
   return thisFPKStress;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor FiniteDeformationElastic3D::getCauchyStressTensor(void)
{
 //  F = this->getF();
   double J = this->getJ();
   stresstensor thisSPKStress = this->getStressTensor();
   stresstensor thisCauchyStress = F("ij") * thisSPKStress("jk") * (F.transpose11())("kl") * (1.0/J);
   return thisCauchyStress;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::commitState (void)
{
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::revertToLastCommit (void)
{
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::revertToStart (void)
{
   Tensor F0("I", 2, def_dim_2);
   F = F0;

   Tensor ss_zero(2,def_dim_2,0.0);
   thisPK2Stress = ss_zero;
   thisGreenStrain = ss_zero;
   Stiffness = getInitialTangentTensor();

   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NDMaterial * FiniteDeformationElastic3D::getCopy (void)
{
    FiniteDeformationElastic3D   *theCopy =
    new FiniteDeformationElastic3D (this->getTag(), this->getWEnergy(), this->getRho());

    theCopy->F = F;
    theCopy->Stiffness = Stiffness;
    theCopy->thisGreenStrain = thisGreenStrain;
    theCopy->thisPK2Stress = thisPK2Stress;

    return theCopy;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NDMaterial * FiniteDeformationElastic3D::getCopy (const char *type)
{

  opserr << "FiniteDeformationElastic3D::getCopy(const char *) - not yet implemented\n";

    return 0;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const char* FiniteDeformationElastic3D::getType (void) const
{
   return "FiniteDeformationElastic3D";
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::getOrder (void) const
{
   return 6;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::sendSelf (int commitTag, Channel &theChannel)
{
   int res = 0;

   static Vector data(4);

   data(0) = this->getTag();
   data(1) = rho;
   data(2) = E;
   data(3) = nu;

   res += theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0)
    {
      opserr << "FiniteDeformationElastic3D::sendSelf -- could not send Vector\n";
      return res;
    }

   return res;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::recvSelf (int commitTag,
                                          Channel &theChannel,
                                          FEM_ObjectBroker &theBroker)
{
   int res = 0;

   static Vector data(4);

   res += theChannel.recvVector(this->getDbTag(), commitTag, data);
   if (res < 0)
    {
      opserr << "FiniteDeformationElastic3D::recvSelf -- could not recv Vector\n";
      return res;
    }

   this->setTag((int)data(0));
   rho = data(1);
   E = data(2);
   nu = data(3);

   return res;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void FiniteDeformationElastic3D::Print (OPS_Stream &s, int flag)
{
   s << "Finite Deformation Elastic 3D model" << endln;
   s << "\trho: " << rho << endln;
   s << "\tE: " << E << endln;
   s << "\tnu: " << nu << endln;
   return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::setParameter(char **argv, int argc, Information &info)
{
   return -1;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::updateParameter(int parameterID, Information &info)
{
   return -1;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::ComputeTrials()
{
   Tensor I_ij("I", 2, def_dim_2);
   Tensor ss_zero(2,def_dim_2,0.0);

   if (F == I_ij)
   {
    thisPK2Stress = ss_zero;
    thisGreenStrain = ss_zero;
    Stiffness = getInitialTangentTensor();
    return 1;
   }

   Vector lambda(3);
   lambda = this->getlambda();
   Vector lambda_wave(3);
   lambda_wave = this->getlambda_wave();
   tensor C = this->getC( );
   double J = this->getJ( );
   Tensor Cinv = C.inverse( );
   Cinv.symmetrize11();
   double lambda1 = lambda(0);
   double lambda2 = lambda(1);
   double lambda3 = lambda(2);
   double I1 = lambda1*lambda1+lambda2*lambda2+lambda3*lambda3;
   double d1 = 2.0 * lambda1*lambda1*lambda1*lambda1
                      - I1 *  lambda1*lambda1 + J*J/(lambda1*lambda1);
   double d2 = 2.0 * lambda2*lambda2*lambda2*lambda2
                      - I1 *  lambda2*lambda2 + J*J/(lambda2*lambda2);
   double d3 = 2.0 * lambda3*lambda3*lambda3*lambda3
                      - I1 *  lambda3*lambda3 + J*J/(lambda3*lambda3);
   Tensor M1 = ( C - I_ij*(I1-lambda1*lambda1) + Cinv *(J*J/(lambda1*lambda1)) ) * (1.0/d1);
   Tensor M2 = ( C - I_ij*(I1-lambda2*lambda2) + Cinv *(J*J/(lambda2*lambda2)) ) * (1.0/d2);
   Tensor M3 = ( C - I_ij*(I1-lambda3*lambda3) + Cinv *(J*J/(lambda3*lambda3)) ) * (1.0/d3);
   double dWdJ = W->dvolwOdJ(J);
   Vector Wa = this->wa();
   Tensor volPK2Stress = Cinv * J * dWdJ;
   Tensor isoPK2Stress = M1*Wa(0) + M2*Wa(1) + M3*Wa(2);
   thisPK2Stress = volPK2Stress + isoPK2Stress;   // This is PK2Stress

   thisGreenStrain = (C - I_ij) * 0.5; // This is Green Strain

   Tensor L_iso = this->FDisoStiffness();
//   L_iso.print("Li","\n");
   Tensor L_vol = this->FDvolStiffness();
//   L_vol.print("Lv","\n");
   Stiffness = L_iso + L_vol;

   return -1;
}



