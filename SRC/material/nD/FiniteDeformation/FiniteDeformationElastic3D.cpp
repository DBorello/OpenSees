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
  Vector  d11 = W->d2isowOdlambda2(lambda_wave);
  d2.val(1,1) = d11(0);
  d2.val(2,2) = d11(1);
  d2.val(3,3) = d11(2);
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
  Tensor CinvCinv = Cinv("ij") * Cinv("kl") ;
  Cinv.null_indices(); CinvCinv.null_indices();

  Tensor ICinv = ( CinvCinv.transpose0110() + CinvCinv.transpose0111() ) * (0.5);

  Tensor CinvCinv_ICinv = CinvCinv - ICinv;

  double lambda1 = lambda(0);
  double lambda2 = lambda(1);
  double lambda3 = lambda(2);

  double I1 = lambda1*lambda1 + lambda2*lambda2 + lambda3*lambda3;

  int caseIndex = this->getCaseIndex();
  Tensor L_iso(2,def_dim_2,0.0);

  if(caseIndex == 0)
  {
    double d1 = 2.0*lambda1*lambda1*lambda1*lambda1
                      - I1*lambda1*lambda1 + J*J/(lambda1*lambda1);
    double d2 = 2.0*lambda2*lambda2*lambda2*lambda2
                      - I1*lambda2*lambda2 + J*J/(lambda2*lambda2);
    double d3 = 2.0*lambda3*lambda3*lambda3*lambda3
                      - I1*lambda3*lambda3 + J*J/(lambda3*lambda3);

    Tensor M1 = ( C - I_ij*(I1-lambda1*lambda1) + Cinv*(J*J/(lambda1*lambda1)) ) * (1.0/d1);
    Tensor M2 = ( C - I_ij*(I1-lambda2*lambda2) + Cinv*(J*J/(lambda2*lambda2)) ) * (1.0/d2);
    Tensor M3 = ( C - I_ij*(I1-lambda3*lambda3) + Cinv*(J*J/(lambda3*lambda3)) ) * (1.0/d3);

    double d1p = 4.0 *lambda1*lambda1*lambda1*lambda1 - I1*lambda1*lambda1
                   - J*J /(lambda1*lambda1);
    double d2p = 4.0 *lambda2*lambda2*lambda2*lambda2 - I1*lambda2*lambda2
                   - J*J /(lambda2*lambda2);
    double d3p = 4.0 *lambda3*lambda3*lambda3*lambda3 - I1*lambda3*lambda3
                   - J*J /(lambda3*lambda3);

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
    Tensor calM2 = ( tempI + (CinvCinv_ICinv -Cm1M2M2Cm1)*(J*J/(lambda2*lambda2)) + dM2M2d*(lambda2*lambda2) - M2M2*d2p ) *(1.0/d2);
    Tensor calM3 = ( tempI + (CinvCinv_ICinv -Cm1M3M3Cm1)*(J*J/(lambda3*lambda3)) + dM3M3d*(lambda3*lambda3) - M3M3*d3p ) *(1.0/d3);

    Vector Wa = this->wa();
    Tensor yab = this->Yab();
    Tensor L_iso_1 = ( calM1*Wa(0) + calM2*Wa(1) + calM3*Wa(2) ) * 2.0;
    Tensor L_iso_2 =  M1("ij") * M1("kl") * yab.cval(1,1)  + M1("ij") * M2("kl") * yab.cval(1,2)  + M1("ij") * M3("kl") * yab.cval(1,3)  +
                      M2("ij") * M1("kl") * yab.cval(2,1)  + M2("ij") * M2("kl") * yab.cval(2,2)  + M2("ij") * M3("kl") * yab.cval(2,3)  +
                      M3("ij") * M1("kl") * yab.cval(3,1)  + M3("ij") * M2("kl") * yab.cval(3,2)  + M3("ij") * M3("kl") * yab.cval(3,3);
    L_iso = L_iso_1 + L_iso_2 ;
  }

  if(caseIndex == 11)
  {
    double d1 = 2.0*lambda1*lambda1*lambda1*lambda1
                      - I1*lambda1*lambda1 + J*J/(lambda1*lambda1);
    Tensor M1 = ( C - I_ij*(I1-lambda1*lambda1) + Cinv*(J*J/(lambda1*lambda1)) ) * (1.0/d1);
    Tensor Mr = Cinv - M1;
    double d1p = 4.0 *lambda1*lambda1*lambda1*lambda1 - I1*lambda1*lambda1
                   - J*J /(lambda1*lambda1);
    Tensor Cm1M1M1Cm1 = Cinv("ij")*M1("kl") + M1("ij")*Cinv("kl");
    Cinv.null_indices(); M1.null_indices(); Cm1M1M1Cm1.null_indices();
    Tensor dM1M1d = I_ij("ij")*M1("kl") + M1("ij")*I_ij("kl");
    I_ij.null_indices(); M1.null_indices(); dM1M1d.null_indices();
    Tensor M1M1 = M1("ij") * M1("kl");
    M1.null_indices(); M1M1.null_indices();
    Tensor calM1 = ( tempI + (CinvCinv_ICinv -Cm1M1M1Cm1)*(J*J/(lambda1*lambda1)) + dM1M1d*(lambda1*lambda1) - M1M1*d1p ) *(1.0/d1);
    Tensor calMr = (ICinv + calM1) * (-1.0);
    Vector Wa = this->wa();
    Tensor yab = this->Yab();
    Tensor L_iso_1 = ( calM1*Wa(0) + calMr*Wa(1) ) * 2.0;
    Tensor L_iso_2 =  M1("ij") * M1("kl") * yab.cval(1,1)  + M1("ij") * Mr("kl") * yab.cval(1,2)  +
                      Mr("ij") * M1("kl") * yab.cval(2,1)  + Mr("ij") * Mr("kl") * yab.cval(2,2);
    L_iso = L_iso_1 + L_iso_2 ;
  }

  if(caseIndex == 12)
  {
    double d2 = 2.0*lambda2*lambda2*lambda2*lambda2
                      - I1*lambda2*lambda2 + J*J/(lambda2*lambda2);
    Tensor M2 = ( C - I_ij*(I1-lambda2*lambda2) + Cinv*(J*J/(lambda2*lambda2)) ) * (1.0/d2);
    Tensor Mr = Cinv - M2;
    double d2p = 4.0 *lambda2*lambda2*lambda2*lambda2 - I1*lambda2*lambda2;
    Tensor Cm1M2M2Cm1 = Cinv("ij")*M2("kl") + M2("ij")*Cinv("kl");
    Cinv.null_indices(); M2.null_indices(); Cm1M2M2Cm1.null_indices();
    Tensor dM2M2d = I_ij("ij")*M2("kl") + M2("ij")*I_ij("kl");
    I_ij.null_indices(); M2.null_indices(); dM2M2d.null_indices();
    Tensor M2M2 = M2("ij") * M2("kl");
    M2.null_indices(); M2M2.null_indices();
    Tensor calM2 = ( tempI + (CinvCinv_ICinv -Cm1M2M2Cm1)*(J*J/(lambda2*lambda2)) + dM2M2d*(lambda2*lambda2) - M2M2*d2p ) *(1.0/d2);
    Tensor calMr = (ICinv + calM2) * (-1.0);
    Vector Wa = this->wa();
    Tensor yab = this->Yab();
    Tensor L_iso_1 = ( calM2*Wa(1) + calMr*Wa(2) ) * 2.0;
    Tensor L_iso_2 =  M2("ij") * M2("kl") * yab.cval(2,2)  + M2("ij") * Mr("kl") * yab.cval(2,3)  +
                      Mr("ij") * M2("kl") * yab.cval(3,2)  + Mr("ij") * Mr("kl") * yab.cval(3,3);
    L_iso = L_iso_1 + L_iso_2 ;
  }

  if(caseIndex == 13)
  {
    double d3 = 2.0*lambda3*lambda3*lambda3*lambda3
                      - I1*lambda3*lambda3 + J*J/(lambda3*lambda3);
    Tensor M3 = ( C - I_ij*(I1-lambda3*lambda3) + Cinv*(J*J/(lambda3*lambda3)) ) * (1.0/d3);
    Tensor Mr = Cinv - M3;
    double d3p = 4.0 *lambda3*lambda3*lambda3*lambda3 - I1*lambda3*lambda3
                   - J*J /(lambda3*lambda3);
    Tensor Cm1M3M3Cm1 = Cinv("ij")*M3("kl") + M3("ij")*Cinv("kl");
    Cinv.null_indices(); M3.null_indices(); Cm1M3M3Cm1.null_indices();
    Tensor dM3M3d = I_ij("ij")*M3("kl") + M3("ij")*I_ij("kl");
    I_ij.null_indices(); M3.null_indices(); dM3M3d.null_indices();
    Tensor M3M3 = M3("ij") * M3("kl");
    M3.null_indices(); M3M3.null_indices();
    Tensor calM3 = ( tempI + (CinvCinv_ICinv -Cm1M3M3Cm1)*(J*J/(lambda3*lambda3)) + dM3M3d*(lambda3*lambda3) - M3M3*d3p ) *(1.0/d3);
    Tensor calMr = (ICinv + calM3) * (-1.0);
    Vector Wa = this->wa();
    Tensor yab = this->Yab();
    Tensor L_iso_1 = ( calM3*Wa(2) + calMr*Wa(0) ) * 2.0;
    Tensor L_iso_2 =  M3("ij") * M3("kl") * yab.cval(3,3)  + M3("ij") * Mr("kl") * yab.cval(3,1)  +
                      Mr("ij") * M3("kl") * yab.cval(1,3)  + Mr("ij") * Mr("kl") * yab.cval(1,1);
    L_iso = L_iso_1 + L_iso_2 ;
  }

    if(caseIndex == 2)
  {
    Vector lambda_wave(3);
    lambda_wave = this->getlambda_wave();
    Vector lambda(3);
    lambda = this->getlambda();
    Vector  d11 = W->d2isowOdlambda2(lambda_wave);
    Vector  d1 = W->disowOdlambda(lambda_wave);
    double G2linear = d11(0)*lambda_wave(0)*lambda_wave(0) + d1(0)*lambda_wave(0);
    if (G2linear == 0.0) G2linear = E / (1.0+nu);
    L_iso = ( ICinv - CinvCinv * (1.0/3.0) ) * G2linear;
  }

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
   return this->getTangentTensor();
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
   stresstensor thisSPKStress = this->getStressTensor();
   stresstensor thisFPKStress = thisSPKStress("ij") * (F.transpose11())("jk") ;
   return thisFPKStress;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor FiniteDeformationElastic3D::getCauchyStressTensor(void)
{
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

   double J = this->getJ( );

   Vector lambda(3);
   lambda = this->getlambda();
   Vector lambda_wave(3);
   lambda_wave = this->getlambda_wave();
   tensor C = this->getC( );
   Tensor Cinv = C.inverse( );
   Cinv.symmetrize11();
   double lambda1 = lambda(0);
   double lambda2 = lambda(1);
   double lambda3 = lambda(2);
   double I1 = lambda1*lambda1+lambda2*lambda2+lambda3*lambda3;

   Tensor isoPK2Stress(2, def_dim_2, 0.0);

   int caseIndex = this->getCaseIndex();

   if (caseIndex == 0)
   {
     double d1 = 2.0 * lambda1*lambda1*lambda1*lambda1
                      - I1 *  lambda1*lambda1 + J*J/(lambda1*lambda1);
     double d2 = 2.0 * lambda2*lambda2*lambda2*lambda2
                      - I1 *  lambda2*lambda2 + J*J/(lambda2*lambda2);
     double d3 = 2.0 * lambda3*lambda3*lambda3*lambda3
                      - I1 *  lambda3*lambda3 + J*J/(lambda3*lambda3);
     Tensor M1 = ( C - I_ij*(I1-lambda1*lambda1) + Cinv *(J*J/(lambda1*lambda1)) ) * (1.0/d1);
     Tensor M2 = ( C - I_ij*(I1-lambda2*lambda2) + Cinv *(J*J/(lambda2*lambda2)) ) * (1.0/d2);
     Tensor M3 = ( C - I_ij*(I1-lambda3*lambda3) + Cinv *(J*J/(lambda3*lambda3)) ) * (1.0/d3);
     Vector Wa = this->wa();
     isoPK2Stress = M1*Wa(0) + M2*Wa(1) + M3*Wa(2);
   }

   if (caseIndex == 11)
   {
     double d1 = 2.0 * lambda1*lambda1*lambda1*lambda1
                      - I1 *  lambda1*lambda1 + J*J/(lambda1*lambda1);
     Tensor M1 = ( C - I_ij*(I1-lambda1*lambda1) + Cinv *(J*J/(lambda1*lambda1)) ) * (1.0/d1);
     Tensor Mr = Cinv - M1;
     Vector Wa = this->wa();
     isoPK2Stress = Mr*Wa(2) + M1*Wa(0);
   }

   if (caseIndex == 12)
   {
     double d2 = 2.0 * lambda2*lambda2*lambda2*lambda2
                      - I1 *  lambda2*lambda2 + J*J/(lambda2*lambda2);
     Tensor M2 = ( C - I_ij*(I1-lambda2*lambda2) + Cinv *(J*J/(lambda2*lambda2)) ) * (1.0/d2);
     Tensor Mr = Cinv - M2;
     Vector Wa = this->wa();
     isoPK2Stress = Mr*Wa(0) + M2*Wa(1);
   }

   if (caseIndex == 13)
   {
     double d3 = 2.0 * lambda3*lambda3*lambda3*lambda3
                      - I1 *  lambda3*lambda3 + J*J/(lambda3*lambda3);
     Tensor M3 = ( C - I_ij*(I1-lambda3*lambda3) + Cinv *(J*J/(lambda3*lambda3)) ) * (1.0/d3);
     Tensor Mr = Cinv - M3;
     Vector Wa = this->wa();
     isoPK2Stress = Mr*Wa(1) + M3*Wa(2);
   }

   if (caseIndex == 2)
   {
     Tensor zeroT(2, def_dim_2, 0.0);
     isoPK2Stress = zeroT;
   }

   double dWdJ = W->dvolwOdJ(J);
   Tensor volPK2Stress = Cinv * J * dWdJ;

   thisPK2Stress = volPK2Stress + isoPK2Stress; // This is PK2Stress

   thisGreenStrain = (C - I_ij) * 0.5; // This is Green Strain

   Tensor L_iso = this->FDisoStiffness();
   Tensor L_vol = this->FDvolStiffness();
   Stiffness = L_iso + L_vol; // This is Langrangian Tangent Stiffness

   return 0;
}
//--------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::getCaseIndex()
{
   Vector lambda(3);
   lambda = this->getlambda();
   double lambda1 = lambda(0);
   double lambda2 = lambda(1);
   double lambda3 = lambda(2);
   double diff12 = fabs(lambda1-lambda2);
   double diff23 = fabs(lambda2-lambda3);
   double diff31 = fabs(lambda3-lambda1);

   double perturbation = pow( d_macheps(), (0.5) );
   static int caseIdx;

   if (diff12 >= perturbation && diff23 >= perturbation && diff31 >= perturbation )
  caseIdx = 0;
   if (diff12 < perturbation && (diff23 >= perturbation || diff31 >= perturbation) )
  caseIdx = 13;
   if (diff23 < perturbation && (diff31 >= perturbation || diff12 >= perturbation) )
  caseIdx = 11;
   if (diff31 < perturbation && (diff12 >= perturbation || diff23 >= perturbation) )
  caseIdx = 12;
   if (diff12 < perturbation &&  diff23 < perturbation && diff31 < perturbation )
  caseIdx = 2;

   if (caseIdx !=  0 &&
       caseIdx != 11 &&
       caseIdx != 12 &&
       caseIdx != 13 &&
       caseIdx !=  2 )
   opserr << "FiniteDeformationElastic3D::getCaseIndex -- unknown case! \n";

   return caseIdx;
}


