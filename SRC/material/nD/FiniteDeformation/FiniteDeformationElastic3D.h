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

#ifndef FiniteDeformationElastic3D_h
#define FiniteDeformationElastic3D_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <Tensor.h>
#include <stresst.h>
#include <straint.h>

#include <ID.h>
#include <Channel.h>
#include <G3Globals.h>


#include <stdlib.h>
#include <iostream.h>

#include <ConsoleErrorHandler.h>


#include <W.h>  // for W Strain Energy Functions

class FiniteDeformationElastic3D : public NDMaterial
{
  public:
    FiniteDeformationElastic3D(int tag, int classTag, WEnergy * , double );
    FiniteDeformationElastic3D(int tag, WEnergy * , double );
    FiniteDeformationElastic3D(int tag, WEnergy * );
    FiniteDeformationElastic3D();
    virtual ~FiniteDeformationElastic3D();

    virtual double getRho(void);
//    virtual double getE(void);
//    virtual double getnu(void);

    virtual int setTrialF(const Tensor &f);
    virtual int setTrialF(const Tensor &f, const Tensor &d);
    virtual int setTrialFIncr(const Tensor &f);
    virtual int setTrialFIncr(const Tensor &f, const Tensor &d);

    virtual const Tensor& getTangentTensor(void) ;
//    virtual const Tensor& getInitailTangentTensor(void) ;

    virtual const  straintensor getStrainTensor(void) ;  // Default Green Strain
    virtual const  stresstensor getStressTensor(void) ; // Default 2nd Piola Kirchhoff Stress

//    virtual const Vector &getStress(void);
//    virtual const Vector &getStrain(void);

//    virtual const stresstensor getCommittedStress(void);
//    virtual const straintensor getCommittedStrain(void);

//    virtual const straintensor getPlasticStrainTensor(void);

    virtual int commitState(void) ;
    virtual int revertToLastCommit(void) ;
    virtual int revertToStart(void) ;

    virtual FiniteDeformationElastic3D *getCopy (void);
    virtual FiniteDeformationElastic3D *getCopy (const char *type){};

    virtual const char *getType (void) const;
    virtual int getOrder (void) const;

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel,
     FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    virtual int setParameter(char **argv, int argc, Information &info);
    virtual int updateParameter(int parameterID, Information &info);

  private:

    WEnergy *getWEnergy(void);

    const Tensor getF(void);
    const Tensor getC(void);
    const double getJ(void) ;
    const Vector getlambda(void) ;
    const Vector getlambda_wave(void) ;

    const Vector wa(void) ;
    const Tensor Yab(void) ;
    const Tensor FDisoStiffness(void) ;
    const Tensor FDvolStiffness(void) ;

    const  stresstensor getPK1StressTensor(void) ;
    const  stresstensor getCauchyStressTensor(void) ;

  protected:

     WEnergy * W;

     double rho;
     double J;

     Tensor F;
     Tensor C;

     Vector lambda;
     Vector lambda_wave;
     Tensor Stiffness;

     straintensor thisGreenStrain;  // from NDMaterial
     stresstensor this2PKStress;   // from NDmaterial


};

#endif

