#ifndef NDMaterial_h
#define NDMaterial_h

// File: ~/material/NDMaterial.h
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for NDMaterial.
// NDMaterial is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) NDMaterial.h, revA"

#include <Material.h>

class Matrix;
class ID;
class Tensor;
class Vector;

class NDMaterial : public Material
{
  public:
    NDMaterial(int tag, int classTag);
    NDMaterial();
    virtual ~NDMaterial();

    // methods to set state and terrieve state using Matrix and Vector classes
    virtual int setTrialStrain(const Vector &v) = 0;
    virtual int setTrialStrain(const Vector &v, const Vector &r) = 0;
    virtual const Matrix &getTangent(void) = 0;
    virtual const Vector &getStress(void) = 0;

    // methods to set and retrieve state using the Tensor class    
    virtual int setTrialStrain(const Tensor &v) = 0;
    virtual int setTrialStrain(const Tensor &v, const Tensor &r) = 0;    
    virtual const Tensor &getTangentTensor(void) = 0;
    virtual const Tensor &getStressTensor(void) = 0;    

    virtual int commitState(void) = 0;
    virtual int revertToLastCommit(void) = 0;

    virtual NDMaterial *getCopy(void) = 0;
    virtual NDMaterial *getCopy(const char *code) = 0;

    virtual const char *getType(void) const = 0;
    virtual int getOrder(void) const = 0;

  protected:
    
  private:
};


#endif

