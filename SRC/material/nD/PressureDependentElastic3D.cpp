/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.4 $                                                              
// $Date: 2002-06-10 22:24:08 $                                     
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PressureDependentElastic3D.cpp,v $

//Boris Jeremic and Zhaohui Yang ___ 07-07-2001
//Pressure dependent elastic isotropic material
                                                                       
                                                                        
#include <PressureDependentElastic3D.h>
#include <Channel.h>
#include <Tensor.h>

PressureDependentElastic3D::PressureDependentElastic3D
(int tag, double E, double nu, double rhop, double expp, double pr, double pop):
 ElasticIsotropicMaterial (tag, ND_TAG_PressureDependentElastic3D, E, nu, rhop),
 sigma(6), D(6,6), epsilon(6), exp(expp), p_ref(pr), po(pop) 
{
  // Set up the elastic constant matrix for 3D elastic isotropic 
  Dt = tensor( 4, def_dim_4, 0.0 ); 
  setInitElasticStiffness();

}

PressureDependentElastic3D::PressureDependentElastic3D():
 ElasticIsotropicMaterial (0, ND_TAG_PressureDependentElastic3D, 0.0, 0.0, 0.0),
 sigma(6), D(6,6), epsilon(6)
{
  Dt = tensor( 4, def_dim_4, 0.0 );
}

PressureDependentElastic3D::~PressureDependentElastic3D ()
{

}

int
PressureDependentElastic3D::setTrialStrain (const Vector &v)
{
  epsilon = v;
  return 0;
}

int
PressureDependentElastic3D::setTrialStrain (const Vector &v, const Vector &r)
{
  epsilon = v;
  return 0;
}

int
PressureDependentElastic3D::setTrialStrainIncr (const Vector &v)
{
  epsilon += v;
  return 0;
}

int
PressureDependentElastic3D::setTrialStrainIncr (const Vector &v, const Vector &r)
{
  epsilon += v;
  return 0;
}

const Matrix&
PressureDependentElastic3D::getTangent (void)
{
  return D;
}

const Vector&
PressureDependentElastic3D::getStress (void)
{
  sigma = D*epsilon;
  return sigma;
}

const Vector&
PressureDependentElastic3D::getStrain (void)
{
  return epsilon;
}

int
PressureDependentElastic3D::setTrialStrain (const Tensor &v)
{
    Strain = v;
    return 0;
}

int
PressureDependentElastic3D::setTrialStrain (const Tensor &v, const Tensor &r)
{
    Strain = v;
    return 0;
}

int
PressureDependentElastic3D::setTrialStrainIncr (const Tensor &v)
{
    //cerr << " before set Tri St Incr " << Strain;
    //cerr << " Strain Incr " << v << endln;
    Strain = Strain + v;
    //cerr << " after setTrialStrainIncr  " << Strain << endln;
    return 0;
}

int
PressureDependentElastic3D::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    Strain = Strain + v;
    return 0;
}

const Tensor&
PressureDependentElastic3D::getTangentTensor (void)
{
    //setElasticStiffness();
    //return Dt;
    return Dt;
}

const stresstensor
PressureDependentElastic3D::getStressTensor (void)
{
    Stress = Dt("ijkl") * Strain("kl");
    //setElasticStiffness();
    return Stress;
}

const straintensor
PressureDependentElastic3D::getStrainTensor (void)
{
    return Strain;
}

const straintensor
PressureDependentElastic3D::getPlasticStrainTensor (void)
{
    //Return zero straintensor
    straintensor t;
    return t;
}

int
PressureDependentElastic3D::commitState (void)
{
    //Set the new Elastic constants
    tensor ret( 4, def_dim_4, 0.0 );
    				       
    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);

    tensor I_ijkl = I2("ij")*I2("kl");


    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    
    //Update E according to 
    Stress = getStressTensor();
    //Dt("ijkl") * Strain("kl");
    double p = Stress.p_hydrostatic();
    //cerr << " p = " <<  p;

    //Cut-off pressure
    if (p <= po) 
      p = po;

    double Ec = E * pow(p/p_ref, exp);
    //cerr << " Eo = " << E << " Ec = " << Ec << " Exp:" << exp<< " p_ref:" << p_ref << " po: " << po<< endl;

    // Building elasticity tensor
    ret = I_ijkl*( Ec*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( Ec / (1.0 + v) );
    
    //ret.print();
    Dt = ret;
    D = Dt;

    return 0;

}

int
PressureDependentElastic3D::revertToLastCommit (void)
{
	// Nothing was previously committed
	return 0;
}

int
PressureDependentElastic3D::revertToStart (void)
{
	// Nothing was previously committed
	return 0;
}

NDMaterial*
PressureDependentElastic3D::getCopy (void)
{
	PressureDependentElastic3D *theCopy =
		new PressureDependentElastic3D (this->getTag(), E, v, rho, exp, p_ref, po);
	//cerr << "In Get copy" <<  *theCopy << endl;
	theCopy->epsilon = this->epsilon;
	theCopy->sigma = this->sigma;
	theCopy->Strain = this->Strain;
	theCopy->Stress = this->Stress;
	// D and Dt are created in the constructor call

	return theCopy;
}

const char*
PressureDependentElastic3D::getType (void) const
{
	return "PressureDependentElastic3D";
}

int
PressureDependentElastic3D::getOrder (void) const
{
	return 6;
}

int 
PressureDependentElastic3D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(6);

	data(0) = this->getTag();
	data(1) = E;
	data(2) = v;
	data(3) = exp;
	data(4) = p_ref;
	data(5) = po;

    	res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"PressureDependentElastic3D::sendSelf");
		return res;
	}

	return res;
}

int
PressureDependentElastic3D::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
	int res = 0;

	static Vector data(6);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"PressureDependentElastic3D::recvSelf");
		return res;
	}
    
	this->setTag((int)data(0));
    	E = data(1);
	v = data(2);
	exp = data(3);
	p_ref = data(4);
	po = data(5);

	// Set up the elastic constant matrix for 3D elastic isotropic
	this->setInitElasticStiffness();
	
	return res;
}
    
void
PressureDependentElastic3D::Print(ostream &s, int flag)
{
	s << "PressureDependentElastic3D" << endl;
	s << "\ttag: " << this->getTag() << endl;
	s << "\tE: " << E << endl;
	s << "\tv: " << v << endl;
	s << "\texp: " << exp << endl;
	s << "\tp_ref: " << p_ref << endl;
	s << "\tpo: " << po << endl;
	//s << "\tD: " << D << endl;
}


//================================================================================
void PressureDependentElastic3D::setInitElasticStiffness(void)
{    
    tensor ret( 4, def_dim_4, 0.0 );
    				       
    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);

    tensor I_ijkl = I2("ij")*I2("kl");


    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    
    //Initialize E according to initial pressure in the gauss point
    //Stress = getStressTensor();
    //Dt("ijkl") * Strain("kl");
    //double po = Stress.p_hydrostatic();

    //cerr << " p_ref = " <<  p_ref << " po = " << po << endln;
    if (po <= 0.0) 
      po = 0.5;
    double Eo = E * pow(po/p_ref, exp);

    //cerr << " E@ref = " << E << " Eo = " << Eo << endln;

    // Building elasticity tensor
    ret = I_ijkl*( Eo*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( Eo / (1.0 + v) );
    
    //ret.print();
    Dt = ret;
    D = Dt;

    return;

}


