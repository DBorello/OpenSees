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
                                                                        
// $Revision: 1.7 $                                                              
// $Date: 2001-05-26 05:33:05 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropic3D.cpp,v $                                                                

//Boris Jeremic and Zhaohui Yang ___ 02-10-2000
                                                                       
                                                                        
#include <ElasticIsotropic3D.h>
#include <Channel.h>
#include <Tensor.h>

Matrix ElasticIsotropic3D::D(6,6);
Vector ElasticIsotropic3D::sigma(6);

ElasticIsotropic3D::ElasticIsotropic3D
(int tag, double E, double nu, double expp, double pr, double pop):
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropic3D, E, nu),
 epsilon(6), exp(expp), p_ref(pr), po(pop) 
{
	// Set up the elastic constant matrix for 3D elastic isotropic 
	Dt = tensor( 4, def_dim_4, 0.0 ); 
	setInitElasticStiffness();
}

ElasticIsotropic3D::ElasticIsotropic3D():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropic3D, 0.0, 0.0),
 epsilon(6)
{
       Dt = tensor( 4, def_dim_4, 0.0 );
}

ElasticIsotropic3D::~ElasticIsotropic3D ()
{

}

int
ElasticIsotropic3D::setTrialStrain (const Vector &v)
{
	epsilon = v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrain (const Vector &v, const Vector &r)
{
	epsilon = v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Vector &v)
{
	epsilon += v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Vector &v, const Vector &r)
{
	epsilon += v;

	return 0;
}

const Matrix&
ElasticIsotropic3D::getTangent (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
			D(j,i) = lam;
		D(i,i) += mu2;
		D(i+3,i+3) = mu;
	}

	return D;
}

const Vector&
ElasticIsotropic3D::getStress (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	double tmp1 = lam+mu2;
	
	double eps0 = epsilon(0);
	double eps1 = epsilon(1);
	double eps2 = epsilon(2);

	sigma(0) = tmp1*eps0 + lam*(eps1+eps2);
	sigma(1) = tmp1*eps1 + lam*(eps2+eps0);
	sigma(2) = tmp1*eps2 + lam*(eps0+eps1);

	sigma(3) = mu*epsilon(3);
	sigma(4) = mu*epsilon(4);
	sigma(5) = mu*epsilon(5);

	return sigma;
}

const Vector&
ElasticIsotropic3D::getStrain (void)
{
	return epsilon;
}

int
ElasticIsotropic3D::setTrialStrain (const Tensor &v)
{
    Strain = v;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrain (const Tensor &v, const Tensor &r)
{
    Strain = v;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Tensor &v)
{
    //cerr << " before set Tri St Incr " << Strain;
    //cerr << " Strain Incr " << v << endln;
    Strain = Strain + v;
    //cerr << " after setTrialStrainIncr  " << Strain << endln;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    Strain = Strain + v;
    return 0;
}

const Tensor&
ElasticIsotropic3D::getTangentTensor (void)
{
    //setElasticStiffness();
    //return Dt;
    return Dt_commit;
}

const stresstensor ElasticIsotropic3D::getStressTensor (void)
{
    Stress = Dt("ijkl") * Strain("kl");
    //setElasticStiffness();
    return Stress;
}

const Tensor&
ElasticIsotropic3D::getStrainTensor (void)
{
    return Strain;
}

int
ElasticIsotropic3D::commitState (void)
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

    if (p <= 0.5) 
      p = 0.5;

    double Ec = E * pow(p/p_ref, exp);
    //cerr << " Eo = " << E << " Ec = " << Ec << endln;

    // Building elasticity tensor
    ret = I_ijkl*( Ec*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( Ec / (1.0 + v) );
    
    //ret.print();
    Dt_commit = ret;
    Dt = ret;
    
    return 0;

}

int
ElasticIsotropic3D::revertToLastCommit (void)
{
	// Nothing was previously committed
	return 0;
}

int
ElasticIsotropic3D::revertToStart (void)
{
	// Nothing was previously committed
	return 0;
}

NDMaterial*
ElasticIsotropic3D::getCopy (void)
{
	ElasticIsotropic3D *theCopy =
		new ElasticIsotropic3D (this->getTag(), E, v);
	theCopy->epsilon = this->epsilon;
	theCopy->Strain = this->Strain;
	theCopy->Stress = this->Stress;
	// D and Dt are created in the constructor call

	return theCopy;
}

const char*
ElasticIsotropic3D::getType (void) const
{
	return "ThreeDimensional";
}

int
ElasticIsotropic3D::getOrder (void) const
{
	return 6;
}

int 
ElasticIsotropic3D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(3);

	data(0) = this->getTag();
	data(1) = E;
	data(2) = v;

    	res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"ElasticIsotropic3D::sendSelf");
		return res;
	}

	return res;
}

int
ElasticIsotropic3D::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
	int res = 0;

	static Vector data(6);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"ElasticIsotropic3D::recvSelf");
		return res;
	}
    
	this->setTag((int)data(0));
    	E = data(1);
	v = data(2);

	// Set up the elastic constant matrix for 3D elastic isotropic
	D.Zero();
	//setElasticStiffness();
	
	return res;
}
    
void
ElasticIsotropic3D::Print(ostream &s, int flag)
{
	s << "ElasticIsotropic3D" << endl;
	s << "\ttag: " << this->getTag() << endl;
	s << "\tE: " << E << endl;
	s << "\tv: " << v << endl;
	s << "\texp: " << exp << endl;
	s << "\tp_ref: " << p_ref << endl;
	//s << "\tD: " << D << endl;
}


//================================================================================
void ElasticIsotropic3D::setInitElasticStiffness(void)
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
    //double po = 100.0; //kPa
    if (po <= 0.5) 
      po = 0.5;
    double Eo = E * pow(po/p_ref, exp);

    //cerr << " E@ref = " << E << " Eo = " << Eo << endln;

    // Building elasticity tensor
    ret = I_ijkl*( Eo*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( Eo / (1.0 + v) );
    
    //ret.print();
    Dt = ret;
    Dt_commit = ret;

    //D = Dt;

    return;

}


