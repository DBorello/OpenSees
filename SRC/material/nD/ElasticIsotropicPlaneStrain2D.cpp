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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-07-16 22:14:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStrain2D.cpp,v $
                                                                        
                                                                        
#include <ElasticIsotropicPlaneStrain2D.h>                                                                        
#include <Channel.h>
#include <Tensor.h>

Vector ElasticIsotropicPlaneStrain2D::sigma(3);
Matrix ElasticIsotropicPlaneStrain2D::D(3,3);

ElasticIsotropicPlaneStrain2D::ElasticIsotropicPlaneStrain2D
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicPlaneStrain2d, E, nu, rho),
 Tepsilon(3), Cepsilon(3)
{

}

ElasticIsotropicPlaneStrain2D::ElasticIsotropicPlaneStrain2D():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicPlaneStrain2d, 0.0, 0.0),
 Tepsilon(3), Cepsilon(3)
{

}

ElasticIsotropicPlaneStrain2D::~ElasticIsotropicPlaneStrain2D ()
{

}

int
ElasticIsotropicPlaneStrain2D::setTrialStrain (const Vector &strain)
{
	Tepsilon = strain;

	return 0;
}

int
ElasticIsotropicPlaneStrain2D::setTrialStrain (const Vector &strain, const Vector &rate)
{
	Tepsilon = strain;

	return 0;
}

int
ElasticIsotropicPlaneStrain2D::setTrialStrainIncr (const Vector &strain)
{
	Tepsilon = Cepsilon;
	Tepsilon.addVector(1.0, strain, 1.0);

	return 0;
}

int
ElasticIsotropicPlaneStrain2D::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
	Tepsilon = Cepsilon;
	Tepsilon.addVector(1.0, strain, 1.0);

	return 0;
}

const Matrix&
ElasticIsotropicPlaneStrain2D::getTangent (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	D(0,0) = D(1,1) = mu2+lam;
	D(0,1) = D(1,0) = lam;
	D(2,2) = mu;

	return D;
}

const Vector&
ElasticIsotropicPlaneStrain2D::getStress (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

    double eps0 = Tepsilon(0);
    double eps1 = Tepsilon(1);

    mu2 += lam;

    //sigma = D*epsilon;
	sigma(0) = mu2*eps0 + lam*eps1;
	sigma(1) = lam*eps0 + mu2*eps1;
	sigma(2) = mu*Tepsilon(2);
	
	return sigma;
}

const Vector&
ElasticIsotropicPlaneStrain2D::getStrain (void)
{
	return Tepsilon;
}

int
ElasticIsotropicPlaneStrain2D::commitState (void)
{
	Cepsilon = Tepsilon;

	return 0;
}

int
ElasticIsotropicPlaneStrain2D::revertToLastCommit (void)
{
	Tepsilon = Cepsilon;

	return 0;
}

int
ElasticIsotropicPlaneStrain2D::revertToStart (void)
{
	Cepsilon.Zero();

	return 0;
}

NDMaterial*
ElasticIsotropicPlaneStrain2D::getCopy (void)
{
	ElasticIsotropicPlaneStrain2D *theCopy =
		new ElasticIsotropicPlaneStrain2D (this->getTag(), E, v, rho);

	theCopy->Cepsilon = Cepsilon;

	return theCopy;
}

const char*
ElasticIsotropicPlaneStrain2D::getType (void) const
{
	return "PlaneStrain";
}

int
ElasticIsotropicPlaneStrain2D::getOrder (void) const
{
	return 3;
}

int 
ElasticIsotropicPlaneStrain2D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(7);

	data(0) = this->getTag();
	data(1) = E;
	data(2) = v;
    data(3) = rho;
	data(4) = Cepsilon(0);
	data(5) = Cepsilon(1);
	data(6) = Cepsilon(2);

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"ElasticIsotropicPlaneStrain2D::sendSelf");
		return res;
	}

	return res;
}

int
ElasticIsotropicPlaneStrain2D::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
	int res = 0;

    static Vector data(7);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"ElasticIsotropicPlaneStrain2D::recvSelf");
		return res;
	}
    
	this->setTag((int)data(0));
    E = data(1);
	v = data(2);
    rho = data(3);
	Cepsilon(0) = data(4);
	Cepsilon(1) = data(5);
	Cepsilon(2) = data(6);

	return res;
}
