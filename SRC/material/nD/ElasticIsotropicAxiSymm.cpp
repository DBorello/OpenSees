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
                                                                        
// $Revision: 1.1 $
// $Date: 2001-07-16 22:59:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicAxiSymm.cpp,v $
                                                                 
#include <ElasticIsotropicAxiSymm.h>                                                                        
#include <Channel.h>
#include <Tensor.h>

Vector ElasticIsotropicAxiSymm::sigma(4);
Matrix ElasticIsotropicAxiSymm::D(4,4);

ElasticIsotropicAxiSymm::ElasticIsotropicAxiSymm
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicAxiSymm, E, nu, rho),
 Tepsilon(4), Cepsilon(4)
{

}

ElasticIsotropicAxiSymm::ElasticIsotropicAxiSymm():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicAxiSymm, 0.0, 0.0),
 Tepsilon(4), Cepsilon(4)
{

}

ElasticIsotropicAxiSymm::~ElasticIsotropicAxiSymm ()
{

}

int
ElasticIsotropicAxiSymm::setTrialStrain (const Vector &strain)
{
	Tepsilon = strain;

	return 0;
}

int
ElasticIsotropicAxiSymm::setTrialStrain (const Vector &strain, const Vector &rate)
{
	Tepsilon = strain;

	return 0;
}

int
ElasticIsotropicAxiSymm::setTrialStrainIncr (const Vector &strain)
{
	Tepsilon = Cepsilon;
	Tepsilon.addVector(1.0, strain, 1.0);

	return 0;
}

int
ElasticIsotropicAxiSymm::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
	Tepsilon = Cepsilon;
	Tepsilon.addVector(1.0, strain, 1.0);

	return 0;
}

const Matrix&
ElasticIsotropicAxiSymm::getTangent (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	D(0,0) = D(1,1) = D(2,2) = mu2+lam;
	D(0,1) = D(1,0) = lam;
    D(0,2) = D(2,0) = lam;
    D(1,2) = D(2,1) = lam;
	D(3,3) = mu;

	return D;
}

const Vector&
ElasticIsotropicAxiSymm::getStress (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

    double eps0 = Tepsilon(0);
    double eps1 = Tepsilon(1);
    double eps2 = Tepsilon(2);

    mu2 += lam;

    //sigma = D*epsilon;
	sigma(0) = mu2*eps0 + lam*(eps1+eps2);
	sigma(1) = mu2*eps1 + lam*(eps0+eps2);
    sigma(2) = mu2*eps2 + lam*(eps0+eps1);
	sigma(3) = mu*Tepsilon(3);
	
	return sigma;
}

const Vector&
ElasticIsotropicAxiSymm::getStrain (void)
{
	return Tepsilon;
}

int
ElasticIsotropicAxiSymm::commitState (void)
{
	Cepsilon = Tepsilon;

	return 0;
}

int
ElasticIsotropicAxiSymm::revertToLastCommit (void)
{
	Tepsilon = Cepsilon;

	return 0;
}

int
ElasticIsotropicAxiSymm::revertToStart (void)
{
	Cepsilon.Zero();

	return 0;
}

NDMaterial*
ElasticIsotropicAxiSymm::getCopy (void)
{
	ElasticIsotropicAxiSymm *theCopy =
		new ElasticIsotropicAxiSymm (this->getTag(), E, v, rho);

	theCopy->Cepsilon = Cepsilon;

	return theCopy;
}

const char*
ElasticIsotropicAxiSymm::getType (void) const
{
	return "AxiSymmetric";
}

int
ElasticIsotropicAxiSymm::getOrder (void) const
{
	return 4;
}

int 
ElasticIsotropicAxiSymm::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(8);

	data(0) = this->getTag();
	data(1) = E;
	data(2) = v;
    data(3) = rho;
	data(4) = Cepsilon(0);
	data(5) = Cepsilon(1);
	data(6) = Cepsilon(2);
    data(7) = Cepsilon(3);

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"ElasticIsotropicAxiSymm::sendSelf");
		return res;
	}

	return res;
}

int
ElasticIsotropicAxiSymm::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
	int res = 0;

    static Vector data(8);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"ElasticIsotropicAxiSymm::recvSelf");
		return res;
	}
    
	this->setTag((int)data(0));
    E = data(1);
	v = data(2);
    rho = data(3);
	Cepsilon(0) = data(4);
	Cepsilon(1) = data(5);
	Cepsilon(2) = data(6);
    Cepsilon(3) = data(7);

	return res;
}
