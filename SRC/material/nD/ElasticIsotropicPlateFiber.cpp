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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlateFiber.cpp,v $
                                                                        
                                                                        

#include <ElasticIsotropicPlateFiber.h>           
#include <Channel.h>
#include <Tensor.h>

Vector ElasticIsotropicPlateFiber::sigma(5);
Matrix ElasticIsotropicPlateFiber::D(5,5);

ElasticIsotropicPlateFiber::ElasticIsotropicPlateFiber
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicPlateFiber, E, nu, rho),
 Tepsilon(5), Cepsilon(5)
{

}

ElasticIsotropicPlateFiber::ElasticIsotropicPlateFiber():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicPlateFiber, 0.0, 0.0),
 Tepsilon(5), Cepsilon(5)
{

}

ElasticIsotropicPlateFiber::~ElasticIsotropicPlateFiber ()
{

}

int
ElasticIsotropicPlateFiber::setTrialStrain (const Vector &strain)
{
	Tepsilon = strain;

	return 0;
}

int
ElasticIsotropicPlateFiber::setTrialStrain (const Vector &strain, const Vector &rate)
{
	Tepsilon = strain;

	return 0;
}

int
ElasticIsotropicPlateFiber::setTrialStrainIncr (const Vector &strain)
{
	Tepsilon = Cepsilon;
	Tepsilon.addVector(1.0, strain, 1.0);

	return 0;
}

int
ElasticIsotropicPlateFiber::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
	Tepsilon = Cepsilon;
	Tepsilon.addVector(1.0, strain, 1.0);

	return 0;
}

const Matrix&
ElasticIsotropicPlateFiber::getTangent (void)
{
    double d00 = E/(1.0-v*v);
    double d01 = v*d00;
    double d22 = 0.5*(d00-d01);

	D(0,0) = D(1,1) = d00;
	D(0,1) = D(1,0) = d01;
	D(2,2) = d22;
    D(2,3) = d22;
    D(4,4) = d22;

	return D;
}

const Vector&
ElasticIsotropicPlateFiber::getStress (void)
{
    double d00 = E/(1.0-v*v);
    double d01 = v*d00;
    double d22 = 0.5*(d00-d01);

    double eps0 = Tepsilon(0);
    double eps1 = Tepsilon(1);

	//sigma = D*epsilon;
	sigma(0) = d00*eps0 + d01*eps1;
	sigma(1) = d01*eps0 + d00*eps1;

	sigma(2) = d22*Tepsilon(2);
	sigma(3) = d22*Tepsilon(3);
	sigma(4) = d22*Tepsilon(4);
	
	return sigma;
}

const Vector&
ElasticIsotropicPlateFiber::getStrain (void)
{
	return Tepsilon;
}

int
ElasticIsotropicPlateFiber::commitState (void)
{
	Cepsilon = Tepsilon;

	return 0;
}

int
ElasticIsotropicPlateFiber::revertToLastCommit (void)
{
	Tepsilon = Cepsilon;

	return 0;
}

int
ElasticIsotropicPlateFiber::revertToStart (void)
{
	Cepsilon.Zero();

	return 0;
}

NDMaterial*
ElasticIsotropicPlateFiber::getCopy (void)
{
	ElasticIsotropicPlateFiber *theCopy =
		new ElasticIsotropicPlateFiber (this->getTag(), E, v, rho);

	theCopy->Cepsilon = Cepsilon;

	return theCopy;
}

const char*
ElasticIsotropicPlateFiber::getType (void) const
{
	return "PlateFiber";
}

int
ElasticIsotropicPlateFiber::getOrder (void) const
{
	return 5;
}

int 
ElasticIsotropicPlateFiber::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(9);

	data(0) = this->getTag();
	data(1) = E;
	data(2) = v;
    data(3) = rho;
	data(4) = Cepsilon(0);
	data(5) = Cepsilon(1);
	data(6) = Cepsilon(2);
    data(7) = Cepsilon(3);
    data(8) = Cepsilon(4);

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"ElasticIsotropicPlateFiber::sendSelf");
		return res;
	}

	return res;
}

int
ElasticIsotropicPlateFiber::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
	int res = 0;

    static Vector data(9);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"ElasticIsotropicPlateFiber::recvSelf");
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
	Cepsilon(4) = data(8);
	
	return res;
}
