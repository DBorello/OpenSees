/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2001-06-14 08:06:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/transformation/NatafXuTransformation.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <XuTransformation.h>
#include <NatafXuTransformation.h>
#include <RandomVariable.h>
#include <CorrelationCoefficient.h>
#include <NormalRV.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixOperations.h>
#include <math.h>
#include <string.h>


NatafXuTransformation::NatafXuTransformation(ReliabilityDomain *passedReliabilityDomain)
:XuTransformation()
{
	theReliabilityDomain = passedReliabilityDomain;


	// Find and set problem size (number of random variables)
	nrv = theReliabilityDomain->getNumberOfRandomVariables();


	// Create/initialize vectors and matrices
	x = new Vector(nrv);
	u = new Vector(nrv);
	jacobian_x_u = new Matrix(nrv,nrv);
	jacobian_u_x = new Matrix(nrv,nrv);
	lowerCholesky = new Matrix(nrv,nrv);
	inverseLowerCholesky = new Matrix(nrv,nrv);


	// Establish correlation matrix according to the Nataf assumption
	setCorrelationMatrix();


	// Create object to do matrix operations on the correlation matrix
	theMatrixOperations = 0;
	theMatrixOperations = new MatrixOperations(*correlationMatrix);
	if (theMatrixOperations == 0) {
		cerr << "NatafXuTransformation::NatafXuTransformation() - could " << endl
			<< " not create the object to perform matrix operations." << endl;
	}


	// Cholesky decomposition of correlation matrix
	int result = theMatrixOperations->computeCholeskyAndItsInverse();
	if (result < 0) {
		cerr << "NatafXuTransformation::NatafXuTransformation() - could not" << endl
			<< " compute the Cholesky decomposition and its inverse " << endl
			<< " for the correlation matrix." << endl;
	}
	(*lowerCholesky) = theMatrixOperations->getLowerCholesky();
	(*inverseLowerCholesky) = theMatrixOperations->getInverseLowerCholesky();

}

NatafXuTransformation::~NatafXuTransformation()
{
	delete correlationMatrix;
	delete lowerCholesky;
	delete inverseLowerCholesky;
	delete jacobian_x_u;
	delete jacobian_u_x;
	delete x;
	delete u;
	delete theMatrixOperations;
}



int 
NatafXuTransformation::set_x(Vector passedx)
{
	(*x) = passedx; // (later: check size of vector, etc.)
	return 0;
}

int 
NatafXuTransformation::set_u(Vector passedu)
{
	(*u) = passedu; // (later: check size of vector, etc.)
	return 0;
}

int 
NatafXuTransformation::transform_x_to_u()
{
	Vector z = x_to_z(*x);
	(*u) = (*inverseLowerCholesky) * z;
	
	return 0;
}

int 
NatafXuTransformation::transform_u_to_x()
{
	Vector z = (*lowerCholesky) * (*u);
	(*x) = z_to_x(z);

	return 0;
}


int 
NatafXuTransformation::transform_u_to_x_andComputeJacobian()
{
	Vector z = (*lowerCholesky) * (*u);
	(*x) = z_to_x(z);


	Matrix Jzx = getJacobian_z_x((*x),z);
	(*jacobian_u_x) = (*inverseLowerCholesky) * Jzx;


	int result = theMatrixOperations->setMatrix((*jacobian_u_x));
	if (result < 0) {
		cerr << "NatafXuTransformation::transform_u_to_x() - could not set " << endl
			<< " the matrix in the object to perform matrix operations." << endl;
		return -1;
	}

	
	result = theMatrixOperations->computeInverse();
	if (result < 0) {
		cerr << "NatafXuTransformation::transform_u_to_x() - could not " << endl
			<< " invert Jacobian_u_x." << endl;
		return -1;
	}
	(*jacobian_x_u) = theMatrixOperations->getInverse();

	return 0;
}







Vector 
NatafXuTransformation::get_x()
{
	return (*x);
}



Vector 
NatafXuTransformation::get_u()
{
	return (*u);
}



Matrix 
NatafXuTransformation::getJacobian_x_u()
{
	return (*jacobian_x_u);
}



Matrix 
NatafXuTransformation::getJacobian_u_x()
{
	return (*jacobian_u_x);
}







Matrix
NatafXuTransformation::getJacobian_z_x(Vector x, Vector z)
{	
	RandomVariable *theRV;
	NormalRV aStandardNormalRV(1, 0.0, 1.0, 0.0);
	Matrix jacobianMatrix_z_x(nrv,nrv);
	for ( int i=0 ; i<nrv ; i++ )
	{
		theRV = theReliabilityDomain->getRandomVariablePtr(i+1);
		if (strcmp(theRV->getType(),"NORMAL")==0) {
			double sigma = theRV->getParameter2();
			jacobianMatrix_z_x(i,i) = 1.0 / sigma;
		}
		else if (strcmp(theRV->getType(),"LOGNORMAL")==0) {
			double zeta = theRV->getParameter2();
			jacobianMatrix_z_x(i,i) = 1.0 / ( zeta * x(i)  );
		}
		else {
			jacobianMatrix_z_x(i,i) = theRV->getPDFvalue(x(i)) / aStandardNormalRV.getPDFvalue(z(i));
			if (jacobianMatrix_z_x(i,i)==0.0) {
			}
		}
	}

	return jacobianMatrix_z_x;
}




Vector
NatafXuTransformation::z_to_x(Vector z)
{
	RandomVariable *theRV;
	NormalRV aStandardNormalRV(1, 0.0, 1.0, 0.0);
	Vector x(nrv);
	for ( int i=0 ; i<nrv ; i++ )
	{
		theRV = theReliabilityDomain->getRandomVariablePtr(i+1);
		if (strcmp(theRV->getType(),"NORMAL")==0) {
			double mju = theRV->getParameter1();
			double sigma = theRV->getParameter2();
			x(i) = z(i) * sigma + mju;
		}
		else if (strcmp(theRV->getType(),"LOGNORMAL")==0) {
			double lambda = theRV->getParameter1();
			double zeta = theRV->getParameter2();
			x(i) = exp ( z(i) * zeta + lambda );
		}
		else {
			x(i) = theRV->getInverseCDFvalue(aStandardNormalRV.getCDFvalue(z(i)));
		}
	}
	return x;
}





Vector
NatafXuTransformation::x_to_z(Vector x)
{
	RandomVariable *theRV;
	NormalRV aStandardNormalRV(1, 0.0, 1.0, 0.0);
	Vector z(nrv);
	for ( int i=0 ; i<nrv ; i++ )
	{
		theRV = theReliabilityDomain->getRandomVariablePtr(i+1);
		if (strcmp(theRV->getType(),"NORMAL")==0) {
			double mju = theRV->getParameter1();
			double sigma = theRV->getParameter2();
			z(i) =   ( x(i) - mju ) / sigma;
		}
		else if (strcmp(theRV->getType(),"LOGNORMAL")==0) {
			double lambda = theRV->getParameter1();
			double zeta = theRV->getParameter2();
			z(i) = ( log ( x(i) ) - lambda ) / zeta;
		}
		else {
			z(i) = aStandardNormalRV.getInverseCDFvalue(theRV->getCDFvalue(x(i)));
		}
	}
	return z;
}













void
NatafXuTransformation::setCorrelationMatrix()
{
	// Initial declarations
	char* typeRv1;
	char* typeRv2;
	int rv1;
	int rv2;
	RandomVariable *rv1Ptr;
	RandomVariable *rv2Ptr;
	double correlation;
	double newCorrelation;

	// Initialize correlation matrix
	correlationMatrix = new Matrix(nrv,nrv);

	// Put 'ones' on the diagonal
	for ( int i=0 ; i<nrv ; i++ )
	{
		(*correlationMatrix)(i,i) = 1.0;
	}

	// Get number of correlation coefficients
	int numberOfCorrelationCoefficients = 
		theReliabilityDomain->getNumberOfCorrelationCoefficients();

	// Modify each coefficient at a time and put it into the correlation matrix
	for ( int j=1 ; j<=numberOfCorrelationCoefficients ; j++ )
	{
		// Get a pointer to the correlation coefficient with tag 'j'
		CorrelationCoefficient *theCorrelationCoefficient = 
			theReliabilityDomain->getCorrelationCoefficientPtr(j);

		// Get value of the correlation
		correlation = theCorrelationCoefficient->getCorrelation();

		// Get tags for the two involved random variables
		rv1 = theCorrelationCoefficient->getRv1();
		rv2 = theCorrelationCoefficient->getRv2();

		// Get pointers to the two random variables
		rv1Ptr = theReliabilityDomain->getRandomVariablePtr(rv1);
		rv2Ptr = theReliabilityDomain->getRandomVariablePtr(rv2);

		// Get the types of the two random variables
		typeRv1 = rv1Ptr->getType();
		typeRv2 = rv2Ptr->getType();

		// Compute the coefficient of variation of the random variables
		double cov1 = rv1Ptr->getStdv() / rv1Ptr->getMean();
		double cov2 = rv2Ptr->getStdv() / rv2Ptr->getMean();

		// Modify the correlation coefficient according to 
		// the type of the the two involved random variables

		/////////////////////////////////////////////////////////////////////////////////
		if ( strcmp(typeRv1,"		NORMAL") == 0  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			double zeta2 = sqrt ( log ( 1.0 + cov2 * cov2 ) );
			newCorrelation = correlation * cov2 / zeta2;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			newCorrelation = (1.001 - 0.007 * cov2 + 0.118 * cov2 * cov2) * correlation;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = 1.107 * correlation;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = 1.014 * correlation;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = 1.023 * correlation;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			newCorrelation = (1.026 + 0.001*correlation - 0.178*u2
				+ 0.268*s2 - 0.001*correlation*correlation
				+ 0.178*u2*u2 - 0.679*s2*s2 - 0.003*s2*correlation)  *  correlation;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = 1.031 * correlation;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = 1.031 * correlation;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			newCorrelation = (1.030 + 0.238*cov2 + 0.364*cov2*cov2) * correlation;
		}
		else if ( strcmp(typeRv1,"NORMAL") == 0  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			newCorrelation = (1.031 - 0.195*cov2 + 0.328*cov2*cov2) * correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			double zeta1 = sqrt ( log ( 1 + pow ( cov1 , 2 ) ) );
			newCorrelation = correlation * cov1 / zeta1;
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			double zeta1 = sqrt ( log ( 1.0 + cov1*cov1 ) );
			double zeta2 = sqrt ( log ( 1.0 + cov2*cov2 ) );
			newCorrelation = 1.0 / ( zeta1 * zeta2 ) * log ( 1 + correlation * cov1 * cov2 );
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			double temp = 1.001 + 0.033*correlation + 0.004*cov1 - 0.016*cov2 
				+ 0.002*correlation*correlation + 0.223*cov1*cov1 + 0.0130*cov2*cov2;
			newCorrelation = correlation * (temp - 0.104*correlation*cov1 + 0.029*cov1*cov2 - 0.119*correlation*cov2);
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.098 + 0.003*correlation + 0.019*cov1 + 0.025*correlation*correlation
				+ 0.303*cov1*cov1 - 0.437*correlation*cov1) * correlation;
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.011 + 0.001*correlation + 0.014*cov1 + 0.004*correlation*correlation 
				+ 0.231*cov1*cov1 - 0.130*correlation*cov1) * correlation;
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.019 + 0.014*cov1 + 0.010*correlation*correlation + 0.249*cov1*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			double uu = u2*u2;
			double ss = s2*s2;
			double xu = correlation * u2;
			double xs = correlation * s2;
			double us = u2*s2;
			double up = u2*cov1;
			double sp = s2 * cov1;
			double temp = 0.979 + 0.053*correlation + 0.181*u2 + 0.293*s2 + 0.002*cov1
				- 0.004*correlation*correlation - 0.181*uu + 5.796*ss + 0.277*cov1*cov1 - 0.107*xu
				- 0.619*xs - 0.190*correlation*cov1 - 3.976*us - 0.097*up + 0.133*sp
				- 14.402*ss*s2 - 0.069*cov1*cov1*cov1
				+ 0.031*correlation*xs + 0.015*correlation*correlation*cov1
				+ 3.976*uu*s2 + 0.097*uu*cov1 - 0.430*ss*cov1
				+ 0.113*sp*cov1 + 1.239*correlation*us + 0.380*correlation*up;
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.029 + 0.001*correlation + 0.014*cov1 + 0.004*correlation*correlation
				+ 0.233*cov1*cov1 - 0.197*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.029 - 0.001*correlation + 0.014*cov1 + 0.004*correlation*correlation
				+ 0.233*cov1*cov1 + 0.197*correlation*cov1) * correlation;
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			double temp = 1.026 + 0.082*correlation - 0.019*cov1 + 0.222*cov2 + 0.018*correlation*correlation
				+ 0.288*cov1*cov1 + 0.379*cov2*cov2;
			newCorrelation = (temp - 0.441*correlation*cov1 + 0.126*cov1*cov2 - 0.277*correlation*cov2);
		}
		else if ( strcmp(typeRv1,"LOGNORMAL") == 0  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			double temp = 1.031+0.052*correlation+0.011*cov1-0.210*cov2+0.002*correlation*correlation+0.220*cov1*cov1+0.350*cov2*cov2;
			newCorrelation = (temp+0.005*correlation*cov1+0.009*cov1*cov2-0.174*correlation*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			newCorrelation = (1.001 - 0.007 * cov1 + 0.118 * cov1 * cov1) * correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			double temp = 1.001 + 0.033*correlation + 0.004*cov2 - 0.016*cov1 
				+ 0.002*correlation*correlation + 0.223*cov2*cov2 + 0.0130*cov1*cov1;
			newCorrelation = correlation * (temp - 0.104*correlation*cov2 + 0.029*cov2*cov1 - 0.119*correlation*cov1);
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			double temp = 1.002+0.022*correlation-0.012*(cov1+cov2)+0.001*correlation*correlation+0.125*(cov1*cov1+cov2*cov2);
			newCorrelation = (temp-0.077*(correlation*cov1+correlation*cov2)+0.014*cov1*cov2)*correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.104+0.003*correlation-0.008*cov1+0.014*correlation*correlation+0.173*cov1*cov1-0.296*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.014+0.001*correlation-0.007*cov1+0.002*correlation*correlation+0.126*cov1*cov1-0.090*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.023+0.000*correlation-0.007*cov1+0.002*correlation*correlation+0.127*cov1*cov1-0.000*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			double uu=u2*u2;
			double ss=s2*s2;
			double xu=correlation*u2;
			double xs=correlation*s2;
			double us=u2*s2;
			double up=u2*cov1;
			double sp=s2*cov1;
			double temp;
			double xp=correlation*cov1;
			if(correlation > 0.0) {
				temp = 0.931+0.050*correlation+0.366*u2+0.549*s2+0.181*cov1
				-0.055*correlation*correlation-0.515*uu+4.804*ss-0.484*cov1*cov1-0.064*xu
				-0.637*xs+0.032*xp-4.059*us-0.319*up-0.211*sp
				+0.052*correlation*correlation*correlation+0.227*uu*u2-10.220*ss*s2+0.559*cov1*cov1*cov1-0.042*correlation*xu
				+0.223*correlation*xs-0.172*correlation*xp+0.028*correlation*uu+0.695*correlation*ss+0.126*correlation*cov1*cov1
				+3.845*uu*s2+0.019*uu*cov1-1.244*us*s2+0.008*up*cov1-2.075*ss*cov1
				+0.167*sp*cov1+0.666*correlation*us+0.386*correlation*up-0.517*correlation*sp+2.125*us*cov1;
			}
			else {
				temp = 1.025+0.050*correlation-0.029*u2+0.047*s2-0.136*cov1
				+0.069*correlation*correlation+0.178*uu+6.281*ss+0.548*cov1*cov1-0.027*xu
				-0.686*xs+0.046*xp-3.513*us+0.231*up+0.299*sp
				+0.063*correlation*correlation*correlation-0.226*uu*u2-17.507*ss*s2-0.366*cov1*cov1*cov1+0.051*correlation*xu
				-0.246*correlation*xs+0.186*correlation*xp-0.001*correlation*uu+0.984*correlation*ss+0.121*correlation*cov1*cov1
				+3.700*uu*s2+0.081*uu*cov1+1.356*us*s2+0.002*up*cov1+1.654*ss*cov1
				-0.135*sp*cov1+0.619*correlation*us+0.410*correlation*up-0.686*correlation*sp-2.205*us*cov1;
			}
			newCorrelation = temp*correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.031+0.001*correlation-0.007*cov1+0.003*correlation*correlation+0.131*cov1*cov1-0.132*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.031-0.001*correlation-0.007*cov1+0.003*correlation*correlation+0.131*cov1*cov1+0.132*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			double temp = 1.029+0.056*correlation-0.030*cov1+0.225*cov2+0.012*correlation*correlation+0.174*cov1*cov1+0.379*cov2*cov2;
			newCorrelation = (temp-0.313*correlation*cov1+0.075*cov1*cov2-0.182*correlation*cov2)*correlation;
		}
		else if ( strcmp(typeRv1,"GAMMA") == 0  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			double temp = 1.032+0.034*correlation-0.007*cov1-0.202*cov2+0.000*correlation*correlation+0.121*cov1*cov1+0.339*cov2*cov2;
			newCorrelation = (temp-0.006*correlation*cov1+0.003*cov1*cov2-0.111*correlation*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			newCorrelation = 1.107 * correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			newCorrelation = (1.098 + 0.003*correlation + 0.019*cov2 + 0.025*correlation*correlation
				+ 0.303*cov2*cov2 - 0.437*correlation*cov2) * correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			newCorrelation = (1.104+0.003*correlation-0.008*cov2+0.014*correlation*correlation+0.173*cov2*cov2-0.296*correlation*cov2)*correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.229-0.367*correlation+0.153*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.123-0.100*correlation+0.021*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.133+0.029*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			double uu=u2*u2;
			double ss=s2*s2;
			double xu=correlation*u2;
			double xs=correlation*s2;
			double us=u2*s2;
			double temp = 1.082-0.004*correlation+0.204*u2+0.432*s2-0.001*correlation*correlation
				- 0.204*uu+7.728*ss+0.008*xu-1.699*xs-5.338*us
				- 19.741*ss*s2+0.135*correlation*xs+5.338*uu*s2+3.397*correlation*us;
			newCorrelation = temp*correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.142-0.154*correlation+0.031*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.142+0.154*correlation+0.031*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			newCorrelation = (1.109-0.152*correlation+0.361*cov2+0.130*correlation*correlation+0.455*cov2*cov2-0.728*correlation*cov2)*correlation;
		}
		else if ( (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 )  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			newCorrelation = (1.147+0.145*correlation-0.271*cov2+0.010*correlation*correlation+0.459*cov2*cov2-0.467*correlation*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			newCorrelation = 1.014 * correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			newCorrelation = (1.011 + 0.001*correlation + 0.014*cov2 + 0.004*correlation*correlation 
				+ 0.231*cov2*cov2 - 0.130*correlation*cov2) * correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			newCorrelation = (1.014+0.001*correlation-0.007*cov2+0.002*correlation*correlation+0.126*cov2*cov2-0.090*correlation*cov2)*correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.123-0.100*correlation+0.021*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.028-0.029*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.038-0.008*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			double temp = 1.037-0.042*correlation-0.182*u2+0.369*s2-0.001*correlation*correlation
				+0.182*u2*u2-1.150*s2*s2+0.084*correlation*u2;
			newCorrelation = temp*correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.046-0.045*correlation+0.006*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.046+0.045*correlation+0.006*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			newCorrelation = (1.036-0.038*correlation+0.266*cov2+0.028*correlation*correlation+0.383*cov2*cov2-0.229*correlation*cov2)*correlation;
		}
		else if ( (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 )  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			newCorrelation = (1.047+0.042*correlation-0.212*cov2+0.353*cov2*cov2-0.136*correlation*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			newCorrelation = 1.023 * correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			newCorrelation = (1.019 + 0.014*cov2 + 0.010*correlation*correlation + 0.249*cov2*cov2)*correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			newCorrelation = (1.023+0.000*correlation-0.007*cov2+0.002*correlation*correlation+0.127*cov2*cov2-0.000*correlation*cov2)*correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.133+0.029*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.038-0.008*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.047-0.047*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			double temp = 1.040+0.015*correlation-0.176*u2+0.432*s2-0.008*correlation*correlation
				+0.176*u2*u2-1.286*s2*s2-0.137*correlation*s2;
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.055+0.015*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.055+0.015*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			newCorrelation = (1.033+0.305*cov2+0.074*correlation*correlation+0.405*cov2*cov2)*correlation;
		}
		else if ( strcmp(typeRv1,"UNIFORM") == 0  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			newCorrelation = (1.061-0.237*cov2-0.005*correlation*correlation+0.379*cov2*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( strcmp(typeRv1,"BETA") == 0  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			newCorrelation = (1.026 + 0.001*correlation - 0.178*u1
				+ 0.268*s1 - 0.001*correlation*correlation
				+ 0.178*u1*u1 - 0.679*s1*s1 - 0.003*s1*correlation)  *  correlation;

		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			double uu = u1*u1;
			double ss = s1*s1;
			double xu = correlation * u1;
			double xs = correlation * s1;
			double us = u1*s1;
			double up = u1*cov2;
			double sp = s1 * cov2;
			double temp = 0.979 + 0.053*correlation + 0.181*u1 + 0.293*s1 + 0.002*cov2
				- 0.004*correlation*correlation - 0.181*uu + 5.796*ss + 0.277*cov2*cov2 - 0.107*xu
				- 0.619*xs - 0.190*correlation*cov2 - 3.976*us - 0.097*up + 0.133*sp
				- 14.402*ss*s1 - 0.069*cov2*cov2*cov2
				+ 0.031*correlation*xs + 0.015*correlation*correlation*cov2
				+ 3.976*uu*s1 + 0.097*uu*cov2 - 0.430*ss*cov2
				+ 0.113*sp*cov2 + 1.239*correlation*us + 0.380*correlation*up;
			newCorrelation = temp * correlation;

		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			double ba2 = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba2; 
			double s2 = rv2Ptr->getStdv() / ba2;
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			double uu=u1*u1;
			double ss=s1*s1;
			double xu=correlation*u1;
			double xs=correlation*s1;
			double us=u1*s1;
			double up=u1*cov2;
			double sp=s1*cov2;
			double temp;
			double xp=correlation*cov2;
			if(correlation > 0.0) {
				temp = 0.931+0.050*correlation+0.366*u1+0.549*s1+0.181*cov2
				-0.055*correlation*correlation-0.515*uu+4.804*ss-0.484*cov2*cov2-0.064*xu
				-0.637*xs+0.032*xp-4.059*us-0.319*up-0.211*sp
				+0.052*correlation*correlation*correlation+0.227*uu*u2-10.220*ss*s2+0.559*cov2*cov2*cov2-0.042*correlation*xu
				+0.223*correlation*xs-0.172*correlation*xp+0.028*correlation*uu+0.695*correlation*ss+0.126*correlation*cov2*cov2
				+3.845*uu*s2+0.019*uu*cov2-1.244*us*s1+0.008*up*cov2-2.075*ss*cov2
				+0.167*sp*cov2+0.666*correlation*us+0.386*correlation*up-0.517*correlation*sp+2.125*us*cov2;
			}
			else {
				temp = 1.025+0.050*correlation-0.029*u1+0.047*s1-0.136*cov2
				+0.069*correlation*correlation+0.178*uu+6.281*ss+0.548*cov2*cov2-0.027*xu
				-0.686*xs+0.046*xp-3.513*us+0.231*up+0.299*sp
				+0.063*correlation*correlation*correlation-0.226*uu*u1-17.507*ss*s2-0.366*cov2*cov2*cov2+0.051*correlation*xu
				-0.246*correlation*xs+0.186*correlation*xp-0.001*correlation*uu+0.984*correlation*ss+0.121*correlation*cov2*cov2
				+3.700*uu*s1+0.081*uu*cov2+1.356*us*s1+0.002*up*cov2+1.654*ss*cov2
				-0.135*sp*cov2+0.619*correlation*us+0.410*correlation*up-0.686*correlation*sp-2.205*us*cov2;
			}
			newCorrelation = temp*correlation;

		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			double ba2 = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba2; 
			double s2 = rv2Ptr->getStdv() / ba2;
			double uu=u1*u1;
			double ss=s1*s1;
			double xu=correlation*u1;
			double xs=correlation*s1;
			double us=u1*s1;
			double temp = 1.082-0.004*correlation+0.204*u1+0.432*s1-0.001*correlation*correlation
				- 0.204*uu+7.728*ss+0.008*xu-1.699*xs-5.338*us
				- 19.741*ss*s2+0.135*correlation*xs+5.338*uu*s1+3.397*correlation*us;
			newCorrelation = temp*correlation;

		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			double temp = 1.037-0.042*correlation-0.182*u1+0.369*s1-0.001*correlation*correlation
				+0.182*u1*u1-1.150*s1*s1+0.084*correlation*u1;
			newCorrelation = temp*correlation;
		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			double temp = 1.040+0.015*correlation-0.176*u1+0.432*s1-0.008*correlation*correlation
				+0.176*u1*u1-1.286*s1*s1-0.137*correlation*s1;
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba1 = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba1; 
			double s1 = rv1Ptr->getStdv() / ba1;
			double ba2 = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba2; 
			double s2 = rv2Ptr->getStdv() / ba2;
			double o=u1;
			double p=s1;
			double q=u2;
			double r=s2;
			double u=o+q;
			double s=p+r;
			double oo=o*o;
			double pp=p*cov1;
			double qq=q*q;
			double rr=r*r;
			double oq=o*q;
			double pr=p*r;
			double us=oo+cov2*cov2;
			double ss=cov1*cov1+rr;
			double uc=oo*o+cov2*cov2*q;
			double sc=cov1*cov1*cov1+rr*r;
			double x=correlation;
			double xx=correlation*correlation;
			double temp;
			if(correlation > 0.0) {
				temp =1.030-0.050*x-0.056*u+0.094*s+0.009*xx-0.084*us+2.583*ss
					+0.100*x*u+0.010*x*s-0.485*u*s+0.270*oq-0.688*pr-0.011*xx*x
					+0.024*uc-10.786*sc+0.013*xx*u-0.035*xx*s+0.001*x*us-0.069*x*ss
					+1.174*us*s+0.004*oq*u+0.227*ss*u+2.783*pr*s+0.058*x*s*u
					-0.260*x*oq-0.352*x*pr-1.609*oq*s+0.194*pr*u;
			}
			else {
				temp=0.939-0.023*x+0.147*u+0.668*s+0.035*xx-0.008*us+3.146*ss
					+0.103*x*u-0.126*x*s-1.866*u*s-0.268*oq-0.304*pr+0.011*xx*x
					-0.024*uc-10.836*sc-0.013*xx*u-0.035*xx*s-0.001*x*us+0.069*x*ss
					+1.175*us*s-0.005*oq*u-0.270*ss*u+2.781*pr*s+0.058*x*u*s
					-0.259*x*oq+0.352*x*pr+1.608*oq*s-0.189*pr*u;
			}
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			double temp = 1.055-0.066*correlation-0.194*u1+0.391*s1+0.003*correlation*correlation
				+0.194*u1*u1-1.134*s1*s1+0.130*correlation*u1+0.003*correlation*s1;
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			double temp = 1.055+0.066*correlation-0.194*u1+0.391*s1+0.003*correlation*correlation
				+0.194*u1*u1-1.134*s1*s1-0.130*correlation*u1-0.003*correlation*s1;
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			double ba1 = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba1; 
			double s1 = rv1Ptr->getStdv() / ba1;
			double uu=u1*u1;
			double ss=s1*s1;
			double xu=correlation*u1;
			double xs=correlation*s1;
			double us=u1*s1;
			double uq=u1*cov2;
			double sq=s1*cov2;
			double xq=correlation*cov2;
			double temp=1.005 + 0.091*correlation + 0.285*u1+ 0.260*s1+ 0.199*cov2
				- 0.023*correlation*correlation - 0.285*uu + 8.180*ss + 0.543*cov2*cov2 - 0.181*xu
				- 1.744*xs - 0.336*xq - 5.450*us - 0.265*uq + 0.514*sq
				-19.661*ss*s1- 0.178*cov2*cov2*cov2
				+ 0.244*correlation*xs + 0.066*correlation*correlation*cov2 - 0.001*correlation*ss
				+ 5.450*uu*s1+ 0.265*uu*cov2 - 0.986*ss*cov2
				+ 0.133*sq*cov2 + 3.488*correlation*us + 0.671*correlation*uq;
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"BETA") == 0  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			double ba = rv1Ptr->getParameter4() - rv1Ptr->getParameter3();
			double u1 = ( rv1Ptr->getMean() - rv1Ptr->getParameter3() ) / ba; 
			double s1 = rv1Ptr->getStdv() / ba;
			double temp = 1.054+0.002*correlation-0.176*u1+0.366*s1-0.201*cov2
				-0.002*correlation*correlation+0.176*u1*u1-1.098*s1*s1+0.340*cov2*cov2
				-0.004*correlation*u1-0.029*s1*cov2;
			newCorrelation = temp * correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			newCorrelation = 1.031 * correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			newCorrelation = (1.029 + 0.001*correlation + 0.014*cov2 + 0.004*correlation*correlation
				+ 0.233*cov2*cov2 - 0.197*correlation*cov2)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			newCorrelation = (1.031+0.001*correlation-0.007*cov2+0.003*correlation*correlation+0.131*cov2*cov2-0.132*correlation*cov2)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.142-0.154*correlation+0.031*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.046-0.045*correlation+0.006*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.055+0.015*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			double temp = 1.055-0.066*correlation-0.194*u2+0.391*s2+0.003*correlation*correlation
				+0.194*u2*u2-1.134*s2*s2+0.130*correlation*u2+0.003*correlation*s2;
			newCorrelation = temp * correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.064-0.069*correlation+0.005*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.064+0.069*correlation+0.005*correlation*correlation)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			newCorrelation = (1.056-0.060*correlation+0.263*cov2+0.020*correlation*correlation+0.383*cov2*cov2-0.332*correlation*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			newCorrelation = (1.064+0.065*correlation-0.210*cov2+0.003*correlation*correlation+0.356*cov2*cov2-0.211*correlation*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			newCorrelation = 1.031 * correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			newCorrelation = (1.029 - 0.001*correlation + 0.014*cov2 + 0.004*correlation*correlation
				+ 0.233*cov2*cov2 + 0.197*correlation*cov2) * correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			newCorrelation = (1.031-0.001*correlation-0.007*cov2+0.003*correlation*correlation+0.131*cov2*cov2+0.132*correlation*cov2)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.142+0.154*correlation+0.031*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.046+0.045*correlation+0.006*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.055+0.015*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			double temp = 1.055+0.066*correlation-0.194*u2+0.391*s2+0.003*correlation*correlation
				+0.194*u2*u2-1.134*s2*s2-0.130*correlation*u2-0.003*correlation*s2;
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.064+0.069*correlation+0.005*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.064-0.069*correlation+0.005*correlation*correlation)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			newCorrelation = (1.056+0.060*correlation+0.263*cov2+0.020*correlation*correlation+0.383*cov2*cov2+0.332*correlation*cov2)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE1SMALLESTVALAUE") == 0  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			newCorrelation = (1.064-0.065*correlation-0.210*cov2+0.003*correlation*correlation+0.356*cov2*cov2+0.211*correlation*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			newCorrelation = (1.030 + 0.238*cov1 + 0.364*cov1*cov1) * correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			double temp = 1.026 + 0.082*correlation - 0.019*cov2 + 0.222*cov1 + 0.018*correlation*correlation
				+ 0.288*cov2*cov2 + 0.379*cov1*cov1;
			newCorrelation = (temp - 0.441*correlation*cov2 + 0.126*cov2*cov1 - 0.277*correlation*cov1);		
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			double temp = 1.029+0.056*correlation-0.030*cov2+0.225*cov1+0.012*correlation*correlation+0.174*cov2*cov2+0.379*cov1*cov1;
			newCorrelation = (temp-0.313*correlation*cov2+0.075*cov1*cov2-0.182*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.109-0.152*correlation+0.361*cov1+0.130*correlation*correlation+0.455*cov1*cov1-0.728*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.036-0.038*correlation+0.266*cov1+0.028*correlation*correlation+0.383*cov1*cov1-0.229*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.033+0.305*cov1+0.074*correlation*correlation+0.405*cov1*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba2 = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba2; 
			double s2 = rv2Ptr->getStdv() / ba2;
			double uu=u2*u2;
			double ss=s2*s2;
			double xu=correlation*u2;
			double xs=correlation*s2;
			double us=u2*s2;
			double uq=u2*cov1;
			double sq=s2*cov1;
			double xq=correlation*cov1;
			double temp=1.005 + 0.091*correlation + 0.285*u2+ 0.260*s2+ 0.199*cov1
				- 0.023*correlation*correlation - 0.285*uu + 8.180*ss + 0.543*cov1*cov1 - 0.181*xu
				- 1.744*xs - 0.336*xq - 5.450*us - 0.265*uq + 0.514*sq
				-19.661*ss*s2- 0.178*cov1*cov1*cov1
				+ 0.244*correlation*xs + 0.066*correlation*correlation*cov1 - 0.001*correlation*ss
				+ 5.450*uu*s2+ 0.265*uu*cov1 - 0.986*ss*cov1
				+ 0.133*sq*cov1 + 3.488*correlation*us + 0.671*correlation*uq;
			newCorrelation = temp * correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.056-0.060*correlation+0.263*cov1+0.020*correlation*correlation+0.383*cov1*cov1-0.332*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.056+0.060*correlation+0.263*cov1+0.020*correlation*correlation+0.383*cov1*cov1+0.332*correlation*cov1)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			double rs=cov1*cov1+cov2*cov2;
			double rc=cov1*cov1*cov1+cov2*cov2*cov2;
			double r=cov1+cov2;
			double pq=cov1*cov2;
			double temp=1.086+0.054*correlation+0.104*r-0.055*correlation*correlation+0.662*rs-0.570*correlation*r+0.203*pq;
			newCorrelation = (temp-0.020*correlation*correlation*correlation-0.218*rc-0.371*correlation*rs+0.257*correlation*correlation*r+0.141*pq*r)*correlation;
		}
		else if ( strcmp(typeRv1,"TYPE2LARGESTVALUE") == 0  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			double temp = 1.065+0.146*correlation+0.241*cov1-0.259*cov2+0.013*correlation*correlation+0.372*cov1*cov1+0.435*cov2*cov2;
			newCorrelation = (temp+0.005*correlation*cov1+0.034*cov1*cov2-0.481*correlation*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  strcmp(typeRv2,"NORMAL") == 0  ) {
			newCorrelation = (1.031 - 0.195*cov1 + 0.328*cov1*cov1) * correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  strcmp(typeRv2,"LOGNORMAL") == 0  ) {
			double temp = 1.031+0.052*correlation+0.011*cov2-0.210*cov1+0.002*correlation*correlation+0.220*cov2*cov2+0.350*cov1*cov1;
			newCorrelation = (temp+0.005*correlation*cov2+0.009*cov2*cov1-0.174*correlation*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  strcmp(typeRv2,"GAMMA") == 0  ) {
			double temp = 1.032+0.034*correlation-0.007*cov2-0.202*cov1+0.000*correlation*correlation+0.121*cov2*cov2+0.339*cov1*cov1;
			newCorrelation = (temp-0.006*correlation*cov2+0.003*cov2*cov1-0.111*correlation*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  (strcmp(typeRv2,"EXPONENTIAL") == 0 || strcmp(typeRv2,"SHIFTEDEXPONENTIAL") == 0 ) ) {
			newCorrelation = (1.147+0.145*correlation-0.271*cov1+0.010*correlation*correlation+0.459*cov1*cov1-0.467*correlation*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  (strcmp(typeRv2,"RAYLEIGH") == 0 || strcmp(typeRv2,"SHIFTEDRAYLEIGH") == 0 ) ) {
			newCorrelation = (1.047+0.042*correlation-0.212*cov1+0.353*cov1*cov1-0.136*correlation*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  strcmp(typeRv2,"UNIFORM") == 0  ) {
			newCorrelation = (1.061-0.237*cov1-0.005*correlation*correlation+0.379*cov1*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  strcmp(typeRv2,"BETA") == 0  ) {
			double ba = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba; 
			double s2 = rv2Ptr->getStdv() / ba;
			double temp = 1.054+0.002*correlation-0.176*u2+0.366*s2-0.201*cov1
				-0.002*correlation*correlation+0.176*u2*u2-1.098*s2*s2+0.340*cov1*cov1
				-0.004*correlation*u2-0.029*s2*cov1;
			newCorrelation = temp * correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  (strcmp(typeRv2,"TYPE1LARGESTVALUE") == 0 || strcmp(typeRv2,"GUMBEL") == 0)  ) {
			newCorrelation = (1.064+0.065*correlation-0.210*cov1+0.003*correlation*correlation+0.356*cov1*cov1-0.211*correlation*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  strcmp(typeRv2,"TYPE1SMALLESTVALAUE") == 0  ) {
			newCorrelation = (1.064-0.065*correlation-0.210*cov1+0.003*correlation*correlation+0.356*cov1*cov1+0.211*correlation*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  strcmp(typeRv2,"TYPE2LARGESTVALUE") == 0  ) {
			double temp = 1.065+0.146*correlation+0.241*cov2-0.259*cov1+0.013*correlation*correlation+0.372*cov2*cov2+0.435*cov1*cov1;
			newCorrelation = (temp+0.005*correlation*cov2+0.034*cov1*cov2-0.481*correlation*cov1)*correlation;
		}
		else if ( (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0)  &&  (strcmp(typeRv2,"TYPE3SMALLESTVALUE") == 0 || strcmp(typeRv2,"WEIBULL") == 0) ) {
			double temp = 1.063-0.004*correlation-0.200*(cov1+cov2)-0.001*correlation*correlation+0.337*(cov1*cov1+cov2*cov2);
			newCorrelation = (temp+0.007*(correlation*cov1+correlation*cov2)-0.007*cov1*cov2)*correlation;
		}
		/////////////////////////////////////////////////////////////////////////////////
		else {
			cerr << "Did not find the given combination of distributions in CorrelationModifier" << endl;
		}

		if(newCorrelation > 1.0) {
			newCorrelation = 0.999999999;
		}
		if(newCorrelation < -1.0) {
			newCorrelation = -0.999999999;
		}
		
		// Put the coefficient into the correlation matrix
		(*correlationMatrix)( ( rv1-1 ) , ( rv2-1 ) ) = newCorrelation;
		(*correlationMatrix)( ( rv2-1 ) , ( rv1-1 ) ) = newCorrelation;
	}

	// Here the correlation matrix should be checked for validity
	// (Whether it is close to singular or not)
	

}
