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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-08-01 00:25:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ExponentialRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <ExponentialRV.h>
#include <math.h>
#include <string.h>

ExponentialRV::ExponentialRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, passedMean, passedStdv, passedStartValue)
{
	// Note: this constructor is void.
	cerr << "WARNING: This type of random variable is not uniquely defined by mean and stdv." << endl;
}
ExponentialRV::ExponentialRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4, passedStartValue)
{
	type = new char[100];
	strcpy(type,"EXPONENTIAL");
	tag = passedTag ;
	lambda = passedParameter1;
	startValue = passedStartValue;
}
ExponentialRV::ExponentialRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, passedMean, passedStdv)
{
	// Note: this constructor is void.
	cerr << "WARNING: This type of random variable is not uniquely defined by mean and stdv." << endl;
}
ExponentialRV::ExponentialRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4)
{
	type = new char[100];
	strcpy(type,"EXPONENTIAL");
	tag = passedTag ;
	lambda = passedParameter1;
	startValue = getMean();
}


ExponentialRV::~ExponentialRV()
{
  if (type != 0)
    delete [] type;
}


void
ExponentialRV::Print(ostream &s, int flag)
{
}


double
ExponentialRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = lambda * exp(-lambda * rvValue);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ExponentialRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = 1 - exp(-lambda*rvValue);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ExponentialRV::getInverseCDFvalue(double probValue)
{
	return -log(1.0-probValue)/lambda;
}


char *
ExponentialRV::getType()
{
	return type;
}


double 
ExponentialRV::getMean()
{
	return 1/lambda;
}



double 
ExponentialRV::getStdv()
{
	return 1/lambda;
}


double 
ExponentialRV::getStartValue()
{
	return startValue;
}


double ExponentialRV::getParameter1()  {return lambda;}
double ExponentialRV::getParameter2()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
double ExponentialRV::getParameter3()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
double ExponentialRV::getParameter4()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
