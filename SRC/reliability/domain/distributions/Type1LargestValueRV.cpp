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
                                                                        
// $Revision: 1.1 $
// $Date: 2001-06-13 05:06:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type1LargestValueRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <Type1LargestValueRV.h>
#include <math.h>
#include <string.h>

Type1LargestValueRV::Type1LargestValueRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, passedMean, passedStdv, passedStartValue)
{
	type = new char[100];
	strcpy(type,"TYPE1LARGESTVALUE");
	tag = passedTag ;
	double pi = 3.14159265358979;
	double gamma = 0.5772156649;
	u = passedMean - (gamma*sqrt(6.0)*passedStdv)/pi;
	alpha = pi / ( sqrt(6.0) * passedStdv );
	startValue = passedStartValue;
}
Type1LargestValueRV::Type1LargestValueRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4, passedStartValue)
{
	type = new char[100];
	strcpy(type,"TYPE1LARGESTVALUE");
	tag = passedTag ;
	u = passedParameter1;
	alpha = passedParameter2;
	startValue = passedStartValue;
}
Type1LargestValueRV::Type1LargestValueRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, passedMean, passedStdv)
{
	type = new char[100];
	strcpy(type,"TYPE1LARGESTVALUE");
	tag = passedTag ;
	double pi = 3.14159265358979;
	double gamma = 0.5772156649;
	u = passedMean - (gamma*sqrt(6.0)*passedStdv)/pi;
	alpha = pi / ( sqrt(6.0) * passedStdv );
	startValue = getMean();
}
Type1LargestValueRV::Type1LargestValueRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4)
{
	type = new char[100];
	strcpy(type,"TYPE1LARGESTVALUE");
	tag = passedTag ;
	u = passedParameter1;
	alpha = passedParameter2;
	startValue = getMean();
}


Type1LargestValueRV::~Type1LargestValueRV()
{
}


void
Type1LargestValueRV::Print(ostream &s, int flag)
{
}


double
Type1LargestValueRV::getPDFvalue(double rvValue)
{
	return alpha * exp(-alpha*(rvValue-u)-exp(-alpha*(rvValue-u)));
}


double
Type1LargestValueRV::getCDFvalue(double rvValue)
{
	return exp(-exp(-alpha*(rvValue-u)));
}


double
Type1LargestValueRV::getInverseCDFvalue(double probValue)
{
	return (alpha*u - log(-log(probValue))) / alpha;
}


char *
Type1LargestValueRV::getType()
{
	return type;
}


double 
Type1LargestValueRV::getMean()
{
	double gamma = 0.5772156649;
	return u+gamma/alpha;
}



double 
Type1LargestValueRV::getStdv()
{
	double pi = 3.14159265358979;
	return pi/(sqrt(6)*alpha);
}


double 
Type1LargestValueRV::getStartValue()
{
	return startValue;
}


double Type1LargestValueRV::getParameter1()  {return u;}
double Type1LargestValueRV::getParameter2()  {return alpha;}
double Type1LargestValueRV::getParameter3()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
double Type1LargestValueRV::getParameter4()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
