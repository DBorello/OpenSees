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
// $Date: 2001-06-14 08:06:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/LaplaceRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <LaplaceRV.h>
#include <math.h>
#include <string.h>

LaplaceRV::LaplaceRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, passedMean, passedStdv, passedStartValue)
{
	type = new char[100];
	strcpy(type,"LAPLACE");
	tag = passedTag ;
	alpha = passedMean;
	beta = sqrt(2.0)/passedStdv;
	startValue = passedStartValue;
}
LaplaceRV::LaplaceRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4, passedStartValue)
{
	type = new char[100];
	strcpy(type,"LAPLACE");
	tag = passedTag ;
	alpha = passedParameter1;
	beta = passedParameter2;
	startValue = passedStartValue;
}
LaplaceRV::LaplaceRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, passedMean, passedStdv)
{
	type = new char[100];
	strcpy(type,"LAPLACE");
	tag = passedTag ;
	alpha = passedMean;
	beta = sqrt(2.0)/passedStdv;
	startValue = getMean();
}
LaplaceRV::LaplaceRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4)
{
	type = new char[100];
	strcpy(type,"LAPLACE");
	tag = passedTag ;
	alpha = passedParameter1;
	beta = passedParameter2;
	startValue = getMean();
}


LaplaceRV::~LaplaceRV()
{
}


void
LaplaceRV::Print(ostream &s, int flag)
{
}


double
LaplaceRV::getPDFvalue(double rvValue)
{
	return 0.5*beta*exp(-beta*fabs(rvValue-alpha));
}


double
LaplaceRV::getCDFvalue(double rvValue)
{
	double result;
	if (rvValue < alpha)  {
		result = 0.5*exp(-beta*fabs(rvValue-alpha));
	}
	else  {
		result = 1 - 0.5*exp(-beta*fabs(rvValue-alpha));
	}
	return result;
}


double
LaplaceRV::getInverseCDFvalue(double rvValue)
{
	return 0.0;
}


char *
LaplaceRV::getType()
{
	return type;
}


double 
LaplaceRV::getMean()
{
	return alpha;
}



double 
LaplaceRV::getStdv()
{
	return sqrt(2.0)/beta;
}


double 
LaplaceRV::getStartValue()
{
	return startValue;
}

double LaplaceRV::getParameter1()  {return alpha;}
double LaplaceRV::getParameter2()  {return beta;}
double LaplaceRV::getParameter3()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
double LaplaceRV::getParameter4()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
