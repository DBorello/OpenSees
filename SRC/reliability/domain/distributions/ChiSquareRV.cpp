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
                                                                        
// $Revision: 1.4 $
// $Date: 2001-08-01 18:06:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ChiSquareRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <ChiSquareRV.h>
#include <GammaRV.h>
#include <math.h>
#include <string.h>

ChiSquareRV::ChiSquareRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, passedMean, passedStdv, passedStartValue)
{
	tag = passedTag ;
	nu = 0.5*passedMean;
	startValue = passedStartValue;
}
ChiSquareRV::ChiSquareRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4, passedStartValue)
{
	tag = passedTag ;
	nu = passedParameter1;
	startValue = passedStartValue;
}
ChiSquareRV::ChiSquareRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, passedMean, passedStdv)
{
	tag = passedTag ;
	nu = 0.5*passedMean;
	startValue = getMean();
}
ChiSquareRV::ChiSquareRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4)
{
	tag = passedTag ;
	nu = passedParameter1;
	startValue = getMean();
}


ChiSquareRV::~ChiSquareRV()
{
}


void
ChiSquareRV::Print(ostream &s, int flag)
{
}


double
ChiSquareRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
		double a = aGammaRV->gammaFunction(0.5*nu);
		result = pow(0.5*rvValue,0.5*nu-1.0)*exp(-0.5*rvValue)/(2.0*a);
		delete aGammaRV;
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ChiSquareRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
		double a = aGammaRV->incompleteGammaFunction(0.5*nu,0.5*rvValue);
		double b = aGammaRV->gammaFunction(0.5*nu);
		result = a/b;
		delete aGammaRV;
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ChiSquareRV::getInverseCDFvalue(double rvValue)
{
	return 0.0;
}


const char *
ChiSquareRV::getType()
{
	return "CHISQUARE";
}


double 
ChiSquareRV::getMean()
{
	return 2*nu;
}



double 
ChiSquareRV::getStdv()
{
	return sqrt(2*nu);
}


double 
ChiSquareRV::getStartValue()
{
	return startValue;
}

double ChiSquareRV::getParameter1()  {return nu;}
double ChiSquareRV::getParameter2()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
double ChiSquareRV::getParameter3()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
double ChiSquareRV::getParameter4()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
