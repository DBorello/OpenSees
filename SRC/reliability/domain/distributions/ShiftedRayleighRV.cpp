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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ShiftedRayleighRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <ShiftedRayleighRV.h>
#include <math.h>
#include <string.h>

ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, passedMean, passedStdv, passedStartValue)
{
	type = new char[100];
	strcpy(type,"SHIFTEDRAYLEIGH");
	tag = passedTag;
	double pi = 3.14159265358979;
	u = 2.0 * passedStdv / sqrt(3.0*pi+4.0);
	x0 = passedMean + passedStdv*sqrt(pi) / sqrt(3.0*pi+4.0);
	startValue = passedStartValue;
}
ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4, passedStartValue)
{
	type = new char[100];
	strcpy(type,"SHIFTEDRAYLEIGH");
	tag = passedTag ;
	u = passedParameter1;
	x0 = passedParameter2;
	startValue = passedStartValue;
}
ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, passedMean, passedStdv)
{
	type = new char[100];
	strcpy(type,"SHIFTEDRAYLEIGH");
	tag = passedTag;
	double pi = 3.14159265358979;
	u = 2.0 * passedStdv / sqrt(3.0*pi+4.0);
	x0 = passedMean + passedStdv*sqrt(pi) / sqrt(3.0*pi+4.0);
	startValue = getMean();
}
ShiftedRayleighRV::ShiftedRayleighRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4)
{
	type = new char[100];
	strcpy(type,"SHIFTEDRAYLEIGH");
	tag = passedTag ;
	u = passedParameter1;
	x0 = passedParameter2;
	startValue = getMean();
}


ShiftedRayleighRV::~ShiftedRayleighRV()
{
}


void
ShiftedRayleighRV::Print(ostream &s, int flag)
{
}


double
ShiftedRayleighRV::getPDFvalue(double rvValue)
{
	double result;
	if ( x0 < rvValue ) {
		result = 2.0*(rvValue-x0)/(u*u) * exp(-pow(((rvValue-x0)/u),2));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ShiftedRayleighRV::getCDFvalue(double rvValue)
{
	double result;
	if ( x0 < rvValue ) {
		result = 1 - exp(-pow(((rvValue-x0)/u),2));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ShiftedRayleighRV::getInverseCDFvalue(double probValue)
{
	return x0 + u * sqrt(-log(1-probValue));
}


char *
ShiftedRayleighRV::getType()
{
	return type;
}


double 
ShiftedRayleighRV::getMean()
{
	double pi = 3.14159265358979;
	return x0 + 0.5 * u * sqrt(pi);
}



double 
ShiftedRayleighRV::getStdv()
{
	double pi = 3.14159265358979;
	return 0.5 * u * sqrt(4.0-pi);
}


double 
ShiftedRayleighRV::getStartValue()
{
	return startValue;
}


double ShiftedRayleighRV::getParameter1()  {return u;}
double ShiftedRayleighRV::getParameter2()  {return x0;}
double ShiftedRayleighRV::getParameter3()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
double ShiftedRayleighRV::getParameter4()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
