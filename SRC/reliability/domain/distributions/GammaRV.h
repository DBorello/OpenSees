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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/GammaRV.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#ifndef GammaRV_h
#define GammaRV_h

#include <RandomVariable.h>

class GammaRV : public RandomVariable
{

public:
	GammaRV(int tag, 
			 double mean,
			 double stdv,
			 double startValue);
	GammaRV(int tag, 
			 double parameter1,
			 double parameter2,
			 double parameter3,
			 double parameter4,
			 double startValue);
	GammaRV(int tag, 
			 double mean,
			 double stdv);
	GammaRV(int tag, 
			 double parameter1,
			 double parameter2,
			 double parameter3,
			 double parameter4);
	~GammaRV();
	void Print(ostream &s, int flag =0);
	double getPDFvalue(double rvValue);
	double getCDFvalue(double rvValue);
	double getInverseCDFvalue(double probValue);
	char * getType();
	double getMean();
	double getStdv();
	double getParameter1();
	double getParameter2();
	double getParameter3();
	double getParameter4();
	double getStartValue();
	double gammaFunction(double x);
	double incompleteGammaFunction(double x, double a);


protected:

private:
	double k;
	double lambda;
	char * type;


};

#endif








