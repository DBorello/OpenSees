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
// $Date: 2001-06-13 05:06:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/CorrelationCoefficient.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#ifndef CorrelationCoefficient_h
#define CorrelationCoefficient_h

#include <ReliabilityDomainComponent.h>

class CorrelationCoefficient : public ReliabilityDomainComponent
{

public:
	CorrelationCoefficient(	int tag, int rv1, int rv2, double correlation);
	~CorrelationCoefficient();

	int		getRv1();
	int		getRv2();
	double	getCorrelation();

	void Print(ostream &s, int flag =0);


protected:

private:
	int tag;
	int rv1;
	int rv2;
	double correlation;


};

#endif
