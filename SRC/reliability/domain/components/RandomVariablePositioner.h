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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariablePositioner.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#ifndef RandomVariablePositioner_h
#define RandomVariablePositioner_h

#include <ReliabilityDomainComponent.h>
#include <Information.h>
#include <DomainComponent.h>

class RandomVariablePositioner : public ReliabilityDomainComponent
{

public:

	RandomVariablePositioner(int tag,
					int RVnumber,
					DomainComponent *theObject,
					char **argv, int argc);

	~RandomVariablePositioner();
	void Print(ostream &s, int flag =0);
	int getRvNumber(void);
	int getTypeOfObject(void);
	int getTagOfObject(void);
	int getTypeOfParameterInObject(void);

	int update (double newValue); 

protected:

private:
	int tag;
	int rvNumber;
	int typeOfObject;
	int tagOfObject;
	int typeOfParameterInObject;

	Information theInfo;
	DomainComponent *theObject;
	int parameterID;


};

#endif

