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
// $Date: 2001-06-13 05:06:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariablePositioner.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <RandomVariablePositioner.h>


RandomVariablePositioner::RandomVariablePositioner (int passedTag,
		int passedRVnumber,
		DomainComponent *object,
		char **argv, int argc)
:ReliabilityDomainComponent(passedTag, 14728)
{
	tag = passedTag;
	rvNumber = passedRVnumber;
	theObject = object;

	if (theObject)
		parameterID = theObject->setParameter (argv, argc, theInfo);

	if (parameterID < 0)
		cerr << "RandomVariablePositioner::RandomVariablePositioner "<< tag <<" -- unable to set parameter" << endl;
}



RandomVariablePositioner::~RandomVariablePositioner()
{
}

int
RandomVariablePositioner::update (double newValue)
{
	theInfo.theDouble = newValue;

	if (parameterID >= 0)
		return theObject->updateParameter (parameterID, theInfo);
	else
		return -1;
}

void
RandomVariablePositioner::Print(ostream &s, int flag)  
{
}


int 
RandomVariablePositioner::getRvNumber(void)
{
	return rvNumber;
}


int 
RandomVariablePositioner::getTypeOfObject(void)
{
	return typeOfObject;
}


int 
RandomVariablePositioner::getTagOfObject(void)
{
	return tagOfObject;
}


int
RandomVariablePositioner::getTypeOfParameterInObject(void)
{
	return typeOfParameterInObject;
}
