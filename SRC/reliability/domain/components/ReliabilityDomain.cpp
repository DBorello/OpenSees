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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ReliabilityDomain.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <ReliabilityDomain.h>

#include <CorrelationCoefficient.h>
#include <RandomVariable.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
#include <ArrayOfTaggedObjects.h>


ReliabilityDomain::ReliabilityDomain()
{
	theRandomVariablesPtr = new ArrayOfTaggedObjects (256);
	theCorrelationCoefficientsPtr = new ArrayOfTaggedObjects (256);
	theLimitStateFunctionsPtr = new ArrayOfTaggedObjects (256);
	theRandomVariablePositionersPtr = new ArrayOfTaggedObjects (256);
	tagOfActiveLimitStateFunction = 1;

}

ReliabilityDomain::~ReliabilityDomain()
{
	if (!theRandomVariablesPtr)
		delete theRandomVariablesPtr;
	if (!theCorrelationCoefficientsPtr)
		delete theCorrelationCoefficientsPtr;
	if (!theLimitStateFunctionsPtr)
		delete theLimitStateFunctionsPtr;
	if (!theRandomVariablePositionersPtr)
		delete theRandomVariablePositionersPtr;
}


bool
ReliabilityDomain::addRandomVariable(RandomVariable *theRandomVariable)
{
	bool result = theRandomVariablesPtr->addComponent(theRandomVariable);
	return result;
}

bool
ReliabilityDomain::addCorrelationCoefficient(CorrelationCoefficient *theCorrelationCoefficient)
{
	bool result = theCorrelationCoefficientsPtr->addComponent(theCorrelationCoefficient);
	return result;
}

bool
ReliabilityDomain::addLimitStateFunction(LimitStateFunction *theLimitStateFunction)
{
	bool result = theLimitStateFunctionsPtr->addComponent(theLimitStateFunction);
	return result;
}

bool
ReliabilityDomain::addRandomVariablePositioner(RandomVariablePositioner *theRandomVariablePositioner)
{
	bool result = theRandomVariablePositionersPtr->addComponent(theRandomVariablePositioner);
	return result;
}




RandomVariable *
ReliabilityDomain::getRandomVariablePtr(int tag)
{
	TaggedObject *theComponent = theRandomVariablesPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	RandomVariable *result = (RandomVariable *) theComponent;
	return result;
}


CorrelationCoefficient * 
ReliabilityDomain::getCorrelationCoefficientPtr(int tag)
{
	TaggedObject *theComponent = theCorrelationCoefficientsPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	CorrelationCoefficient *result = (CorrelationCoefficient *) theComponent;
	return result;
}


LimitStateFunction *
ReliabilityDomain::getLimitStateFunctionPtr(int tag)
{
	TaggedObject *theComponent = theLimitStateFunctionsPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	LimitStateFunction *result = (LimitStateFunction *) theComponent;
	return result;
}


RandomVariablePositioner *
ReliabilityDomain::getRandomVariablePositionerPtr(int tag)
{
	TaggedObject *theComponent = theRandomVariablePositionersPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	RandomVariablePositioner *result = (RandomVariablePositioner *) theComponent;
	return result;
}


int
ReliabilityDomain::getNumberOfRandomVariables()
{
	return theRandomVariablesPtr->getNumComponents();
}


int
ReliabilityDomain::getNumberOfCorrelationCoefficients()
{
	return theCorrelationCoefficientsPtr->getNumComponents();
}


int
ReliabilityDomain::getNumberOfLimitStateFunctions()
{
	return theLimitStateFunctionsPtr->getNumComponents();
}



int
ReliabilityDomain::getNumberOfRandomVariablePositioners()
{
	return theRandomVariablePositionersPtr->getNumComponents();
}





int
ReliabilityDomain::getTagOfActiveLimitStateFunction()
{
	return tagOfActiveLimitStateFunction;
}


void
ReliabilityDomain::setTagOfActiveLimitStateFunction(int passedTag)
{
	tagOfActiveLimitStateFunction = passedTag;
}

