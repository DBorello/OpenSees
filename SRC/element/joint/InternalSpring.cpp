/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-14 23:01:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/InternalSpring.cpp,v $

// Written: A. Altoontash & G. Deierlein 03/02
// Revised:
//
// Purpose: This file contains the implementation for the InternalSpring class.

#include "InternalSpring.h"
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <math.h>
#include <stdlib.h>

Node *InternalSpring::theNodes[1];

InternalSpring::InternalSpring()
   :Element(0,ELE_TAG_InternalSpring), Spring(0),
   connectedExternalNodes(1), nodPtr(0), NDOF(0), DOF1(0), DOF2(0),
   Kd(0) , load(0), rForce(0)
{

}


InternalSpring::InternalSpring(int tag, int nod, int ndof, int dof1, int dof2, UniaxialMaterial &spring )
   :Element(tag,ELE_TAG_InternalSpring), Spring(0),
    connectedExternalNodes(1), nodPtr(0), NDOF(ndof), DOF1(dof1), DOF2(dof2),
	Kd(0) , load(0), rForce(0)
{
	connectedExternalNodes(0) = nod;


	Spring = spring.getCopy();


	if ( Spring == NULL )
	{
		opserr << "ERROR InternalSpring::InternalSpring(): No uniaxial material introduced ";
		exit(-1);
	}
}



// ~InternalSpring():
// 	destructor

InternalSpring::~InternalSpring()
{
	if ( Kd != NULL ) delete Kd;
	if ( load != NULL ) delete load;
	if ( rForce != NULL ) delete rForce;
	if ( Spring != NULL ) delete Spring;
}


int
InternalSpring::getNumExternalNodes(void) const
{
    return 1;
}

const ID &
InternalSpring::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
InternalSpring::getNodePtrs(void) 
{
  theNodes[0] = nodPtr;
  return theNodes;
}

int
InternalSpring::getNumDOF(void) {
  return NDOF;
}


void
InternalSpring::setDomain(Domain *theDomain)
{
  
  
  // first get the node pointers
  int Nd1 = connectedExternalNodes(0);
  nodPtr = theDomain->getNode(Nd1);
  if (nodPtr == 0) {
	    opserr << "WARNING spring2d::setDomain(): Node: ";
	    opserr << Nd1 << "does not exist in model for spring element \n" << *this;
	    return;
  }
  
  
  // now verify the number of dof at the node
  int dofNode = nodPtr->getNumberDOF();
  
  if ( dofNode != NDOF ) {
    opserr << "ERROR InternalSpring::setDomain(): node " << Nd1;
    opserr << " number of degrees of freedom does not match specified value\n";
    return;
  }	
  
  if ( DOF1 < 0 || DOF2 < 0 || DOF1 >= NDOF || DOF2 >= NDOF )
    {
      opserr << "ERROR InternalSpring::setDomain(): Incorrect degrees of freedom specified.";
      return;
    }
  
  // initialize stiffness and force vectors
  Kd = new Matrix(NDOF,NDOF);
  
  load = new Vector(NDOF);
  rForce = new Vector(NDOF);
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  this->revertToStart();
}


int
InternalSpring::commitState()
{	
	return Spring->commitState();
}

int
InternalSpring::revertToLastCommit()
{
	return Spring->revertToLastCommit();    
}

int
InternalSpring::revertToStart()
{
	return Spring->revertToStart();  
}

const Matrix &
InternalSpring::getTangentStiff(void)
{
  double S;
  double delta;
  
  const Vector &disp = nodPtr->getTrialDisp();
  
  delta = disp(DOF2) - disp(DOF1);
  
  int dummy = Spring->setTrialStrain(delta);
  S = Spring->getTangent();

  Kd->Zero();
  (*Kd)(DOF1,DOF1) = S;
  (*Kd)(DOF2,DOF2) = S;
  (*Kd)(DOF1,DOF2) = -S;
  (*Kd)(DOF2,DOF1) = -S;
  
  return (*Kd);
}

const Matrix &
InternalSpring::getInitialStiff(void)
{
  double S;
  
  S = Spring->getInitialTangent();

  Kd->Zero();
  (*Kd)(DOF1,DOF1) = S;
  (*Kd)(DOF2,DOF2) = S;
  (*Kd)(DOF1,DOF2) = -S;
  (*Kd)(DOF2,DOF1) = -S;
  
  return (*Kd);
}

void 
InternalSpring::zeroLoad(void)
{
    load->Zero();
}


int 
InternalSpring::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "InternalSpring::addLoad() - spring " << this->getTag() << ", does not handle ele loads\n"; 
  return -1;
}


int
InternalSpring::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}


const Vector &
InternalSpring::getResistingForce()
{	
    double F;
    double delta;
    
    const Vector &disp = nodPtr->getTrialDisp();
    
    delta = disp(DOF2) - disp(DOF1);
    
    Spring->setTrialStrain(delta);
    
    F = Spring->getStress();
    
    rForce->Zero();
    (*rForce)(DOF1) = -F;
    (*rForce)(DOF2) = F;
    
    return (*rForce);
}


const Vector &
InternalSpring::getResistingForceIncInertia()
{	
    this->getResistingForce();    
    
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      *rForce += this->getRayleighDampingForces();
    
    return (*rForce);
}



int
InternalSpring::sendSelf(int commitTag, Channel &theChannel)
{
	int res;
	
	int dataTag = this->getDbTag();

	
	ID data(7);
	data(0) = this->getTag();
	data(1) = connectedExternalNodes(0);
	data(2) = NDOF;
	data(3) = DOF1;
	data(4) = DOF2;
	data(5) = Spring->getClassTag();

	int matDbTag = Spring->getDbTag();
	if (matDbTag == 0) {
		matDbTag = theChannel.getDbTag();
		if (matDbTag != 0)
			Spring->setDbTag(matDbTag);
	}
	data(6) = matDbTag;
	
	res = theChannel.sendID(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING InternalSpring::sendSelf() - " << this->getTag() << 
		  " failed to send ID\n";
		return -1;
	}
	
	res = Spring->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr << "WARNING InternalSpring::sendSelf() - " << this->getTag() << 
		  " failed to send its Material\n";
		return -3;
	}
	
	return 0;
}

int
InternalSpring::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	static ID data(7);
	res = theChannel.recvID(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING InternalSpring::recvSelf() - failed to receive Vector\n";
		return -1;
	}
	
	this->setTag((int)data(0));


	connectedExternalNodes(0) = data(1);
	NDOF = data(2);
	DOF1 = data(3);
	DOF2 = data(4);
	int matClassTag = data(5);
	int matDbTag = data(6);
	
	if ((Spring == 0) || (Spring->getClassTag() != matClassTag)) {
		if (Spring != 0)
		  delete Spring;
		
		// create a new material object
		Spring = theBroker.getNewUniaxialMaterial(matClassTag);
		if (Spring == 0) {
		  opserr << "WARNING InternalSpring::recvSelf() - " << this->getTag() << 
		    " failed to get a blank Material of type " << matClassTag << endln;
		  return -3;
		}
	}
	
	Spring->setDbTag(matDbTag); // note: we set the dbTag before we receive the material
	res = Spring->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
	  opserr << "WARNING InternalSpring::recvSelf() - " << this->getTag() << 
	    " failed to receive its Material\n";
	  return -3;
	}
	
	return 0;
}


int
InternalSpring::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	return 0;
}

void
InternalSpring::Print(OPS_Stream &s, int flag)
{
    // compute current state
    this->getResistingForce();
    
    s << "Element: " << this->getTag(); 
    s << " type: spring2d  Node: " << connectedExternalNodes(0);
	s << " DOF1 : " << DOF1 << " DOF2 : " << DOF2;
    s << " \nresisting Force: " << rForce;
}


Response *InternalSpring::setResponse(char **argv, int argc, Information &eleInformation)
{
  if (strcmp(argv[0],"moment") == 0 || strcmp(argv[0],"force") == 0 || strcmp(argv[0],"panelMoment") == 0 )
    return new ElementResponse(this, 1, 0.0);

  else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	   strcmp(argv[0],"deformation") == 0)
    return new ElementResponse(this, 2, 0.0);

  // a material quantity    
  else if (strcmp(argv[0],"material") == 0)
    return Spring->setResponse(&argv[1], argc-1, eleInformation);
  
  else
    return 0;

}

int InternalSpring::getResponse(int responseID, Information &eleInformation)
{
	switch (responseID) {
    case 1:
		return eleInformation.setDouble( Spring->getStress() );
      
    case 2:
		return eleInformation.setDouble( Spring->getStrain() );

    default:
		return 0;
	}
}




