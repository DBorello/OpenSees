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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:01:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/MP_Joint2D.cpp,v $

// Written: Arash
// Created: 08/01

// Purpose: This file contains the implementation of class MP_TimeVary.


#include <MP_Joint2D.h>

#include <stdlib.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
 
// constructor for FEM_ObjectBroker
MP_Joint2D::MP_Joint2D()
:MP_Constraint( 0 , CNSTRNT_TAG_MP_Joint2D ),thisDomain(0),
 nodeRetained(0),nodeConstrained(0), MainDOF(0), constraint(0),
 constrDOF(0),retainDOF(0),dbTag1(0), dbTag2(0), RetainedNode(0),
 ConstrainedNode(0), LargeDisplacement(0), Length0(0.0)
{
    
}


// general constructor for ModelBuilder
MP_Joint2D::MP_Joint2D(Domain *theDomain, int tag, int nodeRetain, int nodeConstr,
		int Maindof, int LrgDsp )
:MP_Constraint( tag , CNSTRNT_TAG_MP_Joint2D ), thisDomain(theDomain),
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), MainDOF(Maindof),
 constraint(0), constrDOF(0), retainDOF(0), dbTag1(0), dbTag2(0),
 RetainedNode(0), ConstrainedNode(0), LargeDisplacement( LrgDsp ), Length0(0.0)
{

  this->setTag(tag);

  // get node pointers of constrainted and retained nodes
  ConstrainedNode = theDomain->getNode(nodeConstrained);
  if (ConstrainedNode == NULL)
    {
      opserr << "MP_Joint2D::MP_Joint2D: nodeConstrained: ";
      opserr << nodeConstrained << "does not exist in model\n";
      exit(0);}
  
  RetainedNode = theDomain->getNode(nodeRetained);
  if (RetainedNode == NULL)
    {
      opserr << "MP_Joint2D::MP_Joint2D: nodeRetained: ";
      opserr << nodeRetained << "does not exist in model\n";
      exit(0);}
  
  // check for proper degrees of freedom
  int RnumDOF = RetainedNode->getNumberDOF();
  int CnumDOF = ConstrainedNode->getNumberDOF();
  if (RnumDOF != 4 || CnumDOF != 3 ){
    opserr << "MP_Joint2D::MP_Joint2D - mismatch in numDOF\n DOF not supported by this type of constraint\n";
    return;
  }
  
  // check the main degree of freedom. It should be aone of the rotational degrees of freedom
  if ( MainDOF != 2 && MainDOF!=3 ) {
    opserr << "MP_Joint2D::MP_Joint2D - Wrong main degree of freedom\n";
    return;
  }
  
  
  // check for proper dimensions of coordinate space
  const Vector &crdR = RetainedNode->getCrds();
  int dimR = crdR.Size();
  const Vector &crdC = ConstrainedNode->getCrds();
  int dimC = crdC.Size();
  
  if (dimR != 2 || dimC != 2 ){
    opserr << "MP_Joint2D::MP_Joint2D - mismatch in dimnesion\n dimension not supported by this type of constraint\n";
    return;
  }

  
  // allocate the constranted and retained id's
    constrDOF = new ID(CnumDOF-1);
    retainDOF = new ID(RnumDOF-1); 
 
    if (constrDOF == NULL || retainDOF == NULL ) { 
      opserr << "MP_Joint2D::MP_Joint2D - ran out of memory \ncan not generate ID for nodes\n";
      exit(-1);
    }
    
    (*constrDOF)(0) = 0;
    (*constrDOF)(1) = 1;
    
    (*retainDOF)(0) = 0;
    (*retainDOF)(1) = 1;
    (*retainDOF)(2) = MainDOF;
    
    
    // allocate the constraint matrix
    constraint = new Matrix( CnumDOF-1 , RnumDOF-1 );
    if (constraint == NULL ) { 
	opserr << "MP_Joint2D::MP_Joint2D - ran out of memory 2\n";
	exit(-1);
    }  

    // calculate constraint matrix
    double deltaX = crdC(0) - crdR(0);
    double deltaY = crdC(1) - crdR(1);
    
    Length0 = sqrt( deltaX*deltaX + deltaY*deltaY );
    if ( Length0 <= 1.0e-12 ) { 
      opserr << "MP_Joint2D::MP_Joint2D - The constraint length is zero\n";
    }  
    
    (*constraint) (0,0) = 1.0 ;
    (*constraint) (1,1) = 1.0 ;
    (*constraint) (0,2) = -deltaY ;
    (*constraint) (1,2) = deltaX ;
}



MP_Joint2D::~MP_Joint2D()
{
    // invoke the destructor on the matrix and the two ID objects
    if (constraint != NULL)
	delete constraint;
    if (constrDOF != NULL)
	delete constrDOF;
    if (retainDOF != NULL)
	delete retainDOF;    
}


int
MP_Joint2D::getNodeRetained(void) const
{
    // return id of retained node
    return nodeRetained;
}

int
MP_Joint2D::getNodeConstrained(void) const
{
    // return id of constrained node    
    return nodeConstrained;
}


const ID &
MP_Joint2D::getConstrainedDOFs(void) const
{
    if (constrDOF == NULL) {
	opserr << "MP_Joint2D::getConstrainedDOF - no ID was set, ";
	opserr << "was recvSelf() ever called? or subclass incorrect?\n";	
	exit(-1);
    }

    // return the ID corresponding to constrained DOF of Ccr
    return (*constrDOF);    
}


const ID &
MP_Joint2D::getRetainedDOFs(void) const
{
    if (retainDOF == NULL) {
	opserr << "MP_Joint2D::getRetainedDOFs - no ID was set\n ";
	opserr << "was recvSelf() ever called? or subclass incorrect?\n";		
	exit(-1);
    }

    // return the ID corresponding to retained DOF of Ccr
    return (*retainDOF);    
}


int 
MP_Joint2D::applyConstraint(double timeStamp)
{
    if ( LargeDisplacement != 0 )
	{
		// calculate the constraint at this moment

		// get the coordinates of the two nodes - check dimensions are the same FOR THE MOMENT
		const Vector &crdR = RetainedNode->getCrds();
		const Vector &crdC = ConstrainedNode->getCrds();

		// get commited displacements of nodes to get updated coordinates
		const Vector &dispR = RetainedNode->getDisp();
		const Vector &dispC = ConstrainedNode->getDisp();

		double deltaX = dispC(0) + crdC(0) - dispR(0) - crdR(0);
		double deltaY = dispC(1) + crdC(1) - dispR(1) - crdR(1);

		constraint->Zero();

		(*constraint) (0,0) = 1.0 ;
		(*constraint) (1,1) = 1.0 ;
		(*constraint) (0,2) = -deltaY ;
		(*constraint) (1,2) = deltaX ;
	}

	return 0;
}



bool
MP_Joint2D::isTimeVarying(void) const
{
    if ( LargeDisplacement != 0 ) return true;

	return false;
}


int MP_Joint2D::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}

int MP_Joint2D::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
	return 0;
}


const Matrix &MP_Joint2D::getConstraint(void)
{
    if (constraint == 0) {
	opserr << "MP_Joint2D::getConstraint - no Matrix was set\n";
	exit(-1);
    }    

    if ( LargeDisplacement == 2 )
	{
		// Length correction
		// to correct the trial displacement

		// get the coordinates of the two nodes - check dimensions are the same FOR THE MOMENT
		const Vector &crdR = RetainedNode->getCrds();
		const Vector &crdC = ConstrainedNode->getCrds();

		// get commited displacements of nodes to get updated coordinates
		const Vector &dispR = RetainedNode->getTrialDisp();
		const Vector &dispC = ConstrainedNode->getTrialDisp();

		double deltaX = dispC(0) + crdC(0) - dispR(0) - crdR(0);
		double deltaY = dispC(1) + crdC(1) - dispR(1) - crdR(1);


		Vector Direction(2);
		Direction(0) = deltaX;
		Direction(1) = deltaY;
		double NewLength = Direction.Norm();
		if ( NewLength < 1e-12 ) opserr << "MP_Joint2D::applyConstraint : length of rigid link is too small or zero"; 
		Direction = Direction * (Length0/NewLength);		// correct the length
		// find new displacements of the constrainted node
	
		Vector NewLocation(3);
		NewLocation(0) = Direction(0) + dispR(0) + crdR(0) - crdC(0);
		NewLocation(1) = Direction(1) + dispR(1) + crdR(1) - crdC(1);
		NewLocation(2) = dispC(2);
		int dummy = ConstrainedNode->setTrialDisp( NewLocation );
		// end of length correction procedure

	}

    // return the constraint matrix Ccr
    return (*constraint);
}
    
void MP_Joint2D::Print(OPS_Stream &s, int flag )
{
    s << "MP_Joint2D: " << this->getTag() << "\n";
    s << "\tNode Constrained: " << nodeConstrained;
    s << " node Retained: " << nodeRetained ;
    if (constrDOF != 0)
	s << " constrained dof: " << *constrDOF;    
    if (retainDOF != 0)
	s << " retained dof: " << *retainDOF;        
    if (constraint != 0)
	s << " constraint matrix: " << *constraint << "\n";

}
