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
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/SimpleJoint2D.cpp,v $

// Written: AAA 03/02
// Revised:
// SimpleJoint2D.cpp: implementation of the SimpleJoint2D class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <MP_Constraint.h>
#include <MP_SimpleJoint2D.h>
#include <ElasticMaterial.h>
#include <ElementResponse.h>
#include <SimpleJoint2D.h>
#include <InternalSpring.h>

Matrix SimpleJoint2D::K(12,12);
Vector SimpleJoint2D::V(12);
Node * SimpleJoint2D::theNodes[5];

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SimpleJoint2D::SimpleJoint2D()
  :Element(0, ELE_TAG_SimpleJoint2D ), 
  ExternalNodes(4), InternalConstraints(4), 
  end1Ptr(0), end2Ptr(0), end3Ptr(0), end4Ptr(0),IntNodePtr(0),
  IntNode(0), TheDomain(0), RotElemtag(0), RotElemPtr(0),  
  numDof(0), nodeDbTag(0), dofDbTag(0)
{

}

SimpleJoint2D::SimpleJoint2D(int tag, int nd1, int nd2, int nd3, int nd4,
			     UniaxialMaterial &Spring , Domain *theDomain, 
			     int IntNodeTag, int RotEleTag, int LrgDisp)
  :Element(tag, ELE_TAG_SimpleJoint2D ), 
  ExternalNodes(4), InternalConstraints(4), 
  end1Ptr(0), end2Ptr(0), end3Ptr(0), end4Ptr(0),IntNodePtr(0),
  IntNode(IntNodeTag), TheDomain(0), RotElemtag(RotEleTag), RotElemPtr(0),  
  numDof(0), nodeDbTag(0), dofDbTag(0)
{
  numDof  = 12;

  K.Zero();
  V.Zero();

  TheDomain = theDomain;
  if( TheDomain==NULL ) {
    opserr << "WARNING SimpleJoint2D(): Specified domain does not exist";
    opserr << "Domain = 0\n";
    return;
  }

  // Save external node id's
  ExternalNodes(0) = nd1;
  ExternalNodes(1) = nd2;
  ExternalNodes(2) = nd3;
  ExternalNodes(3) = nd4;
  
  // Set the external node pointers
  end1Ptr = TheDomain->getNode(nd1);
  end2Ptr = TheDomain->getNode(nd2);	
  end3Ptr = TheDomain->getNode(nd3);
  end4Ptr = TheDomain->getNode(nd4);
  
  // check domain for existence of external nodes
  if (end1Ptr == 0) {
    opserr << "WARNING SimpleJoint2D::setDomain(): Nd1: ";
    opserr << nd1 << "does not exist in model for element \n" << *this;
    return;
  }
    
  if (end2Ptr == 0) {
    opserr << "WARNING SimpleJoint2D::setDomain(): Nd2: ";
    opserr << nd2 << "does not exist in model for element\n" << *this;
    return;
  }	
  
  if (end3Ptr == 0) {
    opserr << "WARNING SimpleJoint2D::setDomain(): Nd3: ";
    opserr << nd3 << "does not exist in model for element \n" << *this;
    return;
  }
    
  if (end4Ptr == 0) {
    opserr << "WARNING SimpleJoint2D::setDomain(): Nd4: ";
    opserr << nd4 << "does not exist in model for element\n" << *this;
    return;
  }	

  // 
  const Vector &end1Crd = end1Ptr->getCrds();
  const Vector &end2Crd = end2Ptr->getCrds();	
  const Vector &end3Crd = end3Ptr->getCrds();
  const Vector &end4Crd = end4Ptr->getCrds();

  // check for a two dimensional domain, since this element supports only two dimensions
  int dimNd1 = end1Crd.Size();
  int dimNd2 = end2Crd.Size();
  int dimNd3 = end3Crd.Size();
  int dimNd4 = end4Crd.Size();

  if (dimNd1 != 2 && dimNd2 != 2 && dimNd3 != 2 && dimNd4 != 2 ) {
    opserr << "WARNING SimpleJoint2D::setDomain(): has incorrect space dimension \n";
    opserr << "                                    space dimension not supported by SimpleJoint2D";
    return;
  }
	
  // now verify the number of dof at node ends
  int dofNd1 = end1Ptr->getNumberDOF();
  int dofNd2 = end2Ptr->getNumberDOF();	
  int dofNd3 = end3Ptr->getNumberDOF();
  int dofNd4 = end4Ptr->getNumberDOF();

  if (dofNd1 != 3 && dofNd2 != 3 && dofNd3 != 3 && dofNd4 != 3 ) {
    opserr << "WARNING SimpleJoint2D::SimpleJoint2D: has incorrect degrees of freedom \n";
    opserr << "                                    DOF not supported by SimpleJoint2D";
    return;
  }
  
  // check if joint diagonals lengths are not zero
  Vector Center1(end1Crd);
  Vector Center2(end2Crd);
  Center1 = Center1 - end3Crd;
  Center2 = Center2 - end4Crd;
	
  double L1 = Center1.Norm();
  double L2 = Center2.Norm();
  
  if( Center1.Norm()<1e-6  || Center2.Norm()<1e-6 ) {
    opserr << "WARNING SimpleJoint2D::(): zero length\n";
    return;	
  }
	
  // check if nodes are not located at the same coordinate and construct
  // a parallelogram
  Center1 = end1Crd + end3Crd;
  Center2 = end2Crd + end4Crd;
  
  Center1 = 0.5 * Center1;
  Center2 = 0.5 * Center2;
  
  Vector Center3(Center2);
  Center3 = Center3 - Center1;

  if ( Center3.Norm() > 1e-6 ) {
    opserr << "WARNING SimpleJoint2D::(): can not construct a paralelogram over external nodes\n";
    return;	
  }
	
  // Get the largest node number and compare it with begining MAX_NUMS
  //  int MaxNodeInDomain = getMaxNodeNum( TheDomain ) +1;
  // IntNode = ( MaxNodeInDomain > MAX_NUMS ) ? MaxNodeInDomain : MAX_NUMS;
  
  // Generate internal node
  // I am not using the node copy for future purposes, specially when 
  // the new node coordinates are different	
  IntNodePtr = new Node ( IntNode , 4, Center1(0) , Center1(1) );
  if ( IntNodePtr == NULL ) 
    {
      opserr << "SimpleJoint2D::SimpleJoint2D - Unable to generate new nodes , out of memory\n";
      exit(-1);
    } else {
      if( TheDomain->addNode( IntNodePtr ) == false )		// add intenal nodes to domain
	opserr << "SimpleJoint2D::SimpleJoint2D - unable to add internal nodeto domain\n";
    }
  
  // Generate and add constraints to domain

  // get the constraint numbers
  int startMPtag = theDomain->getNumMPs();
  for ( int i=0 ; i<4 ; i++ ) InternalConstraints(i) = startMPtag + i ;
  
  // create MP_SimpleJoint constraint node 1
  if ( addMP_SimpleJoint( TheDomain, InternalConstraints(0), IntNode, ExternalNodes(0), 2, LrgDisp ) != 0) {
    opserr << "WARNING SimpleJoint2D::SimpleJoint2D(): can not generate ForJoint MP at node 1\n";
    return;
  }

  // create MP_SimpleJoint constraint node 2
  if ( addMP_SimpleJoint( TheDomain, InternalConstraints(1), IntNode, ExternalNodes(1), 3, LrgDisp ) != 0) {
    opserr << "WARNING SimpleJoint2D::SimpleJoint2D(): can not generate ForJoint MP at node 2\n";
    return;
  }

  // create MP_SimpleJoint constraint node 3
  if ( addMP_SimpleJoint( TheDomain, InternalConstraints(2), IntNode, ExternalNodes(2), 2, LrgDisp ) != 0) {
    opserr << "WARNING SimpleJoint2D::SimpleJoint2D(): can not generate ForJoint MP at node 3\n";
    return;
  }
  
  // create MP_SimpleJoint constraint node 4
  if ( addMP_SimpleJoint( TheDomain, InternalConstraints(3), IntNode, ExternalNodes(3), 3, LrgDisp ) != 0) {
    opserr << "WARNING SimpleJoint2D::SimpleJoint2D(): can not generate ForJoint MP at node 4\n";
    return;
  }

  // add rotational spring
  //   RotElemtag = getMaxEleNum( TheDomain ) + 1;
  // RotElemtag = ( RotElemtag > MAX_NUMS ) ? RotElemtag : MAX_NUMS;
  
  RotElemPtr = new InternalSpring(RotElemtag,IntNode,4,2,3, Spring);
  
  if ( RotElemPtr == NULL ) {
    opserr << "SimpleJoint2D::SimpleJoint2D - failed to genrate Rotational spring , out of memory\n";
  } else {
    // add the zero length spring to the domain
    if ( theDomain->addElement( RotElemPtr ) == false) {
      opserr << "SimpleJoint2D::SimpleJoint2D - Unable to add Rotational spring to domain\n";
      delete RotElemPtr;
    }
  }
}


SimpleJoint2D::~SimpleJoint2D()
{
	if ( TheDomain != NULL)
	{
		if ( RotElemPtr != NULL )
		{
			int rotelemtag = RotElemPtr->getTag();
			TheDomain->removeElement( rotelemtag );
			delete RotElemPtr;
		}
		MP_Constraint *Temp_MP;
		for ( int i=0 ; i < 4 ; i++ )
		{
			Temp_MP = TheDomain->getMP_Constraint( InternalConstraints(i) );
			
			if ( Temp_MP != NULL )
			{
				TheDomain->removeMP_Constraint( InternalConstraints(i) );
				delete Temp_MP;
			}
		}
		if ( IntNodePtr != NULL )
		{
			int intnodetag = IntNodePtr->getTag();
			TheDomain->removeNode( intnodetag );
			delete IntNodePtr;
		}
	}
}


void SimpleJoint2D::setDomain(Domain *theDomain)
{
  //Ckeck domain not null - invoked when object removed from a domain
  if (theDomain == 0) {
    end1Ptr = 0;
    end2Ptr = 0;
    end3Ptr = 0;
    end4Ptr = 0;
  }

  this->DomainComponent::setDomain(theDomain);
	
}//setDomain


int SimpleJoint2D::addMP_SimpleJoint(Domain *theDomain, int mpNum, 
				  int RnodeID, int CnodeID, 
				  int AuxDOF, int LrgDispFlag )
{
  MP_Constraint *Temp_MP;

  // create MP_SimpleJoint constraint
  Temp_MP = new MP_SimpleJoint2D( theDomain, mpNum, RnodeID, CnodeID, AuxDOF, LrgDispFlag );
  
  if (Temp_MP == NULL) {
    opserr << "SimpleJoint2D::addMP_SimpleJoint - WARNING ran out of memory for ForJoint MP_Constraint ";
    return -1;
  }
  // Add the multi-point constraint to the domain
  if (theDomain->addMP_Constraint (Temp_MP) == false) {
    opserr << "SimpleJoint2D::addMP_SimpleJoint - WARNING could not add equalDOF MP_Constraint to domain ";
    delete Temp_MP;
    return -2;
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////
// Public methods called, taken care of for 2D element subclasses
//////////////////////////////////////////////////////////////////////

int SimpleJoint2D::update(void)
{
  return 0;
}

int SimpleJoint2D::commitState()
{
  return 0;
}

int SimpleJoint2D::revertToLastCommit()
{
  return 0;
}

int SimpleJoint2D::revertToStart(void)
{
  return 0;
}

int SimpleJoint2D::getNumExternalNodes(void) const
{
  return 4;
}

const ID &SimpleJoint2D::getExternalNodes(void)
{

  return ExternalNodes;
}

Node **
SimpleJoint2D::getNodePtrs(void)
{

  theNodes[0] = end1Ptr;
  theNodes[1] = end2Ptr;
  theNodes[2] = end3Ptr;
  theNodes[3] = end4Ptr;
  theNodes[4] = IntNodePtr;

  return theNodes;
}

int SimpleJoint2D::getNumDOF(void)
{
  return numDof;
}

const Matrix &SimpleJoint2D::getTangentStiff(void)
{

  return K;
}

const Matrix &SimpleJoint2D::getInitialStiff(void)
{
  return K;
}

const Matrix &SimpleJoint2D::getDamp(void)
{	
  return K;
}

const Matrix &SimpleJoint2D::getMass(void)
{
  return K;
}

void SimpleJoint2D::Print(OPS_Stream &s, int flag )
{
  s << "\nElement: " << getTag() << " type: SimpleJoint2D iNode: "
    << ExternalNodes(0) << " jNode: " << ExternalNodes(1) << "\n"
    << " kNode: " << ExternalNodes(2) << " lNode: " << ExternalNodes(3) << "\n";
}

/////////////////////////////////////////////////////////////////////
// methods for applying and returning loads
//////////////////////////////////////////////////////////////////////

void SimpleJoint2D::zeroLoad(void)
{

}

int SimpleJoint2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}

int SimpleJoint2D::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}



const Vector &SimpleJoint2D::getResistingForce()
{
  return V;
}

const Vector &
SimpleJoint2D::getResistingForceIncInertia()
{
  return V;
}


int SimpleJoint2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	// first determine the four corner points of the element based on
	// the display factor (a measure of the distorted image)
	// store this information in 2 3d vectors v1 and v2

	const Vector &node1Crd = end1Ptr->getCrds();
	const Vector &node2Crd = end2Ptr->getCrds();	
	const Vector &node3Crd = end3Ptr->getCrds();
	const Vector &node4Crd = end4Ptr->getCrds();

	const Vector &node1Disp = end1Ptr->getDisp();
	const Vector &node2Disp = end2Ptr->getDisp();    
	const Vector &node3Disp = end3Ptr->getDisp();
	const Vector &node4Disp = end4Ptr->getDisp();  

	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	
	// calculate the current coordinates of four external nodes
	for (int i=0; i<2; i++) 
    {
		v1(i) = node1Crd(i)+node1Disp(i)*fact;
		v2(i) = node2Crd(i)+node2Disp(i)*fact;
		v3(i) = node3Crd(i)+node3Disp(i)*fact;
		v4(i) = node4Crd(i)+node4Disp(i)*fact;
	}

	// draw the center lines
	int dummy;
	dummy = theViewer.drawLine(v1, v3, 1.0, 1.0);
	dummy = theViewer.drawLine(v2, v4, 1.0, 1.0);
	
	// calculate four corners of the element
	Vector vb(3);
	Vector vc(3);

	vb = v1 - v3;
	vc = v2 - v4;

	v1 = v3 - 0.5 * vc;
	v2 = v1 + vb;
	v3 = v2 + vc;
	v4 = v1 + vc;
	
	dummy = theViewer.drawLine(v1, v2, 1.0, 1.0);
	dummy = theViewer.drawLine(v2, v3, 1.0, 1.0);
	dummy = theViewer.drawLine(v3, v4, 1.0, 1.0);
	dummy = theViewer.drawLine(v4, v1, 1.0, 1.0);

	return 0;

}


//most-probably requires to be overridden
Response* SimpleJoint2D::setResponse(char **argv, int argc, Information &eleInformation)
{
//

	if (strcmp(argv[0],"node") == 0 || strcmp(argv[0],"internalNode") == 0 )
    return new ElementResponse(this, 1, Vector(4));

	else if (strcmp(argv[0],"size") == 0 || strcmp(argv[0],"jointSize") == 0 )
    return new ElementResponse(this, 2, Vector(2));

	else 
		return RotElemPtr->setResponse( argv , argc , eleInformation);
  	
}

int SimpleJoint2D::getResponse(int responseID, Information &eleInformation)
{
	switch (responseID) {
	case -1:
		return -1;
	
	case 1:
		if(eleInformation.theVector!=0)
		{
			const Vector& disp = IntNodePtr->getTrialDisp();
			*(eleInformation.theVector) = disp;
		}
		return 0;

	case 2:
		if(eleInformation.theVector!=0)
		{
			
			
			const Vector &node1Crd = end1Ptr->getCrds();
			const Vector &node2Crd = end2Ptr->getCrds();	
			const Vector &node3Crd = end3Ptr->getCrds();
			const Vector &node4Crd = end4Ptr->getCrds();

			const Vector &node1Disp = end1Ptr->getDisp();
			const Vector &node2Disp = end2Ptr->getDisp();    
			const Vector &node3Disp = end3Ptr->getDisp();
			const Vector &node4Disp = end4Ptr->getDisp();  

			Vector v1(2);
			Vector v2(2);
			Vector v3(2);
			Vector v4(2);
	
			// calculate the current coordinates of four external nodes
			for (int i=0; i<2; i++) 
		    {
				v1(i) = node1Crd(i)+node1Disp(i);
				v2(i) = node2Crd(i)+node2Disp(i);
				v3(i) = node3Crd(i)+node3Disp(i);
				v4(i) = node4Crd(i)+node4Disp(i);
			}
			
			v3 = v3 - v1;
			v4 = v4 - v2;

			(*(eleInformation.theVector))(0) = v3.Norm();
			(*(eleInformation.theVector))(1) = v4.Norm();
		}
		return 0;

	default:
		return -1;
	}
}



int SimpleJoint2D::sendSelf(int commitTag, Channel &theChannel)
{
	int res;
	
	int dataTag = this->getDbTag();
	
	static ID data(8);
	data(0) = this->getTag();
	data(1) = IntNode;
	data(2) = RotElemtag;
	data(3) = numDof;
	
	if (ExternalNodes.Size() != 0 && nodeDbTag == 0) nodeDbTag = theChannel.getDbTag();
	if (InternalConstraints.Size() != 0 && dofDbTag == 0) dofDbTag = theChannel.getDbTag();
	data(4) = nodeDbTag;
	data(5) = dofDbTag;

	data(6) = RotElemPtr->getClassTag();
	int elmDbTag = RotElemPtr->getDbTag();
	if (elmDbTag == 0) {
		elmDbTag = theChannel.getDbTag();
		if (elmDbTag != 0)
			RotElemPtr->setDbTag(elmDbTag);
	}
	data(7) = elmDbTag;
	
	res = theChannel.sendID(dataTag, commitTag, data);
	if (res < 0) {
	  opserr << "WARNING SimpleJoint2D::sendSelf() - " << this->getTag() << 
	    " failed to send ID\n";
	  return -1;
	}
	
	// sends the tags of it's external nodes
	res = theChannel.sendID(nodeDbTag, commitTag, ExternalNodes);
	if (res < 0) {
	  opserr << "WARNING SimpleJoint2D::sendSelf() - " << this->getTag() << 
	    " failed to send Vector\n";
	  return -2;
	}
	
	
	// sends the tags of it's internal constraints
	res = theChannel.sendID(dofDbTag, commitTag, InternalConstraints);
	if (res < 0) {
	  opserr << "WARNING SimpleJoint2D::sendSelf() - " << this->getTag() << 
	    " failed to send Vector\n";
	  return -2;
	}
	
	return 0;
}
	  
int SimpleJoint2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	
	int res;
	int dataTag = this->getDbTag();
	
	static ID data(8);
	res = theChannel.recvID(dataTag, commitTag, data);
	if (res < 0) {
	  opserr << "WARNING Truss::recvSelf() - failed to receive Vector\n";
	  return -1;
	}
	
	this->setTag((int)data(0));
	IntNode = data(1);
	RotElemtag = data(2);
	numDof = data(3);
	nodeDbTag = data(4);
	dofDbTag = data(5);
	
	int elmClassTag = data(6);
	int elmDbTag = data(7);


	// receives the tags of it's external nodes
	res = theChannel.recvID(nodeDbTag, commitTag, ExternalNodes);
	if (res < 0) {
	  opserr << "WARNING SimpleJoint2D::recvSelf() -  failed to receive external nodes\n";
	  return -2;
	}

	
	// receives the tags of it's constraint tags
	res = theChannel.recvID(dofDbTag, commitTag, InternalConstraints);
	if (res < 0) {
	  opserr << "WARNING SimpleJoint2D::recvSelf() - " << this->getTag() << " failed to receive internal constraints\n";
	  return -2;
	}
	
	return 0;
}
