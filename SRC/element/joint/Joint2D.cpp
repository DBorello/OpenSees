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
// $Date: 2003-02-25 23:32:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/Joint2D.cpp,v $

// Written: AAA 03/02
// Revised:
// Joint2D.cpp: implementation of the Joint2D class.
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
#include <MP_Joint2D.h>
#include <ElasticMaterial.h>
#include <ElementResponse.h>
#include <Joint2D.h>
#include <InternalSpring.h>

Matrix Joint2D::K(16,16);
Vector Joint2D::V(16);
Node *Joint2D::theNodes[5];

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Joint2D::Joint2D()
  :Element(0, ELE_TAG_Joint2D ), 
  ExternalNodes(5), InternalConstraints(0), 
  Spring1(0), Spring2(0), Spring3(0), Spring4(0), SpringC(0),
  end1Ptr(0), end2Ptr(0), end3Ptr(0), end4Ptr(0),IntNodePtr(0),
  IntNode(0), TheDomain(0), numDof(0), nodeRecord(0), dofRecord(0)
{

}

Joint2D::Joint2D(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
			     UniaxialMaterial &spring1,	UniaxialMaterial &spring2,
			     UniaxialMaterial &spring3, UniaxialMaterial &spring4,
			     UniaxialMaterial &springC, Domain *theDomain, int LrgDisp)
  :Element(tag, ELE_TAG_Joint2D ), 
  ExternalNodes(5), InternalConstraints(4), 
  Spring1(0), Spring2(0), Spring3(0), Spring4(0), SpringC(0),
  end1Ptr(0), end2Ptr(0), end3Ptr(0), end4Ptr(0),IntNodePtr(0),
  IntNode(IntNodeTag), TheDomain(0), numDof(0), nodeRecord(0), dofRecord(0)
{
  numDof  = 16;

  K.Zero();
  V.Zero();

  TheDomain = theDomain;
  if( TheDomain==NULL ) {
    opserr << "WARNING Joint2D(): Specified domain does not exist";
    opserr << "Domain = 0\n";
    return;
  }

  // Save external node id's
  ExternalNodes(0) = nd1;
  ExternalNodes(1) = nd2;
  ExternalNodes(2) = nd3;
  ExternalNodes(3) = nd4;
  ExternalNodes(4) = IntNodeTag;
  
  // Set the external node pointers
  end1Ptr = TheDomain->getNode(nd1);
  end2Ptr = TheDomain->getNode(nd2);	
  end3Ptr = TheDomain->getNode(nd3);
  end4Ptr = TheDomain->getNode(nd4);
  
  // check domain for existence of external nodes
  if (end1Ptr == 0) {
    opserr << "WARNING Joint2D::setDomain(): Nd1: ";
    opserr << nd1 << "does not exist in model for element \n" << *this;
    return;
  }
    
  if (end2Ptr == 0) {
    opserr << "WARNING Joint2D::setDomain(): Nd2: ";
    opserr << nd2 << "does not exist in model for element\n" << *this;
    return;
  }	
  
  if (end3Ptr == 0) {
    opserr << "WARNING Joint2D::setDomain(): Nd3: ";
    opserr << nd3 << "does not exist in model for element \n" << *this;
    return;
  }
    
  if (end4Ptr == 0) {
    opserr << "WARNING Joint2D::setDomain(): Nd4: ";
    opserr << nd4 << "does not exist in model for element\n" << *this;
    return;
  }	

  // check for a two dimensional domain, since this element supports only two dimensions 
  const Vector &end1Crd = end1Ptr->getCrds();
  const Vector &end2Crd = end2Ptr->getCrds();	
  const Vector &end3Crd = end3Ptr->getCrds();
  const Vector &end4Crd = end4Ptr->getCrds();

  int dimNd1 = end1Crd.Size();
  int dimNd2 = end2Crd.Size();
  int dimNd3 = end3Crd.Size();
  int dimNd4 = end4Crd.Size();

  if (dimNd1 != 2 && dimNd2 != 2 && dimNd3 != 2 && dimNd4 != 2 ) {
    opserr << "WARNING Joint2D::setDomain(): has incorrect space dimension \n";
    opserr << "                                    space dimension not supported by Joint2D";
    return;
  }
	
  // now verify the number of dof at node ends
  int dofNd1 = end1Ptr->getNumberDOF();
  int dofNd2 = end2Ptr->getNumberDOF();	
  int dofNd3 = end3Ptr->getNumberDOF();
  int dofNd4 = end4Ptr->getNumberDOF();

  if (dofNd1 != 3 && dofNd2 != 3 && dofNd3 != 3 && dofNd4 != 3 ) {
    opserr << "WARNING Joint2D::Joint2D: has incorrect degrees of freedom \n";
    opserr << "                                    DOF not supported by Joint2D";
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
    opserr << "WARNING Joint2D::(): zero length\n";
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
    opserr << "WARNING Joint2D::(): can not construct a paralelogram over external nodes\n";
    return;	
  }
	
  // Generate internal node and add it up to domain
  IntNodePtr = new Node ( IntNode , 4, Center1(0) , Center1(1) );
  if ( IntNodePtr == NULL ) 
    {
     opserr << "Joint2D::Joint2D - Unable to generate new nodes , out of memory\n";
     exit(-1);
    } else {
      if( TheDomain->addNode( IntNodePtr ) == false )		// add intenal nodes to domain
	opserr << "Joint2D::Joint2D - unable to add internal nodeto domain\n";
    }

  // make copy of the uniaxial materials for the element

  Spring1 = spring1.getCopy();
  Spring2 = spring2.getCopy();
  Spring3 = spring3.getCopy();
  Spring4 = spring4.getCopy();
  SpringC = springC.getCopy();
  
  if ( Spring1 == NULL || Spring2 == NULL || Spring3 == NULL || Spring4 == NULL || SpringC == NULL )
  {
    opserr << "ERROR Joint2D::Joint2D(): Can not make copy of uniaxial materials, out of memory ";
    exit(-1);
  }

  // Generate and add constraints to domain

  // get the constraint numbers
  int startMPtag = theDomain->getNumMPs();
  for ( int i=0 ; i<4 ; i++ ) InternalConstraints(i) = startMPtag + i ;
  
  // create MP_Joint constraint node 1
  if ( addMP_Joint( TheDomain, InternalConstraints(0), IntNode, ExternalNodes(0), 2, LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 1\n";
    return;
  }

  // create MP_Joint constraint node 2
  if ( addMP_Joint( TheDomain, InternalConstraints(1), IntNode, ExternalNodes(1), 3, LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 2\n";
    return;
  }

  // create MP_Joint constraint node 3
  if ( addMP_Joint( TheDomain, InternalConstraints(2), IntNode, ExternalNodes(2), 2, LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 3\n";
    return;
  }
  
  // create MP_Joint constraint node 4
  if ( addMP_Joint( TheDomain, InternalConstraints(3), IntNode, ExternalNodes(3), 3, LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 4\n";
    return;
  }
}


Joint2D::~Joint2D()
{
	if ( Spring1 != NULL ) delete Spring1;
	if ( Spring2 != NULL ) delete Spring2;
	if ( Spring3 != NULL ) delete Spring3;
	if ( Spring4 != NULL ) delete Spring4;
	if ( SpringC != NULL ) delete SpringC;
}


void Joint2D::setDomain(Domain *theDomain)
{
	//Ckeck domain not null - invoked when object removed from a domain
	if (theDomain == 0)
	{
		end1Ptr = 0;
		end2Ptr = 0;
		end3Ptr = 0;
		end4Ptr = 0;
	}

	this->DomainComponent::setDomain(theDomain);
	
}//setDomain


int Joint2D::addMP_Joint(Domain *theDomain, int mpNum, 
				  int RnodeID, int CnodeID, 
				  int MainDOF, int LrgDispFlag )
{
	MP_Constraint *Temp_MP;

	// create MP_ForJoint constraint
	Temp_MP = new MP_Joint2D( theDomain, mpNum, RnodeID, CnodeID, MainDOF, LrgDispFlag );
  
	if (Temp_MP == NULL)
	{
		opserr << "Joint2D::addMP_Joint - WARNING ran out of memory for ForJoint MP_Constraint ";
		return -1;
	}
	// Add the multi-point constraint to the domain
	if (theDomain->addMP_Constraint (Temp_MP) == false)
	{
		opserr << "Joint2D::addMP_Joint - WARNING could not add equalDOF MP_Constraint to domain ";
		delete Temp_MP;
		return -2;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////
// Public methods called, taken care of for 2D element subclasses
//////////////////////////////////////////////////////////////////////

int Joint2D::update(void)
{
  return 0;
}

int Joint2D::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "Joint2D::commitState () - failed in base class";
    }    

    int CS1 = Spring1->commitState();
    int CS2 = Spring2->commitState();
    int CS3 = Spring3->commitState();
    int CS4 = Spring4->commitState();
    int CSC = SpringC->commitState();
    
    if ( CS1!=0 || CS2!=0 || CS3!=0 || CS4!=0 || retVal != 0 ) return -1;

    return CS1;
}

int Joint2D::revertToLastCommit()
{
	int RLS1 = Spring1->revertToLastCommit();
	int RLS2 = Spring2->revertToLastCommit();
	int RLS3 = Spring3->revertToLastCommit();
	int RLS4 = Spring4->revertToLastCommit();
	int RLSC = SpringC->revertToLastCommit();

	if ( RLS1!=RLS2 || RLS2!=RLS3 || RLS3!=RLS4 || RLS4!=RLSC ) return -1;
	return RLS1;
}

int Joint2D::revertToStart(void)
{
	int RS1 = Spring1->revertToStart();
	int RS2 = Spring2->revertToStart();
	int RS3 = Spring3->revertToStart();
	int RS4 = Spring4->revertToStart();
	int RSC = SpringC->revertToStart();

	if ( RS1!=RS2 || RS2!=RS3 || RS3!=RS4 || RS4!=RSC ) return -1;
	return RS1;
}


int Joint2D::getNumExternalNodes(void) const
{
	return 5;
}

const ID &Joint2D::getExternalNodes(void)
{
	return ExternalNodes;
}

Node **
Joint2D::getNodePtrs(void)
{
  theNodes[0] = end1Ptr;
  theNodes[1] = end2Ptr;
  theNodes[2] = end3Ptr;
  theNodes[3] = end4Ptr;
  theNodes[4] = IntNodePtr;

  return theNodes;
}

int Joint2D::getNumDOF(void)
{
  return numDof;
}

const Matrix &Joint2D::getTangentStiff(void)
{
	const Vector &disp1 = end1Ptr->getTrialDisp();
	const Vector &disp2 = end2Ptr->getTrialDisp();
	const Vector &disp3 = end3Ptr->getTrialDisp();
	const Vector &disp4 = end4Ptr->getTrialDisp();
	const Vector &dispC = IntNodePtr->getTrialDisp();

	double Delta1 = disp1(2) - dispC(3);
	double Delta2 = disp2(2) - dispC(2);
	double Delta3 = disp3(2) - dispC(3);
	double Delta4 = disp4(2) - dispC(2);
	double DeltaC = dispC(3) - dispC(2);

	int dummy;

	dummy = Spring1->setTrialStrain(Delta1);
	dummy = Spring2->setTrialStrain(Delta2);
	dummy = Spring3->setTrialStrain(Delta3);
	dummy = Spring4->setTrialStrain(Delta4);
	dummy = SpringC->setTrialStrain(DeltaC);

	double K1 = Spring1->getTangent();
	double K2 = Spring2->getTangent();
	double K3 = Spring3->getTangent();
	double K4 = Spring4->getTangent();
	double KC = SpringC->getTangent();

	K.Zero();

	K(2,2)  =  K1;
	K(2,15) = -K1;
	K(5,5)  =  K2;
	K(5,14) = -K2;
	K(8,8)  =  K3;
	K(8,15) = -K3;
	K(11,11)=  K4;
	K(11,14)= -K4;
	K(14,5) = -K2;
	K(14,11)= -K4;
	K(14,14)=  K2 + K4 + KC;
	K(14,15)= -KC;
	K(15,2) = -K1;
	K(15,8) = -K3;
	K(15,14)= -KC;
	K(15,15)=  K1 + K3 + KC;

	return K;
}

const Matrix &Joint2D::getInitialStiff(void)
{
  double K1 = Spring1->getInitialTangent();
  double K2 = Spring2->getInitialTangent();
  double K3 = Spring3->getInitialTangent();
  double K4 = Spring4->getInitialTangent();
  double KC = SpringC->getInitialTangent();
  
  K.Zero();
  
  K(2,2)  =  K1;
  K(2,15) = -K1;
  K(5,5)  =  K2;
  K(5,14) = -K2;
  K(8,8)  =  K3;
  K(8,15) = -K3;
  K(11,11)=  K4;
  K(11,14)= -K4;
  K(14,5) = -K2;
  K(14,11)= -K4;
  K(14,14)=  K2 + K4 + KC;
  K(14,15)= -KC;
  K(15,2) = -K1;
  K(15,8) = -K3;
  K(15,14)= -KC;
  K(15,15)=  K1 + K3 + KC;
  
  return K;
}


void Joint2D::Print(OPS_Stream &s, int flag )
{
  s << "\nElement: " << getTag() << " type: Joint2D iNode: "
    << ExternalNodes(0) << " jNode: " << ExternalNodes(1) << "\n"
    << " kNode: " << ExternalNodes(2) << " lNode: " << ExternalNodes(3) << "\n"
	<< " Internal node: " << ExternalNodes(4) << "\n";
}

/////////////////////////////////////////////////////////////////////
// methods for applying and returning loads
//////////////////////////////////////////////////////////////////////

void Joint2D::zeroLoad(void)
{

}

int Joint2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}

int Joint2D::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}



const Vector &Joint2D::getResistingForce()
{
	const Vector &disp1 = end1Ptr->getTrialDisp();
	const Vector &disp2 = end2Ptr->getTrialDisp();
	const Vector &disp3 = end3Ptr->getTrialDisp();
	const Vector &disp4 = end4Ptr->getTrialDisp();
	const Vector &dispC = IntNodePtr->getTrialDisp();

	double Delta1 = disp1(2) - dispC(3);
	double Delta2 = disp2(2) - dispC(2);
	double Delta3 = disp3(2) - dispC(3);
	double Delta4 = disp4(2) - dispC(2);
	double DeltaC = dispC(3) - dispC(2);

	int dummy;

	dummy = Spring1->setTrialStrain(Delta1);
	dummy = Spring2->setTrialStrain(Delta2);
	dummy = Spring3->setTrialStrain(Delta3);
	dummy = Spring4->setTrialStrain(Delta4);
	dummy = SpringC->setTrialStrain(DeltaC);

	double F1 = Spring1->getStress();
	double F2 = Spring2->getStress();
	double F3 = Spring3->getStress();
	double F4 = Spring4->getStress();
	double FC = SpringC->getStress();

	V.Zero();

	V(2) = F1;
	V(5) = F2;
	V(8) = F3;
	V(11)= F4;
	V(14)= -FC - F2 - F4;
	V(15)= FC - F1 - F3;

	return V;
}

const Vector &
Joint2D::getResistingForceIncInertia()
{
  this->getResistingForce();

  // add rayleigh damping forces if factors have been set
  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    V += this->getRayleighDampingForces();

  return V;
}



int Joint2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
Response* Joint2D::setResponse(const char **argv, int argc, Information &eleInformation)
{
//
// we compare argv[0] for known response types for the Truss
//

	if (strcmp(argv[0],"node") == 0 || strcmp(argv[0],"internalNode") == 0 )
    return new ElementResponse(this, 1, Vector(4));

	else if (strcmp(argv[0],"size") == 0 || strcmp(argv[0],"jointSize") == 0 )
    return new ElementResponse(this, 2, Vector(2));

	else if (strcmp(argv[0],"moment") == 0 || strcmp(argv[0],"force") == 0 )
    return new ElementResponse(this, 3, Vector(5));

	else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	   strcmp(argv[0],"deformation") == 0 )
    return new ElementResponse(this, 4, Vector(5));

	else if (strcmp(argv[0],"defoANDforce") == 0 || strcmp(argv[0],"deformationANDforce") == 0 ||
	   strcmp(argv[0],"deformationsANDforces") == 0 )
    return new ElementResponse(this, 5, Vector(10));

	else if ( strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0 )
    return new ElementResponse(this, 6, Matrix(16,16) );

	else 
		return 0;
  	
}

int Joint2D::getResponse(int responseID, Information &eleInformation)
{
	switch (responseID) {
	case -1:
		return -1;
	
	case 1:
		if(eleInformation.theVector!=0)
		{
			const Vector& disp = IntNodePtr->getTrialDisp();
			(*(eleInformation.theVector))(0) = disp(0);
			(*(eleInformation.theVector))(1) = disp(1);
			(*(eleInformation.theVector))(2) = disp(2);
			(*(eleInformation.theVector))(3) = disp(3);
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

			v1(0) = v3.Norm();
			v1(1) = v4.Norm();

			*(eleInformation.theVector) = v1;
		}
		return 0;

	case 3:
		if(eleInformation.theVector!=0)
		{
			(*(eleInformation.theVector))(0) = Spring1->getStress();
			(*(eleInformation.theVector))(1) = Spring2->getStress();
			(*(eleInformation.theVector))(2) = Spring3->getStress();
			(*(eleInformation.theVector))(3) = Spring4->getStress();
			(*(eleInformation.theVector))(4) = SpringC->getStress();
		}
		return 0;

	case 4:
		if(eleInformation.theVector!=0)
		{
			(*(eleInformation.theVector))(0) = Spring1->getStrain();
			(*(eleInformation.theVector))(1) = Spring2->getStrain();
			(*(eleInformation.theVector))(2) = Spring3->getStrain();
			(*(eleInformation.theVector))(3) = Spring4->getStrain();
			(*(eleInformation.theVector))(4) = SpringC->getStrain();
		}
		return 0;

	case 5:
		if(eleInformation.theVector!=0)
		{
			(*(eleInformation.theVector))(0) = Spring1->getStrain();
			(*(eleInformation.theVector))(1) = Spring1->getStress();
			(*(eleInformation.theVector))(2) = Spring2->getStrain();
			(*(eleInformation.theVector))(3) = Spring2->getStress();
			(*(eleInformation.theVector))(4) = Spring3->getStrain();
			(*(eleInformation.theVector))(5) = Spring3->getStress();
			(*(eleInformation.theVector))(6) = Spring4->getStrain();
			(*(eleInformation.theVector))(7) = Spring4->getStress();
			(*(eleInformation.theVector))(8) = SpringC->getStrain();
			(*(eleInformation.theVector))(9) = SpringC->getStress();
		}
		return 0;

	case 6:
		return eleInformation.setMatrix(this->getTangentStiff());

	default:
		return -1;
	}
}
