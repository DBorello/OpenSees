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
                                                                        
// $Revision: 1.27 $
// $Date: 2010-02-04 01:17:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/CoupledZeroLengthoupled.cpp,v $

// Written: GLF
// Created: 12/99
// Revision: A
//
// Description: This file contains the implementation for the CoupledZeroLengthoupled class.
//
// What: "@(#) CoupledZeroLengthoupled.C, revA"

#include <CoupledZeroLength.h>
#include <G3Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static Matrix CoupledZeroLengthM2(2,2);   // class wide matrix for 2*2
static Matrix CoupledZeroLengthM4(4,4);   // class wide matrix for 4*4
static Matrix CoupledZeroLengthM6(6,6);   // class wide matrix for 6*6
static Matrix CoupledZeroLengthM12(12,12);  // class wide matrix for 12*12
static Vector CoupledZeroLengthV2(2);   // class wide Vector for size 2
static Vector CoupledZeroLengthV4(4);   // class wide Vector for size 4
static Vector CoupledZeroLengthV6(6);   // class wide Vector for size 6
static Vector CoupledZeroLengthV12(12);  // class wide Vector for size 12


#include <elementAPI.h>

void *
OPS_CoupledZeroLength()
{

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new CoupledZeroLength();
    return theEle;
  }

  if (numRemainingArgs != 6 && numRemainingArgs != 7) {
    opserr << "ERROR - CoupledZeroLength not enough args provided, want: element CoupledZeroLength tag? iNode? jNode? dirn1? dirn2? matTag? <useRayleigh?>\n";
  }

  // get the id and end nodes 
  int iData[7];
  int numData;

  iData[6] = 0; // turn off rayleigh damping by default

  numData = numRemainingArgs;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];
  int matID = iData[5];
  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matID);
  
  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matID << "not found for element " << eleTag << endln;
    return 0;
  }

  // now create the truss and add it to the Domain

  theEle = new CoupledZeroLength(eleTag, iData[1], iData[2], *theMaterial, iData[3]-1, iData[4]-1, iData[6]);

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    delete theMaterial;
    return 0;
  }

  return theEle;
}



//  Constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the CoupledZeroLength end nodes.

//  Construct element with one unidirectional material (numMaterials1d=1)
CoupledZeroLength::CoupledZeroLength(int tag,
			 int Nd1, int Nd2, 
			 UniaxialMaterial &theMat,
			 int direction1, int direction2,
			 int doRayleigh)
  :Element(tag,ELE_TAG_CoupledZeroLength),     
   connectedExternalNodes(2),
   dimension(0), numDOF(0), transformation(3,3), useRayleighDamping(doRayleigh),
   theMatrix(0), theVector(0),
   theMaterial(0), dirn1(direction1), dirn2(direction2), d0(0), v0(0)
{
  // allocate memory for numMaterials1d uniaxial material models
  theMaterial = theMat.getCopy();
  if ( theMaterial == 0) {
    opserr << "FATAL CoupledZeroLength::CoupledZeroLength - failed to create a 1d  material\n";
    exit(-1);
  }

  // initialize uniaxial materials and directions and check for valid values
  if (direction1 < 0 || direction1 > 5 || direction2 < 0 || direction2 > 5) {
    opserr << "FATAL: CoupledZeroLength::CoupledZeroLength - invalid diection\n";
    exit(-1);
  }

  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  dX = 0.0;
  dY = 0.0;
  fX = 0.0;
  fY = 0.0;
}


//   Constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
CoupledZeroLength::CoupledZeroLength(void)
  :Element(0,ELE_TAG_CoupledZeroLength),     
   connectedExternalNodes(2),
   dimension(0), numDOF(0), transformation(3,3),
   theMatrix(0), theVector(0),
   theMaterial(0), dirn1(0), dirn2(0), d0(0), v0(0)
{
  // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2)
    opserr << "FATAL CoupledZeroLength::CoupledZeroLength - failed to create an ID of correct size\n";

  dX = 0.0;
  dY = 0.0;
  fX = 0.0;
  fY = 0.0;
}


//  Destructor:
//  delete must be invoked on any objects created by the object
//  and on the matertial object.
CoupledZeroLength::~CoupledZeroLength()
{
  if (theMaterial != 0)
    delete theMaterial;

  if (d0 != 0)
    delete d0;
  
  if (v0 != 0)
    delete v0;
}


int
CoupledZeroLength::getNumExternalNodes(void) const
{
    return 2;
}


const ID &
CoupledZeroLength::getExternalNodes(void) 
{
    return connectedExternalNodes;
}



Node **
CoupledZeroLength::getNodePtrs(void) 
{
  return theNodes;
}

int
CoupledZeroLength::getNumDOF(void) 
{
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the CoupledZeroLength element, we set matrix and vector pointers,
//    allocate space for t matrix and define it as the basic deformation-
//    displacement transformation matrix.
void
CoupledZeroLength::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	return;
    }

    // set default values for error conditions
    numDOF = 2;
    theMatrix = &CoupledZeroLengthM2;
    theVector = &CoupledZeroLengthV2;
    
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	

    // if can't find both - send a warning message
    if ( theNodes[0] == 0 || theNodes[1] == 0 ) {
      if (theNodes[0] == 0) 
        opserr << "WARNING CoupledZeroLength::setDomain() - Nd1: " << Nd1 << " does not exist in ";
      else
        opserr << "WARNING CoupledZeroLength::setDomain() - Nd2: " << Nd2 << " does not exist in ";

      opserr << "model for CoupledZeroLength ele: " << this->getTag() << endln;

      return;
    }

    // now determine the number of dof and the dimension    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if ( dofNd1 != dofNd2 ) {
      opserr << "WARNING CoupledZeroLength::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
	"have differing dof at ends for CoupledZeroLength " << this->getTag() << endln;
      return;
    }	

    // Check that length is zero within tolerance
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    dimension = end1Crd.Size();

    Vector diff = end1Crd - end2Crd;
    double L  = diff.Norm();
    double v1 = end1Crd.Norm();
    double v2 = end2Crd.Norm();
    double vm;
    
    vm = (v1<v2) ? v2 : v1;


    if (L > LENTOL*vm)
      opserr << "WARNING CoupledZeroLength::setDomain(): Element " << this->getTag() << " has L= " << L << 
	", which is greater than the tolerance\n";
        
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set the number of dof for element and set matrix and vector pointer
    if (dimension == 1 && dofNd1 == 1) {
	numDOF = 2;    
	theMatrix = &CoupledZeroLengthM2;
	theVector = &CoupledZeroLengthV2;
	elemType  = D1N2;
    }
    else if (dimension == 2 && dofNd1 == 2) {
	numDOF = 4;
	theMatrix = &CoupledZeroLengthM4;
	theVector = &CoupledZeroLengthV4;
	elemType  = D2N4;
    }
    else if (dimension == 2 && dofNd1 == 3) {
	numDOF = 6;	
	theMatrix = &CoupledZeroLengthM6;
	theVector = &CoupledZeroLengthV6;
	elemType  = D2N6;
    }
    else if (dimension == 3 && dofNd1 == 3) {
	numDOF = 6;	
	theMatrix = &CoupledZeroLengthM6;
	theVector = &CoupledZeroLengthV6;
	elemType  = D3N6;
    }
    else if (dimension == 3 && dofNd1 == 6) {
	numDOF = 12;	    
	theMatrix = &CoupledZeroLengthM12;
	theVector = &CoupledZeroLengthV12;
	elemType  = D3N12;
    }
    else {
      opserr << "WARNING CoupledZeroLength::setDomain cannot handle " << dimension << 
	"dofs at nodes in " << dofNd1 << " d problem\n"; 
      return;
    }

    // get trial displacements and take difference
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    Vector  diffD  = disp2-disp1;
    const Vector& vel1  = theNodes[0]->getTrialVel();
    const Vector& vel2  = theNodes[1]->getTrialVel();
    Vector  diffV = vel2-vel1;

    if (diffD != 0.0)
      d0 = new Vector(diffD);
    
    if (diffV != 0)
      v0 = new Vector(diffV);
}   	 



int
CoupledZeroLength::commitState()
{
    int code=0;

    // call element commitState to do any base class stuff
    if ((code = this->Element::commitState()) != 0) {
      opserr << "CoupledZeroLength::commitState () - failed in base class";
    }    

    // commit 1d materials
    code += theMaterial->commitState();

    double strain = theMaterial->getStrain();
    double force = theMaterial->getStress();

    if (strain != 0.0 && force != 0.0) {
      fX = force*dX/strain;
      fY = force*dY/strain;
    }

    return code;
}

int
CoupledZeroLength::revertToLastCommit()
{
    int code=0;
    
    // revert state for 1d materials
    code += theMaterial->revertToLastCommit();
    
    return code;
}


int
CoupledZeroLength::revertToStart()
{   
    int code=0;
    
    // revert to start for 1d materials
    code += theMaterial->revertToStart();
    
    return code;
}


int
CoupledZeroLength::update(void)
{
    double strain;
    double strainRate;

    // get trial displacements and take difference
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    Vector  diff  = disp2-disp1;
    const Vector& vel1  = theNodes[0]->getTrialVel();
    const Vector& vel2  = theNodes[1]->getTrialVel();
    Vector  diffv = vel2-vel1;


    if (d0 != 0)
      diff -= *d0;

    if (v0 != 0)
      diffv -= *v0;

    dX = diffv(dirn1);
    dY = diffv(dirn2);

    strainRate = sqrt(dX*dX + dY*dY);


    dX = diff(dirn1);
    dY = diff(dirn2);
    strain = sqrt(dX*dX + dY*dY);

    return theMaterial->setTrialStrain(strain,strainRate);
}

const Matrix &
CoupledZeroLength::getTangentStiff(void)
{
    double E;

    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // zero stiffness matrix
    stiff.Zero();

    E = theMaterial->getTangent();    

    int numNodeDof = numDOF/2;
    int dirn1b = dirn1+numNodeDof;
    int dirn2b = dirn2+numNodeDof;

    stiff(dirn1,dirn1)   = E;
    stiff(dirn1b,dirn1b) = E;      
    stiff(dirn1,dirn1b)  = -E;
    stiff(dirn1b,dirn1)  = -E;      

    stiff(dirn2,dirn2)   = E;
    stiff(dirn2b,dirn2b) = E;      
    stiff(dirn2,dirn2b)  = -E;
    stiff(dirn2b,dirn2)  = -E;      

    return stiff;
}


const Matrix &
CoupledZeroLength::getInitialStiff(void)
{
    double E;

    E = theMaterial->getInitialTangent();

    Matrix& stiff = *theMatrix;
    stiff.Zero();

    int numNodeDof = numDOF/2;
    int dirn1b = dirn1+numNodeDof;
    int dirn2b = dirn2+numNodeDof;

    stiff(dirn1,dirn1)   = E;
    stiff(dirn1b,dirn1b) = E;      
    stiff(dirn1,dirn1b)  = -E;
    stiff(dirn1b,dirn1)  = -E;      

    stiff(dirn2,dirn2)   = E;
    stiff(dirn2b,dirn2b) = E;      
    stiff(dirn2,dirn2b)  = -E;
    stiff(dirn2b,dirn2)  = -E;      

    return stiff;
}
    

const Matrix &
CoupledZeroLength::getDamp(void)
{
  if (useRayleighDamping == 1)
    return this->Element::getDamp();

    double E;

    E = theMaterial->getDampTangent();

    Matrix& stiff = *theMatrix;
    stiff.Zero();

    int numNodeDof = numDOF/2;
    int dirn1b = dirn1+numNodeDof;
    int dirn2b = dirn2+numNodeDof;

    stiff(dirn1,dirn1)   = E;
    stiff(dirn1b,dirn1b) = E;      
    stiff(dirn1,dirn1b)  = -E;
    stiff(dirn1b,dirn1)  = -E;      

    stiff(dirn2,dirn2)   = E;
    stiff(dirn2b,dirn2b) = E;      
    stiff(dirn2,dirn2b)  = -E;
    stiff(dirn2b,dirn2)  = -E;      
  
    return stiff;
}


const Matrix &
CoupledZeroLength::getMass(void)
{
  // no mass 
  theMatrix->Zero();    
  return *theMatrix; 
}


void 
CoupledZeroLength::zeroLoad(void)
{
  // does nothing now
}

int 
CoupledZeroLength::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "CoupledZeroLength::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  return -1;
}

int 
CoupledZeroLength::addInertiaLoadToUnbalance(const Vector &accel)
{
  // does nothing as element has no mass yet!
  return 0;
}


const Vector &
CoupledZeroLength::getResistingForce()
{
  double force, strain;
    
  // zero the residual
  theVector->Zero();

  // get resisting force for material
  force = theMaterial->getStress();
  strain = theMaterial->getStrain();

  double Fx = force;
  double Fy = force;

  if (strain != 0.0) {
    Fx *= dX/strain;
    Fy *= dY/strain;
  } else {
    double oldF = sqrt(fX*fX+fY*fY);
    if (oldF != 0.0) {
      Fx *= fX/oldF;
      Fy *= fY/oldF;
    }
  }

  int numNodeDof = numDOF/2;
  int dirn1b = dirn1+numNodeDof;
  int dirn2b = dirn2+numNodeDof;
  
  (*theVector)(dirn1)   = -Fx;
  (*theVector)(dirn1b)  =  Fx;      
  (*theVector)(dirn2)   = -Fy;
  (*theVector)(dirn2b)  =  Fy;      
  
  //  opserr << "CoupledZeroLength::getResistingForce() " << force << " forces: " << *theVector;

  return *theVector;
}


const Vector &
CoupledZeroLength::getResistingForceIncInertia()
{	
    // there is no mass, so return
    
    this->getResistingForce();

    if (useRayleighDamping == 1)
      if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	*theVector += this->getRayleighDampingForces();

    return *theVector;
}


int
CoupledZeroLength::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// CoupledZeroLength packs its data into an ID and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments

	// Make one size bigger so not a multiple of 3, otherwise will conflict
	// with classTags ID
	static ID idData(10);

	idData(0) = this->getTag();
	idData(1) = dimension;
	idData(2) = numDOF;
	idData(3) = connectedExternalNodes(0);
	idData(4) = connectedExternalNodes(1);
	idData(5) = useRayleighDamping;
	idData(6) = dirn1;
	idData(7) = dirn2;

	int matDbTag = theMaterial->getDbTag();
	if (matDbTag == 0) {
	  matDbTag = theChannel.getDbTag();
	  if (matDbTag != 0)
	    theMaterial->setDbTag(matDbTag);
	}
	idData(8)= matDbTag;
	idData(9)= theMaterial->getClassTag();

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
	  opserr << "CoupledZeroLength::sendSelf -- failed to send ID data\n";
	  return res;
	}

	res += theMaterial->sendSelf(commitTag, theChannel);

	return res;
}

int
CoupledZeroLength::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // CoupledZeroLength creates an ID, receives the ID and then sets the 
  // internal data with the data in the ID

  static ID idData(10);

  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "CoupledZeroLength::recvSelf -- failed to receive ID data\n";
			    
    return res;
  }

  res += theChannel.recvMatrix(dataTag, commitTag, transformation);
  if (res < 0) {
    opserr << "CoupledZeroLength::recvSelf -- failed to receive transformation Matrix\n";
			    
    return res;
  }

  this->setTag(idData(0));
  dimension = idData(1);
  numDOF = idData(2);
  connectedExternalNodes(0) = idData(3);
  connectedExternalNodes(1) = idData(4);
  useRayleighDamping = idData(5);
  dirn1 = idData(6);
  dirn1 = idData(7);

  int matDbTag = idData(8);
  int matClassTag = idData(9);

  // If null, get a new one from the broker
  if (theMaterial == 0 || theMaterial->getClassTag() != matClassTag) {
    if (theMaterial != 0)
      delete theMaterial;

    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
	opserr << "CoupledZeroLength::recvSelf  -- failed to allocate new Material " << endln;
	return -1;
    }
  }

  // Receive the materials
  theMaterial->setDbTag(matDbTag);
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "CoupledZeroLength::recvSelf  -- failed to receive new Material1d " << endln;
  }
  
  return res;
}


int
CoupledZeroLength::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // ensure setDomain() worked
    if (theNodes[0] == 0 || theNodes[1] == 0 )
       return 0;

    // first determine the two end points of the CoupledZeroLength based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();    

    if (displayMode == 1 || displayMode == 2) {
	Vector v1(3);
	Vector v2(3);
	for (int i=0; i<dimension; i++) {
	    v1(i) = end1Crd(i)+end1Disp(i)*fact;
	    v2(i) = end2Crd(i)+end2Disp(i)*fact;    
	}
	
	// don't display strain or force
	double strain = 0.0;
	double force  = 0.0;
    
	if (displayMode == 2) // use the strain as the drawing measure
	    return theViewer.drawLine(v1, v2, (float)strain, (float)strain);	
	else { // otherwise use the axial force as measure
	    return theViewer.drawLine(v1,v2, (float)force, (float)force);
	}
    }
    return 0;
}


void
CoupledZeroLength::Print(OPS_Stream &s, int flag)
{
    // compute the strain and axial force in the member
    if (flag == 0) { // print everything
	s << "Element: " << this->getTag(); 
	s << " type: CoupledZeroLength  iNode: " << connectedExternalNodes(0);
	s << " jNode: " << connectedExternalNodes(1) << endln;
	s << "\tMaterial1d, tag: " << theMaterial->getTag() ;
	  s << *(theMaterial);
    } else if (flag == 1) {
      s << this->getTag() << "  " << theMaterial->getStrain() << "  ";
    }
}

Response*
CoupledZeroLength::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","CoupledZeroLength");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForces") == 0) || (strcmp(argv[0],"globalforces") == 0)) {

            char outputData[10];
            int numDOFperNode = numDOF/2;
            for (int i=0; i<numDOFperNode; i++) {
                sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
            for (int j=0; j<numDOFperNode; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse = new ElementResponse(this, 1, Vector(numDOF));

    // a material quantity
    } else if (strcmp(argv[0],"material") == 0) {
      if (argc > 1) {
	theResponse =  theMaterial->setResponse(&argv[1], argc-1, output);
      }
    }

    output.endTag();

    return theResponse;
}

int 
CoupledZeroLength::getResponse(int responseID, Information &eleInformation)
{
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    const Vector  diff  = disp2-disp1;

    switch (responseID) {
    case -1:
        return -1;

    case 1:
        return eleInformation.setVector(this->getResistingForce());

    default:
        return -1;
    }
}

int
CoupledZeroLength::setParameter(const char **argv, int argc, Parameter &param)
{
  int result = -1;  

  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "material") == 0) {
      if (argc > 1) {
	return theMaterial->setParameter(&argv[1], argc-1, param);
      } else {
	return -1;
      }
  }

  int res = theMaterial->setParameter(argv, argc, param);
  if (res != -1) {
    result = res;
  }

  return result;
}

int
CoupledZeroLength::updateParameter (int parameterID, Information &info)
{
  return 0;
}

int
CoupledZeroLength::activateParameter(int passedParameterID)
{
  
  return 0;
}


