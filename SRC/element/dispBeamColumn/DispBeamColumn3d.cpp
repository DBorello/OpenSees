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
                                                                        
// $Revision: 1.6 $
// $Date: 2002-03-22 20:47:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/DispBeamColumn3d.cpp,v $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumn3d.

#include <DispBeamColumn3d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>

#include <G3Globals.h>

Matrix DispBeamColumn3d::K(12,12);
Vector DispBeamColumn3d::P(12);
double DispBeamColumn3d::workArea[200];
GaussQuadRule1d01 DispBeamColumn3d::quadRule;

DispBeamColumn3d::DispBeamColumn3d(int tag, int nd1, int nd2,
		int numSec, SectionForceDeformation &s,
		CrdTransf3d &coordTransf, double r)
:Element (tag, ELE_TAG_DispBeamColumn3d),
numSections(numSec), theSections(0), crdTransf(0),
connectedExternalNodes(2), L(0.0), nd1Ptr(0), nd2Ptr(0),
Q(12), q(6), rho(r)
{
    // Allocate arrays of pointers to SectionForceDeformations
    theSections = new SectionForceDeformation *[numSections];
    
	if (theSections == 0)
	    g3ErrorHandler->fatal("%s - failed to allocate section model pointer",
			"DispBeamColumn3d::DispBeamColumn3d");

    for (int i = 0; i < numSections; i++) {

		// Get copies of the material model for each integration point
		theSections[i] = s.getCopy();
			
		// Check allocation
		if (theSections[i] == 0)
			g3ErrorHandler->fatal("%s -- failed to get a copy of section model",
				"DispBeamColumn3d::DispBeamColumn3d");
	}

	crdTransf = coordTransf.getCopy();

	if (crdTransf == 0)
	    g3ErrorHandler->fatal("%s - failed to copy coordinate transformation",
			"DispBeamColumn3d::DispBeamColumn3d");

	// Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
}

DispBeamColumn3d::DispBeamColumn3d()
:Element (0, ELE_TAG_DispBeamColumn3d),
numSections(0), theSections(0), crdTransf(0),
connectedExternalNodes(2), L(0.0), nd1Ptr(0), nd2Ptr(0),
Q(12), q(6), rho(0.0)
{

}

DispBeamColumn3d::~DispBeamColumn3d()
{    
    for (int i = 0; i < numSections; i++) {
		if (theSections[i])
			delete theSections[i];
	}

    // Delete the array of pointers to SectionForceDeformation pointer arrays
    if (theSections)
		delete [] theSections;

	if (crdTransf)
		delete crdTransf;
}

int
DispBeamColumn3d::getNumExternalNodes() const
{
    return 2;
}

const ID&
DispBeamColumn3d::getExternalNodes()
{
    return connectedExternalNodes;
}

int
DispBeamColumn3d::getNumDOF()
{
    return 12;
}

void
DispBeamColumn3d::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nd1Ptr = 0;
	nd2Ptr = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);

    nd1Ptr = theDomain->getNode(Nd1);
    nd2Ptr = theDomain->getNode(Nd2);

    if (nd1Ptr == 0 || nd2Ptr == 0) {
	//g3ErrorHandler->fatal("FATAL ERROR DispBeamColumn3d (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = nd1Ptr->getNumberDOF();
    int dofNd2 = nd2Ptr->getNumberDOF();
    
    if (dofNd1 != 6 || dofNd2 != 6) {
	//g3ErrorHandler->fatal("FATAL ERROR DispBeamColumn3d (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());
	
	return;
    }

	if (crdTransf->initialize(nd1Ptr, nd2Ptr)) {
		// Add some error check
	}

	L = crdTransf->getInitialLength();

	if (L == 0.0) {
		// Add some error check
	}

    this->DomainComponent::setDomain(theDomain);

	this->update();
}

int
DispBeamColumn3d::commitState()
{
    int retVal = 0;

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int
DispBeamColumn3d::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
DispBeamColumn3d::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
DispBeamColumn3d::update(void)
{
	// Update the transformation
	crdTransf->update();

    // Get basic deformations
    static Vector v(6);
	v = crdTransf->getBasicTrialDisp();

	double oneOverL = 1.0/L;
	const Matrix &pts = quadRule.getIntegrPointCoords(numSections);

	// Assuming member is prismatic ... have to move inside
	// the loop if it is not prismatic
	int order = theSections[0]->getOrder();
	const ID &code = theSections[0]->getType();
	Vector e(workArea, order);

	// Loop over the integration points
	for (int i = 0; i < numSections; i++) {

		double xi6 = 6.0*pts(i,0);

		int j;
		for (j = 0; j < order; j++) {
			switch(code(j)) {
			case SECTION_RESPONSE_P:
				e(j) = oneOverL*v(0); break;
			case SECTION_RESPONSE_MZ:
				e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)); break;
			case SECTION_RESPONSE_MY:
				e(j) = oneOverL*((xi6-4.0)*v(3) + (xi6-2.0)*v(4)); break;
			case SECTION_RESPONSE_T:
				e(j) = oneOverL*v(5); break;
			default:
				e(j) = 0.0; break;
			}
		}

		// Set the section deformations
		theSections[i]->setTrialSectionDeformation(e);
	}

	return 0;
}

const Matrix&
DispBeamColumn3d::getTangentStiff()
{
	static Matrix kb(6,6);

	// Zero for integral
	kb.Zero();
	q.Zero();

	const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
	const Vector &wts = quadRule.getIntegrPointWeights(numSections);

	// Assuming member is prismatic ... have to move inside
	// the loop if it is not prismatic
	int order = theSections[0]->getOrder();
	const ID &code = theSections[0]->getType();

	double oneOverL = 1.0/L;
	Matrix ka(workArea, order, 6);

	// Loop over the integration points
	for (int i = 0; i < numSections; i++) {

		// Get the section tangent stiffness and stress resultant
		const Matrix &ks = theSections[i]->getSectionTangent();
		const Vector &s = theSections[i]->getStressResultant();

		double xi6 = 6.0*pts(i,0);
		ka.Zero();

		// Perform numerical integration
		//kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
		double wti = wts(i)*oneOverL;
		double tmp;
		int j, k;
		for (j = 0; j < order; j++) {
			switch(code(j)) {
			case SECTION_RESPONSE_P:
				for (k = 0; k < order; k++)
					ka(k,0) += ks(k,j)*wti;
				break;
			case SECTION_RESPONSE_MZ:
				for (k = 0; k < order; k++) {
					tmp = ks(k,j)*wti;
					ka(k,1) += (xi6-4.0)*tmp;
					ka(k,2) += (xi6-2.0)*tmp;
				}
                break;
			case SECTION_RESPONSE_MY:
				for (k = 0; k < order; k++) {
					tmp = ks(k,j)*wti;
					ka(k,3) += (xi6-4.0)*tmp;
					ka(k,4) += (xi6-2.0)*tmp;
				}
				break;
			case SECTION_RESPONSE_T:
				for (k = 0; k < order; k++)
					ka(k,5) += ks(k,j)*wti;
				break;
			default:
				break;
			}
		}
		for (j = 0; j < order; j++) {
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				for (k = 0; k < 6; k++)
					kb(0,k) += ka(j,k);
				break;
			case SECTION_RESPONSE_MZ:
				for (k = 0; k < 6; k++) {
					tmp = ka(j,k);
					kb(1,k) += (xi6-4.0)*tmp;
					kb(2,k) += (xi6-2.0)*tmp;
				}
				break;
			case SECTION_RESPONSE_MY:
				for (k = 0; k < 6; k++) {
					tmp = ka(j,k);
					kb(3,k) += (xi6-4.0)*tmp;
					kb(4,k) += (xi6-2.0)*tmp;
				}
				break;
			case SECTION_RESPONSE_T:
				for (k = 0; k < 6; k++)
					kb(5,k) += ka(j,k);
				break;
			default:
				break;
			}
		}

		//q.addMatrixTransposeVector(1.0, *B, s, wts(i));
		double si;
		for (j = 0; j < order; j++) {
			si = s(j)*wts(i);
			switch(code(j)) {
			case SECTION_RESPONSE_P:
				q(0) += si; break;
			case SECTION_RESPONSE_MZ:
				q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si; break;
			case SECTION_RESPONSE_MY:
				q(3) += (xi6-4.0)*si; q(4) += (xi6-2.0)*si; break;
			case SECTION_RESPONSE_T:
				q(5) += si; break;
			default:
				break;
			}
		}

	}

	// Transform to global stiffness
	K = crdTransf->getGlobalStiffMatrix(kb, q);

	return K;
}

const Matrix&
DispBeamColumn3d::getSecantStiff()
{
	return this->getTangentStiff();
}

const Matrix&
DispBeamColumn3d::getDamp()
{
	K.Zero();
	
	return K;
}

const Matrix&
DispBeamColumn3d::getMass()
{
	K.Zero();

	if (rho == 0.0)
		return K;

	double m = 0.5*rho*L;

	K(0,0) = K(1,1) = K(2,2) = K(6,6) = K(7,7) = K(8,8) = m;

	return K;
}

void
DispBeamColumn3d::zeroLoad(void)
{
	Q.Zero();

	return;
}

int 
DispBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  g3ErrorHandler->warning("DispBeamColumn3d::addLoad - load type unknown for truss with tag: %d\n",
			  this->getTag());

  return -1;
}

int 
DispBeamColumn3d::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0) 
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = nd1Ptr->getRV(accel);
	const Vector &Raccel2 = nd2Ptr->getRV(accel);

    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
		g3ErrorHandler->warning("DispBeamColumn3d::addInertiaLoadToUnbalance %s\n",
				"matrix and vector sizes are incompatable");
		return -1;
    }

	double m = 0.5*rho*L;

    // Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	Q(0) -= m*Raccel1(0);
	Q(1) -= m*Raccel1(1);
    Q(2) -= m*Raccel1(2);
	Q(6) -= m*Raccel2(0);
	Q(7) -= m*Raccel2(1);
    Q(8) -= m*Raccel2(2);

    return 0;
}

const Vector&
DispBeamColumn3d::getResistingForce()
{
	const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
	const Vector &wts = quadRule.getIntegrPointWeights(numSections);

	// Assuming member is prismatic ... have to move inside
	// the loop if it is not prismatic
	int order = theSections[0]->getOrder();
	const ID &code = theSections[0]->getType();

	// Zero for integration
	q.Zero();

	// Loop over the integration points
	for (int i = 0; i < numSections; i++) {

		double xi6 = 6.0*pts(i,0);

		// Get section stress resultant
		const Vector &s = theSections[i]->getStressResultant();

		// Perform numerical integration on internal force
		//q.addMatrixTransposeVector(1.0, *B, s, wts(i));

		double si;
		for (int j = 0; j < order; j++) {
			si = s(j)*wts(i);
			switch(code(j)) {
			case SECTION_RESPONSE_P:
				q(0) += si; break;
			case SECTION_RESPONSE_MZ:
				q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si; break;
			case SECTION_RESPONSE_MY:
				q(3) += (xi6-4.0)*si; q(4) += (xi6-2.0)*si; break;
			case SECTION_RESPONSE_T:
				q(5) += si; break;
			default:
				break;
			}
		}

	}

	// Transform forces
	static Vector dummy(3);		// No distributed loads
	P = crdTransf->getGlobalResistingForce(q,dummy);

	// Subtract other external nodal loads ... P_res = P_int - P_ext
	P.addVector(1.0, Q, -1.0);

	return P;
}

const Vector&
DispBeamColumn3d::getResistingForceIncInertia()
{
	// Check for a quick return
	if (rho == 0.0)
		return this->getResistingForce();

	const Vector &accel1 = nd1Ptr->getTrialAccel();
	const Vector &accel2 = nd2Ptr->getTrialAccel();
	
	// Compute the current resisting force
	this->getResistingForce();

	double m = 0.5*rho*L;

	P(0) += m*accel1(0);
	P(1) += m*accel1(1);
    P(2) += m*accel1(2);
	P(6) += m*accel2(0);
	P(7) += m*accel2(1);
    P(8) += m*accel2(2);

	return P;
}

int
DispBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
DispBeamColumn3d::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
	return -1;
}

void
DispBeamColumn3d::Print(ostream &s, int flag)
{
	s << "\nDispBeamColumn3d, element id:  " << this->getTag() << endl;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tmass density:  " << rho << endl;
	theSections[0]->Print(s,flag);
}


int
DispBeamColumn3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = nd1Ptr->getCrds();
    const Vector &end2Crd = nd2Ptr->getCrds();	

    const Vector &end1Disp = nd1Ptr->getDisp();
    const Vector &end2Disp = nd2Ptr->getDisp();

	static Vector v1(3);
	static Vector v2(3);

	for (int i = 0; i < 3; i++) {
		v1(i) = end1Crd(i) + end1Disp(i)*fact;
		v2(i) = end2Crd(i) + end2Disp(i)*fact;    
	}
	
	return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
DispBeamColumn3d::setResponse(char **argv, int argc, Information &eleInfo)
{
    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
		|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
		return new ElementResponse(this, 1, P);

    // local force -
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
		return new ElementResponse(this, 2, P);
    
    // section response -
    else if (strcmp(argv[0],"section") ==0) {
		if (argc <= 2)
			return 0;
	
		int sectionNum = atoi(argv[1]);
		if (sectionNum > 0 && sectionNum <= numSections)
			return theSections[sectionNum-1]->setResponse(&argv[2], argc-2, eleInfo);
		else
			return 0;
	}
    
	else
		return 0;
}

int 
DispBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  double V;

  switch (responseID) {
    case 1:  // global forces
      return eleInfo.setVector(this->getResistingForce());

    case 2:
		// Axial
		P(6) = q(0);
		P(0) = -q(0);

		// Torsion
		P(11) = q(5);
		P(5)  = -q(5);

		// Moments about z and shears along y
		P(2) = q(1);
		P(8) = q(2);
		V = (q(1)+q(2))/L;
		P(1) = V;
		P(7) = -V;

		// Moments about y and shears along z
		P(4)  = q(3);
		P(10) = q(4);
		V = (q(3)+q(4))/L;
		P(3) = -V;
		P(9) = V;
      return eleInfo.setVector(P);
      
    default: 
	  return -1;
  }
}
