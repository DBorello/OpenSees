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
                                                                        
// $Revision: 1.2 $
// $Date: 2001-08-20 00:37:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/DispBeamColumn2d.cpp,v $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumn2d.

#include <DispBeamColumn2d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <GaussQuadRule1d01.h>
#include <CrdTransf2d.h>
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

Matrix DispBeamColumn2d::K(6,6);
Vector DispBeamColumn2d::P(6);
double DispBeamColumn2d::workArea[100];

DispBeamColumn2d::DispBeamColumn2d(int tag, int nd1, int nd2,
		int numSec, SectionForceDeformation &s,
		CrdTransf2d &coordTransf, double r)
:Element (tag, ELE_TAG_DispBeamColumn2d), rho(r),
 Q(6), q(3), connectedExternalNodes(2), L(0.0),
 numSections(numSec), theSections(0), crdTransf(0)
{
    // Allocate arrays of pointers to SectionForceDeformations
    theSections = new SectionForceDeformation *[numSections];
    
	if (theSections == 0)
	    g3ErrorHandler->fatal("%s - failed to allocate section model pointer",
			"DispBeamColumn2d::DispBeamColumn2d");

    for (int i = 0; i < numSections; i++) {

		// Get copies of the material model for each integration point
		theSections[i] = s.getCopy();
			
		// Check allocation
		if (theSections[i] == 0)
			g3ErrorHandler->fatal("%s -- failed to get a copy of section model",
				"DispBeamColumn2d::DispBeamColumn2d");
	}

	crdTransf = coordTransf.getCopy();

	if (crdTransf == 0)
	    g3ErrorHandler->fatal("%s - failed to copy coordinate transformation",
			"DispBeamColumn2d::DispBeamColumn2d");

	// Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;

// AddingSensitivity:BEGIN /////////////////////////////////////
	gradientIdentifier = 0;
	gradientSectionTag = 0;
	gradientMaterialTag = 0;
// AddingSensitivity:END //////////////////////////////////////
}

DispBeamColumn2d::DispBeamColumn2d()
:Element (0, ELE_TAG_DispBeamColumn2d), rho(0.0),
 Q(6), q(3), connectedExternalNodes(2), L(0.0),
 numSections(0), theSections(0), crdTransf(0)
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	gradientIdentifier = 0;
	gradientSectionTag = 0;
	gradientMaterialTag = 0;
// AddingSensitivity:END //////////////////////////////////////
}

DispBeamColumn2d::~DispBeamColumn2d()
{    
    for (int i = 0; i < numSections; i++) {
		if (theSections[i])
			delete theSections[i];
	}

    // Delete the array of pointers to SectionForceDeformation pointer arrays
    if (theSections)
		delete [] theSections;

}

int
DispBeamColumn2d::getNumExternalNodes() const
{
    return 2;
}

const ID&
DispBeamColumn2d::getExternalNodes()
{
    return connectedExternalNodes;
}

int
DispBeamColumn2d::getNumDOF()
{
    return 6;
}

void
DispBeamColumn2d::setDomain(Domain *theDomain)
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
	//g3ErrorHandler->fatal("FATAL ERROR DispBeamColumn2d (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = nd1Ptr->getNumberDOF();
    int dofNd2 = nd2Ptr->getNumberDOF();
    
    if (dofNd1 != 3 || dofNd2 != 3) {
	//g3ErrorHandler->fatal("FATAL ERROR DispBeamColumn2d (tag: %d), has differing number of DOFs at its nodes",
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
DispBeamColumn2d::commitState()
{
    int retVal = 0;

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int
DispBeamColumn2d::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
DispBeamColumn2d::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
DispBeamColumn2d::update(void)
{
	// Update the transformation
	crdTransf->update();

    // Get basic deformations
    static Vector v(3);
	v = crdTransf->getBasicTrialDisp();

	double oneOverL = 1.0/L;
	GaussQuadRule1d01 quadrat(numSections);
	const Matrix &pts = quadrat.getIntegrPointCoords();
	const Vector &wts = quadrat.getIntegrPointWeights();

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
DispBeamColumn2d::getTangentStiff()
{
	static Matrix kb(3,3);

	// Zero for integral
	kb.Zero();
	q.Zero();

	GaussQuadRule1d01 quadrat(numSections);
	const Matrix &pts = quadrat.getIntegrPointCoords();
	const Vector &wts = quadrat.getIntegrPointWeights();

	// Assuming member is prismatic ... have to move inside
	// the loop if it is not prismatic
	int order = theSections[0]->getOrder();
	const ID &code = theSections[0]->getType();

	double oneOverL = 1.0/L;
	Matrix ka(workArea, order, 3);

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
			default:
				break;
			}
		}
		for (j = 0; j < order; j++) {
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				for (k = 0; k < 3; k++)
					kb(0,k) += ka(j,k);
				break;
			case SECTION_RESPONSE_MZ:
				for (k = 0; k < 3; k++) {
					tmp = ka(j,k);
					kb(1,k) += (xi6-4.0)*tmp;
					kb(2,k) += (xi6-2.0)*tmp;
				}
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
DispBeamColumn2d::getSecantStiff()
{
	return this->getTangentStiff();
}

const Matrix&
DispBeamColumn2d::getDamp()
{
	K.Zero();
	
	return K;
}

const Matrix&
DispBeamColumn2d::getMass()
{
	K.Zero();

	if (rho == 0.0)
		return K;

	double m = 0.5*rho*L;

	K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;

	return K;
}

void
DispBeamColumn2d::zeroLoad(void)
{
	Q(0) = 0.0;
	Q(1) = 0.0;
	Q(2) = 0.0;
	Q(3) = 0.0;
	Q(4) = 0.0;
	Q(5) = 0.0;

	return;
}

int 
DispBeamColumn2d::addLoad(const Vector &addLoad)
{
	if (addLoad.Size() != 6) {
		g3ErrorHandler->warning("DispBeamColumn2d::addLoad %s\n",
				"Vector not of correct size");
		return -1;
	}

	// Add to the external nodal loads
	//Q.addVector(1.0, addLoad, 1.0);
	Q(0) += addLoad(0);
	Q(1) += addLoad(1);
	Q(2) += addLoad(2);
	Q(3) += addLoad(3);
	Q(4) += addLoad(4);
	Q(5) += addLoad(5);

	return 0;
}

int 
DispBeamColumn2d::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0) 
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = nd1Ptr->getRV(accel);
	const Vector &Raccel2 = nd2Ptr->getRV(accel);

    if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
		g3ErrorHandler->warning("DispBeamColumn2d::addInertiaLoadToUnbalance %s\n",
				"matrix and vector sizes are incompatable");
		return -1;
    }

	double m = 0.5*rho*L;

    // Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	Q(0) -= m*Raccel1(0);
	Q(1) -= m*Raccel1(1);
	Q(3) -= m*Raccel2(0);
	Q(4) -= m*Raccel2(1);

    return 0;
}

const Vector&
DispBeamColumn2d::getResistingForce()
{
	GaussQuadRule1d01 quadrat(numSections);
	const Matrix &pts = quadrat.getIntegrPointCoords();
	const Vector &wts = quadrat.getIntegrPointWeights();

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
			default:
				break;
			}
		}

	}

	// Transform forces
	static Vector dummy(2);		// No distributed loads
	P = crdTransf->getGlobalResistingForce(q,dummy);

	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P.addVector(1.0, Q, -1.0);
	P(0) -= Q(0);
	P(1) -= Q(1);
	P(2) -= Q(2);
	P(3) -= Q(3);
	P(4) -= Q(4);
	P(5) -= Q(5);

	return P;
}

const Vector&
DispBeamColumn2d::getResistingForceIncInertia()
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
	P(3) += m*accel2(0);
	P(4) += m*accel2(1);

	return P;
}

int
DispBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
DispBeamColumn2d::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
	return -1;
}

void
DispBeamColumn2d::Print(ostream &s, int flag)
{
	s << "\nDispBeamColumn2d, element id:  " << this->getTag() << endl;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tmass density:  " << rho << endl;
	theSections[0]->Print(s,flag);
}


int
DispBeamColumn2d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = nd1Ptr->getCrds();
    const Vector &end2Crd = nd2Ptr->getCrds();	

    const Vector &end1Disp = nd1Ptr->getDisp();
    const Vector &end2Disp = nd2Ptr->getDisp();

	static Vector v1(3);
	static Vector v2(3);

	for (int i = 0; i < 2; i++) {
		v1(i) = end1Crd(i) + end1Disp(i)*fact;
		v2(i) = end2Crd(i) + end2Disp(i)*fact;    
	}
	
	return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
DispBeamColumn2d::setResponse(char **argv, int argc, Information &eleInfo)
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
DispBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
  double V;

  switch (responseID) {
    case 1:  // global forces
      return eleInfo.setVector(P);

    case 2:
      P(3) = q(0);
      P(0) = -q(0);
      P(2) = q(1);
      P(5) = q(2);
      V = (q(1)+q(2))/L;
      P(1) = V;
      P(4) = -V;
      return eleInfo.setVector(P);
      
    default: 
	  return -1;
  }
}



int
DispBeamColumn2d::setParameter (char **argv, int argc, Information &info)
{
	//
	// From the parameterID value it should be possible to extract
	// information about:
	//  1) Which parameter is in question. The parameter could
	//     be at element, section, or material level. 
	//  2) Which section and material number (tag) it belongs to. 
	//
	// To accomplish this the parameterID is given the following value:
	//     parameterID = type + 1000*matrTag + 100000*sectionTag
	// ...where 'type' is an integer in the range (1-99) and added 100
	// for each level (from material to section to element). 
	//
	// Example:
	//    If 'E0' (case 2) is random in material #3 of section #5
	//    the value of the parameterID at this (element) level would be:
	//    parameterID = 2 + 1000*3 + 100000*5 = 503002
	//    As seen, all given information can be extracted from this number. 
	//

	// Initial declarations
	int parameterID;

	// If the parameter belongs to the element itself
	if (strcmp(argv[0],"rho") == 0) {
		info.theType = DoubleType;
		return 1;
	}

	// If the parameter is belonging to a section or lower
	else if (strcmp(argv[0],"section") == 0) {

		// For now, no parameters of the section itself:
		if (argc<5) {
			cerr << "For now: cannot handle parameters of the section itself." << endl;
			return -1;
		}

		// Get section and material tag numbers from user input
		int paramSectionTag = atoi(argv[1]);

		// Find the right section and call its setParameter method
		for (int i=0; i<numSections; i++) {
			if (paramSectionTag == theSections[i]->getTag()) {
				parameterID = theSections[i]->setParameter(&argv[2], argc-2, info);
			}
		}
		
		// Check if the parameterID is valid
		if (parameterID < 0) {
			cerr << "DispBeamColumn2d::setParameter() - could not set parameter. " << endl;
			return -1;
		}
		else {
			// Return the parameterID value (according to the above comments)
			return parameterID;
		}
	}
    
	// Otherwise parameter is unknown for this class
	else {
		return -1;
	}
}

int
DispBeamColumn2d::updateParameter (int parameterID, Information &info)
{
	// If the parameterID value is not equal to 1 it belongs 
	// to section or material further down in the hierarchy. 

	if (parameterID == 1) {

		this->rho = info.theDouble;
		return 0;

	}
	else if (parameterID > 0 ) {

		// Extract the section number
		int sectionNumber = (int)( floor((double)parameterID) / (100000) );

		int ok = -1;
		for (int i=0; i<numSections; i++) {
			if (sectionNumber == theSections[i]->getTag()) {
				ok = theSections[i]->updateParameter(parameterID, info);
			}
		}

		if (ok < 0) {
			cerr << "DispBeamColumn2d::updateParameter() - could not update parameter. " << endl;
			return ok;
		}
		else {
			return ok;
		}
	}
	else {
		cerr << "DispBeamColumn2d::updateParameter() - could not update parameter. " << endl;
		return -1;
	}       
}

// AddingSensitivity:BEGIN ///////////////////////////////////
const Vector &
DispBeamColumn2d::gradient(bool compute, int identifier)
{
/*	The gradient method can be called with four different purposes:
	1) To clear the sensitivity flag so that the object does not contribute:
			gradient(false, 0)
	2) To set the sensitivity flag so that the object contributes
	   (the sensitivity flag is stored as the value of parameterID):
			gradient(false, parameterID)
	3) To obtain the gradient vector from the object (like for the residual):
			gradient(true, 0)
	4) To commit unconditional sensitivities for path-dependent problems:
			gradient(true, gradNumber)
*/

	if (compute) { // If yes: compute or commit gradients

		if ( gradientIdentifier != 0 ) { // Check if this element contributes has a parameter
			
			if (identifier == 0) { // "Phase 1": return gradient vector

				GaussQuadRule1d01 quadrat(numSections);
				const Matrix &pts = quadrat.getIntegrPointCoords();
				const Vector &wts = quadrat.getIntegrPointWeights();

				// Assuming member is prismatic ... have to move inside
				// the loop if it is not prismatic
				int order = theSections[0]->getOrder();
				const ID &code = theSections[0]->getType();

				// Zero for integration
				q.Zero();

				// Loop over the integration points
				for (int i = 0; i < numSections; i++) {

					double xi6 = 6.0*pts(i,0);

					// Get section stress resultant gradient
					Vector s(order);
					int ok = theSections[i]->gradient(true,identifier,s);

					// Perform numerical integration on internal force gradient
					//q.addMatrixTransposeVector(1.0, *B, s, wts(i));

					double si;
					for (int j = 0; j < order; j++) {
						si = s(j)*wts(i);
						switch(code(j)) {
						case SECTION_RESPONSE_P:
							q(0) += si; break;
						case SECTION_RESPONSE_MZ:
							q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si; break;
						default:
							break;
						}
					}

				}

				// Transform forces
				static Vector dummy(2);		// No distributed loads
				P = crdTransf->getGlobalResistingForce(q,dummy);

			}

			else { // "Phase 2": commit unconditional sensitivities

				// Path-dependent problems: Not treated here yet
			}

		}

		else {

			// Return zero if gradientIdentifier is zero
			P.Zero();
		}
	}
	
	else { // Just set private data flags

		// Set gradient flag
		gradientIdentifier = identifier;

		if (identifier == 0 ) {

			// "Zero out" all flags downwards through sections/materials 
			Vector dummy(2);
			int ok = -1;
			for (int i=0; i<numSections; i++) {
				ok = theSections[i]->gradient(false, identifier, dummy);
			}
		}

		else if (gradientIdentifier == 1) {
			// Don't treat the 'rho' for now
		}

		else {

			// Extract section and material tags
			gradientSectionTag = (int)( floor((double)identifier) / (100000) );
			int tempIdentifier = identifier - gradientSectionTag*100000;
			gradientMaterialTag = (int)( floor((double)tempIdentifier) / (1000) );

			// Go down to the sections and set appropriate flags
			Vector dummy(2);
			int ok = -1;
			for (int i=0; i<numSections; i++) {
				if (gradientSectionTag == theSections[i]->getTag()) {
					ok = theSections[i]->gradient(false, identifier, dummy);
				}
			}
		}
	}

	return P;
}
// AddingSensitivity:END /////////////////////////////////////////////

