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
                                                                        
// $Revision: 1.9 $
// $Date: 2002-05-16 00:07:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beamWithHinges/BeamWithHinges2d.cpp,v $

#include <BeamWithHinges2d.h>
#include <Element.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Matrix.h>
#include <Vector.h>						
#include <Node.h>
#include <MatrixUtil.h>
#include <math.h>
#include <stdlib.h>
#include <iostream.h>
#include <string.h>
#include <stdio.h>

#include <SectionForceDeformation.h>
#include <CrdTransf2d.h>

#include <Information.h>
#include <ElementResponse.h>
#include <Renderer.h>

Matrix BeamWithHinges2d::theMatrix(6,6);
Vector BeamWithHinges2d::theVector(6);
double BeamWithHinges2d::workArea[100];

BeamWithHinges2d::BeamWithHinges2d(void)
  :Element(0, ELE_TAG_BeamWithHinges2d),
   E(0.0), A(0.0), I(0.0),
   lp1(0.0), lp2(0.0), rho(0.0),
   theCoordTransf(0),
   connectedExternalNodes(2),
   node1Ptr(0), node2Ptr(0),
   kb(3,3), q(3), load(6),
   kbCommit(3,3), qCommit(3),
   initialFlag(0), maxIter(0), tolerance(0.0)
{
  section[0] = 0;
  section[1] = 0;
}

BeamWithHinges2d::BeamWithHinges2d(int tag, int nodeI, int nodeJ,
				   double e, double a, double i,
				   SectionForceDeformation &sectionRefI, double lpi,
				   SectionForceDeformation &sectionRefJ, double lpj,
				   CrdTransf2d &coordTransf,
				   double r, int max, double tol)
  :Element(tag, ELE_TAG_BeamWithHinges2d),
   E(e), A(a), I(i),
   lp1(lpi), lp2(lpj), rho(r),
   theCoordTransf(0),
   connectedExternalNodes(2),
   node1Ptr(0), node2Ptr(0),
   kb(3,3), q(3), load(6),
   kbCommit(3,3), qCommit(3),
   initialFlag(0), maxIter(max), tolerance(tol)
{
  if (E <= 0.0)  {
    g3ErrorHandler->fatal("%s -- input parameter E is <= 0.0",
			  "BeamWithHinges2d::BeamWithHinges2d");
  }
  
  if (I <= 0.0)  {
    g3ErrorHandler->fatal("%s -- input parameter I is <= 0.0",
			  "BeamWithHinges2d::BeamWithHinges2d");
  }
  
  if (A <= 0.0)  {
    g3ErrorHandler->fatal("%s -- input parameter A is <= 0.0",
			  "BeamWithHinges2d::BeamWithHinges2d");
  }
  
  // Get copies of sections
  section[0] = sectionRefI.getCopy();
  
  if (section[0] == 0)
    g3ErrorHandler->fatal("%s -- failed to get copy of section I",
			  "BeamWithHinges2d::BeamWithHinges2d");
  
  section[1] = sectionRefJ.getCopy();
  
  if (section[1] == 0)
    g3ErrorHandler->fatal("%s -- failed to get copy of section J",
			  "BeamWithHinges2d::BeamWithHinges2d");
  
  theCoordTransf = coordTransf.getCopy();
  
  if (theCoordTransf == 0)
    g3ErrorHandler->fatal("%s -- failed to get copy of coordinate transformation",
			  "BeamWithHinges2d::BeamWithHinges2d");
  
  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;
}

BeamWithHinges2d::~BeamWithHinges2d(void)
{
  for (int i = 0; i < 2; i++)
    if (section[i] != 0)
      delete section[i];
  
  if (theCoordTransf)
    delete theCoordTransf;
}

int 
BeamWithHinges2d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
BeamWithHinges2d::getExternalNodes(void)
{
  return connectedExternalNodes;
}

int 
BeamWithHinges2d::getNumDOF(void)
{
  return 6;
}

void
BeamWithHinges2d::setDomain(Domain *theDomain)
{
  //This function may be called after a beam is constructed, so
  //geometry may change.  Therefore calculate all beam geometry here.
  
  if(theDomain == 0) {
    node1Ptr = 0;
    node2Ptr = 0;
    return;
  }
  
  // set the node pointers and verify them
  this->setNodePtrs(theDomain);
  
  // call the DomainComponent version of the function
  this->DomainComponent::setDomain(theDomain);
  
  if (theCoordTransf->initialize(node1Ptr, node2Ptr) != 0)
    g3ErrorHandler->fatal("%s -- failed to initialize coordinate transformation",
			  "BeamWithHinges2d::setDomain()");
  
  // get element length
  double L = theCoordTransf->getInitialLength();
  if (L == 0.0)
    g3ErrorHandler->fatal("%s -- element has zero length",
			  "BeamWithHinges2d::setDomain()");

  // Set up section interpolation and hinge lengths
  this->setHinges();

  if (initialFlag == 2)
    theCoordTransf->update();
  else
    this->update();
}

int
BeamWithHinges2d::commitState(void)
{
  int err = 0;
  
  for (int i = 0; i < 2; i++) {
    if (section[i] != 0)
      err += section[i]->commitState();
  }
  
  err += theCoordTransf->commitState();
  
  kbCommit = kb;
  qCommit = q;
  
  //initialFlag = 0;
  
  return err;
}

int
BeamWithHinges2d::revertToLastCommit(void)
{
  int err = 0;
  
  // Revert the sections and then get their last commited
  // deformations, stress resultants, and flexibilities
  for (int i = 0; i < 2; i++) {
    if (section[i] != 0) {
      err += section[i]->revertToLastCommit();
      e[i] = section[i]->getSectionDeformation();
      sr[i] = section[i]->getStressResultant();
      fs[i] = section[i]->getSectionFlexibility();
    }
  }
  
  // Commit the coordinate transformation
  err += theCoordTransf->revertToLastCommit();
  
  kb = kbCommit;
  q = qCommit;
  
  initialFlag = 0;
  this->update();

  return err;
}

int
BeamWithHinges2d::revertToStart(void)
{
  int err = 0;
  
  for (int i = 0; i < 2; i++) {
    if (section[i] != 0) {
      err += section[i]->revertToStart();
      fs[i].Zero();
      e[i].Zero();
      sr[i].Zero();
    }
  }
  
  err += theCoordTransf->revertToStart();
  
  kb.Zero();
  q.Zero();
  
  initialFlag = 0;
  this->update();

  return err;
}

const Matrix &
BeamWithHinges2d::getTangentStiff(void)
{
  // Will remove once we clean up the corotational 2d transformation -- MHS
  theCoordTransf->update();

  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &
BeamWithHinges2d::getDamp(void)
{
  theMatrix.Zero();

  return theMatrix;
}

const Matrix &
BeamWithHinges2d::getMass(void)
{
  theMatrix.Zero();

  if (rho != 0.0) {
    double L = theCoordTransf->getInitialLength();  
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(3,3) = theMatrix(4,4) = 0.5*L*rho;
  }
  
  return theMatrix;
}

void 
BeamWithHinges2d::zeroLoad(void)
{
  load.Zero();
}

int
BeamWithHinges2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  g3ErrorHandler->warning("%s -- load type unknown for ele with tag: %d",
			  "BeamWithHinges2d::addLoad", this->getTag());
  
  return -1;
}

int
BeamWithHinges2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;
  
  const Vector &Raccel1 = node1Ptr->getRV(accel);
  const Vector &Raccel2 = node2Ptr->getRV(accel);
  
  double L = theCoordTransf->getInitialLength();
  double mass = 0.5*L*rho;
  
  int i,j;
  for (i = 0, j = 3; i < 2; i++, j++) {
    load(i) += mass*Raccel1(i);
    load(j) += mass*Raccel2(i);	// Yes, this should be 'i'
  }
  
  return 0;
}

const Vector &
BeamWithHinges2d::getResistingForce(void)
{
  static Vector dummy(2);
  
  // Will remove once we clean up the corotational 2d transformation -- MHS
  theCoordTransf->update();

  return theCoordTransf->getGlobalResistingForce (q, dummy);
}

const Vector &
BeamWithHinges2d::getResistingForceIncInertia(void)
{
  if (rho == 0.0)
    return this->getResistingForce();
  
  double ag[6];
  
  const Vector &accel1 = node1Ptr->getTrialAccel();
  const Vector &accel2 = node2Ptr->getTrialAccel();
  
  ag[0] = accel1(0);
  ag[1] = accel1(1);
  //ag[2] = accel1(2); // no rotational element mass
  ag[3] = accel2(0);
  ag[4] = accel2(1);
  //ag[5] = accel2(2); // no rotational element mass
  
  theVector = this->getResistingForce();
  
  double L = theCoordTransf->getInitialLength();
  double mass = 0.5*L*rho;
  
  int i,j;
  for (i = 0, j = 3; i < 2; i++, j++) {
    theVector(i) += mass*ag[i];
    theVector(j) += mass*ag[j];
  }
  
  return theVector;
}

int
BeamWithHinges2d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
BeamWithHinges2d::recvSelf(int commitTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  initialFlag = 2;

  return -1;
}

void 
BeamWithHinges2d::Print(ostream &s, int flag)
{
  s << "\nBeamWithHinges2d, tag: " << this->getTag() << endl;
  s << "\tConnected Nodes: " << connectedExternalNodes;
  s << "\tE: " << E << endl;
  s << "\tA: " << A << endl;
  s << "\tI: " << I << endl;
  
  double P, V, M1, M2;
  double L = theCoordTransf->getInitialLength();
  P = qCommit(0);
  M1 = qCommit(1);
  M2 = qCommit(2);
  V = (M1+M2)/L;

  s << "\tEnd 1 Forces (P V M): "
    << -P << ' ' <<  V << ' ' << M1 << endl;
  s << "\tEnd 2 Forces (P V M): "
    <<  P << ' ' << -V << ' ' << M2 << endl;
  
  if (section[0] != 0) {
    s << "Hinge 1, section tag: " << section[0]->getTag() << 
      ", length: " << lp1 << endl;
    section[0]->Print(s,flag);
  }
  
  if (section[1] != 0) {
    s << "Hinge 2, section tag: " << section[2]->getTag() << 
      ", length: " << lp2 << endl;
    section[1]->Print(s,flag);
  }
}

//////////////////////////////
//Private Function Definitions

void 
BeamWithHinges2d::setNodePtrs(Domain *theDomain)
{
  node1Ptr = theDomain->getNode(connectedExternalNodes(0));
  node2Ptr = theDomain->getNode(connectedExternalNodes(1));
  
  if(node1Ptr == 0) {
    g3ErrorHandler->fatal("%s -- node 1 does not exist",
			  "BeamWithHinges2d::setNodePtrs()");
  }
  
  if(node2Ptr == 0) {
    g3ErrorHandler->fatal("%s -- node 2 does not exist",
			  "BeamWithHinges2d::setNodePtrs()");
  }
  
  // check for correct # of DOF's
  int dofNd1 = node1Ptr->getNumberDOF();
  int dofNd2 = node2Ptr->getNumberDOF();
  if ((dofNd1 != 3) || (dofNd2 != 3))  {
    g3ErrorHandler->fatal("%s -- nodal dof is not three",
			  "BeamWithHinges2d::setNodePtrs()");
  }
}

int
BeamWithHinges2d::update(void)
{
  // Update the coordinate transformation
  theCoordTransf->update();
  
  // Convert to basic system from local coord's (eliminate rb-modes)
  static Vector v(3);				// basic system deformations
  v = theCoordTransf->getBasicTrialDisp();
  
  static Vector dv(3);
  dv = theCoordTransf->getBasicIncrDeltaDisp();
  
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;

  // Section locations along element length ...
  double xi[2];
  // and their integration weights
  double lp[2];
  
  xi[0] = 0.5*lp1;
  xi[1] = L-0.5*lp2;
  
  lp[0] = lp1;
  lp[1] = lp2;

  // element properties
  static Matrix f(3,3);	// element flexibility
  static Vector vr(3);	// Residual element deformations
  
  static Matrix Iden(3,3);   // an identity matrix for matrix inverse
  Iden.Zero();
  for (int i = 0; i < 3; i++)
    Iden(i,i) = 1.0;

  // Length of elastic interior
  double Le = L-lp1-lp2;
  double LoverEA  = Le/(E*A);
  double Lover3EI = Le/(3*E*I);
  double Lover6EI = 0.5*Lover3EI;
  
  // Elastic flexibility of element interior
  static Matrix fe(2,2);
  fe(0,0) = fe(1,1) =  Lover3EI;
  fe(0,1) = fe(1,0) = -Lover6EI;
  
  // Equilibrium transformation matrix
  static Matrix B(2,2);
  double beta1 = lp1*oneOverL;
  double beta2 = lp2*oneOverL;
  B(0,0) = 1.0 - beta1;
  B(1,1) = 1.0 - beta2;
  B(0,1) = -beta1;
  B(1,0) = -beta2;
  
  // Transform the elastic flexibility of the element
  // interior to the basic system
  static Matrix fElastic(2,2);
  fElastic.addMatrixTripleProduct(0.0, B, fe, 1.0);

  // calculate nodal force increments and update nodal forces
  static Vector dq(3);
  //dq = kb * dv;   // using previous stiff matrix k,i
  dq.addMatrixVector(0.0, kb, dv, 1.0);
  
  for (int j = 0; j < maxIter; j++) {
    
    // q += dq;
    q.addVector(1.0, dq, 1.0);
    
    // Set element flexibility to flexibility of elastic region
    f(0,0) = LoverEA;
    f(1,1) = fElastic(0,0);
    f(2,2) = fElastic(1,1);
    f(1,2) = fElastic(0,1);
    f(2,1) = fElastic(1,0);
    f(0,1) = f(1,0) = f(0,2) = f(2,0) = 0.0;
    
    // vr = fElastic * q;
    vr(0) = LoverEA*q(0);
    vr(1) = fElastic(0,0)*q(1) + fElastic(0,1)*q(2);
    vr(2) = fElastic(1,0)*q(1) + fElastic(1,1)*q(2);
    
    for (int i = 0; i < 2; i++) {
      
      if (section[i] == 0 || lp[i] <= 0.0)
	continue;
      
      // Get section information
      int order = section[i]->getOrder();
      const ID &code = section[i]->getType();
      
      Vector s(workArea, order);
      Vector ds(&workArea[order], order);
      Vector de(&workArea[2*order], order);
      
      Matrix fb(&workArea[3*order], order, 3);
      
      double xL = xi[i]*oneOverL;
      double xL1 = xL-1.0;
      
      int ii;
      // Section forces
      // s = b*q + bp*currDistrLoad;
      //this->getForceInterpMatrix(b, xi[i], code);
      //s.addMatrixVector(0.0, b, q, 1.0);
      for (ii = 0; ii < order; ii++) {
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  s(ii) = q(0); break;
	case SECTION_RESPONSE_MZ:
	  s(ii) = xL1*q(1) + xL*q(2); break;
	case SECTION_RESPONSE_VY:
	  s(ii) = oneOverL*(q(1)+q(2)); break;
	default:
	  s(ii) = 0.0; break;
	}
      }
      
      // UNCOMMENT WHEN DISTRIBUTED LOADS ARE ADDED TO INTERFACE
      // this->getDistrLoadInterpMatrix(bp, xi[i], code);
      // s.addMatrixVector(1.0, bp, currDistrLoad, 1.0);
      
      // Increment in section forces
      // ds = s - sr
      ds = s;
      ds.addVector(1.0, sr[i], -1.0);
      
      // compute section deformation increments and update current deformations
      // e += fs * ds;
      de.addMatrixVector(0.0, fs[i], ds, 1.0);
      if (initialFlag != 0)
	e[i].addVector(1.0, de, 1.0);
      
      // set section deformations
      section[i]->setTrialSectionDeformation(e[i]);
      
      // get section resisting forces
      sr[i] = section[i]->getStressResultant();
      
      // get section flexibility matrix
      fs[i] = section[i]->getSectionFlexibility();
      
      // ds = s - sr;
      ds = s;
      ds.addVector(1.0, sr[i], -1.0);
      
      de.addMatrixVector(0.0, fs[i], ds, 1.0);
      
      // integrate section flexibility matrix
      // f += (b^ fs * b) * lp[i];
      //f.addMatrixTripleProduct(1.0, b, fSec, lp[i]);
      int jj;
      fb.Zero();
      double tmp;
      const Matrix &fSec = fs[i];
      for (ii = 0; ii < order; ii++) {
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  for (jj = 0; jj < order; jj++)
	    fb(jj,0) += fSec(jj,ii)*lp[i];
	  break;
	case SECTION_RESPONSE_MZ:
	  for (jj = 0; jj < order; jj++) {
	    tmp = fSec(jj,ii)*lp[i];
	    fb(jj,1) += xL1*tmp;
	    fb(jj,2) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VY:
	  for (jj = 0; jj < order; jj++) {
	    tmp = oneOverL*fSec(jj,ii)*lp[i];
	    fb(jj,1) += tmp;
	    fb(jj,2) += tmp;
	  }
	  break;
	default:
	  break;
	}
      }
      for (ii = 0; ii < order; ii++) {
	switch (code(ii)) {
	case SECTION_RESPONSE_P:
	  for (jj = 0; jj < 3; jj++)
	    f(0,jj) += fb(ii,jj);
	  break;
	case SECTION_RESPONSE_MZ:
	  for (jj = 0; jj < 3; jj++) {
	    tmp = fb(ii,jj);
	    f(1,jj) += xL1*tmp;
	    f(2,jj) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VY:
	  for (jj = 0; jj < 3; jj++) {
	    tmp = oneOverL*fb(ii,jj);
	    f(1,jj) += tmp;
	    f(2,jj) += tmp;
	  }
	  break;
	default:
	  break;
	}
      }
      
      // UNCOMMENT WHEN DISTRIBUTED LOADS ARE ADDED TO INTERFACE
      // vr.addMatrixVector(1.0, vElastic, currDistrLoad, 1.0);
      
      // vr += (b^ (e+de)) * lp[i];
      de.addVector(1.0, e[i], 1.0);
      //vr.addMatrixTransposeVector(1.0, b, de, lp[i]);
      double dei;
      for (ii = 0; ii < order; ii++) {
	dei = de(ii)*lp[i];
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  vr(0) += dei; break;
	case SECTION_RESPONSE_MZ:
	  vr(1) += xL1*dei; vr(2) += xL*dei; break;
	case SECTION_RESPONSE_VY:
	  tmp = oneOverL*dei;
	  vr(1) += tmp; vr(2) += tmp; break;
	default:
	  break;
	}
      }
    }
    
    // calculate element stiffness matrix
    //invert3by3Matrix(f, kb);
    if (f.Solve(Iden,kb) < 0)
      g3ErrorHandler->warning("%s -- could not invert flexibility",
			      "BeamWithHinges2d::update()");    

    // dv = v - vr;
    dv = v;
    dv.addVector(1.0, vr, -1.0);
    
    // determine resisting forces
    // dq = kb * dv;
    dq.addMatrixVector(0.0, kb, dv, 1.0);
    
    double dW = dv^ dq;
    
    if (fabs(dW) < tolerance)
      break;
  }
  
  // q += dq;
  q.addVector(1.0, dq, 1.0);
  
  initialFlag = 1;
  
  return 0;
}

void
BeamWithHinges2d::setHinges(void)
{
  for (int i = 0; i < 2; i++) {
    if (section[i] == 0)
      continue;
    
    // Get the number of section response quantities
    int order = section[i]->getOrder();
    
    fs[i] = Matrix(order,order);
    e[i]  = Vector(order);
    sr[i] = Vector(order);
  }
}

void
BeamWithHinges2d::getForceInterpMatrix(Matrix &b, double x, const ID &code)
{			
  b.Zero();
  
  double L = theCoordTransf->getInitialLength();
  double xi = x/L;
  
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
      b(i,1) = xi - 1.0;
      b(i,2) = xi;
      break;
    case SECTION_RESPONSE_P:		// Axial, P, interpolation
      b(i,0) = 1.0;
      break;
    case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
      b(i,1) = b(i,2) = 1.0/L;
      break;
    default:
      break;
    }
  }
}

void
BeamWithHinges2d::getDistrLoadInterpMatrix(Matrix &bp, double x, const ID & code)
{
  bp.Zero();

  double L = theCoordTransf->getInitialLength();
  double xi = x/L;
  
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
      bp(i,1) = 0.5*xi*(xi-1);
      break;
    case SECTION_RESPONSE_P:		// Axial, P, interpolation
      bp(i,0) = 1.0-xi;
      break;
    case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
      bp(i,1) = xi-0.5;
      break;
    default:
      break;
    }
  }
}

Response*
BeamWithHinges2d::setResponse(char **argv, int argc, Information &info)
{
  // hinge rotations
  if (strcmp(argv[0],"plasticDeformation") == 0 ||
      strcmp(argv[0],"plasticRotation") == 0)
    return new ElementResponse(this, 1, Vector(3));
  
  // global forces
  else if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
	   strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    return new ElementResponse(this, 2, theVector);
  
  // stiffness
  else if (strcmp(argv[0],"stiffness") == 0)
    return new ElementResponse(this, 3, theMatrix);
  
  // local forces
  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    return new ElementResponse(this, 4, theVector);
  
  // section response
  else if (strcmp(argv[0],"section") == 0) {
    int sectionNum = atoi(argv[1]) - 1;
    
    if (sectionNum >= 0 && sectionNum < 2)
      if (section[sectionNum] != 0)
	return section[sectionNum]->setResponse(&argv[2], argc-2, info);
    return 0;
  }
  else
    return 0;
}

int
BeamWithHinges2d::getResponse(int responseID, Information &eleInfo)
{
  double V;
  double L = theCoordTransf->getInitialLength();
  static Vector force(6);
  static Vector def(3);
  
  switch (responseID) {
  case 1: {
    const Vector &v = theCoordTransf->getBasicTrialDisp();
    double LoverEA  = L/(E*A);
    double Lover3EI = L/(3*E*I);
    double Lover6EI = 0.5*Lover3EI;

    double q1 = qCommit(1);
    double q2 = qCommit(2);

    def(0) = v(0) - LoverEA*qCommit(0);
    def(1) = v(1) - Lover3EI*q1 + Lover6EI*q2;
    def(2) = v(2) + Lover6EI*q1 - Lover3EI*q2;
    
    return eleInfo.setVector(def);
  }
  
  case 2: // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 3: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());
    
  case 4: // local forces
    // Axial
    force(3) = q(0);
    force(0) = -q(0);
    // Moment
    force(2) = q(1);
    force(5) = q(2);
    // Shear
    V = (q(1)+q(2))/L;
    force(1) = V;
    force(4) = -V;
    return eleInfo.setVector(force);
    
  default:
    return -1;
  }
}

int
BeamWithHinges2d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  // first determine the end points of the quad based on
  // the display factor (a measure of the distorted image)
  const Vector &end1Crd = node1Ptr->getCrds();
  const Vector &end2Crd = node2Ptr->getCrds();	
  
  const Vector &end1Disp = node1Ptr->getDisp();
  const Vector &end2Disp = node2Ptr->getDisp();
  
  static Vector v1(3);
  static Vector v2(3);
  
  for (int i = 0; i < 2; i++) {
    v1(i) = end1Crd(i) + end1Disp(i)*fact;
    v2(i) = end2Crd(i) + end2Disp(i)*fact;    
  }
  
  return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

int
BeamWithHinges2d::setParameter(char **argv, int argc, Information &info)
{
  // E of the beam interior
  if (strcmp(argv[0],"E") == 0) {
    info.theType = DoubleType;
    return 1;
  }
  
  // A of the beam interior
  else if (strcmp(argv[0],"A") == 0) {
    info.theType = DoubleType;
    return 3;
  }
  
  // I of the beam interior
  else if (strcmp(argv[0],"I") == 0) {
    info.theType = DoubleType;
    return 4;
  }
  
  // Section parameter
  else if (strcmp(argv[0],"section") ==0) {
    if (argc <= 2)
      return -1;
    
    int sectionNum = atoi(argv[1]);
    
    int ok = -1;
    
    if (sectionNum == 1)
      ok = section[0]->setParameter (&argv[2], argc-2, info);
    if (sectionNum == 2)
      ok = section[1]->setParameter (&argv[2], argc-2, info);
    
    if (ok < 0)
      return -1;
    else if (ok < 100)
      return sectionNum*100 + ok;
    else 
      return -1;
  }
  
  // Unknown parameter
  else
    return -1;
}	

int
BeamWithHinges2d::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    this->E = info.theDouble;
    return 0;
  case 3:
    this->A = info.theDouble;
    return 0;
  case 4:
    this->I = info.theDouble;
    return 0;
  default:
    if (parameterID >= 100) { // section quantity
      int sectionNum = parameterID/100; 
      if (sectionNum == 1)
	return section[0]->updateParameter (parameterID-100, info);
      if (sectionNum == 2)
	return section[1]->updateParameter (parameterID-2*100, info);
      else
	return -1;
    }
    else // unknown
      return -1;
  }	
}	
