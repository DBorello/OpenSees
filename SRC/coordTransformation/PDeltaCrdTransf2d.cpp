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
// $Date: 2001-11-26 22:59:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/PDeltaCrdTransf2d.cpp,v $
                                                                        
                                                                        
// File: ~/crdTransf/PDeltaCrdTransf2d.C
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
// 
// Purpose: This file contains the implementation for the 
// PDeltaCrdTransf2d class. PDeltaCrdTransf2d is a linear
// transformation for a planar frame between the global 
// and basic coordinate systems


#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>

#include <iomanip.h>

#include <PDeltaCrdTransf2d.h>

// constructor:
PDeltaCrdTransf2d::PDeltaCrdTransf2d(int tag)
  :CrdTransf2d(tag, CRDTR_TAG_PDeltaCrdTransf2d),
   nodeIPtr(0), nodeJPtr(0),
   nodeIOffset(0), nodeJOffset(0),
   cosTheta(0), sinTheta(0), L(0), ul14(0)
{
	// Does nothing
}

// constructor:
PDeltaCrdTransf2d::PDeltaCrdTransf2d(int tag,
				     const Vector &rigJntOffset1,
				     const Vector &rigJntOffset2)
  :CrdTransf2d(tag, CRDTR_TAG_PDeltaCrdTransf2d),
   nodeIPtr(0), nodeJPtr(0),
   nodeIOffset(0), nodeJOffset(0),
   cosTheta(0), sinTheta(0), L(0), ul14(0)
{
	// check rigid joint offset for node I
	if (&rigJntOffset1 == 0 || rigJntOffset1.Size() != 2 ) {
		cerr << "PDeltaCrdTransf2d::PDeltaCrdTransf2d:  Invalid rigid joint offset vector for node I\n";
		cerr << "Size must be 2\n";      
	}
	else {
		nodeIOffset = new double[2];
		nodeIOffset[0] = rigJntOffset1(0);
		nodeIOffset[1] = rigJntOffset1(1);
	}
   
   // check rigid joint offset for node J
	if (&rigJntOffset2 == 0 || rigJntOffset2.Size() != 2 ) {
		cerr << "PDeltaCrdTransf2d::PDeltaCrdTransf2d:  Invalid rigid joint offset vector for node J\n";
		cerr << "Size must be 2\n";      
	}
	else {
		nodeJOffset = new double[2];
		nodeJOffset[0] = rigJntOffset2(0);
		nodeJOffset[1] = rigJntOffset2(1);
	}
}



 
// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
PDeltaCrdTransf2d::PDeltaCrdTransf2d()
 :CrdTransf2d(0, CRDTR_TAG_PDeltaCrdTransf2d),
  nodeIPtr(0), nodeJPtr(0),
  nodeIOffset(0), nodeJOffset(0),
  cosTheta(0), sinTheta(0), L(0), ul14(0)
{

}



// destructor:
PDeltaCrdTransf2d::~PDeltaCrdTransf2d() 
{
	if (nodeIOffset)
		delete [] nodeIOffset;
	if (nodeJOffset)
		delete [] nodeJOffset;
}


int
PDeltaCrdTransf2d::commitState(void)
{
   return 0;
}


int
PDeltaCrdTransf2d::revertToLastCommit(void)
{
   return 0;
}


int
PDeltaCrdTransf2d::revertToStart(void)
{
   return 0;
}


int 
PDeltaCrdTransf2d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
   int error;

   nodeIPtr = nodeIPointer;
   nodeJPtr = nodeJPointer;

   if ((!nodeIPtr) || (!nodeJPtr))
   {
      cerr << "\nPDeltaCrdTransf2d::initialize";
      cerr << "\ninvalid pointers to the element nodes\n";
      return -1;
   }
       
   // get element length and orientation
   if ((error = this->computeElemtLengthAndOrient()))
      return error;
      
   return 0;
}


int
PDeltaCrdTransf2d::update(void)
{
	const Vector &nodeIDisp = nodeIPtr->getTrialDisp();
	const Vector &nodeJDisp = nodeJPtr->getTrialDisp();

	double ul1;
	double ul4;

	if (nodeIOffset == 0) {
		ul1 = -sinTheta*nodeIDisp(0) + cosTheta*nodeIDisp(1);
		ul4 = -sinTheta*nodeJDisp(0) + cosTheta*nodeJDisp(1);
	}
	else {
		double t12 = sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		double t45 = sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];

		ul1 = -sinTheta*nodeIDisp(0) + cosTheta*nodeIDisp(1) + t12*nodeIDisp(2);
		ul4 = -sinTheta*nodeJDisp(0) + cosTheta*nodeJDisp(1) + t45*nodeJDisp(2);
	}

	ul14 = ul1-ul4;

	return 0;
}


int 
PDeltaCrdTransf2d::computeElemtLengthAndOrient()
{
   // element projection
   static Vector dx(2);

   const Vector &ndICoords = nodeIPtr->getCrds();
   const Vector &ndJCoords = nodeJPtr->getCrds();

   if (nodeIOffset == 0) {
	   dx(0) = ndJCoords(0) - ndICoords(0);
	   dx(1) = ndJCoords(1) - ndICoords(1);
   }
   else {
	   dx(0) = ndJCoords(0) + nodeJOffset[0] - ndICoords(0) - nodeIOffset[0];
	   dx(1) = ndJCoords(1) + nodeJOffset[1] - ndICoords(1) - nodeIOffset[1];
   }
   
   // calculate the element length
   L = dx.Norm();

   if (L == 0.0) 
   {
      cerr << "\nPDeltaCrdTransf2d::computeElemtLengthAndOrien: 0 length\n";
      return -2;  
   }

   // calculate the element local x axis components (direction cossines)
   // wrt to the global coordinates 
   cosTheta = dx(0)/L;
   sinTheta = dx(1)/L;
   
   return 0;
}





double 
PDeltaCrdTransf2d::getInitialLength(void)
{
   return L;
}


double 
PDeltaCrdTransf2d::getDeformedLength(void)
{
   return L;
}


const Vector &
PDeltaCrdTransf2d::getBasicTrialDisp (void)
{
	// determine global displacements
	const Vector &disp1 = nodeIPtr->getTrialDisp();
	const Vector &disp2 = nodeJPtr->getTrialDisp();

	static double ug[6];
	for (int i = 0; i < 3; i++) {
		ug[i]   = disp1(i);
		ug[i+3] = disp2(i);
	}

	static Vector ub(3);

	double oneOverL = 1.0/L;
	double sl = sinTheta*oneOverL;
	double cl = cosTheta*oneOverL;

	if (nodeIOffset == 0) {
		ub(0) = -cosTheta*ug[0] - sinTheta*ug[1] +
			cosTheta*ug[3] + sinTheta*ug[4];

		ub(1) = -sl*ug[0] + cl*ug[1] + ug[2] +
			sl*ug[3] - cl*ug[4];

		//ub(2) = -sl*ug[0] + cl*ug[1] +
		//	sl*ug[3] - cl*ug[4] + ug[5];
		ub(2) = ub(1) + ug[5] - ug[2];
	}
	else {
		double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
		double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
		double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];

		ub(0) = -cosTheta*ug[0] - sinTheta*ug[1] - t02*ug[2] +
			cosTheta*ug[3] + sinTheta*ug[4] + t35*ug[5];

		ub(1) = -sl*ug[0] + cl*ug[1] + (1.0+oneOverL*t12)*ug[2] +
			sl*ug[3] - cl*ug[4] - oneOverL*t45*ug[5];

		//ub(2) = -sl*ug[0] + cl*ug[1] + oneOverL*t12*ug[2] +
		//	sl*ug[3] - cl*ug[4] + (1.0-oneOverL*t45)*ug[5];
		ub(2) = ub(1) + ug[5] - ug[2];
	}

	return ub;
}


const Vector &
PDeltaCrdTransf2d::getBasicIncrDisp (void)
{
	// determine global displacements
	const Vector &disp1 = nodeIPtr->getIncrDisp();
	const Vector &disp2 = nodeJPtr->getIncrDisp();

	static double dug[6];
	for (int i = 0; i < 3; i++) {
		dug[i]   = disp1(i);
		dug[i+3] = disp2(i);
	}

	static Vector dub(3);

	double oneOverL = 1.0/L;
	double sl = sinTheta*oneOverL;
	double cl = cosTheta*oneOverL;

	if (nodeIOffset == 0) {
		dub(0) = -cosTheta*dug[0] - sinTheta*dug[1] +
			cosTheta*dug[3] + sinTheta*dug[4];

		dub(1) = -sl*dug[0] + cl*dug[1] + dug[2] +
			sl*dug[3] - cl*dug[4];

		//dub(2) = -sl*dug[0] + cl*dug[1] +
		//	sl*dug[3] - cl*dug[4] + dug[5];
		dub(2) = dub(1) + dug[5] - dug[2];
	}
	else {
		double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
		double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
		double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];

		dub(0) = -cosTheta*dug[0] - sinTheta*dug[1] - t02*dug[2] +
			cosTheta*dug[3] + sinTheta*dug[4] + t35*dug[5];

		dub(1) = -sl*dug[0] + cl*dug[1] + (1.0+oneOverL*t12)*dug[2] +
			sl*dug[3] - cl*dug[4] - oneOverL*t45*dug[5];

		//dub(2) = -sl*dug[0] + cl*dug[1] + oneOverL*t12*dug[2] +
		//	sl*dug[3] - cl*dug[4] + (1.0-oneOverL*t45)*dug[5];
		dub(2) = dub(1) + dug[5] - dug[2];
	}

	return dub;
}


const Vector &
PDeltaCrdTransf2d::getBasicIncrDeltaDisp(void)
{
	// determine global displacements
	const Vector &disp1 = nodeIPtr->getIncrDeltaDisp();
	const Vector &disp2 = nodeJPtr->getIncrDeltaDisp();

	static double Dug[6];
	for (int i = 0; i < 3; i++) {
		Dug[i]   = disp1(i);
		Dug[i+3] = disp2(i);
	}

	static Vector Dub(3);

	double oneOverL = 1.0/L;
	double sl = sinTheta*oneOverL;
	double cl = cosTheta*oneOverL;

	if (nodeIOffset == 0) {
		Dub(0) = -cosTheta*Dug[0] - sinTheta*Dug[1] +
			cosTheta*Dug[3] + sinTheta*Dug[4];

		Dub(1) = -sl*Dug[0] + cl*Dug[1] + Dug[2] +
			sl*Dug[3] - cl*Dug[4];

		//Dub(2) = -sl*Dug[0] + cl*Dug[1] +
		//	sl*Dug[3] - cl*Dug[4] + Dug[5];
		Dub(2) = Dub(1) + Dug[5] - Dug[2];
	}
	else {
		double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
		double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
		double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];

		Dub(0) = -cosTheta*Dug[0] - sinTheta*Dug[1] - t02*Dug[2] +
			cosTheta*Dug[3] + sinTheta*Dug[4] + t35*Dug[5];

		Dub(1) = -sl*Dug[0] + cl*Dug[1] + (1.0+oneOverL*t12)*Dug[2] +
			sl*Dug[3] - cl*Dug[4] - oneOverL*t45*Dug[5];

		//Dub(2) = -sl*Dug[0] + cl*Dug[1] + oneOverL*t12*Dug[2] +
		//	sl*Dug[3] - cl*Dug[4] + (1.0-oneOverL*t45)*Dug[5];
		Dub(2) = Dub(1) + Dug[5] - Dug[2];
	}

	return Dub;
}


const Vector &
PDeltaCrdTransf2d::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
	// transform resisting forces from the basic system to local coordinates
	static double pl[6];

	double q0 = pb(0);
	double q1 = pb(1);
	double q2 = pb(2);

	double oneOverL = 1.0/L;

	double V = oneOverL*(q1+q2);
	pl[0] = -q0;
	pl[1] =  V;
	pl[2] =  q1;
	pl[3] =  q0;
	pl[4] = -V;
	pl[5] =  q2;

	// add end forces due to element p0 loads
	pl[0] += p0(0);
	pl[1] += p0(1);
	pl[4] += p0(2);
     
	// Include leaning column effects (P-Delta)
	double NoverL = ul14*q0*oneOverL;             
	pl[1] += NoverL;
	pl[4] -= NoverL;

	// transform resisting forces  from local to global coordinates
	static Vector pg(6);

	pg(0) = cosTheta*pl[0] - sinTheta*pl[1];
	pg(1) = sinTheta*pl[0] + cosTheta*pl[1];
   
	pg(3) = cosTheta*pl[3] - sinTheta*pl[4];
	pg(4) = sinTheta*pl[3] + cosTheta*pl[4];
	
	if (nodeIOffset == 0) {
		pg(2) = pl[2];
		pg(5) = pl[5];
	}
	else {
		double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
		double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
		double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];

		pg(2) = t02*pl[0] + t12*pl[1] + pl[2];
		pg(5) = t35*pl[3] + t45*pl[4] + pl[5];
	}

	return pg;
}

const Matrix &
PDeltaCrdTransf2d::getGlobalStiffMatrix (const Matrix &kb, const Vector &pb)
{
	static Matrix kg(6,6);
	static double kl[6][6];
	static double tmp[6][6];

	double oneOverL = 1.0/L;

	// Basic stiffness
	double kb00, kb01, kb02, kb10, kb11, kb12, kb20, kb21, kb22;
	kb00 = kb(0,0);		kb01 = kb(0,1);		kb02 = kb(0,2);
	kb10 = kb(1,0);		kb11 = kb(1,1);		kb12 = kb(1,2);
	kb20 = kb(2,0);		kb21 = kb(2,1);		kb22 = kb(2,2);

	// Transform basic stiffness to local system
	kl[0][0] =  kb00;
	kl[1][0] = -oneOverL*(kb10+kb20);
	kl[2][0] = -kb10;
	kl[3][0] = -kb00;
	kl[4][0] = -kl[1][0];
	kl[5][0] = -kb20;

	kl[0][1] = -oneOverL*(kb01+kb02);
	kl[1][1] =  oneOverL*oneOverL*(kb11+kb12+kb21+kb22);
	kl[2][1] =  oneOverL*(kb11+kb12);
	kl[3][1] = -kl[0][1];
	kl[4][1] = -kl[1][1];
	kl[5][1] =  oneOverL*(kb21+kb22);

	kl[0][2] = -kb01;
	kl[1][2] =  oneOverL*(kb11+kb21);
	kl[2][2] =  kb11;
	kl[3][2] =  kb01;
	kl[4][2] = -kl[1][2];
	kl[5][2] =  kb21;

	kl[0][3] = -kl[0][0];
	kl[1][3] = -kl[1][0];
	kl[2][3] = -kl[2][0];
	kl[3][3] = -kl[3][0];
	kl[4][3] = -kl[4][0];
	kl[5][3] = -kl[5][0];

	kl[0][4] = -kl[0][1];
	kl[1][4] = -kl[1][1];
	kl[2][4] = -kl[2][1];
	kl[3][4] = -kl[3][1];
	kl[4][4] = -kl[4][1];
	kl[5][4] = -kl[5][1];

	kl[0][5] = -kb02;
	kl[1][5] =  oneOverL*(kb12+kb22);
	kl[2][5] =  kb12;
	kl[3][5] =  kb02;
	kl[4][5] = -kl[1][5];
	kl[5][5] =  kb22;

	// Include geometric stiffness effects in local system
	double NoverL = pb(0)*oneOverL;
	kl[1][1] += NoverL;
	kl[4][4] += NoverL;
	kl[1][4] -= NoverL;
	kl[4][1] -= NoverL;

	double t02 = 0.0;
	double t12 = 0.0;
	double t35 = 0.0;
	double t45 = 0.0;

	if (nodeIOffset) {
		t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
		t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
		t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
	}

	// Now transform from local to global ... compute kl*T
	tmp[0][0] = kl[0][0]*cosTheta - kl[0][1]*sinTheta;
	tmp[1][0] = kl[1][0]*cosTheta - kl[1][1]*sinTheta;
	tmp[2][0] = kl[2][0]*cosTheta - kl[2][1]*sinTheta;
	tmp[3][0] = kl[3][0]*cosTheta - kl[3][1]*sinTheta;
	tmp[4][0] = kl[4][0]*cosTheta - kl[4][1]*sinTheta;
	tmp[5][0] = kl[5][0]*cosTheta - kl[5][1]*sinTheta;

	tmp[0][1] = kl[0][0]*sinTheta + kl[0][1]*cosTheta;
	tmp[1][1] = kl[1][0]*sinTheta + kl[1][1]*cosTheta;
	tmp[2][1] = kl[2][0]*sinTheta + kl[2][1]*cosTheta;
	tmp[3][1] = kl[3][0]*sinTheta + kl[3][1]*cosTheta;
	tmp[4][1] = kl[4][0]*sinTheta + kl[4][1]*cosTheta;
	tmp[5][1] = kl[5][0]*sinTheta + kl[5][1]*cosTheta;

	if (nodeIOffset) {
		tmp[0][2] = kl[0][0]*t02 + kl[0][1]*t12 + kl[0][2];
		tmp[1][2] = kl[1][0]*t02 + kl[1][1]*t12 + kl[1][2];
		tmp[2][2] = kl[2][0]*t02 + kl[2][1]*t12 + kl[2][2];
		tmp[3][2] = kl[3][0]*t02 + kl[3][1]*t12 + kl[3][2];
		tmp[4][2] = kl[4][0]*t02 + kl[4][1]*t12 + kl[4][2];
		tmp[5][2] = kl[5][0]*t02 + kl[5][1]*t12 + kl[5][2];
	}
	else {
		tmp[0][2] = kl[0][2];
		tmp[1][2] = kl[1][2];
		tmp[2][2] = kl[2][2];
		tmp[3][2] = kl[3][2];
		tmp[4][2] = kl[4][2];
		tmp[5][2] = kl[5][2];
	}
	
	tmp[0][3] = kl[0][3]*cosTheta - kl[0][4]*sinTheta;
	tmp[1][3] = kl[1][3]*cosTheta - kl[1][4]*sinTheta;
	tmp[2][3] = kl[2][3]*cosTheta - kl[2][4]*sinTheta;
	tmp[3][3] = kl[3][3]*cosTheta - kl[3][4]*sinTheta;
	tmp[4][3] = kl[4][3]*cosTheta - kl[4][4]*sinTheta;
	tmp[5][3] = kl[5][3]*cosTheta - kl[5][4]*sinTheta;

	tmp[0][4] = kl[0][3]*sinTheta + kl[0][4]*cosTheta;
	tmp[1][4] = kl[1][3]*sinTheta + kl[1][4]*cosTheta;
	tmp[2][4] = kl[2][3]*sinTheta + kl[2][4]*cosTheta;
	tmp[3][4] = kl[3][3]*sinTheta + kl[3][4]*cosTheta;
	tmp[4][4] = kl[4][3]*sinTheta + kl[4][4]*cosTheta;
	tmp[5][4] = kl[5][3]*sinTheta + kl[5][4]*cosTheta;

	if (nodeJOffset) {
		tmp[0][5] = kl[0][3]*t35 + kl[0][4]*t45 + kl[0][5];
		tmp[1][5] = kl[1][3]*t35 + kl[1][4]*t45 + kl[1][5];
		tmp[2][5] = kl[2][3]*t35 + kl[2][4]*t45 + kl[2][5];
		tmp[3][5] = kl[3][3]*t35 + kl[3][4]*t45 + kl[3][5];
		tmp[4][5] = kl[4][3]*t35 + kl[4][4]*t45 + kl[4][5];
		tmp[5][5] = kl[5][3]*t35 + kl[5][4]*t45 + kl[5][5];
	}
	else {
		tmp[0][5] = kl[0][5];
		tmp[1][5] = kl[1][5];
		tmp[2][5] = kl[2][5];
		tmp[3][5] = kl[3][5];
		tmp[4][5] = kl[4][5];
		tmp[5][5] = kl[5][5];
	}

	// Now compute T'*(kl*T)
	kg(0,0) = cosTheta*tmp[0][0] - sinTheta*tmp[1][0];
	kg(0,1) = cosTheta*tmp[0][1] - sinTheta*tmp[1][1];
	kg(0,2) = cosTheta*tmp[0][2] - sinTheta*tmp[1][2];
	kg(0,3) = cosTheta*tmp[0][3] - sinTheta*tmp[1][3];
	kg(0,4) = cosTheta*tmp[0][4] - sinTheta*tmp[1][4];
	kg(0,5) = cosTheta*tmp[0][5] - sinTheta*tmp[1][5];

	kg(1,0) = sinTheta*tmp[0][0] + cosTheta*tmp[1][0];
	kg(1,1) = sinTheta*tmp[0][1] + cosTheta*tmp[1][1];
	kg(1,2) = sinTheta*tmp[0][2] + cosTheta*tmp[1][2];
	kg(1,3) = sinTheta*tmp[0][3] + cosTheta*tmp[1][3];
	kg(1,4) = sinTheta*tmp[0][4] + cosTheta*tmp[1][4];
	kg(1,5) = sinTheta*tmp[0][5] + cosTheta*tmp[1][5];

	if (nodeIOffset) {
		kg(2,0) = t02*tmp[0][0] + t12*tmp[1][0] + tmp[2][0];
		kg(2,1) = t02*tmp[0][1] + t12*tmp[1][1] + tmp[2][1];
		kg(2,2) = t02*tmp[0][2] + t12*tmp[1][2] + tmp[2][2];
		kg(2,3) = t02*tmp[0][3] + t12*tmp[1][3] + tmp[2][3];
		kg(2,4) = t02*tmp[0][4] + t12*tmp[1][4] + tmp[2][4];
		kg(2,5) = t02*tmp[0][5] + t12*tmp[1][5] + tmp[2][5];
	}
	else {
		kg(2,0) = tmp[2][0];
		kg(2,1) = tmp[2][1];
		kg(2,2) = tmp[2][2];
		kg(2,3) = tmp[2][3];
		kg(2,4) = tmp[2][4];
		kg(2,5) = tmp[2][5];
	}

	kg(3,0) = cosTheta*tmp[3][0] - sinTheta*tmp[4][0];
	kg(3,1) = cosTheta*tmp[3][1] - sinTheta*tmp[4][1];
	kg(3,2) = cosTheta*tmp[3][2] - sinTheta*tmp[4][2];
	kg(3,3) = cosTheta*tmp[3][3] - sinTheta*tmp[4][3];
	kg(3,4) = cosTheta*tmp[3][4] - sinTheta*tmp[4][4];
	kg(3,5) = cosTheta*tmp[3][5] - sinTheta*tmp[4][5];

	kg(4,0) = sinTheta*tmp[3][0] + cosTheta*tmp[4][0];
	kg(4,1) = sinTheta*tmp[3][1] + cosTheta*tmp[4][1];
	kg(4,2) = sinTheta*tmp[3][2] + cosTheta*tmp[4][2];
	kg(4,3) = sinTheta*tmp[3][3] + cosTheta*tmp[4][3];
	kg(4,4) = sinTheta*tmp[3][4] + cosTheta*tmp[4][4];
	kg(4,5) = sinTheta*tmp[3][5] + cosTheta*tmp[4][5];
	
	if (nodeJOffset) {
		kg(5,0) = t35*tmp[3][0] + t45*tmp[4][0] + tmp[5][0];
		kg(5,1) = t35*tmp[3][1] + t45*tmp[4][1] + tmp[5][1];
		kg(5,2) = t35*tmp[3][2] + t45*tmp[4][2] + tmp[5][2];
		kg(5,3) = t35*tmp[3][3] + t45*tmp[4][3] + tmp[5][3];
		kg(5,4) = t35*tmp[3][4] + t45*tmp[4][4] + tmp[5][4];
		kg(5,5) = t35*tmp[3][5] + t45*tmp[4][5] + tmp[5][5];
	}
	else {
		kg(5,0) = tmp[5][0];
		kg(5,1) = tmp[5][1];
		kg(5,2) = tmp[5][2];
		kg(5,3) = tmp[5][3];
		kg(5,4) = tmp[5][4];
		kg(5,5) = tmp[5][5];
	}

	return kg;
}
  



CrdTransf2d *
PDeltaCrdTransf2d::getCopy(void)
{
  // create a new instance of PDeltaCrdTransf2d 

  PDeltaCrdTransf2d *theCopy;

  if (nodeIOffset) {
	  Vector offsetI(nodeIOffset, 2);
	  Vector offsetJ(nodeJOffset, 2);

	  theCopy = new PDeltaCrdTransf2d(this->getTag(), offsetI, offsetJ);
  }
  else
	  theCopy = new PDeltaCrdTransf2d(this->getTag());

  theCopy->nodeIPtr = nodeIPtr;
  theCopy->nodeJPtr = nodeJPtr;
  theCopy->cosTheta = cosTheta;
  theCopy->sinTheta = sinTheta;
  theCopy->L = L;
  theCopy->ul14 = ul14;
  
  return theCopy;
}


int 
PDeltaCrdTransf2d::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(9);

	data(0) = this->getTag();
	data(6) = L;
	data(7) = cosTheta;
	data(8) = sinTheta;

	res += theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s - failed to send Vector",
			"PDeltaCrdTransf2d::sendSelf");
		return res;
	}

    return res;
}

    

int 
PDeltaCrdTransf2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	static Vector data(9);

	res += theChannel.recvVector(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s - failed to receive Vector",
			"PDeltaCrdTransf2d::recvSelf");
		return res;
	}

	this->setTag((int)data(0));
	L = data(6);
	cosTheta = data(7);
	sinTheta = data(8);

    return res;
}
 	

const Vector &
PDeltaCrdTransf2d::getPointGlobalCoordFromLocal(const Vector &xl)
{
   static Vector xg(2);

   const Vector &nodeICoords = nodeIPtr->getCrds();
   xg(0) = nodeICoords(0);
   xg(1) = nodeICoords(1);

   if (nodeIOffset) {
	   xg(0) += nodeIOffset[0];
	   xg(1) += nodeIOffset[1];
   }

   // xg = xg + Rlj'*xl
   xg(0) += cosTheta*xl(0) - sinTheta*xl(1);
   xg(1) += sinTheta*xl(0) + cosTheta*xl(1);
     
   return xg;  
}

    
const Vector &
PDeltaCrdTransf2d::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
   // determine global displacements
   const Vector &disp1 = nodeIPtr->getTrialDisp();
   const Vector &disp2 = nodeJPtr->getTrialDisp();

   static Vector ug(6);
   for (int i = 0; i < 3; i++)
   {
      ug(i)   = disp1(i);
      ug(i+3) = disp2(i);
   }

   // transform global end displacements to local coordinates
   static Vector ul(6);      // total displacements

	ul(0) =  cosTheta*ug(0) + sinTheta*ug(1);
	ul(1) = -sinTheta*ug(0) + cosTheta*ug(1);
	ul(2) =  ug(2);
	ul(3) =  cosTheta*ug(3) + sinTheta*ug(4);
	ul(4) = -sinTheta*ug(3) + cosTheta*ug(4);
	ul(5) =  ug(5);

   	if (nodeIOffset != 0) {
		double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
		double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
		double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];

		ul(0) += t02*ug(2);
		ul(1) += t12*ug(2);
		ul(3) += t35*ug(5);
		ul(4) += t45*ug(5);
	}
   
   // compute displacements at point xi, in local coordinates
   static Vector uxl(2),  uxg(2);

   uxl(0) = uxb(0) +        ul(0);
   uxl(1) = uxb(1) + (1-xi)*ul(1) + xi*ul(4);

   // rotate displacements to global coordinates
   // uxg = RljT*uxl
   uxg(0) = cosTheta*uxl(0) - sinTheta*uxl(1);
   uxg(1) = sinTheta*uxl(0) + cosTheta*uxl(1);
     
   return uxg;  
}




void
PDeltaCrdTransf2d::Print(ostream &s, int flag)
{
   s << "\nCrdTransf: " << this->getTag() << " Type: PDeltaCrdTransf2d";
   s << "\tnodeI Offset: " << nodeIOffset;
   s << "\tnodeJ Offset: " << nodeJOffset;
}








