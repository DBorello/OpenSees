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
                                                                        
// $Revision: 1.1 $
// $Date: 2007-08-07 16:37:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CorotCrdTransf2d.cpp,v $

// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 05/2000
// Revision: rms 06/2000 (using Assemble, and AssembleTranspose)
// 
// Purpose: This file contains the implementation for the 
// CorotCrdTransf2d class. CorotCrdTransf2d is a Corot
// transformation for a planar frame between the global 
// and basic coordinate systems.


#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>

#include <CorotCrdTransf2d.h>

// initialize static variables
Matrix CorotCrdTransf2d::Tlg(6,6);
Matrix CorotCrdTransf2d::Tbl(3,6);
Vector CorotCrdTransf2d::uxg(3); 
Vector CorotCrdTransf2d::pg(6); 
Vector CorotCrdTransf2d::dub(3); 
Vector CorotCrdTransf2d::Dub(3); 
Matrix CorotCrdTransf2d::kg(6,6);

 
// constructor:
CorotCrdTransf2d::CorotCrdTransf2d(int tag, 
				   const Vector &rigJntOffsetI,
                                   const Vector &rigJntOffsetJ):
  CrdTransf2d(tag, CRDTR_TAG_CorotCrdTransf2d),
  nodeIOffset(2), nodeJOffset(2), cosTheta(0), sinTheta(0),
  cosAlpha(0), sinAlpha(0),
  nodeIPtr(0), nodeJPtr(0), L(0), Ln(0), ub(3), ubcommit(3), ubpr(3)
{
   // check rigid joint offset for node I
   if (&rigJntOffsetI == 0 || rigJntOffsetI.Size() != 2 )
   {
      opserr << "CorotCrdTransf2d::CorotCrdTransf2d:  Invalid rigid joint offset vector for node I\n";
      opserr << "Size must be 2\n";      
      nodeIOffset.Zero();      
   }
   else
     nodeIOffset = rigJntOffsetI;
   
   // check rigid joint offset for node J
   if (&rigJntOffsetJ == 0 || rigJntOffsetJ.Size() != 2 )
   {
      opserr << "CorotCrdTransf2d::CorotCrdTransf2d:  Invalid rigid joint offset vector for node J\n";
      opserr << "Size must be 2\n";      
      nodeJOffset.Zero(); 
   }
   else
     nodeJOffset = rigJntOffsetJ;

   // temporary
   if (nodeIOffset.Norm() != 0 || nodeJOffset.Norm() != 0)
   {
      opserr << "CorotCrdTransf2d::CorotCrdTransf2d: rigid joint zones not implemented yet\n";
      opserr << "Using zero values\n"; 
      nodeIOffset.Zero();
      nodeJOffset.Zero();
   }
}  

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
CorotCrdTransf2d::CorotCrdTransf2d():
  CrdTransf2d(0, CRDTR_TAG_CorotCrdTransf2d),
  nodeIOffset(2), nodeJOffset(2), cosTheta(0), sinTheta(0),
  cosAlpha(0), sinAlpha(0),
  nodeIPtr(0), nodeJPtr(0), L(0), Ln(0), ub(3), ubcommit(3), ubpr(3)
{
}

    
// destructor:
CorotCrdTransf2d::~CorotCrdTransf2d() 
{
}


int
CorotCrdTransf2d::commitState(void)
{
   ubcommit = ub;
   return 0;
}


int
CorotCrdTransf2d::revertToLastCommit(void)
{
   ub = ubcommit;

   this->update();
   return 0;
}


int
CorotCrdTransf2d::revertToStart(void)
{
   ub.Zero();
   this->update();
   return 0;
}


int 
CorotCrdTransf2d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
   int error;

   nodeIPtr = nodeIPointer;
   nodeJPtr = nodeJPointer;

   if ((!nodeIPtr) || (!nodeJPtr))
   {
      opserr << "\nCorotCrdTransf2d::initialize";
      opserr << "\ninvalid pointers to the element nodes\n";
      return -1;
   }

   // get element length and orientation
   if ((error = this->compElemtLengthAndOrient()))
      return error;

   return 0;
}



int  
CorotCrdTransf2d::update(void)
{       
   // get matrix that transforms displacements from global to local coordinates
   this->getTransfMatrixLocalGlobal(Tlg);     // OPTIMIZE LATER

   // get global displacements 
   const Vector &dispI = nodeIPtr->getTrialDisp();
   const Vector &dispJ = nodeJPtr->getTrialDisp();
 
   static Vector ug(6);    
   for (int i = 0; i < 3; i++) {
      ug(i  ) = dispI(i);
      ug(i+3) = dispJ(i);
   }
   
   // transform global end displacements to  local coordinates
   static Vector ul(6);

   // ul.addMatrixVector(0.0, Tlg, ug, 1.0);   // ul = Tlg * ug
   ul(0) = cosTheta*ug(0) + sinTheta*ug(1);
   ul(1) = cosTheta*ug(1) - sinTheta*ug(0);
   ul(2) = ug(2);
   ul(3) = cosTheta*ug(3) + sinTheta*ug(4);
   ul(4) = cosTheta*ug(4) - sinTheta*ug(3);
   ul(5) = ug(5);

   // get deformed element length and orientation with respect to the local system
   this->compElemtLengthAndOrientWRTLocalSystem(ul);
   
   // determine displacements in the basic system eliminating rigid body modes 
   ubpr = ub;
   this->transfLocalDisplsToBasic(ul);

   // compute the transformation matrix from local to the basic system
   this->getTransfMatrixBasicLocal(Tbl);

   return 0;
}
   

int 
CorotCrdTransf2d::compElemtLengthAndOrient(void)
{
   // element projection
   static Vector dx(3);

   dx = (nodeJPtr->getCrds() + nodeJOffset) - (nodeIPtr->getCrds() + nodeIOffset);  
	    
   // calculate the element length

   L = dx.Norm();

   if (L == 0.0) 
   {
      opserr << "\nCorotCrdTransf2d::compElemtLengthAndOrien: 0 length\n";
      return -2;  
   }

   // calculate the element local x axis components (direction cossines)
   // wrt to the global coordinates 
   cosTheta = dx(0)/L;
   sinTheta = dx(1)/L;

   return 0;
}


   
int
CorotCrdTransf2d::compElemtLengthAndOrientWRTLocalSystem(const Vector &ul)
{
   // get length and chord orientation of the deformed element with respect 
   // to the local coordinate system 
   // (deformed chord corresponding to the basic system)   
    
   double dulx, duly, Lc;

   dulx = ul(3) - ul(0);           // horizontal relative displacement  
   duly = ul(4) - ul(1);           // vertical relative displacement

   Lc = L + dulx;                  // horizontal component of the deformed member
   
   Ln = sqrt (Lc*Lc + duly*duly);  // chord length of the deformed configuration
   
   if (Ln == 0.0) {
      opserr << "\nCorotCrdTransf2d::compElemtLengthAndOrientWRTLocalSystem: 0 length\n";
      return -2;  
   }
   
   cosAlpha = Lc / Ln;             // horizontal projection of the displaced chord 
   
   sinAlpha = duly / Ln;            // vertical projection of the displaced chord
   
   return 0;
} 





void CorotCrdTransf2d::transfLocalDisplsToBasic(const Vector &ul)
{
    // Eliminate rigid body modes, determining displacements wrt the basic system
   double alpha;
   
   alpha = atan2 (sinAlpha, cosAlpha);

   ub(0) = Ln - L;
   ub(1) = ul(2) - alpha;
   ub(2) = ul(5) - alpha;
}




void
CorotCrdTransf2d::getTransfMatrixLocalGlobal (Matrix &Tlg) 
{
   // setup transformation matrix
   Tlg.Zero();

   Tlg(0,0) = Tlg(3,3) =  cosTheta;           
   Tlg(0,1) = Tlg(3,4) =  sinTheta;
   Tlg(1,0) = Tlg(4,3) = -sinTheta;
   Tlg(1,1) = Tlg(4,4) =  cosTheta;
   Tlg(2,2) = Tlg(5,5) =  1.0;
}



void
CorotCrdTransf2d::getTransfMatrixBasicLocal(Matrix &Tbl)
{
   // set up exact force transformation matrix from basic to local coordinates

   Tbl(0,0) = -cosAlpha;      
   Tbl(1,0) = -sinAlpha/Ln;
   Tbl(2,0) = -sinAlpha/Ln;
   
   Tbl(0,1) = -sinAlpha;
   Tbl(1,1) =  cosAlpha/Ln;
   Tbl(2,1) =  cosAlpha/Ln;
   
   Tbl(0,2) =  0;
   Tbl(1,2) =  1;  
   Tbl(2,2) =  0;  
   
   Tbl(0,3) =  cosAlpha;
   Tbl(1,3) =  sinAlpha/Ln;
   Tbl(2,3) =  sinAlpha/Ln;
   
   Tbl(0,4) =  sinAlpha;
   Tbl(1,4) = -cosAlpha/Ln;
   Tbl(2,4) = -cosAlpha/Ln;
   
   Tbl(0,5) =  0;
   Tbl(1,5) =  0; 
   Tbl(2,5) =  1;  
}

   

const Vector &
CorotCrdTransf2d::getBasicTrialDisp (void)
{
  return ub;    
}



const Vector &
CorotCrdTransf2d::getBasicIncrDeltaDisp (void)
{
  // dub = ub - ubpr;
   dub = ub;
   dub.addVector (1.0, ubpr, -1.0);
   
   return dub;        
}



const Vector &
CorotCrdTransf2d::getBasicIncrDisp(void)
{
  // Dub = ub - ubcommit;
  Dub = ub;
  Dub.addVector(1.0, ubcommit, -1.0);

   return Dub;        
}



const Vector &
CorotCrdTransf2d::getGlobalResistingForce(const Vector &pb, const Vector &unifLoad)
{

   // transform resisting forces from the basic system to local coordinates
   this->getTransfMatrixBasicLocal(Tbl);
   static Vector pl(6);
   pl.addMatrixTransposeVector(0.0, Tbl, pb, 1.0);    // pl = Tbl ^ pb;

   ///opserr << "pl: " << pl;
   // check distributed load is zero (not implemented yet)
   
   if (unifLoad.Norm() != 0.0)
   {
       opserr << "CorotCrdTransf2d::getGlobalResistingForce: affect of Po not implemented yet.";
       opserr << "using zero value";
   }  
       
   // transform resisting forces  from local to global coordinates
   this->getTransfMatrixLocalGlobal(Tlg);     // OPTIMIZE LATER
   pg.addMatrixTransposeVector(0.0, Tlg, pl, 1.0);   // pg = Tlg ^ pl; residual

   return pg;
}



const Matrix &
CorotCrdTransf2d::getGlobalStiffMatrix (const Matrix &kb, const Vector &pb)
{
   // transform tangent stiffness matrix from the basic system to local coordinates
   static Matrix kl(6,6);
   this->getTransfMatrixBasicLocal(Tbl);
   kl.addMatrixTripleProduct(0.0, Tbl, kb, 1.0);      // kl = Tbl ^ kb * Tbl;

   // add geometric stiffness matrix
   kl.addMatrix(1.0, this->getGeomStiffMatrix(pb), 1.0);
  
   // transform tangent  stiffness matrix from local to global coordinates
   
   // kg.addMatrixTripleProduct(0.0, Tlg, kl, 1.0);
   double s2, c2, cs;  

   s2 = sinTheta*sinTheta;
   c2 = cosTheta*cosTheta;
   cs = sinTheta*cosTheta;
   
   double k11, k12, k13, k21, k22, k23, k31, k32, k33;

   k11 = kl(0,0);    k12 = kl(0,1);    k13 = kl(0,2);
   k21 = kl(1,0);    k22 = kl(1,1);    k23 = kl(1,2);
   k31 = kl(2,0);    k32 = kl(2,1);    k33 = kl(2,2);
   
   kg(0,0) = c2*k11+s2*k22-cs*(k21+k12); 
   kg(1,0) = c2*k21-s2*k12+cs*(k11-k22);
   kg(2,0) = cosTheta*k31-sinTheta*k32;

   kg(0,1) = c2*k12-s2*k21+cs*(k11-k22);
   kg(1,1) = c2*k22+s2*k11+cs*(k21+k12);
   kg(2,1) = sinTheta*k31+cosTheta*k32;

   kg(0,2) = cosTheta*k13 - sinTheta*k23;
   kg(1,2) = sinTheta*k13 + cosTheta*k23;
   kg(2,2) = k33;

   k11 = kl(0,3);    k12 = kl(0,4);    k13 = kl(0,5);
   k21 = kl(1,3);    k22 = kl(1,4);    k23 = kl(1,5);
   k31 = kl(2,3);    k32 = kl(2,4);    k33 = kl(2,5);
   
   kg(0,3) = c2*k11+s2*k22-cs*(k21+k12); 
   kg(1,3) = c2*k21-s2*k12+cs*(k11-k22);
   kg(2,3) = cosTheta*k31-sinTheta*k32;

   kg(0,4) = c2*k12-s2*k21+cs*(k11-k22);
   kg(1,4) = c2*k22+s2*k11+cs*(k21+k12);
   kg(2,4) = sinTheta*k31+cosTheta*k32;

   kg(0,5) = cosTheta*k13 - sinTheta*k23;
   kg(1,5) = sinTheta*k13 + cosTheta*k23;
   kg(2,5) = k33;

   k11 = kl(3,0);    k12 = kl(3,1);    k13 = kl(3,2);
   k21 = kl(4,0);    k22 = kl(4,1);    k23 = kl(4,2);
   k31 = kl(5,0);    k32 = kl(5,1);    k33 = kl(5,2);
   
   kg(3,0) = c2*k11+s2*k22-cs*(k21+k12); 
   kg(4,0) = c2*k21-s2*k12+cs*(k11-k22);
   kg(5,0) = cosTheta*k31-sinTheta*k32;

   kg(3,1) = c2*k12-s2*k21+cs*(k11-k22);
   kg(4,1) = c2*k22+s2*k11+cs*(k21+k12);
   kg(5,1) = sinTheta*k31+cosTheta*k32;

   kg(3,2) = cosTheta*k13 - sinTheta*k23;
   kg(4,2) = sinTheta*k13 + cosTheta*k23;
   kg(5,2) = k33;

   k11 = kl(3,3);    k12 = kl(3,4);    k13 = kl(3,5);
   k21 = kl(4,3);    k22 = kl(4,4);    k23 = kl(4,5);
   k31 = kl(5,3);    k32 = kl(5,4);    k33 = kl(5,5);
   
   kg(3,3) = c2*k11+s2*k22-cs*(k21+k12); 
   kg(4,3) = c2*k21-s2*k12+cs*(k11-k22);
   kg(5,3) = cosTheta*k31-sinTheta*k32;

   kg(3,4) = c2*k12-s2*k21+cs*(k11-k22);
   kg(4,4) = c2*k22+s2*k11+cs*(k21+k12);
   kg(5,4) = sinTheta*k31+cosTheta*k32;

   kg(3,5) = cosTheta*k13 - sinTheta*k23;
   kg(4,5) = sinTheta*k13 + cosTheta*k23;
   kg(5,5) = k33;

   return kg;
}



const Matrix &
CorotCrdTransf2d::getInitialGlobalStiffMatrix (const Matrix &kb)
{
   // transform tangent stiffness matrix from the basic system to local coordinates
   static Matrix kl(6,6);
   static Matrix T(3,6);

   T(0,0) = -1.0;
   T(1,0) = 0;
   T(2,0) = 0;
   
   T(0,1) = 0;
   T(1,1) = 1/L;
   T(2,1) = 1/L;
   
   T(0,2) =  0;
   T(1,2) =  1;  
   T(2,2) =  0;  
   
   T(0,3) =  1;
   T(1,3) =  0;
   T(2,3) =  0;
   
   T(0,4) =  0;
   T(1,4) = -1/L;
   T(2,4) = -1/L;
   
   T(0,5) =  0;
   T(1,5) =  0; 
   T(2,5) =  1;  

   kl.addMatrixTripleProduct(0.0, T, kb, 1.0);      // kl = Tbl ^ kb * Tbl;

   // add geometric stiffness matrix
   // kl.addMatrix(1.0, this->getGeomStiffMatrix(pb), 1.0);
  
   // transform tangent  stiffness matrix from local to global coordinates
   
   // kg.addMatrixTripleProduct(0.0, Tlg, kl, 1.0);
   double s2, c2, cs;  

   s2 = sinTheta*sinTheta;
   c2 = cosTheta*cosTheta;
   cs = sinTheta*cosTheta;
   
   double k11, k12, k13, k21, k22, k23, k31, k32, k33;

   k11 = kl(0,0);    k12 = kl(0,1);    k13 = kl(0,2);
   k21 = kl(1,0);    k22 = kl(1,1);    k23 = kl(1,2);
   k31 = kl(2,0);    k32 = kl(2,1);    k33 = kl(2,2);
   
   kg(0,0) = c2*k11+s2*k22-cs*(k21+k12); 
   kg(1,0) = c2*k21-s2*k12+cs*(k11-k22);
   kg(2,0) = cosTheta*k31-sinTheta*k32;

   kg(0,1) = c2*k12-s2*k21+cs*(k11-k22);
   kg(1,1) = c2*k22+s2*k11+cs*(k21+k12);
   kg(2,1) = sinTheta*k31+cosTheta*k32;

   kg(0,2) = cosTheta*k13 - sinTheta*k23;
   kg(1,2) = sinTheta*k13 + cosTheta*k23;
   kg(2,2) = k33;

   k11 = kl(0,3);    k12 = kl(0,4);    k13 = kl(0,5);
   k21 = kl(1,3);    k22 = kl(1,4);    k23 = kl(1,5);
   k31 = kl(2,3);    k32 = kl(2,4);    k33 = kl(2,5);
   
   kg(0,3) = c2*k11+s2*k22-cs*(k21+k12); 
   kg(1,3) = c2*k21-s2*k12+cs*(k11-k22);
   kg(2,3) = cosTheta*k31-sinTheta*k32;

   kg(0,4) = c2*k12-s2*k21+cs*(k11-k22);
   kg(1,4) = c2*k22+s2*k11+cs*(k21+k12);
   kg(2,4) = sinTheta*k31+cosTheta*k32;

   kg(0,5) = cosTheta*k13 - sinTheta*k23;
   kg(1,5) = sinTheta*k13 + cosTheta*k23;
   kg(2,5) = k33;

   k11 = kl(3,0);    k12 = kl(3,1);    k13 = kl(3,2);
   k21 = kl(4,0);    k22 = kl(4,1);    k23 = kl(4,2);
   k31 = kl(5,0);    k32 = kl(5,1);    k33 = kl(5,2);
   
   kg(3,0) = c2*k11+s2*k22-cs*(k21+k12); 
   kg(4,0) = c2*k21-s2*k12+cs*(k11-k22);
   kg(5,0) = cosTheta*k31-sinTheta*k32;

   kg(3,1) = c2*k12-s2*k21+cs*(k11-k22);
   kg(4,1) = c2*k22+s2*k11+cs*(k21+k12);
   kg(5,1) = sinTheta*k31+cosTheta*k32;

   kg(3,2) = cosTheta*k13 - sinTheta*k23;
   kg(4,2) = sinTheta*k13 + cosTheta*k23;
   kg(5,2) = k33;

   k11 = kl(3,3);    k12 = kl(3,4);    k13 = kl(3,5);
   k21 = kl(4,3);    k22 = kl(4,4);    k23 = kl(4,5);
   k31 = kl(5,3);    k32 = kl(5,4);    k33 = kl(5,5);
   
   kg(3,3) = c2*k11+s2*k22-cs*(k21+k12); 
   kg(4,3) = c2*k21-s2*k12+cs*(k11-k22);
   kg(5,3) = cosTheta*k31-sinTheta*k32;

   kg(3,4) = c2*k12-s2*k21+cs*(k11-k22);
   kg(4,4) = c2*k22+s2*k11+cs*(k21+k12);
   kg(5,4) = sinTheta*k31+cosTheta*k32;

   kg(3,5) = cosTheta*k13 - sinTheta*k23;
   kg(4,5) = sinTheta*k13 + cosTheta*k23;
   kg(5,5) = k33;

   return kg;
}
  



const Matrix &
CorotCrdTransf2d::getGeomStiffMatrix(const Vector &pb) const
{
   // get  geometric stiffness matrix present in the transformation 
   // from basic to local system

   double s2, c2, cs;  

   s2 = sinAlpha*sinAlpha;
   c2 = cosAlpha*cosAlpha;
   cs = sinAlpha*cosAlpha;

   static Matrix kg0(6,6), kg12(6,6);
   kg0.Zero();

   kg12.Zero();
   
   kg0(0,0) = kg0(3,3) =  s2;
   kg0(0,1) = kg0(3,4) = -cs;
   kg0(1,0) = kg0(4,3) = -cs;
   kg0(1,1) = kg0(4,4) =  c2;
	  
   kg0(0,3) = kg0(3,0) = -s2;
   kg0(0,4) = kg0(3,1) =  cs;
   kg0(1,3) = kg0(4,0) =  cs;
   kg0(1,4) = kg0(4,1) = -c2;
	  
   kg0 *= pb(0)/Ln;

   kg12(0,0) = kg12(3,3) = -2*cs;
   kg12(0,1) = kg12(3,4) =  c2-s2;
   kg12(1,0) = kg12(4,3) =  c2-s2;
   kg12(1,1) = kg12(4,4) =  2*cs;
	  
   kg12(0,3) = kg12(3,0) =  2*cs;
   kg12(0,4) = kg12(3,1) = -c2+s2;
   kg12(1,3) = kg12(4,0) = -c2+s2;
   kg12(1,4) = kg12(4,1) = -2*cs;
	  
   kg12 *= (pb(1)+pb(2))/(Ln*Ln);

   static Matrix kg(6,6);
   // kg = kg0 + kg12;
   kg = kg0;
   kg.addMatrix(1.0, kg12, 1.0);
   
   return kg;
}
  




double 
CorotCrdTransf2d::getInitialLength(void)
{
   return L;
}


double 
CorotCrdTransf2d::getDeformedLength(void)
{
   return Ln;
}



CrdTransf2d *
CorotCrdTransf2d::getCopy(void)
{
  // create a new instance of CorotCrdTransf2d 

  CorotCrdTransf2d *theCopy = new CorotCrdTransf2d (this->getTag(), nodeIOffset, nodeJOffset);
  
  if (!theCopy)
  {
    opserr << "CorotCrdTransf2d::getCopy() - out of memory creating copy\n";
    exit(-1);
  }    
    
  theCopy->nodeIPtr = nodeIPtr;
  theCopy->nodeJPtr = nodeJPtr;
  theCopy->cosTheta = cosTheta;
  theCopy->sinTheta = sinTheta;
  theCopy->cosAlpha = cosAlpha;
  theCopy->sinAlpha = sinAlpha;
  theCopy->L = L;
  theCopy->Ln = Ln;
  theCopy->ub = ub;
  theCopy->ubcommit = ubcommit;
      
  return theCopy;
}


int 
CorotCrdTransf2d::sendSelf(int cTag, Channel &theChannel)
{
  Vector data(7);
  data(0) = ubcommit(0);
  data(1) = ubcommit(1);
  data(2) = ubcommit(2);
  data(3) = nodeIOffset(0);
  data(4) = nodeIOffset(1);
  data(5) = nodeJOffset(0);
  data(6) = nodeJOffset(1);

  if (theChannel.sendVector(this->getTag(), cTag, data) < 0) {
    opserr << " CorotCrdTransf2d::sendSelf() - data could not be sent\n" ;
    return -1;
  }
  return 0;
}

    

int 
CorotCrdTransf2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(7);
  if (theChannel.recvVector(this->getTag(), cTag, data) < 0) {
    opserr << " CorotCrdTransf2d::recvSelf() - data could not be received\n" ;
    return -1;
  }

  ubcommit(0) = data(0);
  ubcommit(1) = data(1);
  ubcommit(2) = data(2);
  nodeIOffset(0) =data(4);
  nodeIOffset(1) =data(4);
  nodeJOffset(0) =data(5);
  nodeJOffset(1) =data(6);

  ub = ubcommit;

  return 0;
}
 	
const Vector &
CorotCrdTransf2d::getPointGlobalCoordFromLocal(const Vector &xl)
{
   static Vector xg(3);
   opserr << " CorotCrdTransf2d::getPointGlobalCoordFromLocal: not implemented yet" ;
     
   return xg;  
}

    
const Vector &
CorotCrdTransf2d::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
   opserr << " CorotCrdTransf2d::getPointGlobalDisplFromBasic: not implemented yet" ;
     
   return uxg;  
}




void
CorotCrdTransf2d::Print(OPS_Stream &s, int flag)
{
   s << "\nCrdTransf: " << this->getTag() << " Type: LinearCrdTransf2d";
   s << "\tnodeI Offset: " << nodeIOffset;
   s << "\tnodeJ Offset: " << nodeJOffset;
}
