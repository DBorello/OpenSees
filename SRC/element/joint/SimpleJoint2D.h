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
// $Date: 2002-06-10 22:41:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/SimpleJoint2D.h,v $

// Written: AAA 03/02
// Revised:

// SimpleJoint2d.h: interface for the SimpleJoint2d class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SimpleJoint2D_h
#define SimpleJoint2D_h

#include <bool.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Element.h>
#include <ID.h>
#include <Domain.h>

class UniaxialMaterial;
class Response;

class SimpleJoint2D : public Element  
{
public:
  SimpleJoint2D();

  SimpleJoint2D(int tag, int nd1, int nd2, int nd3, int nd4,
		UniaxialMaterial &Spring, 
		Domain *theDomain,
		int newNodeTag, 
		int newRotEleTag,
		int LrgDispFlag =0);
	
  ~SimpleJoint2D();

  // methods dealing with domain
  int	getNumExternalNodes(void) const;
  const	ID &getExternalNodes(void);
  int	getNumDOF(void);
  void	setDomain(Domain *theDomain);  
  bool	isSubdomain(void) { return false; } ;
	
  // methods dealing with committed state and update
  int update(void);
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  // methods to return the current linearized stiffness,
  // damping and mass matrices
  const	Matrix &getTangentStiff(void);
  const	Matrix &getDamp(void);
  const	Matrix &getMass(void);
	
  // methods for returning and applying loads
  //virtual Vector &getUVLoadVector(double q1, double q2);
  void	zeroLoad(void); 
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const	Vector &getResistingForce(void);
  const	Vector &getResistingForceIncInertia(void);     

  // method for graphics
  int	displaySelf(Renderer &theViewer, int displayMode, float fact);  
	
  // method for obtaining information specific to an element
  Response* setResponse(char **argv, int argc, Information &eleInformation);
  int getResponse(int responseID, Information &eleInformation);
  int sendSelf(int commitTag, Channel &theChannel) {return -1;}
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {return -1;}
  void Print(ostream &s, int flag =0);

 protected:
  int 	addMP_ForJoint(Domain *theDomain, int mpNum, int RnodeID, int CnodeID, int AuxDOF, int LrgDispFlag );

 private:
  ID		ExternalNodes, InternalConstraints;	
  Node	*end1Ptr, *end2Ptr, *end3Ptr, *end4Ptr, *IntNodePtr;
  int		IntNode;
  Domain	*TheDomain;
  int		RotElemtag;
  Element *RotElemPtr;
  int		numDof, nodeRecord, dofRecord;
  static	Matrix K;
  static	Vector V;
};

#endif
