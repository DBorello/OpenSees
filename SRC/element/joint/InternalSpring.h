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
// $Date: 2002-04-30 21:33:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/InternalSpring.h,v $
                                                                        
// Written: AAA 03/02
// Revised:
//
// Purpose: This file contains the class definition for InternalSpring.
// InternalSpring is a single node element applying a nonlinear uniaxial material 
// to certain degrees of freedom of a node.

#ifndef InternalSpring_h
#define InternalSpring_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>

#include <UniaxialMaterial.h>

class InternalSpring : public Element
{
  public:
    InternalSpring();    
    InternalSpring(int tag, int nod, int ndof, int dof1, int dof2, UniaxialMaterial &spring );
    ~InternalSpring();

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);            
    
    const Matrix &getTangentStiff(void);
    const Matrix &getSecantStiff(void);    
    const Matrix &getDamp(void);    
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    void Print(ostream &s, int flag =0);
	Response *setResponse(char **argv, int argc, Information &eleInformation);
	int getResponse(int responseID, Information &eleInformation);
	
    
	private:
    const Matrix &getStiff(void);    
    
    // private attributes - a copy for each object of the class        
	UniaxialMaterial *Spring; 

    ID  connectedExternalNodes;    
    Node *nodPtr;
	int NDOF,DOF1,DOF2;

    Matrix *Kd; // the stiffness matrix
    Matrix *m; // the mass or damping matrix	

    Vector *load;
	Vector *rForce;
};

#endif

