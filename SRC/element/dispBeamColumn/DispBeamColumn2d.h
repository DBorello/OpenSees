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
// $Date: 2002-05-16 00:07:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/DispBeamColumn2d.h,v $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumn2d.
// The element displacement field gives rise to constant axial strain and
// linear curvature.

#ifndef DispBeamColumn2d_h
#define DispBeamColumn2d_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <GaussQuadRule1d01.h>

class Node;
class SectionForceDeformation;
class CrdTransf2d;
class Response;

class DispBeamColumn2d : public Element
{
  public:
    DispBeamColumn2d(int tag, int nd1, int nd2,
		  int numSections, SectionForceDeformation &s,
		  CrdTransf2d &coordTransf, double rho = 0.0);
    DispBeamColumn2d();
    virtual ~DispBeamColumn2d();

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);
    void Print(ostream &s, int flag =0);

	Response *setResponse(char **argv, int argc, Information &eleInfo);
    int getResponse(int responseID, Information &eleInfo);

	int setParameter(char **argv, int argc, Information &info);
	int updateParameter(int parameterID, Information &info);

// AddingSensitivity:BEGIN //////////////////////////////////////////
	const Vector & gradient(bool compute, int identifier);
// AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
	int numSections;
    SectionForceDeformation **theSections; // pointer to the ND material objects
	CrdTransf2d *crdTransf;        // pointer to coordinate tranformation object 

    ID connectedExternalNodes; // Tags of quad nodes
	double L;

    Node *nd1Ptr;		// Pointers to quad nodes
    Node *nd2Ptr;

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
	Vector Q;		// Applied nodal loads
	Vector q;		// Basic force
    double rho;			// Mass density per unit length

	static double workArea[];

	static GaussQuadRule1d01 quadRule;

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int gradientIdentifier;
	int gradientSectionTag;
	int gradientMaterialTag;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif

