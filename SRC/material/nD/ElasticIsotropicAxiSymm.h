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
// $Date: 2001-07-16 22:59:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicAxiSymm.h,v $

#ifndef ElasticIsotropicAxiSymm_h
#define ElasticIsotropicAxiSymm_h

// File: ~/material/ElasticIsotropicAxiSymm.h
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: 
//
// What: "@(#) ElasticIsotropicAxiSymm.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <Tensor.h>

class ElasticIsotropicAxiSymm : public ElasticIsotropicMaterial
{
  public:
    ElasticIsotropicAxiSymm (int tag, double E, double nu, double rho);
	ElasticIsotropicAxiSymm ();
    ~ElasticIsotropicAxiSymm ();

    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);
    const Matrix &getTangent (void);
    const Vector &getStress (void);
	const Vector &getStrain (void);
    
	int commitState (void);
	int revertToLastCommit (void);
    int revertToStart (void);
    
	NDMaterial *getCopy (void);
	const char *getType (void) const;
	int getOrder (void) const;

	int sendSelf(int commitTag, Channel &theChannel);  
	int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    

  protected:

  private:
  	static Vector sigma;	// Stress vector ... class-wide for returns
	static Matrix D;		// Elastic constants
	Vector Tepsilon;	// Trial strains
	Vector Cepsilon;	// Committed strains
};


#endif


