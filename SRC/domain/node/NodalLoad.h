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
// $Date: 2001-06-14 06:21:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/NodalLoad.h,v $
                                                                        
                                                                        
#ifndef NodalLoad_h
#define NodalLoad_h

// File: ~/domain/node/NodalLoad.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class interface for NodalLoad.
// NodalLoad is a class for appling nodal loads to the model.

#include <Load.h>
#include <Node.h>
#include <Vector.h>

class NodalLoad : public Load
{
  public:
    NodalLoad(int classTag);
    NodalLoad(int tag, int node, int classTag);
    NodalLoad(int tag, int node, const Vector &load, bool isLoadConstant = false);
    ~NodalLoad();

    virtual void setDomain(Domain *newDomain);
    virtual int getNodeTag(void) const;
    virtual void applyLoad(double loadFactor);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    
    virtual void Print(ostream &s, int flag =0);   
    
    int setParameter (char **argv, int argc, Information &info);
    int updateParameter (int parameterID, Information &info);

  protected:

  private:
    int  myNode;        // tag indicating associated Node objects tag
    Node *myNodePtr;    // pointer to Node object on which load acts
    Vector *load;       // the reference load - pointer to new copy or 0
    bool  konstant;     // true if load is load factor independent

};

#endif

