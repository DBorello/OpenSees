//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):                 
//#                ``This  source  code  is  Copyrighted in
//#                U.S.,  for  an  indefinite  period,  and
//#                anybody  caught  using  it  without  our
//#                permission,  will be mighty good friends
//#                of  ourn,  cause  we  don't give a darn.
//#                Hack  it.  Compile it. Debug it. Run it.
//#                Yodel  it. Enjoy it. We wrote it, that's
//#                all we wanted to do.''																		
//#              
//#  
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Plastic Bowl (aka Domain Reduction) implementation:
//#                    This file contains the class definition for PBowlLoading.
//#                    PBowlLoading is a subclass of loadPattern, 
//#                    which implements the plastic bowl loading 
//#                    (aka Domain Reduction Method) as described 
//#                    by Jacobo Bielak et al.


//# CLASS:             PBowlLoading
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic
//# PROGRAMMER(S):     Jinxou Liao, Zhaohui Yang, Boris Jeremic
//#
//#
//# DATE:              21Oct2002
//# UPDATE HISTORY:    31Oct2002 fixed some memory leaks
//#                    04Nov2002 changed the way plastic bowl elements are 
//#                    input.  
//#
//#
//#
//#
//#
//===============================================================================

#ifndef PBowlLoading_h
#define PBowlLoading_h

#include <LoadPattern.h>
#include <Matrix.h>

#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <stdlib.h>
#include <Channel.h>
#include <ErrorHandler.h>

#include <fstream.h>
#include <iostream.h>


class Vector;
class Matrix;

class PBowlLoading : public LoadPattern
{
  public:
    PBowlLoading();
    PBowlLoading(int tag,
                 char *PBEfName,
                 char *DispfName,
                 char *AccefName,
                 double theTimeIncr=1.0,
                 double theFactor=1.0);
    ~PBowlLoading();

    void setDomain(Domain *theDomain);
    void applyLoad(double time);
    void Print(ostream &s, int flag =0);

    // methods for o/p
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
            FEM_ObjectBroker &theBroker);

    //  method to obtain a blank copy of the LoadPattern
    LoadPattern *getCopy(void);

  protected:
    //void addPBElements(const ID &PBEle);    //Adding plastic bowl elements
    //void addPBNodes(const ID &PBNodes);     //Adding plastic bowl nodes
    //void addPBLoads(const Matrix &PBLoads); //Adding plastic bowl loades
    void CompPBLoads();        //Finding all plastic bowl nodes and compute the equivalent forces from plastic bowl loading
    const Vector & getNodalLoad(int node, double time); //Getting the nodal load computed from plastic bowl loading corresponding to time

  private:
    ID *PBowlElements;   // vector containing the plastic bowling elements
    ID *PBowlNodes;      // vector containing the plastic bowling nodes
    Matrix *PBowlLoads;  // matrix containing the plastic bowling loads

    Matrix *U;           // vector to store input displ. for all nodes and all time steps
    int UnumDataPoints;   // number of data points
    Matrix *Udd;         // vector to store input accel. for all nodes and all time steps
    int UddnumDataPoints;// number of data points

    int thetimeSteps;
    double PBTimeIncr;   // specifies the time increment used in load path vector
    double cFactor;      // additional factor on the returned load factor

    bool LoadComputed;   // flag to indicate whether the equivalent force has been computed
};

#endif
