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
// $Date: 2009-05-11 21:08:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/ArpackSOE.cpp,v $

// Written: fmk
// Created: 05/09
//
// Description: This file contains the class definition for ArpackSOE

#include <ArpackSOE.h>
#include <ArpackSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>

ArpackSOE::ArpackSOE(ArpackSolver &theSolvr,  
		     AnalysisModel &aModel, 
		     LinearSOE &theLinSOE,
		     double theShift)
:EigenSOE(theSolvr, EigenSOE_TAGS_ArpackSOE),
 M(0), Msize(0), mDiagonal(false), shift(theShift), 
 theModel(&aModel), theSOE(&theLinSOE)
{
  theSolvr.setEigenSOE(*this);
}


int
ArpackSOE::getNumEqn(void) const
{
  return theSOE->getNumEqn();
}
    
ArpackSOE::~ArpackSOE()
{
  if (M != 0) delete [] M;

  if (theSOE != 0)
    delete theSOE;
}

int 
ArpackSOE::setSize(Graph &theGraph)
{
  theSOE->setSize(theGraph);

  int result = 0;
  int size = theGraph.getNumVertex();

  if (size != Msize && size > 0) {

    if (M != 0) 
      delete [] M;

    M = new double[size];
    
    if (M == 0) {
      opserr << "WARNING ArpackSOE::ArpackSOE : - out of memory creating memory for M\n";
      Msize = 0;
    } else
      Msize = size;
  }

  
  // invoke setSize() on the Solver

  EigenSolver *theSolvr = this->getSolver();
  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    opserr << "WARNING:ArpackSOE::setSize() -  solver failed setSize()\n";
    return solverOK;
  } 

  return result;    
}

int 
ArpackSOE::addA(const Matrix &m, const ID &id, double fact)
{

  // check for a quick return 
  if (fact == 0.0)  return 0;

  return theSOE->addA(m, id, fact);
}


void 
ArpackSOE::zeroA(void)
{
  return theSOE->zeroA();
}

int 
ArpackSOE::addM(const Matrix &m, const ID &id, double fact)
{
  int res = this->addA(m, id, -shift);

  if (res < 0)
    return res;

  if (mDiagonal == false)
    return  res;

  int idSize = id.Size();
  for (int i=0; i<idSize; i++) {
    int locI = id(i);
    if (locI >= 0 && locI < Msize) {
      for (int j=0; j<idSize; j++) {
	int locJ = id(j);
	if (locJ >= 0 && locJ < Msize) {
	  if (locI == locJ) {
	    M[locI] += m(i,i);
	  } else {
	    if (m(i,j) != 0.0) {
	      mDiagonal = false;
	      return res;
	    }
	  }
	}
      }
    }
  }

  return 0;
}   
 
void 
ArpackSOE::zeroM(void)
{
  mDiagonal = true;

  for (int i=0; i<Msize; i++)
    M[i] = 0;
}


double 
ArpackSOE::getShift(void)
{
    return shift;
}


int 
ArpackSOE::sendSelf(int commitTag, Channel &theChannel)
{
  
    return 0;
}
    
int 
ArpackSOE::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
    return 0;
}

