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
// $Date: 2001-11-27 06:55:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dPointLoad.cpp,v $
                                                                        
// Written: fmk 

// Purpose: This file contains the class implementation Beam2dPointLoad.

#include <Beam2dPointLoad.h>
#include <Vector.h>

Vector Beam2dPointLoad::data(2);

Beam2dPointLoad::Beam2dPointLoad(int tag, double mag, double dist, const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam2dPointLoad, theElementTags), P(mag), x(dist)
{

}

Beam2dPointLoad::Beam2dPointLoad()
  :ElementalLoad(LOAD_TAG_Beam2dPointLoad), P(0.0), x(0.0)
{

}

Beam2dPointLoad::~Beam2dPointLoad()
{

}

const Vector &
Beam2dPointLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dPointLoad;
  data(0) = P;
  data(1) = x;
  return data;
}

int 
Beam2dPointLoad::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
Beam2dPointLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
Beam2dPointLoad::Print(ostream &s, int flag)
{
  s << "Beam2dPointLoad - reference load : " << P << " acting at : " << x << " relative to length\n";
  s << "  elements acted on: " << this->getElementTags();
}
