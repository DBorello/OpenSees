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
// $Date: 2001-11-26 22:55:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dUniformLoad.cpp,v $
                                                                        

// Written: fmk 
//
// Purpose: This file contains the class implementation of Beam2dUniformLoad.

#include <Beam2dUniformLoad.h>
#include <Vector.h>

Vector Beam2dUniformLoad::data(1);

Beam2dUniformLoad::Beam2dUniformLoad(int tag, double value, const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam2dUniformLoad, theElementTags), w(value)
{

}

Beam2dUniformLoad::Beam2dUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam2dUniformLoad), w(0.0)
{

}

Beam2dUniformLoad::~Beam2dUniformLoad()
{

}

const Vector &
Beam2dUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dUniformLoad;
  data(0) = w;
  return data;
}


int 
Beam2dUniformLoad::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
Beam2dUniformLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
Beam2dUniformLoad::Print(ostream &s, int flag =0)
{
  s << "Beam2dUniformLoad - reference load : " << w << endl;
  s << "  elements acted on: " << this->getElementTags();
}
