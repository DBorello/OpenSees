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
                                                                        
// $Revision: 1.3 $
// $Date: 2002-06-06 18:24:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dUniformLoad.cpp,v $
                                                                        

// Written: fmk 
//
// Purpose: This file contains the class implementation of Beam2dUniformLoad.

#include <Beam2dUniformLoad.h>
#include <Vector.h>

Vector Beam2dUniformLoad::data(2);

Beam2dUniformLoad::Beam2dUniformLoad(int tag, double wt, double wa,
				     const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam2dUniformLoad, theElementTags),
   wTrans(wt), wAxial(wa)
{

}

Beam2dUniformLoad::Beam2dUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam2dUniformLoad), wTrans(0.0), wAxial(0.0)
{

}

Beam2dUniformLoad::~Beam2dUniformLoad()
{

}

const Vector &
Beam2dUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dUniformLoad;
  data(0) = wTrans;
  data(1) = wAxial;
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
Beam2dUniformLoad::Print(ostream &s, int flag)
{
  s << "Beam2dUniformLoad - Reference load" << endl;
  s << "  Transverse: " << wTrans << endl;
  s << "  Axial:      " << wAxial << endl;
  s << "  Elements acted on: " << this->getElementTags();
}
