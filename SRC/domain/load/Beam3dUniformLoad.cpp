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
// $Date: 2002-06-07 17:39:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dUniformLoad.cpp,v $
                                                                        

// Written: fmk 
//
// Purpose: This file contains the class implementation of Beam3dUniformLoad.

#include <Beam3dUniformLoad.h>
#include <Vector.h>

Vector Beam3dUniformLoad::data(3);

Beam3dUniformLoad::Beam3dUniformLoad(int tag, double wY, double wZ, double wX,
				     const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam3dUniformLoad, theElementTags),
   wy(wY), wz(wZ), wx(wX)
{

}

Beam3dUniformLoad::Beam3dUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam3dUniformLoad),
   wy(0.0), wz(0.0), wx(0.0)
{

}

Beam3dUniformLoad::~Beam3dUniformLoad()
{

}

const Vector &
Beam3dUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dUniformLoad;
  data(0) = wy;
  data(1) = wz;
  data(2) = wx;
  return data;
}


int 
Beam3dUniformLoad::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
Beam3dUniformLoad::recvSelf(int commitTag, Channel &theChannel,
			    FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
Beam3dUniformLoad::Print(ostream &s, int flag)
{
  s << "Beam3dUniformLoad - Reference load" << endl;
  s << "  Transverse (y): " << wy << endl;
  s << "  Transverse (z): " << wz << endl;
  s << "  Axial (x):      " << wx << endl;
  s << "  Elements acted on: " << this->getElementTags();
}
