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
// $Date: 2002-06-07 17:39:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dPointLoad.cpp,v $
                                                                        
// Written: fmk 

// Purpose: This file contains the class implementation Beam3dPointLoad.

#include <Beam3dPointLoad.h>
#include <Vector.h>

Vector Beam3dPointLoad::data(4);

Beam3dPointLoad::Beam3dPointLoad(int tag, double py, double pz, double dist,
				 const ID &theElementTags, double px)
  :ElementalLoad(tag, LOAD_TAG_Beam3dPointLoad, theElementTags),
   Py(py), Pz(pz), Px(px), x(dist)
{

}

Beam3dPointLoad::Beam3dPointLoad()
  :ElementalLoad(LOAD_TAG_Beam3dPointLoad),
   Py(0.0), Pz(0.0), Px(0.0), x(0.0)
{

}

Beam3dPointLoad::~Beam3dPointLoad()
{

}

const Vector &
Beam3dPointLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dPointLoad;
  data(0) = Py;
  data(1) = Pz;
  data(2) = Px;
  data(3) = x;
  return data;
}

int 
Beam3dPointLoad::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
Beam3dPointLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
Beam3dPointLoad::Print(ostream &s, int flag)
{
  s << "Beam3dPointLoad - Reference load" << endl;
  s << "  Transverse (y): " << Py << endl;
  s << "  Transverse (z): " << Pz << endl;
  s << "  Axial (x):      " << Px << endl;
  s << "  Elements acted on: " << this->getElementTags();
}
