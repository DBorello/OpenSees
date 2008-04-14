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
// $Date: 2008-04-14 21:30:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/DamageModel.cpp,v $                                                                        
                                                                        
// Written: Arash Altoontash, Gregory Deierlein
// Created: 08/02
// Revision: AA
//
// Description: This file contains the class implementation for DamageModel

#include <DamageModel.h>
#include <Response.h>
#include <DamageResponse.h>
#include <string.h>

DamageModel::DamageModel(int tag, int clasTag)
:TaggedObject(tag), MovableObject(clasTag)
{

}


DamageModel::~DamageModel()
{
  // does nothing


}

Response*
DamageModel::setResponse(const char **argv, int argc, Information &info)
{
  return 0;
}

int 
DamageModel::getResponse(int responseID, Information &info)
{
  return -1;
}
