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
// $Date: 2001-08-13 21:27:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasHardeningMaterial.h,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasHardeningMaterial. FedeasHardeningMaterial wraps the FEDEAS
// 1d material subroutine Hard_1.

#ifndef FedeasHardeningMaterial_h
#define FedeasHardeningMaterial_h

#include <FedeasMaterial.h>

class FedeasHardeningMaterial : public FedeasMaterial
{
  public:
    FedeasHardeningMaterial(int tag,
		double E, double sigmaY, double Hiso, double Hkin);
	FedeasHardeningMaterial(void);
    virtual ~FedeasHardeningMaterial();

  protected:

  private:

};


#endif

