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
// $Date: 2001-08-17 16:28:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasConcr3Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasConcr3Material. FedeasConcr3Material wraps the FEDEAS
// 1d material subroutine Concr_3.

#include <FedeasConcr3Material.h>

FedeasConcr3Material::FedeasConcr3Material(int tag,
					 double fc, double ec, double fu, double eu,
					 double rat, double ft, double epst0,
					 double ft0, double beta, double epstu):
// 2 history variables and 10 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasConcrete3, 2, 10)
{
	data[0]  = fc;
	data[1]  = ec;
	data[2]  = fu;
	data[3]  = eu;
	data[4]  = rat;
	data[5]  = ft;
	data[6]  = epst0;
	data[7]  = ft0;
	data[8]  = beta;
	data[9]  = epstu;
}

FedeasConcr3Material::FedeasConcr3Material(void):
FedeasMaterial(0, MAT_TAG_FedeasConcrete3, 2, 10)
{
	// Does nothing
}

FedeasConcr3Material::~FedeasConcr3Material(void)
{
	// Does nothing
}
