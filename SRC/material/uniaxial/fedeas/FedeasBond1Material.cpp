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
// $Date: 2001-08-17 16:28:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasBond1Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasBond1Material. FedeasBond1Material wraps the FEDEAS
// 1d material subroutine Bond_1.

#include <FedeasBond1Material.h>

FedeasBond1Material::FedeasBond1Material(int tag,
	double u1p, double q1p, double u2p, double u3p, double q3p,
	double u1n, double q1n, double u2n, double u3n, double q3n,
	double s0, double bb):
// 6 history variables and 12 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasBond1, 6, 12)
{
	// Fill in material parameters
	data[0]  = u1p;
	data[1]  = q1p;
	data[2]  = u2p;
	data[3]  = u3p;
	data[4]  = q3p;
	data[5]  = u1n;
	data[6]  = q1n;
	data[7]  = u2n;
	data[8]  = u3n;
	data[9]  = q3n;
	data[10] = s0;
	data[11] = bb;
}

FedeasBond1Material::FedeasBond1Material(void):
FedeasMaterial(0, MAT_TAG_FedeasBond1, 6, 12)
{
	// Does nothing
}

FedeasBond1Material::~FedeasBond1Material(void)
{
	// Does nothing
}
