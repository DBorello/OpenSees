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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasBond2Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasBond2Material. FedeasBond2Material wraps the FEDEAS
// 1d material subroutine Bond_2.

#include <FedeasBond2Material.h>

FedeasBond2Material::FedeasBond2Material(int tag,
	double u1p, double q1p, double u2p, double u3p, double q3p,
	double u1n, double q1n, double u2n, double u3n, double q3n,
	double s0, double bb, double alp, double aln, double En0):
// 27 history variables and 15 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasBond2, 27, 15)
{
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
	data[12] = alp;
	data[13] = aln;
	data[14] = En0;
}

FedeasBond2Material::FedeasBond2Material(void):
FedeasMaterial(0, MAT_TAG_FedeasBond2, 27, 15)
{
	// Does nothing
}

FedeasBond2Material::~FedeasBond2Material(void)
{
	// Does nothing
}
