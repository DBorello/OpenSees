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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasHyster1Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasHyster1Material. FedeasHyster1Material wraps the FEDEAS
// 1d material subroutine Hyster_1.

#include <FedeasHyster1Material.h>

FedeasHyster1Material::FedeasHyster1Material(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double pinchX, double pinchY, double damfc1, double damfc2):
// 6 history variables and 12 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasHysteretic1, 6, 12)
{
	data[0]  = mom1p;
	data[1]  = rot1p;
	data[2]  = mom2p;
	data[3]  = rot2p;

	data[4]  = mom1n;
	data[5]  = rot1n;
	data[6]  = mom2n;
	data[7]  = rot2n;

	data[8]  = pinchX;
	data[9]  = pinchY;
	data[10] = damfc1;
	data[11] = damfc2;
}

FedeasHyster1Material::FedeasHyster1Material(void):
FedeasMaterial(0, MAT_TAG_FedeasHysteretic1, 6, 12)
{
	// Does nothing
}

FedeasHyster1Material::~FedeasHyster1Material(void)
{
	// Does nothing
}
