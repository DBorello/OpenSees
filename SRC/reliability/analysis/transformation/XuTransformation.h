/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2001-06-13 05:06:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/transformation/XuTransformation.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#ifndef XuTransformation_h
#define XuTransformation_h

#include <Vector.h>
#include <Matrix.h>

class XuTransformation
{

public:
	XuTransformation();
	virtual ~XuTransformation();

	virtual int set_x(Vector x) =0;
	virtual int set_u(Vector u) =0;

	virtual int transform_x_to_u() =0;
	virtual int transform_u_to_x() =0;
	virtual int transform_u_to_x_andComputeJacobian() =0;

	virtual Vector get_x() =0;
	virtual Vector get_u() =0;
	virtual Matrix getJacobian_x_u() =0;
	virtual Matrix getJacobian_u_x() =0;

protected:

private:

};

#endif
