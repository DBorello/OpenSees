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
                                                                        
// $Revision: 1.4 $
// $Date: 2001-08-20 00:37:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/OpenSeesGFunEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//			haukaas 06/22/01 (analysis commands batch read from file)
//			haukaas 08/19/01 (modifications for Release 1.2 of OpenSees)
//

#ifndef OpenSeesGFunEvaluator_h
#define OpenSeesGFunEvaluator_h

#include <GFunEvaluator.h>
#include <Vector.h>
#include <ReliabilityDomain.h>

#include <tcl.h>


class OpenSeesGFunEvaluator : public GFunEvaluator
{

public:
	OpenSeesGFunEvaluator(Tcl_Interp *passedTclInterp,
						ReliabilityDomain *passedReliabilityDomain,
						char *fileName);
	~OpenSeesGFunEvaluator();

	int		evaluate_g(Vector passed_x);
	double	get_g();

protected:

private:
	double g;
	Tcl_Interp *theTclInterp;
	ReliabilityDomain *theReliabilityDomain;
	char fileName[256];

};

#endif
