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
// $Date: 2007-10-26 15:55:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/system/MVNcdf.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef MVNcdf_h
#define MVNcdf_h

#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <Vector.h>
#include <Matrix.h>

class MVNcdf : public SystemAnalysis
{

public:
	MVNcdf(ReliabilityDomain *passedReliabilityDomain,
				   TCL_Char *fileName, int analysisType);
	~MVNcdf();

	int		analyze(void);

protected:

private:
	double	MVNcdffunc(const Vector&, const Matrix&, double);
	char fileName[256];
    int analysisType;
	
};

#endif

