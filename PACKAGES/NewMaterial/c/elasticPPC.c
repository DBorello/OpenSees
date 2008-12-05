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

/*                                                                        
** $Revision: 1.2 $
** $Date: 2008-12-05 19:14:30 $
** $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/c/elasticPPC.c,v $
                                                                        
** Written: fmk 
**
** Description: This file contains the implementation of elasticPP material
*/

#include "OPS_ProceduralAPI.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define DBL_EPSILON 1e-18

#ifdef _USRDLL
#include <windows.h>
#define DllExport _declspec(dllexport)
#elif _MACOSX
#define DllExport __attribute__((visibility("default")))
#else
#define DllExport
#endif

extern "C" DllExport void
elasticPPC (matObj *thisObj, modelState *model, double *strain, double *tang, double *stress, int *isw, int *result) 
{
  if (*isw == ISW_INIT) {

    double dData[2];
    int    iData[1];

    /* get the input data  - tag? E? eyp? */
    int numData = 1;
    OPS_GetIntInput(&numData, iData);
    numData = 2;
    OPS_GetDoubleInput(&numData, dData); 

    /* Allocate the element state */
    thisObj->tag = iData[0];
    thisObj->nParam = 2;  /* E, eyp */
    thisObj->nState = 2;  /* strain, ep */
    OPS_AllocateMaterial(thisObj);

    double E = dData[0];
    double eyp = dData[1];


    thisObj->theParam[0] = E;
    thisObj->theParam[1] = eyp;

  } else if (*isw == ISW_COMMIT) {

    double trialStrain = thisObj->tState[0];
    double f;		// yield function

    double E = thisObj->theParam[0];
    double eyp = thisObj->theParam[1];
    double ep = thisObj->cState[1];

    double fyp = E*eyp;
    double fyn = -fyp;

    // compute trial stress
    double sigtrial = E * ( trialStrain - ep );

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    double fYieldSurface = - E * DBL_EPSILON;
    if ( f > fYieldSurface ) {
      // plastic
      if ( sigtrial > 0.0 ) {
	ep += f / E;
      } else {
	ep -= f / E;
      }
    }

    thisObj->cState[0] = trialStrain;    
    thisObj->cState[1] = ep;

  } else if (*isw == ISW_REVERT) {
    
    for (int i=0; i<3; i++)
      thisObj->tState[i] = thisObj->cState[i];

  } else if (*isw == ISW_REVERT_TO_START) {

    for (int i=0; i<2; i++) {
      thisObj->cState[i] = 0.0;
      thisObj->tState[i] = 0.0;
    }
    
  } else if (*isw == ISW_FORM_TANG_AND_RESID) {

    double trialStrain = *strain;
    double f;		// yield function
    double trialStress, trialTangent;

    double E = thisObj->theParam[0];
    double eyp = thisObj->theParam[1];
    double ep = thisObj->cState[1];

    double fyp = E*eyp;
    double fyn = -fyp;
    
    // compute trial stress
    double sigtrial = E * ( trialStrain - ep );

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    double fYieldSurface = - E * DBL_EPSILON;
    if ( f <= fYieldSurface ) {

      // elastic
      trialStress = sigtrial;
      trialTangent = E;

    } else {

      // plastic
      if ( sigtrial > 0.0 ) {
	trialStress = fyp;
      } else {
	trialStress = fyn;
      }

      trialTangent = 0.0;
    }

    thisObj->tState[0] = trialStrain;
    *stress = trialStress;
    *tang = trialTangent;
  }

  *result = 0;
}
