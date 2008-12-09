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
** $Date: 2008-12-09 22:36:45 $
** $Source: /usr/local/cvs/OpenSees/PACKAGES/NewElement/c/trussC.c,v $
                                                                        
** Written: fmk 
**
** Description: This file contains the implementation of truss element.
*/

#include "elementAPI.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef _USRDLL
#include <windows.h>
#define DllExport _declspec(dllexport)
#elif _MACOSX
#define DllExport __attribute__((visibility("default")))
#else
#define DllExport
#endif

extern "C" DllExport void
trussC (eleObj *thisObj, modelState *model, double *tang, double *resid, int *isw, int *errFlag) 
{
  double matStrain[1];
  double matTang[1];
  double matStress[1];

  if (*isw == ISW_INIT) {

    double dData[1];
    int    iData[4];

    /* get the input data  - tag? nd1? nd2? A? matTag? */
    int numData = 3;
    OPS_GetIntInput(&numData, iData);
    numData = 1;
    OPS_GetDoubleInput(&numData, dData);
    OPS_GetIntInput(&numData, &iData[3]);

    /* Allocate the element state */
    thisObj->tag = iData[0];
    int n1 = iData[1];
    int n2 = iData[2];
    int matTag = iData[3];

    thisObj->nNode = 2;
    thisObj->nParam = 4;
    thisObj->nDOF = 4;
    thisObj->nState = 0;
    thisObj->nMat = 1;
    iData[0] = matTag;

    int *matData = iData;

    int matType = OPS_UNIAXIAL_MATERIAL_TYPE;
    OPS_AllocateElement(thisObj, matData, &matType);

    /* fill in the element state */
    thisObj->node[0] = n1;
    thisObj->node[1] = n2;
    double nd1Crd[2];
    double nd2Crd[2];
    int numCrd = 2;

    thisObj->param[0] = dData[0];    

    OPS_GetNodeCrd(&n1, &numCrd, nd1Crd);
    OPS_GetNodeCrd(&n2, &numCrd, nd2Crd);

    double dx = nd2Crd[0]-nd1Crd[0];
    double dy = nd2Crd[1]-nd1Crd[1];	
    double L = sqrt(dx*dx + dy*dy);

    thisObj->param[1] = L;    
    if (L == 0.0) {
      OPS_Error("Warning - truss element has zero length\n", 1);
      return;
    }

    double cs = dx/L;
    double sn = dy/L;
    
    thisObj->param[2] = cs;
    thisObj->param[3] = sn;

    /* ******************************
       placed in AllocateElement
       matObj *theMat = OPS_GetMaterial(&(iData[3]));
       if (theMat == 0) {
       //      OPS_Error("Warning - truss element could not find material\n",1);
       return;
       }
       thisObj->mats[0] = theMat;
       ********************************************** */

    *errFlag = 0;

  } else if (*isw == ISW_COMMIT) {

      matObj *theMat = thisObj->mats[0];
      theMat->matFunctPtr(theMat, model, matStrain, matTang, matStress, isw, errFlag); 
      
  } else if (*isw == ISW_REVERT) {
    
    matObj *theMat = thisObj->mats[0];
    theMat->matFunctPtr(theMat, model, matStrain, matTang, matStress, isw, errFlag); 

  } else if (*isw == ISW_REVERT_TO_START) {

    matObj *theMat = thisObj->mats[0];
    theMat->matFunctPtr(theMat, model, matStrain, matTang, matStress, isw, errFlag); 
    
  } else if (*isw == ISW_FORM_TANG_AND_RESID) {

    double L = thisObj->param[1];
    if (L == 0.0) {
      //      OPS_Error("Warning - truss element has zero length\n", 1);
      return;
    }

    double d1[2];
    double d2[2];

    int nd1 = thisObj->node[0];
    int nd2 = thisObj->node[1];

    int numDOF = 2;
    OPS_GetNodeDisp(&nd1, &numDOF, d1);
    OPS_GetNodeDisp(&nd2, &numDOF, d2);    

    double A = thisObj->param[0];
    double cs = thisObj->param[2];
    double sn = thisObj->param[3];

    double tran[4];
    tran[0] = -cs;
    tran[1] = -sn;
    tran[2] = cs;
    tran[3] = sn;

    double dLength = 0.0;
    for (int i=0; i<2; i++){
      dLength -= (d2[i]-d1[i]) * tran[i];
    }

    matStrain[0] = dLength/L;

    matObj *theMat = thisObj->mats[0];

    theMat->matFunctPtr(theMat, model, matStrain, matTang, matStress, isw, errFlag); 

    /* ******************* instead of call material funtion
    int matNum = 0;
    *errFlag = OPS_InvokeMaterial(thisObj, &matNum, model, matStrain, matStress, matTang, isw); 
    fprintf(stderr,"strain, tang, stress: %e %e %e\n",matStrain[0], matTang[0], matStress[0]);


    *errFlag = OPS_InvokeMaterialDirectly(theMat, model, matStrain, matStress, matTang, isw); 
    ******************************************************** */

    if (*errFlag == 0) {

      double force = A*matStress[0];

      for (int i=0; i<4; i++)
	resid[i] = tran[i]*force;
      
      double k = A*matTang[0]/L;

      // tang(j,i)
      for (int i = 0; i<4; i++) 
	for (int j=0; j < 4; j++)
	  tang[i+j*4] = k * tran[i]*tran[j];

    }

  } else if (*isw == ISW_FORM_MASS) {

    double L = thisObj->param[1];
    if (L == 0.0) {
      //      OPS_Error("Warning - truss element has zero length\n", 1);
      return;
    }
    double A = thisObj->param[0];

    double rho = thisObj->param[4];
    for (int i=0; i<16; i++) 
      tang[i] = 0.0;

    if (rho != 0.0) {
      double massV  = rho * A * L/2;
      tang[0] = massV;
      tang[1+1*4] = massV;
      tang[2+2*4] = massV;
      tang[3+3*4] = massV;
    }
    
    *errFlag = 0;
  }
}



