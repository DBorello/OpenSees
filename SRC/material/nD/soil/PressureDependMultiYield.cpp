// $Revision: 1.20 $
// $Date: 2002-05-16 00:07:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureDependMultiYield.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000

//
// PressureDependMultiYield.cpp
// -------------------
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include <PressureDependMultiYield.h>
#include <Information.h>

int PressureDependMultiYield::loadStage = 0;
double PressureDependMultiYield::pAtm = 101.;
Matrix PressureDependMultiYield::theTangent(6,6);
T2Vector PressureDependMultiYield::trialStrain;
T2Vector PressureDependMultiYield::subStrainRate;
Vector PressureDependMultiYield::workV6(6);
T2Vector PressureDependMultiYield::workT2V;

const	double pi = 3.14159265358979;

PressureDependMultiYield::PressureDependMultiYield (int tag, int nd, 
								double r, double refShearModul,
						    double refBulkModul, double frictionAng,
						    double peakShearStra, double refPress, 
						    double pressDependCoe,
						    double phaseTransformAng, 
						    double contractionParam1,
						    double dilationParam1,
						    double dilationParam2,
						    double liquefactionParam1,
						    double liquefactionParam2,
						    double liquefactionParam4,
						    int numberOfYieldSurf, 
						    double ei,
								double volLim1, double volLim2, double volLim3,
								double atm, double cohesi)
 : NDMaterial(tag,ND_TAG_PressureDependMultiYield), currentStress(),
   trialStress(), currentStrain(), strainRate(),
   reversalStress(), PPZPivot(), PPZCenter(), 
   lockStress(), reversalStressCommitted(), 
   PPZPivotCommitted(), PPZCenterCommitted(),
   lockStressCommitted()
{
  if (nd !=2 && nd !=3) {
    cerr << "FATAL:PressureDependMultiYield:: dimension error" << endl;
    cerr << "Dimension has to be 2 or 3, you give nd= " << nd << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (refShearModul <= 0) {
    cerr << "FATAL:PressureDependMultiYield:: refShearModulus <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (refBulkModul <= 0) {
    cerr << "FATAL:PressureDependMultiYield:: refBulkModulus <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (frictionAng <= 0.) {
    cerr << "FATAL:PressureDependMultiYield:: frictionAngle <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (frictionAng >= 90.) {
    cerr << "FATAL:PressureDependMultiYield:: frictionAngle >= 90" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (phaseTransformAng <= 0.) {
    cerr << "FATAL:PressureDependMultiYield:: phaseTransformAng <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (phaseTransformAng > frictionAng) {
    cerr << "WARNING:PressureDependMultiYield:: phaseTransformAng > frictionAng" << endl;
    cerr << "Will set phaseTransformAng = frictionAng." <<endl;
    phaseTransformAng = frictionAng;
  }
  if (cohesi < 0) {
    cerr << "WARNING:PressureDependMultiYield:: cohesion < 0" << endl;
    cerr << "Will reset cohesion to zero." << endl;
    cohesi = 0.;
  } 
  if (peakShearStra <= 0) {
    cerr << "FATAL:PressureDependMultiYield:: peakShearStra <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (refPress <= 0) {
    cerr << "FATAL:PressureDependMultiYield:: refPress <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (pressDependCoe < 0) {
    cerr << "WARNING:PressureDependMultiYield:: pressDependCoe < 0" << endl;
    cerr << "Will reset pressDependCoe to zero." << endl;
    pressDependCoe = 0.;
  }
  if (numberOfYieldSurf <= 0) {
    cerr << "WARNING:PressureDependMultiYield:: numberOfSurfaces <= 0" << endl;
    cerr << "Will use 10 yield surfaces." << endl;
    numberOfYieldSurf = 10;
  }
  if (numberOfYieldSurf > 40) {
    cerr << "WARNING:PressureDependMultiYield::PressureDependMultiYield: numberOfSurfaces > 40" << endl;
    cerr << "Will use 40 yield surfaces." << endl;
    numberOfYieldSurf = 40;
  }
  if (volLim1 < 0) {
    cerr << "WARNING:PressureDependMultiYield:: volLim1 < 0" << endl;
    cerr << "Will reset volLimit to 0.8" << endl;
    volLim1 = 0.8;
  }
  if (r < 0) {
    cerr << "FATAL:PressureDependMultiYield:: rho <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (ei < 0) {
    cerr << "FATAL:PressureDependMultiYield:: e <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }

  ndm = nd;
  loadStage = 0;   //default
  refShearModulus = refShearModul;
  refBulkModulus = refBulkModul;
  frictionAngle = frictionAng;
  peakShearStrain = peakShearStra;
  refPressure = -refPress;  //compression is negative
  cohesion = cohesi;
  pressDependCoeff = pressDependCoe;
  numOfSurfaces = numberOfYieldSurf;
  phaseTransfAngle = phaseTransformAng;
  contractParam1 = contractionParam1;
  dilateParam1 = dilationParam1;
  dilateParam2 = dilationParam2;
  volLimit1 = volLim1;
  volLimit2 = volLim2;
  volLimit3 = volLim3;
  liquefyParam1 = liquefactionParam1;
  liquefyParam2 = liquefactionParam2;
  liquefyParam4 = liquefactionParam4;
  rho = r;
	einit = ei;
  pAtm = atm;

  e2p = committedActiveSurf = activeSurfaceNum = 0; 
  onPPZCommitted = onPPZ = -1 ; 
  PPZSizeCommitted = PPZSize = 0.;
  pressureDCommitted = pressureD = modulusFactor = 0.;
  cumuDilateStrainOctaCommitted    = cumuDilateStrainOcta = 0.;
  maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta = 0.;
  cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta = 0.;
  prePPZStrainOctaCommitted = prePPZStrainOcta = 0.;
  oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta = 0.;

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1]; 

  setUpSurfaces();  // residualPress and stressRatioPT are calculated inside.
}
   

PressureDependMultiYield::PressureDependMultiYield () 
 : NDMaterial(0,ND_TAG_PressureDependMultiYield), 
   currentStress(), trialStress(), currentStrain(), 
  strainRate(), reversalStress(), PPZPivot(), 
  PPZCenter(), lockStress(), reversalStressCommitted(), 
  PPZPivotCommitted(), PPZCenterCommitted(),
  lockStressCommitted()
{
  ndm = 3;
  refShearModulus = refBulkModulus = frictionAngle = 0.;
  peakShearStrain = refPressure = cohesion = pressDependCoeff = 0.;
  numOfSurfaces = 1;
  phaseTransfAngle = contractParam1 = 0.;
  dilateParam1 = dilateParam2 = volLimit1 = volLimit2 = volLimit3 = 0.;
  liquefyParam1 = liquefyParam2 = liquefyParam4 = modulusFactor = 0.;
  rho = einit = 0.;

  theSurfaces = new MultiYieldSurface[1];
  committedSurfaces = new MultiYieldSurface[1];
  activeSurfaceNum = committedActiveSurf = 0; 
}


PressureDependMultiYield::PressureDependMultiYield (const PressureDependMultiYield & a)
 : NDMaterial(a.getTag(),ND_TAG_PressureDependMultiYield), 
   currentStress(a.currentStress), trialStress(a.trialStress), 
  currentStrain(a.currentStrain), strainRate(a.strainRate), 
   reversalStress(a.reversalStress), PPZPivot(a.PPZPivot), PPZCenter(a.PPZCenter), 
  lockStress(a.lockStress), reversalStressCommitted(a.reversalStressCommitted), 
  PPZPivotCommitted(a.PPZPivotCommitted), 
  PPZCenterCommitted(a.PPZCenterCommitted),
  lockStressCommitted(a.lockStressCommitted)
{
  ndm = a.ndm;
  refShearModulus = a.refShearModulus;
  refBulkModulus = a.refBulkModulus;
  frictionAngle = a.frictionAngle;
  peakShearStrain = a.peakShearStrain;
  refPressure = a.refPressure;
  cohesion = a.cohesion;
  pressDependCoeff = a.pressDependCoeff;
  numOfSurfaces = a.numOfSurfaces;
  phaseTransfAngle = a.phaseTransfAngle;
  contractParam1 = a.contractParam1;
  dilateParam1 = a.dilateParam1;
  dilateParam2 = a.dilateParam2;
  volLimit1 = a.volLimit1;
  volLimit2 = a.volLimit2;
  volLimit3 = a.volLimit3;
  liquefyParam1 = a.liquefyParam1;
  liquefyParam2 = a.liquefyParam2;
  liquefyParam4 = a.liquefyParam4;
  rho = a.rho;
	einit = a.einit;

  e2p = a.e2p;
  residualPress = a.residualPress;
  stressRatioPT = a.stressRatioPT;
  strainPTOcta = a.strainPTOcta;
  modulusFactor = a.modulusFactor;
  activeSurfaceNum = a.activeSurfaceNum;
  committedActiveSurf = a.committedActiveSurf;
  pressureDCommitted     = a.pressureDCommitted;
  onPPZCommitted = a.onPPZCommitted; 
  PPZSizeCommitted      = a.PPZSizeCommitted;
  cumuDilateStrainOctaCommitted    = a.cumuDilateStrainOctaCommitted;
  maxCumuDilateStrainOctaCommitted = a.maxCumuDilateStrainOctaCommitted;
  cumuTranslateStrainOctaCommitted = a.cumuTranslateStrainOctaCommitted;
  prePPZStrainOctaCommitted        = a.prePPZStrainOctaCommitted;
  oppoPrePPZStrainOctaCommitted    = a.oppoPrePPZStrainOctaCommitted;
  pressureD     = a.pressureD;
  onPPZ = a.onPPZ; 
  PPZSize      = a.PPZSize;
  cumuDilateStrainOcta    = a.cumuDilateStrainOcta;
  maxCumuDilateStrainOcta = a.maxCumuDilateStrainOcta;
  cumuTranslateStrainOcta = a.cumuTranslateStrainOcta;
  prePPZStrainOcta        = a.prePPZStrainOcta;
  oppoPrePPZStrainOcta    = a.oppoPrePPZStrainOcta;

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1];  //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];  
  for(int i=1; i<=numOfSurfaces; i++) {
    committedSurfaces[i] = a.committedSurfaces[i];  
    theSurfaces[i] = a.theSurfaces[i];  
  }
}


PressureDependMultiYield::~PressureDependMultiYield ()
{
  if (theSurfaces != 0) delete [] theSurfaces;
  if (committedSurfaces != 0) delete [] committedSurfaces;
}


void PressureDependMultiYield::elast2Plast(void)
{
  if (loadStage == 0 || e2p == 1) return;
  e2p = 1;

  if (currentStress.volume() > 0.) {
    //cerr << "WARNING:PressureDependMultiYield::elast2Plast(): material in tension." << endl;
    currentStress.setData(currentStress.deviator(),0);
  }

  // Active surface is 0, return
  if (currentStress.deviatorLength() == 0.) return;

  // Find active surface
  while (yieldFunc(currentStress, committedSurfaces, ++committedActiveSurf) > 0) {
    if (committedActiveSurf == numOfSurfaces) {
      //cerr <<"WARNING:PressureDependMultiYield::elast2Plast(): stress out of failure surface"<<endl;
      deviatorScaling(currentStress, committedSurfaces, numOfSurfaces);
      initSurfaceUpdate();
      return;
    }
  } 

  committedActiveSurf--;
  initSurfaceUpdate();
}

 
int PressureDependMultiYield::setTrialStrain (const Vector &strain)
{
  if (ndm==3 && strain.Size()==6) 
    workV6 = strain;
  else if (ndm==2 && strain.Size()==3) {
    workV6[0] = strain[0];
    workV6[1] = strain[1];
    workV6[2] = 0.0;
    workV6[3] = strain[2];
    workV6[4] = 0.0;
    workV6[5] = 0.0;
  }
  else {
    cerr << "Fatal:PressureDependMultiYield:: Material dimension is: " << ndm << endl;
    cerr << "But strain vector size is: " << strain.Size() << endl;
    g3ErrorHandler->fatal(" ");
  }

  //strainRate.setData(workV6-currentStrain.t2Vector(1),1);
  workV6 -= currentStrain.t2Vector(1);
  strainRate.setData(workV6, 1);

  return 0;
}


int PressureDependMultiYield::setTrialStrain (const Vector &strain, const Vector &rate)
{
  return setTrialStrain (strain);
}


int PressureDependMultiYield::setTrialStrainIncr (const Vector &strain)
{
  if (ndm==3 && strain.Size()==6) 
    workV6 = strain;
  else if (ndm==2 && strain.Size()==3) {
    workV6[0] = strain[0];
    workV6[1] = strain[1];
    workV6[2] = 0.0;
    workV6[3] = strain[2];
    workV6[4] = 0.0;
    workV6[5] = 0.0;
  }
  else {
    cerr << "Fatal:PressureDependMultiYield:: Material dimension is: " << ndm << endl;
    cerr << "But strain vector size is: " << strain.Size() << endl;
    g3ErrorHandler->fatal(" ");
  }

  strainRate.setData(workV6,1);
  return 0;
}


int PressureDependMultiYield::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return setTrialStrainIncr(strain);
}


const Matrix & PressureDependMultiYield::getTangent (void)
{
  if (loadStage != 0 && e2p == 0) elast2Plast();
  
  if (loadStage==0) {  //linear elastic
    for (int i=0;i<6;i++) 
      for (int j=0;j<6;j++) {
	theTangent(i,j) = 0.;
        if (i==j) theTangent(i,j) += refShearModulus;
        if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
      }
  }
  else {
    double coeff1, coeff2;
    double factor = getModulusFactor(currentStress);
    double shearModulus = factor*refShearModulus;
    double bulkModulus = factor*refBulkModulus;		
	
    if (loadStage!=0 && committedActiveSurf > 0) {
      getSurfaceNormal(currentStress, workT2V);
      workV6 = workT2V.deviator();
      double volume = workT2V.volume();
      double Ho = 9.*bulkModulus*volume*volume + 2.*shearModulus*(workV6 && workV6);
      double plastModul = committedSurfaces[committedActiveSurf].modulus();
      coeff1 = 9.*bulkModulus*bulkModulus*volume*volume/(Ho+plastModul);
      coeff2 = 4.*shearModulus*shearModulus/(Ho+plastModul);
    }
    else {
      coeff1 = coeff2 = 0.;
      workV6.Zero();
    }

    for (int i=0;i<6;i++) 
      for (int j=0;j<6;j++) {
	theTangent(i,j) = - coeff2*workV6[i]*workV6[j];
        if (i==j) theTangent(i,j) += shearModulus;
        if (i<3 && j<3 && i==j) theTangent(i,j) += shearModulus;
	if (i<3 && j<3) theTangent(i,j) += (bulkModulus - 2.*shearModulus/3. - coeff1);
      }
  }


  if (ndm==3) 
    return theTangent;
  else {
    static Matrix workM(3,3);
    workM(0,0) = theTangent(0,0);
    workM(0,1) = theTangent(0,1);
    workM(0,2) = 0.; 
    //workM(0,2) = theTangent(0,3);
    workM(1,0) = theTangent(1,0);
    workM(1,1) = theTangent(1,1);
    workM(1,2) = 0.; 
    workM(2,0) = 0.; 
    workM(2,1) = 0.; 
    //workM(1,2) = theTangent(1,3);
    //workM(2,0) = theTangent(3,0);
    //workM(2,1) = theTangent(3,1);
    workM(2,2) = theTangent(3,3);
    return workM;
  }
}


const Vector & PressureDependMultiYield::getStress (void)
{
  int i;
  if (loadStage != 0 && e2p == 0) 
    elast2Plast();

  if (loadStage==0) {  //linear elastic
    trialStrain.setData(currentStrain.t2Vector() + strainRate.t2Vector());
    getTangent();
    //workV6 = theTangent * trialStrain.t2Vector(1);
	  workV6.addMatrixVector(0.0, theTangent, trialStrain.t2Vector(1), 1.0);

    trialStress.setData(workV6);
  }
  else {
    for (i=1; i<=numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
    activeSurfaceNum = committedActiveSurf;
    pressureD = pressureDCommitted;
    reversalStress = reversalStressCommitted;
    onPPZ = onPPZCommitted;
    PPZSize = PPZSizeCommitted;
    cumuDilateStrainOcta = cumuDilateStrainOctaCommitted;
    maxCumuDilateStrainOcta = maxCumuDilateStrainOctaCommitted;
    cumuTranslateStrainOcta = cumuTranslateStrainOctaCommitted;
    prePPZStrainOcta = prePPZStrainOctaCommitted;
    oppoPrePPZStrainOcta = oppoPrePPZStrainOctaCommitted;
    PPZPivot = PPZPivotCommitted;
    PPZCenter = PPZCenterCommitted;
    lockStress = lockStressCommitted;

    subStrainRate = strainRate;
    setTrialStress(currentStress);
    if (activeSurfaceNum>0 && isLoadReversal()) {
      updateInnerSurface();
      activeSurfaceNum = 0;
    }

    if (activeSurfaceNum==0 && !isCrossingNextSurface()) {
 	    workV6 = currentStrain.t2Vector();
	    workV6.addVector(1.0, strainRate.t2Vector(), 1.0);
	    trialStrain.setData(workV6);
		}
		else {
      int numSubIncre = setSubStrainRate();

      for (i=0; i<numSubIncre; i++) {
//      trialStrain.setData(currentStrain.t2Vector() 
//			     + subStrainRate.t2Vector()*(i+1));
 	      workV6 = currentStrain.t2Vector();
	      workV6.addVector(1.0, subStrainRate.t2Vector(), (i+1));
	      trialStrain.setData(workV6);

        if (i==0)  setTrialStress(currentStress);
        else setTrialStress(trialStress); 
        if (activeSurfaceNum>0 && isLoadReversal()) {
          updateInnerSurface();
          activeSurfaceNum = 0;
				}
        if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
        if (activeSurfaceNum==0) activeSurfaceNum++;
        int lock = stressCorrection(0);
        if(lock==0) updateActiveSurface();
			}
    }
  }

  if (ndm==3)
    return trialStress.t2Vector();
  else {
		static Vector workV(3);
    workV[0] = trialStress.t2Vector()[0];
    workV[1] = trialStress.t2Vector()[1];
    workV[2] = trialStress.t2Vector()[3];
    return workV;
  }
}


const Vector & PressureDependMultiYield::getStrain (void)
{
  return getCommittedStrain();
}


int PressureDependMultiYield::commitState (void)
{
  currentStress = trialStress;
  //currentStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
  workV6 = currentStrain.t2Vector();
  workV6 += strainRate.t2Vector();
  currentStrain.setData(workV6);

  workV6.Zero();
  strainRate.setData(workV6);

  if (loadStage) {
    committedActiveSurf = activeSurfaceNum;
    for (int i=1; i<=numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
    pressureDCommitted = pressureD;
    reversalStressCommitted = reversalStress;
    onPPZCommitted = onPPZ;
    PPZSizeCommitted = PPZSize;
    cumuDilateStrainOctaCommitted = cumuDilateStrainOcta;
    maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta;
    cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta;
    prePPZStrainOctaCommitted = prePPZStrainOcta;
    oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta;
    PPZPivotCommitted = PPZPivot;
    PPZCenterCommitted = PPZCenter;
    lockStressCommitted = lockStress;
  }
 
  return 0;
}


int PressureDependMultiYield::revertToLastCommit (void)
{
  return 0;
}


NDMaterial * PressureDependMultiYield::getCopy (void)
{
  PressureDependMultiYield * copy = new PressureDependMultiYield(*this);
  return copy;
}


NDMaterial * PressureDependMultiYield::getCopy (const char *code)
{
  if (strcmp(code,"PressureDependMultiYield") == 0 || strcmp(code,"PlaneStrain") == 0
      || strcmp(code,"ThreeDimensional") == 0) {
    PressureDependMultiYield * copy = new PressureDependMultiYield(*this);
    return copy;
  }

  return 0;
}


const char * PressureDependMultiYield::getType (void) const
{
  return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int PressureDependMultiYield::getOrder (void) const
{
  return (ndm == 2) ? 3 : 6;
}


int PressureDependMultiYield::updateParameter(int responseID, Information &info)
{
  loadStage = responseID;
  return 0;
}


int PressureDependMultiYield::sendSelf(int commitTag, Channel &theChannel)
{
	int i, res = 0;

	static Vector data(393);
	data(0) = this->getTag();
	data(1) = ndm;
	data(2) = loadStage;
  data(3) = rho;
  data(4) = einit;
  data(5) = refShearModulus;
	data(6) = refBulkModulus;
	data(7) = frictionAngle;
  data(8) = peakShearStrain;
  data(9) = refPressure;
	data(10) = cohesion;
	data(11) = pressDependCoeff;
	data(12) = numOfSurfaces;
	data(13) = phaseTransfAngle;
  data(14) = contractParam1;
  data(15) = dilateParam1;
  data(16) = dilateParam2;
  data(17) = volLimit1;
  data(18) = volLimit2;
  data(19) = volLimit3;
  data(20) = pAtm;
  data(21) = liquefyParam1;
  data(22) = liquefyParam2;
  data(23) = liquefyParam4;
  data(24) = residualPress;
  data(25) = stressRatioPT;   
  data(26) = e2p; 
  data(27) = committedActiveSurf;
  data(28) = strainPTOcta;
  data(29) = pressureDCommitted;
  data(30) = onPPZCommitted;
  data(31) = PPZSizeCommitted;
  data(32) = cumuDilateStrainOctaCommitted;
  data(33) = maxCumuDilateStrainOctaCommitted;
  data(34) = cumuTranslateStrainOctaCommitted;
  data(35) = prePPZStrainOctaCommitted;
  data(36) = oppoPrePPZStrainOctaCommitted;

  workV6 = currentStress.t2Vector();
	for(i = 0; i < 6; i++) data(i+37) = workV6[i];

  workV6 = currentStrain.t2Vector();
	for(i = 0; i < 6; i++) data(i+43) = workV6[i];
	  
  workV6 = PPZPivotCommitted.t2Vector();
	for(i = 0; i < 6; i++) data(i+49) = workV6[i];

  workV6 = PPZCenterCommitted.t2Vector();
	for(i = 0; i < 6; i++) data(i+55) = workV6[i];

  workV6 = lockStressCommitted.t2Vector();
	for(i = 0; i < 6; i++) data(i+61) = workV6[i];

  workV6 = reversalStressCommitted.t2Vector();
	for(i = 0; i < 6; i++) data(i+67) = workV6[i];

	for(i = 0; i < numOfSurfaces; i++) {
		int k = 73 + i*8;
		data(k) = committedSurfaces[i+1].size();
		data(k+1) = committedSurfaces[i+1].modulus();
		workV6 = committedSurfaces[i+1].center();
    data(k+2) = workV6(0);
    data(k+3) = workV6(1);
    data(k+4) = workV6(2);
    data(k+5) = workV6(3);
    data(k+6) = workV6(4);
    data(k+7) = workV6(5);
	}

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"PressureDependMultiYield::sendSelf");
		return res;
	}

	return res;
}


int PressureDependMultiYield::recvSelf(int commitTag, Channel &theChannel, 
				       FEM_ObjectBroker &theBroker)    
{
	int i, res = 0;

	static Vector data(393);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"PressureDependMultiYield::recvSelf");
		return res;
	}
    
	this->setTag((int)data(0));
	ndm = data(1);
	loadStage = data(2);
  rho = data(3);
	einit = data(4);
  refShearModulus = data(5);
	refBulkModulus = data(6);
	frictionAngle = data(7);
  peakShearStrain = data(8);
  refPressure = data(9);
	cohesion = data(10);
	pressDependCoeff = data(11);
	numOfSurfaces = data(12);
	phaseTransfAngle = data(13);
  contractParam1 = data(14);
  dilateParam1 = data(15);
  dilateParam2 = data(16);
  volLimit1 = data(17);
  volLimit2 = data(18);
  volLimit3 = data(19);
  pAtm = data(20);
  liquefyParam1 = data(21);
  liquefyParam2 = data(22);
  liquefyParam4 = data(23);
  residualPress = data(24);
  stressRatioPT = data(25);   
  e2p = data(26); 
  committedActiveSurf = data(27);
  strainPTOcta = data(28);
  pressureDCommitted = data(29);
  onPPZCommitted = data(30);
  PPZSizeCommitted = data(31);
  cumuDilateStrainOctaCommitted = data(32);
  maxCumuDilateStrainOctaCommitted = data(33);
  cumuTranslateStrainOctaCommitted = data(34);
  prePPZStrainOctaCommitted = data(35);
  oppoPrePPZStrainOctaCommitted = data(36);

	for(i = 0; i < 6; i++) workV6[i] = data(i+37);
  currentStress.setData(workV6);

	for(i = 0; i < 6; i++) workV6[i] = data(i+43);
  currentStrain.setData(workV6);

	for(i = 0; i < 6; i++) workV6[i] = data(i+49);
  PPZPivotCommitted.setData(workV6);

	for(i = 0; i < 6; i++) workV6[i] = data(i+55);
  PPZCenterCommitted.setData(workV6);

	for(i = 0; i < 6; i++) workV6[i] = data(i+61);
  lockStressCommitted.setData(workV6);

	for(i = 0; i < 6; i++) workV6[i] = data(i+67);
  reversalStressCommitted.setData(workV6);

	for(i = 0; i < numOfSurfaces; i++) {
		int k = 73 + i*8;
    workV6(0) = data(k+2);
    workV6(1) = data(k+3);
    workV6(2) = data(k+4);
    workV6(3) = data(k+5);
    workV6(4) = data(k+6);
    workV6(5) = data(k+7);
		committedSurfaces[i+1].setData(workV6, data(k), data(k+1));
	}

	return res;
}


int PressureDependMultiYield::getResponse (int responseID, Information &matInfo)
{
  switch (responseID) {
  case -1:
    return -1;
  case 1:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = getCommittedStress();
    return 0;
  case 2:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = getCommittedStrain();
    return 0;
  case 3:
    if (matInfo.theMatrix != 0)
      *(matInfo.theMatrix) = getTangent();
    return 0;
  default:
    return -1;
  }
}


void PressureDependMultiYield::Print(ostream &s, int flag )
{
  s << "PressureDependMultiYield" << endl;
}


const Vector & PressureDependMultiYield::getCommittedStress (void)
{
	double scale = currentStress.deviatorRatio(residualPress)/committedSurfaces[numOfSurfaces].size();
  if (ndm==3) {
		static Vector temp7(7);
		workV6 = currentStress.t2Vector();
    temp7[0] = workV6[0];
    temp7[1] = workV6[1];
    temp7[2] = workV6[2];
    temp7[3] = workV6[3];
    temp7[4] = workV6[4];
    temp7[5] = workV6[5];
    temp7[6] = scale;
		return temp7;
	}

  else {
    static Vector temp4(4);  
		workV6 = currentStress.t2Vector();
    temp4[0] = workV6[0];
    temp4[1] = workV6[1];
    temp4[2] = workV6[3];
    temp4[3] = scale;

    /*
		double volume = -currentStress.volume();
    double ratio = (volume+residualPress)/(-refPressure+residualPress);
    ratio = pow(ratio, 1.-pressDependCoeff);
		temp4[4] = onPPZCommitted;
    temp4[5] = PPZSizeCommitted;
    temp4[6] = cumuDilateStrainOctaCommitted;
    temp4[7] = maxCumuDilateStrainOctaCommitted;
    temp4[8] = cumuTranslateStrainOctaCommitted;
    temp4[9] = prePPZStrainOctaCommitted;
    temp4[10] = PPZPivot.t2Vector()[3];
		temp4[11] = PPZCenter.t2Vector()[3];
		temp4[12] = pressureDCommitted;
    temp4[13] = junk; //!!!!!!!!!!!!!!!!!!!!!!!!
		temp4[14] = committedActiveSurf; */

    return temp4;
  }
}


const Vector & PressureDependMultiYield::getCommittedStrain (void)
{
  if (ndm==3)
    return currentStrain.t2Vector(1);
  else {
		static Vector workV(3);
		workV6 = currentStrain.t2Vector(1);
    workV[0] = workV6[0];
    workV[1] = workV6[1];
    workV[2] = workV6[3];
    return workV;
  }
}


// NOTE: surfaces[0] is not used 
void PressureDependMultiYield::setUpSurfaces (void)
{ 
  double refStrain, peakShear, coneHeight;

  double sinPhi = sin(frictionAngle * pi/180.);
  double Mnys = 6.*sinPhi/(3.-sinPhi);
  double sinPhiPT = sin(phaseTransfAngle * pi/180.);
  stressRatioPT = 6.*sinPhiPT/(3.-sinPhiPT);
  residualPress = 3.* cohesion / (sqrt(2.) * Mnys);

  // a small nonzero residualPress for numerical purpose only
  if (residualPress < 0.01) residualPress = 0.01; 
  coneHeight = - (refPressure - residualPress);
  peakShear = sqrt(2.) * coneHeight * Mnys / 3.; 
  refStrain = (peakShearStrain * peakShear) 
			        / (refShearModulus * peakShearStrain - peakShear);

  double stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
  double ratio1, ratio2;
  double stressInc = peakShear / numOfSurfaces;

  for (int ii=1; ii<=numOfSurfaces; ii++){
    stress1 = ii * stressInc; 
    stress2 = stress1 + stressInc;
    ratio1 = 3. * stress1 / sqrt(2.) / coneHeight;
    ratio2 = 3. * stress2 / sqrt(2.) / coneHeight;
    strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
    strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);
    
    if (ratio1 <= stressRatioPT && ratio2 >= stressRatioPT) {
      double ratio = (ratio2 - stressRatioPT)/(ratio2 - ratio1);
      strainPTOcta = strain2 - ratio * (strain2 - strain1);
    }

    size = ratio1;
    elasto_plast_modul = 2.*(stress2 - stress1)/(strain2 - strain1);
    if ( (2.*refShearModulus - elasto_plast_modul) <= 0) 
      plast_modul = UP_LIMIT;
    else 
      plast_modul = (2.*refShearModulus * elasto_plast_modul)/
	(2.*refShearModulus - elasto_plast_modul);
    if (plast_modul < 0) plast_modul = 0;
    if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;
    if (ii==numOfSurfaces) plast_modul = 0;
    workV6.Zero();
    committedSurfaces[ii] = MultiYieldSurface(workV6,size,plast_modul);
  }  // ii   
}


double PressureDependMultiYield::yieldFunc(const T2Vector & stress, 
					   const MultiYieldSurface * surfaces, 
					   int surfaceNum)
{
  double coneHeight = stress.volume() - residualPress;
  //workV6 = stress.deviator() - surfaces[surfaceNum].center()*coneHeight;
  workV6 = stress.deviator();
  workV6.addVector(1.0, surfaces[surfaceNum].center(), -coneHeight);

  double sz = surfaces[surfaceNum].size()*coneHeight;

  return 3./2.*(workV6 && workV6) - sz * sz;
}


void PressureDependMultiYield::deviatorScaling(T2Vector & stress, 
					       const MultiYieldSurface * surfaces, 
					       int surfaceNum)
{
  double diff = yieldFunc(stress, surfaces, surfaceNum);
  double coneHeight = stress.volume() - residualPress;

  if ( surfaceNum < numOfSurfaces && diff < 0. ) {
    double sz = -surfaces[surfaceNum].size()*coneHeight;
    double deviaSz = sqrt(sz*sz + diff);
    static Vector devia(6);
    devia = stress.deviator(); 
    workV6 = devia;
    workV6.addVector(1.0, surfaces[surfaceNum].center(), -coneHeight);
    double coeff = (sz-deviaSz) / deviaSz;
    if (coeff < 1.e-13) coeff = 1.e-13;
    //devia += workV6 * coeff;
    devia.addVector(1.0, workV6, coeff);
    stress.setData(devia, stress.volume());
    deviatorScaling(stress, surfaces, surfaceNum);  // recursive call
  }

  if (surfaceNum==numOfSurfaces && fabs(diff) > LOW_LIMIT) {
    double sz = -surfaces[surfaceNum].size()*coneHeight;
    //workV6 = stress.deviator() * sz/sqrt(diff+sz*sz);
    workV6 = stress.deviator();
    workV6 *= sz/sqrt(diff+sz*sz);
    stress.setData(workV6, stress.volume());
  }
}


void PressureDependMultiYield::initSurfaceUpdate(void)
{
  if (committedActiveSurf == 0) return; 
  
  double coneHeight = - (currentStress.volume() - residualPress);
  static Vector devia(6); 
  devia = currentStress.deviator();
  double Ms = sqrt(3./2.*(devia && devia));

  if (committedActiveSurf < numOfSurfaces) { // failure surface can't move
    //workV6 = devia * (1. - committedSurfaces[committedActiveSurf].size()*coneHeight / Ms);
    workV6.addVector(0.0, devia, (1. - committedSurfaces[committedActiveSurf].size()*coneHeight / Ms));
	 
    //workV6 = workV6 / -coneHeight;
    workV6 /= -coneHeight;
    committedSurfaces[committedActiveSurf].setCenter(workV6);
  }

  for (int i=1; i<committedActiveSurf; i++) {
    //workV6 = devia * (1. - committedSurfaces[i].size()*coneHeight / Ms);
    workV6.addVector(0.0, devia, (1. - committedSurfaces[i].size()*coneHeight / Ms));
    //workV6 = workV6 / -coneHeight;
    workV6 /= -coneHeight;
    committedSurfaces[i].setCenter(workV6); 
    theSurfaces[i] = committedSurfaces[i];
  }
  activeSurfaceNum = committedActiveSurf;
}


void PressureDependMultiYield::initStrainUpdate(void)
{
  // elastic strain state
  double stressRatio = currentStress.deviatorRatio(residualPress);
  double ratio = (-currentStress.volume()+residualPress)/(-refPressure+residualPress);
  ratio = pow(ratio, 1.-pressDependCoeff);
  modulusFactor = getModulusFactor(currentStress);
  double shearCoeff = 1./(2.*refShearModulus*modulusFactor);
  double bulkCoeff = 1./(3.*refBulkModulus*modulusFactor);

  //currentStrain = currentStress.deviator()*shearCoeff
  //              + currentStress.volume()*bulkCoeff;

  // modified fmk as discussed with z.yang
  workV6.addVector(0.0, currentStress.deviator(), shearCoeff);
  currentStrain.setData(workV6, currentStress.volume()*bulkCoeff);
	
  double octalStrain = currentStrain.octahedralShear(1);
  if (octalStrain <= LOW_LIMIT) octalStrain = LOW_LIMIT;

  // plastic strain state (scaled from elastic strain)
  double scale, PPZLimit;
  if (stressRatio >= stressRatioPT) {  //above PT
    onPPZ = 2;
    prePPZStrainOcta = strainPTOcta * ratio;
    PPZLimit = getPPZLimits(1,currentStress);
    scale = sqrt(prePPZStrainOcta+PPZLimit)/octalStrain;
  }
  else {  // below PT
    onPPZ = -1;
    prePPZStrainOcta = octalStrain;
    if (prePPZStrainOcta > strainPTOcta * ratio) prePPZStrainOcta=strainPTOcta*ratio;
    scale = sqrt(prePPZStrainOcta)/octalStrain;
  }
  //currentStrain.setData(currentStrain.deviator()*scale, currentStrain.volume());
  workV6.addVector(0.0, currentStrain.deviator(), scale);
  currentStrain.setData(workV6, currentStrain.volume());
  PPZPivot = currentStrain;
}


double PressureDependMultiYield::getModulusFactor(T2Vector & stress)
{
  double conHeig = stress.volume() - residualPress;
  double scale = conHeig / (refPressure-residualPress);
           
  return pow(scale, pressDependCoeff); 
}


void PressureDependMultiYield::setTrialStress(T2Vector & stress)
{
  modulusFactor = getModulusFactor(stress);
  //workV6 = stress.deviator() 
  //	             + subStrainRate.deviator()*2.*refShearModulus*modulusFactor;
  workV6 = stress.deviator();
  workV6.addVector(1.0, subStrainRate.deviator(), 2*refShearModulus*modulusFactor);
  
  double volume = stress.volume() 
    + subStrainRate.volume()*3.*refBulkModulus*modulusFactor;
  if (volume > 0.) volume = 0.;
  trialStress.setData(workV6, volume);
}


int PressureDependMultiYield::setSubStrainRate(void)
{
  if (activeSurfaceNum==numOfSurfaces) return 1;
  if (strainRate.isZero()) return 0;  

  double elast_plast_modulus;
  double conHeig = -(currentStress.volume() - residualPress);
  double factor = getModulusFactor(currentStress);
  if (activeSurfaceNum==0) 
    elast_plast_modulus = 2*refShearModulus*factor;
  else {
    double plast_modulus = theSurfaces[activeSurfaceNum].modulus()*factor;
    elast_plast_modulus = 2*refShearModulus*factor*plast_modulus 
      / (2*refShearModulus*factor+plast_modulus);
  }
  //workV6 = strainRate.deviator()*elast_plast_modulus;
  workV6.addVector(0.0, strainRate.deviator(), elast_plast_modulus);
  workT2V.setData(workV6,0);

  double singleCross = theSurfaces[numOfSurfaces].size()*conHeig / numOfSurfaces;
  double totalCross = 3.*workT2V.octahedralShear() / sqrt(2.);
  int numOfSub = totalCross/singleCross + 1;
  if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;
	
  int numOfSub1 = strainRate.octahedralShear(1) / 1.0e-4;
  if (numOfSub1 > numOfSub) numOfSub = numOfSub1;

  workV6.addVector(0.0, strainRate.t2Vector(), 1.0/numOfSub);

  subStrainRate.setData(workV6);

  return numOfSub;
}


void
PressureDependMultiYield::getContactStress(T2Vector &contactStress)
{
  double conHeig = trialStress.volume() - residualPress;
  static Vector center(6);
  center = theSurfaces[activeSurfaceNum].center(); 
  //workV6 = trialStress.deviator() - center*conHeig;
  workV6 = trialStress.deviator();
  workV6.addVector(1.0, center, -conHeig);
  double Ms = sqrt(3./2.*(workV6 && workV6));
  //workV6 = workV6 * theSurfaces[activeSurfaceNum].size()*(-conHeig) / Ms + center*conHeig;
  workV6.addVector(theSurfaces[activeSurfaceNum].size()*(-conHeig) / Ms, center, conHeig);
  //return T2Vector(workV6,trialStress.volume()); 
  contactStress.setData(workV6,trialStress.volume()); 
}


int PressureDependMultiYield::isLoadReversal(void)
{
  if(activeSurfaceNum == 0) return 0;

  getSurfaceNormal(currentStress, workT2V);
  //if (((trialStress.t2Vector() - currentStress.t2Vector()) 
  //	&& workT2V.t2Vector()) < 0) return 1;

  workV6 = trialStress.t2Vector();
  workV6 -= currentStress.t2Vector();
  if ((workV6 && workT2V.t2Vector()) < 0) return 1;

  return 0;   
}


void
PressureDependMultiYield::getSurfaceNormal(const T2Vector & stress, T2Vector &normal)
{
  double conHeig = stress.volume() - residualPress;
  workV6 = stress.deviator();
  static Vector center(6);
  center = theSurfaces[activeSurfaceNum].center(); 
  double sz = theSurfaces[activeSurfaceNum].size();
  double volume = conHeig*((center && center) - 2./3.*sz*sz) - (workV6 && center);
  //workT2V.setData((workV6-center*conHeig)*3., volume);
  workV6.addVector(1.0, center, -conHeig);
  workV6 *= 3.0;
  workT2V.setData(workV6, volume);
  
  normal.setData(workT2V.unitT2Vector());
}


double PressureDependMultiYield::getPlasticPotential(const T2Vector & contactStress,
						     const T2Vector & surfaceNormal)
{
  double plasticPotential, contractRule, unloadRule, dilateRule, shearLoading, temp;

  double contactRatio = contactStress.deviatorRatio(residualPress);
  temp = contactRatio/stressRatioPT;
  double factorPT = (temp*temp - 1)/(temp*temp + 1)/3.;
  double volume = contactStress.volume();
  contractRule = factorPT*contractParam1;
  if (contractRule > 0.) contractRule = -contractRule;
	if (contractRule<-5.0e4) contractRule = -5.0e4;
  temp = currentStress.volume() - pressureD;
  if (temp >= 0.) unloadRule = 0.;
  else {
    double conHeiD = pressureD-residualPress;
    double temp1 = sqrt(3./2.)*currentStress.deviatorLength() 
      + stressRatioPT*conHeiD;
    temp = temp1/(-temp);
    if (temp < theSurfaces[numOfSurfaces].size()) 
      temp = theSurfaces[numOfSurfaces].size();
    temp1 = (reversalStress.volume()-residualPress)/conHeiD;
    unloadRule = -sqrt(3./2.)*surfaceNormal.deviatorLength()*temp1/temp;
  }
	
  double currentRatio = currentStress.deviatorRatio(residualPress);
  double trialRatio = trialStress.deviatorRatio(residualPress);
  shearLoading = currentStress.deviator() && trialStress.deviator();

  if (factorPT < 0.) {  //below PT
    if (pressureD == 0.) 
      plasticPotential = contractRule;
    else if (trialStress.volume() >= pressureD) {
      pressureD = 0.;
      plasticPotential = contractRule;
    }
    else if (trialRatio > currentRatio && shearLoading >= 0.)
      plasticPotential = contractRule;
    else  plasticPotential = unloadRule;
  }
  
  else {  //above PT
    if (trialRatio > currentRatio && shearLoading >= 0.) {  //dilation
      if (pressureD == 0.) pressureD = currentStress.volume();
      reversalStress = currentStress;
      updatePPZ(contactStress);  
      if (onPPZ==-1 || onPPZ==1) return LOCK_VALUE; 
      if (isCriticalState(contactStress))  
	      dilateRule = 0;
      else
        dilateRule = factorPT*dilateParam1*exp(dilateParam2*cumuDilateStrainOcta);
			if (dilateRule>5.0e4) dilateRule = 5.0e4;
      return dilateRule;
    }
    else if (pressureD == 0.) plasticPotential = contractRule;
    else if (trialStress.volume() >= pressureD) {
      pressureD = 0.;
      plasticPotential = contractRule;
    }
    else plasticPotential = unloadRule;
  }

  if (onPPZ > 0) onPPZ = 0;
  if (onPPZ != -1) PPZTranslation(contactStress);  
  if (isCriticalState(contactStress)) return 0.;
  return plasticPotential;
}


int PressureDependMultiYield::isCriticalState(const T2Vector & stress)
{
  double vol = trialStrain.volume()*3.0;
	double etria = einit + vol + vol*einit;
	double ecr1 = volLimit1 - volLimit2*pow(-stress.volume()/pAtm, volLimit3);

	vol = currentStrain.volume()*3.0;
	double ecurr = einit + vol + vol*einit;
	double ecr2 = volLimit1 - volLimit2*pow(-currentStress.volume()/pAtm, volLimit3);

	if (ecurr < ecr2 && etria < ecr1) return 0;
	if (ecurr > ecr2 && etria > ecr1) return 0;	
	return 1;
}


void PressureDependMultiYield::updatePPZ(const T2Vector & contactStress)
{
  // PPZ inactive if liquefyParam1==0.
  if (liquefyParam1==0.) {
    if (onPPZ==2) {
		  workT2V.setData(trialStrain.t2Vector() - PPZPivot.t2Vector());
      cumuDilateStrainOcta = workT2V.octahedralShear(1);
    }
    else if (onPPZ != 2) {
      onPPZ = 2;
      PPZPivot = trialStrain;
      cumuDilateStrainOcta = 0.;
    }
    return;
  }

  // dilation: calc. cumulated dilative strain
  if (onPPZ==2) {
    PPZPivot = trialStrain;
    workV6 = PPZPivot.t2Vector();
    workV6 -= PPZCenter.t2Vector();
    workT2V.setData(workV6);
    //cumuDilateStrainOcta = workT2V.octahedralShear(1) - PPZSize;
		cumuDilateStrainOcta += subStrainRate.octahedralShear(1);
    if (cumuDilateStrainOcta > maxCumuDilateStrainOcta) 
      maxCumuDilateStrainOcta = cumuDilateStrainOcta;
    return;
  }

  // calc. PPZ size.
  double PPZLimit = getPPZLimits(1,contactStress);
	double TransLimit = getPPZLimits(2,contactStress);

  if (onPPZ==-1 || onPPZ==0) {
    workV6 = trialStrain.t2Vector();
    workV6 -= PPZPivot.t2Vector();
    workT2V.setData(workV6);
		double temp = workT2V.octahedralShear(1);
		if (temp > cumuDilateStrainOcta) {
      double volume = -contactStress.volume();
      oppoPrePPZStrainOcta = prePPZStrainOcta;
      double ratio = (volume+residualPress)/(-refPressure+residualPress);
      ratio = pow(ratio, 1.-pressDependCoeff);
      prePPZStrainOcta = ratio * strainPTOcta;
      if (oppoPrePPZStrainOcta == 0.) oppoPrePPZStrainOcta = prePPZStrainOcta;
		}
  }
  if (onPPZ > -1) PPZSize = PPZLimit 
    + (prePPZStrainOcta+oppoPrePPZStrainOcta+TransLimit+maxCumuDilateStrainOcta)/2.;
	else PPZSize = PPZLimit 
    + (prePPZStrainOcta+oppoPrePPZStrainOcta+maxCumuDilateStrainOcta)/2.;

  // calc. new PPZ center.
  if (onPPZ==0 || onPPZ==1) { 
    //workT2V.setData(PPZPivot.t2Vector() - PPZCenter.t2Vector());
    workV6 = PPZPivot.t2Vector();
    workV6 -= PPZCenter.t2Vector();
    workT2V.setData(workV6);
    
    double coeff = (PPZSize-cumuTranslateStrainOcta)/workT2V.octahedralShear(1);
    //PPZCenter.setData(PPZPivot.t2Vector() - workT2V.t2Vector()*coeff);
    workV6 = PPZPivot.t2Vector();
    workV6.addVector(1.0, workT2V.t2Vector(), -coeff);
    PPZCenter.setData(workV6);
  }

  //workT2V.setData(trialStrain.t2Vector() - PPZCenter.t2Vector());
  workV6 = trialStrain.t2Vector();
  workV6 -= PPZCenter.t2Vector();
  workT2V.setData(workV6);
  double temp = subStrainRate.t2Vector() && workV6;

  if (workT2V.octahedralShear(1) > PPZSize && temp > 0.) {  //outside PPZ
    workV6 = trialStrain.t2Vector();
    workV6 -= PPZPivot.t2Vector();
    workT2V.setData(workV6);
		double temp1 = workT2V.octahedralShear(1);
		if (temp1 > cumuDilateStrainOcta) {
      cumuDilateStrainOcta = 0.;
      if (PPZLimit == 0.) maxCumuDilateStrainOcta = 0.;
		}
    onPPZ = 2;
    PPZPivot = trialStrain;
    cumuTranslateStrainOcta = 0.;
  }
  else {  //inside PPZ
    if (onPPZ == 0 || onPPZ == 1) PPZTranslation(contactStress);
    if (onPPZ == -1 || onPPZ == 0) lockStress = contactStress;
    if (onPPZ == 0) onPPZ = 1;
  }
}


void PressureDependMultiYield::PPZTranslation(const T2Vector & contactStress)
{
  if (liquefyParam1==0.) return;
  
  double PPZTransLimit = getPPZLimits(2,contactStress);

  //workT2V.setData(trialStrain.deviator()-PPZPivot.deviator());
  workV6 = trialStrain.deviator();
  workV6 -= PPZPivot.deviator();
  workT2V.setData(workV6);
  
  double temp = workT2V.octahedralShear(1);
  if (cumuTranslateStrainOcta < temp) cumuTranslateStrainOcta = temp;
	if (maxCumuDilateStrainOcta == 0.) temp = PPZTransLimit;
	else temp = PPZTransLimit*cumuDilateStrainOcta/maxCumuDilateStrainOcta;
  if (cumuTranslateStrainOcta > temp) cumuTranslateStrainOcta = temp;
}


double PressureDependMultiYield::getPPZLimits(int which, const T2Vector & contactStress)
{
  double PPZLimit, temp;
  double volume = -contactStress.volume();

  if (volume >= liquefyParam1) PPZLimit = 0.;
  else {
    temp = volume*pi/liquefyParam1/2.;
		// liquefyParam3 = 3.0 default
    PPZLimit = liquefyParam2 * pow(cos(temp), 3.);
  }
  
  if (which==1) 
    return PPZLimit;
  else if (which==2) 
    return liquefyParam4 * PPZLimit;
  else {
    cerr << "FATAL:PressureDependMultiYield::getPPZLimits: unknown argument value" << endl;
    g3ErrorHandler->fatal(" ");
    return 0.0;
  }
}


double PressureDependMultiYield::getLoadingFunc(const T2Vector & contactStress, 
						const T2Vector & surfaceNormal,
						double plasticPotential,
						int crossedSurface)
{
  double loadingFunc, limit;
  double modul = theSurfaces[activeSurfaceNum].modulus();
  double temp1 = 2. * refShearModulus * modulusFactor 
    * (surfaceNormal.deviator() && surfaceNormal.deviator()) ;
  double temp2 = 9. * refBulkModulus * modulusFactor 
    * surfaceNormal.volume() * plasticPotential ;
	
  //for the first crossing 
  double temp = temp1 + temp2 + modul * modulusFactor;
  if (activeSurfaceNum == numOfSurfaces) 
    limit = theSurfaces[activeSurfaceNum-1].modulus() * modulusFactor /2.;
  else limit = modul * modulusFactor /2.;
  if (temp < limit) {
    plasticPotential = (temp2 + limit - temp)/(9. * refBulkModulus * modulusFactor 
					       * surfaceNormal.volume());
    temp = limit;
  }
  //loadingFunc = (surfaceNormal.t2Vector() 
  //	             && (trialStress.deviator()-contactStress.deviator()))/temp;
  workV6 = trialStress.deviator();
  workV6 -= contactStress.deviator();  
  loadingFunc = (surfaceNormal.t2Vector() && workV6) /temp;
  
  if (loadingFunc < 0.) loadingFunc = 0;

  //for more than one crossing 
  if(crossedSurface) {
    temp = (theSurfaces[activeSurfaceNum-1].modulus() - modul)
      /theSurfaces[activeSurfaceNum-1].modulus();
    loadingFunc *= temp;
  }

  return loadingFunc;
}


int PressureDependMultiYield::stressCorrection(int crossedSurface)
{
  static T2Vector contactStress;
  getContactStress(contactStress);
  static T2Vector surfNormal;
  getSurfaceNormal(contactStress, surfNormal);
  double plasticPotential = getPlasticPotential(contactStress,surfNormal);
  if (plasticPotential==LOCK_VALUE && (onPPZ == -1 || onPPZ == 1)) {
    trialStress = lockStress;
    return 1;
  }

	double tVolume = trialStress.volume();
  double loadingFunc = getLoadingFunc(contactStress, surfNormal, 
				      plasticPotential, crossedSurface);
  double volume = tVolume - plasticPotential*3*refBulkModulus*modulusFactor*loadingFunc;

  //workV6 = trialStress.deviator() 
  //	         - surfNormal.deviator()*2*refShearModulus*modulusFactor*loadingFunc;
  workV6 = trialStress.deviator();
	
  if (volume > 0. && volume != tVolume) {
		double coeff = tVolume / (tVolume - volume);
    coeff *= -2*refShearModulus*modulusFactor*loadingFunc;
    workV6.addVector(1.0, surfNormal.deviator(), coeff);
		volume = 0.;
  } else if (volume > 0.) {
		volume = 0.;
	} else {
		double coeff = -2*refShearModulus*modulusFactor*loadingFunc;
    workV6.addVector(1.0, surfNormal.deviator(), coeff);
	}
	/*
	  if (volume>0.)volume = 0.;
		double coeff = -2*refShearModulus*modulusFactor*loadingFunc;
    workV6.addVector(1.0, surfNormal.deviator(), coeff);
  */
	trialStress.setData(workV6, volume);
  deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

  if (isCrossingNextSurface()) {
    activeSurfaceNum++;
    return stressCorrection(1);  //recursive call
  }
  return 0;
}


void PressureDependMultiYield::updateActiveSurface(void)
{
  if (activeSurfaceNum == numOfSurfaces) return;

  double A, B, C, X;
  static Vector t1(6);
  static Vector t2(6);
  static Vector center(6);
  static Vector outcenter(6);
  double conHeig = trialStress.volume() - residualPress;
  center = theSurfaces[activeSurfaceNum].center();
  double size = theSurfaces[activeSurfaceNum].size();
  outcenter = theSurfaces[activeSurfaceNum+1].center();
  double outsize = theSurfaces[activeSurfaceNum+1].size();

  //t1 = trialStress.deviator() - center*conHeig;
  //t2 = (center - outcenter)*conHeig;
  t1 = trialStress.deviator();
  t1.addVector(1.0, center, -conHeig);
  t2 = center;
  t2 -= outcenter;
  t2 *=conHeig;
 
  A = t1 && t1;
  B = 2. * (t1 && t2);
  C = (t2 && t2) - 2./3.* outsize * outsize * conHeig * conHeig;
  X = secondOrderEqn(A,B,C,0);
  if ( fabs(X-1.) < LOW_LIMIT ) X = 1.;
  if (X < 1.) return;

  if (X < 1.){
    //t2 = trialStress.deviator() - outcenter*conHeig;
    t2 = trialStress.deviator();
    t2.addVector(1.0, outcenter, -conHeig);
    
    double xx1 = (t2 && t2) - 2./3.* outsize * outsize * conHeig * conHeig;
    double xx2 = (t1 && t1) - 2./3.* size * size * conHeig * conHeig;
    cerr << "FATAL:PressureDependMultiYield::updateActiveSurface(): error in Direction of surface motion.\n" 
	 << "X-1= " << X-1 <<" A= "<<A<<" B= "<<B<<" C= "<<C <<" M= "<<activeSurfaceNum<<" low_limit="<<LOW_LIMIT<<"\n"
	 << "diff1= "<<xx1 <<" diff2= "<<xx2 <<" p= "<<conHeig<<" size= "<<size<<" outs= "<<outsize<<endl; 
    g3ErrorHandler->fatal(" ");
  }

  //workV6 = (t1 * X + center*conHeig) * (1. - size / outsize) 
  //	     - (center - outcenter * size / outsize) * conHeig;
  
  workV6.addVector(0.0, t1, X);
  workV6.addVector(1.0, center, conHeig);
  workV6 *= (1.0 - size / outsize);
  t2 = center;
  t2.addVector(1.0, outcenter, -size/outsize);
  t2 *= conHeig;
  workV6 -= t2;

  workT2V.setData(workV6);
  if (workT2V.deviatorLength() < LOW_LIMIT) return;

  workV6 = workT2V.deviator();  
  A = conHeig * conHeig * (workV6 && workV6);
  B = 2 * conHeig * (t1 && workV6);
  if (fabs(B) < LOW_LIMIT) B = 0.; 
  C = (t1 && t1) - 2./3.* size * size * conHeig * conHeig;
  if ( fabs(C) < LOW_LIMIT || fabs(C)/(t1 && t1) < LOW_LIMIT ) return;
  if (B > 0. || C < 0.) {
    cerr << "FATAL:PressureDependMultiYield::updateActiveSurface(): error in surface motion.\n" 
	 << "A= " <<A <<" B= " <<B <<" C= "<<C <<" (t1&&t1)= "<<(t1&&t1) <<endl; 
    g3ErrorHandler->fatal(" ");
  }
  X = secondOrderEqn(A,B,C,1);  

  center.addVector(1.0, workV6, -X);
  theSurfaces[activeSurfaceNum].setCenter(center);
}      


void PressureDependMultiYield::updateInnerSurface(void)
{
	if (activeSurfaceNum <= 1) return;
	static Vector devia(6);
	static Vector center(6);

	double conHeig = currentStress.volume() - residualPress;
	devia = currentStress.deviator();
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();

	for (int i=1; i<activeSurfaceNum; i++) {
		workV6.addVector(0.0, center, conHeig);
		workV6 -= devia;
		workV6 *= theSurfaces[i].size()/size;
		workV6 += devia;

		//workV6 = devia - (devia - center*conHeig) * theSurfaces[i].size() / size;
		
		workV6 /= conHeig;
		theSurfaces[i].setCenter(workV6);
	}
}


int PressureDependMultiYield:: isCrossingNextSurface(void)
{
  if (activeSurfaceNum == numOfSurfaces) return 0;  

  if(yieldFunc(trialStress, theSurfaces, activeSurfaceNum+1) > 0) return 1;
  
  return 0;
}
 
