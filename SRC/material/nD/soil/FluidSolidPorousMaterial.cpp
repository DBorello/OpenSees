// $Revision: 1.7 $
// $Date: 2002-02-08 19:51:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/FluidSolidPorousMaterial.cpp,v $
                                                                        
// Written: ZHY

//
// FluidSolidPorousMaterial.cpp
// -------------------
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include <FluidSolidPorousMaterial.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>

int FluidSolidPorousMaterial::loadStage = 0;
double FluidSolidPorousMaterial::AtmoPress = 0.;

Vector FluidSolidPorousMaterial::workV3(3);
Vector FluidSolidPorousMaterial::workV6(6);
Matrix FluidSolidPorousMaterial::workM3(3,3);
Matrix FluidSolidPorousMaterial::workM6(6,6);

FluidSolidPorousMaterial::FluidSolidPorousMaterial (int tag, int nd, NDMaterial &soilMat,
                                      double combinedBulkModul, double atm)
 : NDMaterial(tag, ND_TAG_FluidSolidPorousMaterial)
{
	if (combinedBulkModul < 0) {
		cerr << "WARNING:FluidSolidPorousMaterial::FluidSolidPorousMaterial: combinedBulkModulus < 0" << endl;
	  cerr << "Will reset to 0." <<endl;
    combinedBulkModul = 0.;
  }
	ndm = nd;
	loadStage = 0;  //default
  theSoilMaterial = soilMat.getCopy();
	AtmoPress = atm;
  combinedBulkModulus = combinedBulkModul;
  trialExcessPressure = currentExcessPressure = 0.;
	trialVolumeStrain = currentVolumeStrain = 0.;
}
   

FluidSolidPorousMaterial::FluidSolidPorousMaterial () 
 : NDMaterial(0,ND_TAG_FluidSolidPorousMaterial), theSoilMaterial(0)
{
	ndm = 3; 
	combinedBulkModulus = 0.;
  trialExcessPressure = currentExcessPressure = 0.;
	trialVolumeStrain = currentVolumeStrain = 0.;
}


FluidSolidPorousMaterial::FluidSolidPorousMaterial (const FluidSolidPorousMaterial & a)
 : NDMaterial(a.getTag(),ND_TAG_FluidSolidPorousMaterial)
{
	ndm = a.ndm;
	combinedBulkModulus = a.combinedBulkModulus;
  theSoilMaterial = a.theSoilMaterial->getCopy();
  trialExcessPressure = a.trialExcessPressure;
	currentExcessPressure = a.currentExcessPressure;
	trialVolumeStrain = a.trialVolumeStrain;
	currentVolumeStrain = a.currentVolumeStrain;
}


FluidSolidPorousMaterial::~FluidSolidPorousMaterial ()
{
	if (theSoilMaterial != 0)
		delete theSoilMaterial;
}


int FluidSolidPorousMaterial::setTrialStrain (const Vector &strain)
{
	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = strain[0]+strain[1]+strain[2];
	else {
		cerr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal(" ");
	}

  return theSoilMaterial->setTrialStrain(strain);
}


int FluidSolidPorousMaterial::setTrialStrain (const Vector &strain, const Vector &rate)
{
	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = strain[0]+strain[1]+strain[2];
	else {
		cerr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal(" ");
	}

  return theSoilMaterial->setTrialStrain(strain, rate);
}


int FluidSolidPorousMaterial::setTrialStrainIncr (const Vector &strain)
{
	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1]+strain[2];
	else {
		cerr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal(" ");
	}

  return theSoilMaterial->setTrialStrainIncr(strain);
}


int FluidSolidPorousMaterial::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1]+strain[2];
	else {
		cerr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal(" ");
	}

  return theSoilMaterial->setTrialStrainIncr(strain, rate);
}


const Matrix & FluidSolidPorousMaterial::getTangent (void)
{
	Matrix *workM = (ndm == 2) ? &workM3 : &workM6;
  
	*workM = theSoilMaterial->getTangent();

	if (loadStage != 0) { 
	  for (int i=0; i<ndm; i++) 
		  for (int j=0; j<ndm; j++) 
			  (*workM)(i,j) = (*workM)(i,j) + combinedBulkModulus;
  }

	return *workM;
}


const Vector & FluidSolidPorousMaterial::getStress (void)
{
	Vector *workV = (ndm == 2) ? &workV3 : &workV6;
  
	*workV = theSoilMaterial->getStress();

	if (loadStage != 0) { 
		trialExcessPressure = currentExcessPressure;
    trialExcessPressure += 
			       (trialVolumeStrain - currentVolumeStrain) * combinedBulkModulus;
	  for (int i=0; i<ndm; i++) 
      (*workV)[i] += trialExcessPressure;
  }

  return *workV;
}


int FluidSolidPorousMaterial::updateParameter(int responseID, Information &info)
{
	loadStage = responseID;
	return 0;
}


const Vector & FluidSolidPorousMaterial::getCommittedStress (void)
{
	return theSoilMaterial->getCommittedStress();
}


const Vector & FluidSolidPorousMaterial::getCommittedStrain (void)
{
	return theSoilMaterial->getCommittedStrain();
}


const Vector & FluidSolidPorousMaterial::getCommittedPressure (void)
{
	static Vector temp(2);

	temp[0] = currentExcessPressure;
  if (temp[0] == 0.) temp[1] = 0.;
	else {
  	if (ndm == 3) 
	    temp[1] = (theSoilMaterial->getCommittedStress()[0]
		            + theSoilMaterial->getCommittedStress()[1]
			  				+ theSoilMaterial->getCommittedStress()[2])/3.;
	  else 
	    temp[1] = (theSoilMaterial->getCommittedStress()[0]
		            + theSoilMaterial->getCommittedStress()[1])/2.;

		temp[1] = fabs(temp[0])/(fabs(temp[1])+fabs(temp[0]));
	}
	return temp;
}


const Vector & FluidSolidPorousMaterial::getStrain (void)
{
  return theSoilMaterial->getStrain();
}


int FluidSolidPorousMaterial::commitState (void)
{
	currentVolumeStrain = trialVolumeStrain;
	if (loadStage != 0) 
		currentExcessPressure = trialExcessPressure;
	else
    currentExcessPressure = 0.;

	return theSoilMaterial->commitState();
}


int FluidSolidPorousMaterial::revertToLastCommit (void)
{
	return theSoilMaterial->revertToLastCommit();
}

int FluidSolidPorousMaterial::revertToStart (void)
{
	return theSoilMaterial->revertToStart();
}

NDMaterial * FluidSolidPorousMaterial::getCopy (void)
{
  FluidSolidPorousMaterial * copy = new FluidSolidPorousMaterial(*this);
	return copy;
}


NDMaterial * FluidSolidPorousMaterial::getCopy (const char *code)
{
	if (strcmp(code,"FluidSolidPorous") == 0 || strcmp(code,"PlaneStrain") == 0 ||
		strcmp(code,"ThreeDimensional") == 0) {
     FluidSolidPorousMaterial * copy = new FluidSolidPorousMaterial(*this);
	   return copy;
	}

	return 0;
}


const char * FluidSolidPorousMaterial::getType (void) const
{
	return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int FluidSolidPorousMaterial::getOrder (void) const
{
	return (ndm == 2) ? 3 : 6;
}


int FluidSolidPorousMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(7);
	data(0) = this->getTag();
	data(1) = ndm;
	data(2) = loadStage;
  data(3) = AtmoPress;
  data(4) = combinedBulkModulus;
	data(5) = currentExcessPressure;
  data(6) = currentVolumeStrain;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"FluidSolidPorousMaterial::sendSelf");
		return res;
	}

	ID classTags(2);

	classTags(0) = theSoilMaterial->getClassTag();
	int matDbTag = theSoilMaterial->getDbTag();
	// NOTE: we do have to ensure that the material has a database
	// tag if we are sending to a database channel.
	if (matDbTag == 0) {
		matDbTag = theChannel.getDbTag();
		if (matDbTag != 0)
			theSoilMaterial->setDbTag(matDbTag);
	}
	classTags(1) = matDbTag;

	res += theChannel.sendID(this->getDbTag(), commitTag, classTags);
	if (res < 0) {
		g3ErrorHandler->warning("WARNING FluidSolidPorousMaterial::sendSelf() - %d failed to send ID\n",
			this->getTag());
		return res;
	}

	// Finally, asks the material object to send itself
	res += theSoilMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("WARNING FluidSolidPorousMaterial::sendSelf() - %d failed to send its Material\n",this->getTag());
		return res;
	}

	return res;
}


int FluidSolidPorousMaterial::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)    
{
	int res = 0;

	static Vector data(7);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"FluidSolidPorousMaterial::recvSelf");
		return res;
	}
    
	this->setTag((int)data(0));
	ndm = data(1);
	loadStage = data(2);
  AtmoPress = data(3);
  combinedBulkModulus = data(4);
	currentExcessPressure = data(5);
  currentVolumeStrain = data(6);

	// now receives the ids of its material
	ID classTags(2);

	res += theChannel.recvID(this->getDbTag(), commitTag, classTags);
	if (res < 0)  {
		g3ErrorHandler->warning("FluidSolidPorousMaterial::recvSelf() - %s\n",
			    "failed to recv ID data");
		return res;
	}    
	
	int matClassTag = classTags(0);
	int matDbTag = classTags(1);
	// Check that material is of the right type; if not,
	// delete it and create a new one of the right type
	if (theSoilMaterial->getClassTag() != matClassTag) {
		delete theSoilMaterial;
		theSoilMaterial = theBroker.getNewNDMaterial(matClassTag);
		if (theSoilMaterial == 0) {
			g3ErrorHandler->fatal("FluidSolidPorousMaterial::recvSelf() - %s %d\n",
				"Broker could not create NDMaterial of class type",matClassTag);
			return -1;
		}
	}

	// Receive the material
	theSoilMaterial->setDbTag(matDbTag);
	res += theSoilMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("FluidSolidPorousMaterial::recvSelf() - material failed to recv itself");
		return res;
	}

	return res;
}



Response*
FluidSolidPorousMaterial::setResponse (char **argv, int argc, Information &matInfo)
{
    if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getCommittedStress());

    else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getCommittedStrain());
    
	else if (strcmp(argv[0],"tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());

    else if (strcmp(argv[0],"pressure") == 0)
		return new MaterialResponse(this, 4, this->getCommittedPressure());

    else
		return 0;
}


int FluidSolidPorousMaterial::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
		case 1:
			return matInfo.setVector(this->getCommittedStress());

		case 2:
			return matInfo.setVector(this->getCommittedStrain());

		case 3:
			return matInfo.setMatrix(this->getTangent());
			
		case 4:
			return matInfo.setVector(this->getCommittedPressure());

		default:
			return -1;
	}
}


void FluidSolidPorousMaterial::Print(ostream &s, int flag )
{
	s << "FluidSolidPorousMaterial" << endl;
}










