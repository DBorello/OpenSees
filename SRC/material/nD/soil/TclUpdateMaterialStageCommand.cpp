// $Revision: 1.8 $
// $Date: 2003-08-29 00:23:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/TclUpdateMaterialStageCommand.cpp,v $
                                                                        
// Written: ZHY

#include <TclModelBuilder.h>

#include <PressureIndependMultiYield.h>
#include <PressureDependMultiYield.h>
#include <FluidSolidPorousMaterial.h>
#include <Information.h>

#include <string.h>

// by ZHY
int
TclModelBuilderUpdateMaterialStageCommand(ClientData clientData, 
					  Tcl_Interp *interp, 
					  int argc,
					  TCL_Char **argv, 
					  TclModelBuilder *theTclBuilder)
{
  if (argc < 5) {
      opserr << "WARNING insufficient number of UpdateMaterialStage arguments\n";
      opserr << "Want: UpdateMaterialStage material matTag? stage value?" << endln;
      return TCL_ERROR;
  }

  if (strcmp(argv[1],"-material") != 0) {
      opserr << "WARNING UpdateMaterialStage: Only accept parameter '-material' for now" << endln;
      return TCL_ERROR;		
  }		

  int tag, value; 

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING MYSstage: invalid material tag" << endln;
      return TCL_ERROR;		
  }

  NDMaterial * a = theTclBuilder->getNDMaterial(tag);
  if (a==0) {
      opserr << "WARNING UpdateMaterialStage: couldn't get NDmaterial tagged: " << tag << endln;
      return TCL_ERROR;		
  }

  if (strcmp(argv[3],"-stage") != 0) {
      opserr << "WARNING UpdateMaterialStage: Only accept parameter '-stage' for now" << endln;
      return TCL_ERROR;		
  }		

  if (Tcl_GetInt(interp, argv[4], &value) != TCL_OK) {
      opserr << "WARNING UpdateMaterialStage: invalid parameter value" << endln;
      return TCL_ERROR;		
  }	

  const char * c = a->getType();
	if (strcmp(c, "PlaneStrain") == 0 || 
      strcmp(c, "ThreeDimensional") == 0 ) {
      Information info;
      a->updateParameter(value,info); 
  }
  else {
      opserr << "WARNING UpdateMaterialStage: The tagged is not a "<<endln;
      opserr << "PressureDependMultiYield/PressureIndependMultiYield/FluidSolidPorous material. " << endln;
      return TCL_ERROR;		
  }

  return TCL_OK;
}


int
TclModelBuilderUpdateParameterCommand(ClientData clientData, 
					  Tcl_Interp *interp, 
					  int argc,
					  char **argv, 
					  TclModelBuilder *theTclBuilder)
{
  if (argc < 5) {
      opserr << "WARNING insufficient number of updateParameter arguments\n";
      opserr << "Want: updateParameter -material matNum? -param? newValue?" << endln;
      return TCL_ERROR;
  }

  if (strcmp(argv[1],"-material") != 0) {
      opserr << "WARNING UpdateParameter: Only accept parameter '-material' for now" << endln;
      return TCL_ERROR;		
  }		

  int tag, id; double value; 

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING UpdateParameter: invalid material tag" << endln;
      return TCL_ERROR;		
  }

  NDMaterial * a = theTclBuilder->getNDMaterial(tag);
  if (a==0) {
      opserr << "WARNING UpdateParameter: couldn't get NDmaterial tagged: " << tag << endln;
      return TCL_ERROR;		
  }

  if (strcmp(argv[3],"-refG") == 0) 
		id = 10;
	else if (strcmp(argv[3],"-refB") == 0) 
	  id = 11;
	else {
      opserr << "WARNING UpdateParameter: Only accept parameter '-refG' or '-refB' for now" << endln;
      return TCL_ERROR;		
  }		

  if (Tcl_GetDouble(interp, argv[4], &value) != TCL_OK) {
      opserr << "WARNING UpdateParameter: invalid parameter value" << endln;
      return TCL_ERROR;		
  }	

  const char * c = a->getType();
	if (strcmp(c, "PlaneStrain") == 0 || 
      strcmp(c, "ThreeDimensional") == 0 ) {
      Information info;
      info.setDouble(value);
      a->updateParameter(id,info); 
  }
  else {
      opserr << "WARNING UpdateParameter: The tagged is not a "<<endln;
      opserr << "PressureDependMultiYield/PressureIndependMultiYield/FluidSolidPorous material. " << endln;
      return TCL_ERROR;		
  }

  return TCL_OK;
}
