// NewCommand.cpp : Defines the entry point for the DLL application and
//    a function that can be called to set the global pointer variables in 
//    the dll to be the same as those in the existing process address space.

#include <elementAPI.h>

#include <OPS_Stream.h>
#include <Domain.h>
#include <TclModelBuilder.h>

#include <windows.h>

#include <SimulationInformation.h>
SimulationInformation simulationInfo;

#define DllExport _declspec(dllexport)

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    return TRUE;
}

OPS_Stream *opserrPtr =0;

/*
extern "C" DllExport
void setGlobalPointers(OPS_Stream *theErrorStreamPtr)
{
   opserrPtr = theErrorStreamPtr;
}
*/
/*
int __cdecl OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  Tcl_Exit(0);
  return 0;
}
*/

typedef int (*OPS_ErrorPtrType)(char *, int);
typedef int (*OPS_GetIntInputPtrType)(int *, int *);
typedef int (*OPS_GetDoubleInputPtrType)(int *, double *);
typedef int (*OPS_AllocateElementPtrType)(eleObj *, int *matTags, int *maType);
typedef int (*OPS_AllocateMaterialPtrType)(matObj *);
typedef UniaxialMaterial *(*OPS_GetUniaxialMaterialPtrType)(int matTag);
typedef NDMaterial * (*OPS_GetNDMaterialPtrType)(int matTag);

//int    OPS_InvokeMaterial(struct eleObj *, int *,modelState *, double *, double *, double *, int *);

OPS_ErrorPtrType OPS_ErrorPtr =0;
OPS_GetIntInputPtrType   OPS_GetIntInputPtr =0;
OPS_GetDoubleInputPtrType OPS_GetDoubleInputPtr =0;
OPS_AllocateElementPtrType OPS_AllocateElementPtr =0;
OPS_AllocateMaterialPtrType OPS_AllocateMaterialPtr =0;
OPS_GetUniaxialMaterialPtrType OPS_GetUniaxialMaterialPtr = 0;
OPS_GetNDMaterialPtrType OPS_GetNDMaterialPtr = 0;

extern "C" DllExport
void setGlobalPointers(OPS_Stream *theErrorStreamPtr,
				OPS_ErrorPtrType          errorFunct,
				OPS_GetIntInputPtrType    getIntInputFunct,
				OPS_GetDoubleInputPtrType getDoubleInputFunct,
				OPS_AllocateElementPtrType  allocateElementFunct,
				OPS_AllocateMaterialPtrType allocateMaterialFunct,
				OPS_GetUniaxialMaterialPtrType OPS_GetUniaxialMaterialFunct,
				OPS_GetNDMaterialPtrType OPS_GetNDMaterialFunct)
{
	opserrPtr = theErrorStreamPtr;
	OPS_ErrorPtr = errorFunct;
	OPS_GetIntInputPtr = getIntInputFunct;
	OPS_GetDoubleInputPtr =getDoubleInputFunct;
	OPS_AllocateElementPtr = allocateElementFunct;
	OPS_AllocateMaterialPtr =allocateMaterialFunct;
	OPS_GetUniaxialMaterialPtr = OPS_GetUniaxialMaterialFunct;
	OPS_GetNDMaterialPtr = OPS_GetNDMaterialFunct;
}


int OPS_Error(char *data, int length)
{
  return (*OPS_ErrorPtr)(data, length);
}

extern "C" int  OPS_GetIntInput(int *numData, int*data)
{
  opserr << "int OPS_GetIntInput(int *numData, int*data)\n";
  return (*OPS_GetIntInputPtr)(numData, data);
}


extern "C" int OPS_GetDoubleInput(int *numData, double *data)
{
  opserr << "int        OPS_GetDoubleInput(int *numData, int*data)\n";
  return (*OPS_GetDoubleInputPtr)(numData, data);
}

extern "C" int OPS_AllocateMaterial(matObj *mat)
{
	return (*OPS_AllocateMaterialPtr)(mat);
}
extern UniaxialMaterial *OPS_GetUniaxialMaterial(int matTag)
{
	return (*OPS_GetUniaxialMaterialPtr)(matTag);
}
