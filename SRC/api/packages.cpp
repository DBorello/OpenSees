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
**                                                                    **
** ****************************************************************** */

/*                                                                        
** $Revision: 1.6 $
** $Date: 2009-01-10 17:45:55 $
** $Source: /usr/local/cvs/OpenSees/SRC/api/packages.cpp,v $
                                                                        
** Written: fmk 
*/                                                                        

#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>


#ifdef _WIN32

#include <windows.h>
#include <elementAPI.h>
#else
#include <dlfcn.h>
#endif

int 
getLibraryFunction(const char *libName, const char *funcName, void **libHandle, void **funcHandle) {

  int result = 0;
  
  *libHandle = NULL;
  *funcHandle = NULL;
  
#ifdef _WIN32
  
  //
  // first try and open dll
  //
  
  int libNameLength = strlen(libName);
  char *localLibName = new char[libNameLength+5];
  strcpy(localLibName, libName);
  strcpy(&localLibName[libNameLength], ".dll");
  
  HINSTANCE hLib = LoadLibrary(localLibName);
  
  delete [] localLibName;

  if (hLib != NULL) {
    char mod[124];
    GetModuleFileName((HMODULE)hLib, (LPTSTR)mod, 124);
    
    //
    // Now look for function with funcName
    //
    
    (*funcHandle) = (void *)GetProcAddress((HMODULE)hLib, funcName);
    
    
    if (*funcHandle == NULL) {
      FreeLibrary((HMODULE)hLib);
      return -2;
    }	

    //
    // we need to set the OpenSees pointer global variables if function there
    //

    typedef int (_cdecl *LocalInitPtrType)();
    typedef int (_cdecl *OPS_ErrorPtrType)(char *, int);
    typedef int (_cdecl *OPS_GetIntInputPtrType)(int *, int *);
    typedef int (_cdecl *OPS_GetDoubleInputPtrType)(int *, double *);
    typedef int (_cdecl *OPS_AllocateElementPtrType)(eleObj *, int *matTags, int *maType);
    typedef int (_cdecl *OPS_AllocateMaterialPtrType)(matObj *);
    typedef UniaxialMaterial *(*OPS_GetUniaxialMaterialPtrType)(int matTag);
    typedef NDMaterial * (*OPS_GetNDMaterialPtrType)(int matTag);
    typedef int (_cdecl *OPS_GetNodeInfoPtrType)(int *, int *, double *);
    
    //typedef void(_cdecl *setGlobalPointersFunction)(OPS_Stream *);
    typedef void (_cdecl *setGlobalPointersFunction)(OPS_Stream *,
						     OPS_ErrorPtrType,
						     OPS_GetIntInputPtrType,
						     OPS_GetDoubleInputPtrType,
						     OPS_AllocateElementPtrType,
						     OPS_AllocateMaterialPtrType,
						     OPS_GetUniaxialMaterialPtrType,
						     OPS_GetNDMaterialPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNodeInfoPtrType); 

    setGlobalPointersFunction funcPtr;
    
    // look for pointer function
    funcPtr = (setGlobalPointersFunction)GetProcAddress((HMODULE)hLib,"setGlobalPointers");
    if (!funcPtr) {
      FreeLibrary((HMODULE)hLib);
      return -2;
    }
    
    // invoke pointer function
    (funcPtr)(opserrPtr, OPS_Error, OPS_GetIntInput, OPS_GetDoubleInput,
	      OPS_AllocateElement, OPS_AllocateMaterial, OPS_GetUniaxialMaterial, 
	      OPS_GetNDMaterial, OPS_GetNodeCrd, OPS_GetNodeDisp, OPS_GetNodeVel,
	      OPS_GetNodeAcc);
   

   LocalInitPtrType initPtr;
	initPtr = (LocalInitPtrType)GetProcAddress((HMODULE)hLib,"localInit");
	if (!initPtr) {
      initPtr();
    }
    
  } else // no lib exists
    return -1;
  
  libHandle =  (void **)&hLib;

#else

  int libNameLength = strlen(libName);
  char *localLibName = new char[libNameLength+10];
  strcpy(localLibName, libName);

#ifdef _MACOSX
  strcpy(&localLibName[libNameLength], ".dylib");
#else
  strcpy(&localLibName[libNameLength], ".so");
#endif

  char *error;

  *libHandle = dlopen (localLibName, RTLD_NOW);
  
  if (*libHandle != NULL) {    

    void *funcPtr = dlsym(*libHandle, funcName);

    //    *funcHandle = dlsym(*libHandle, funcName);
    error = dlerror();
    if (funcPtr == NULL ) {

      char *underscoreFunctionName = new char[strlen(funcName)+2];
      strcpy(underscoreFunctionName, "_");
      strcpy(&underscoreFunctionName[1], funcName);    
      void *funcPtr = dlsym(*libHandle, underscoreFunctionName);
      delete [] underscoreFunctionName;

      if (funcPtr == NULL)  {
	result =  -1;
	dlclose(*libHandle);
	delete [] localLibName;
	return result;
      }
    } else
      ;

    *funcHandle = funcPtr;
    
    typedef int (*localInitPtrType)();
    localInitPtrType initFunct;
    funcPtr = dlsym(*libHandle, "localInit");
    if (funcPtr != NULL ) {
      initFunct = (localInitPtrType)funcPtr;
      initFunct();
    } 
  } else {
    result =  -1;
  }
  
  delete [] localLibName;

#endif

  return result;
}
