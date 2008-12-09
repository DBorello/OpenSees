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
** $Revision: 1.2 $
** $Date: 2008-12-09 00:09:40 $
** $Source: /usr/local/cvs/OpenSees/SRC/api/packages.cpp,v $
                                                                        
** Written: fmk 
*/                                                                        

#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>


#ifdef _WIN32

#include <windows.h>
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
    
    typedef void(_cdecl *setGlobalPointersFunction)(OPS_Stream *);
    setGlobalPointersFunction funcPtr;
    
    // look for pointer function
    funcPtr = (setGlobalPointersFunction)GetProcAddress((HMODULE)hLib,"setGlobalPointers");
    if (!funcPtr) {
      FreeLibrary((HMODULE)hLib);
      return -2;
    }
    
    // invoke pointer function
    (funcPtr)(opserrPtr);
    
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

  opserr << "packages - START " << localLibName << " " << funcName << endln;

  *libHandle = dlopen (localLibName, RTLD_NOW);

  if (*libHandle != NULL) {    

    opserr << "packages - FOUND LIB " << localLibName << " " << funcName << endln;
    
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
	opserr << "packages - NO FUNCTION " << localLibName << " " << funcName << endln;
	result =  -1;
	dlclose(*libHandle);
      }
    } else
      opserr << "packages - FOUND FUNCTION " << localLibName << " " << funcName << endln;

    *funcHandle = funcPtr;
  } else {
    opserr << "packages - NO LIB " << localLibName << " " << funcName << " result: " << result << endln;
    result =  -1;
  }
  
  delete [] localLibName;
  

#endif
  opserr << "packages - " << libName << " " << funcName << " result: " << result << endln;

  return result;
}
