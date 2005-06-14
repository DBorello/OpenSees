
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>

#ifdef _WIN32

#else

#include <dlfcn.h>
#endif

int 
getLibraryFunction(const char *libName, const char *funcName, void **libHandle, void **funcHandle) {

  int result = 0;
  
  *libHandle = NULL;
  *funcHandle = NULL;
  
#ifdef _WIN32
  
#else
  
  int libNameLength = strlen(libName);
  char *localLibName = new char[libNameLength+4];
  strcpy(localLibName, libName);
  strcpy(&localLibName[libNameLength], ".so");
  
  char *error;
  *libHandle = dlopen (localLibName, RTLD_NOW);
  if (*libHandle != NULL) {    
    
    *funcHandle = dlsym(*libHandle, funcName);
    if ((error = dlerror()) != NULL) {
      opserr << *error;
      result = -2;
    } else
      result = -1;
  }

  delete [] localLibName;
  
#endif
  
  return result;
}
