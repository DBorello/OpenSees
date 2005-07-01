
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
opserr << "FOUND LIB\n";
		char mod[124];
		GetModuleFileName((HMODULE)hLib, (LPTSTR)mod, 124);
	
		//
		// now look for function with funcName
		//

		(*funcHandle) = (void *)GetProcAddress((HMODULE)hLib, funcName);
	
		if (*funcHandle == NULL) {
			FreeLibrary((HMODULE)hLib);
			return -2;
		}	

		opserr << "FOUND FUNC\n";

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
