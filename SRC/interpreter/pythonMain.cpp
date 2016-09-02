
#include "PythonInterpreter.h"

int main(int argc, char **argv) {
  PythonInterpreter* theInterpreter = new PythonInterpreter(argc, argv);
  int res = theInterpreter->run();
  
  delete theInterpreter;
  theInterpreter = 0;
  return res;
}

