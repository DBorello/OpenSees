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
                                                                        
// $Revision: 1.11 $
// $Date: 2007-10-02 20:53:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/SimulationInformation.cpp,v $
//
// Description: This file contains the class definition for SimulationInformation.
// SimulationInformation is a stopwatch.
//
// What: "@(#) SimulationInformation.h, revA"



#include<SimulationInformation.h>
#include <FileIter.h>
#include <OPS_Globals.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <tcl.h>

static int numSimulationInformation = 0;
static SimulationInformation *theLastSimulationInformation = 0;

SimulationInformation::SimulationInformation() 
  :title(0), description(0), contactName(0),
   lengthUnit(0), forceUnit(0), timeUnit(0)
{
  strcpy(startTime," ");
  strcpy(endTime," ");

  if (numSimulationInformation == 0) {
    numSimulationInformation++;
  }

  theLastSimulationInformation = this;

  theFiles = new File("", "", true);

  this->start();
}

#include <XmlFileStream.h>

SimulationInformation::~SimulationInformation()
{ 
  if (strcmp(endTime," ") == 0)
    this->end();

  XmlFileStream b("test.out");
  b << *this;

  if (title != 0)
    delete [] title;
  if (description != 0)
    delete [] description;
  if (contactName != 0)
    delete [] contactName;
  if ( lengthUnit != 0)
    delete [] lengthUnit;
  if ( forceUnit!= 0)
    delete [] forceUnit;
  if ( timeUnit != 0)
    delete [] timeUnit;

  if (theLastSimulationInformation == this)
    theLastSimulationInformation = 0;

  delete theFiles;
}


int
SimulationInformation::start(void)
{
  paramNames.clear();
  paramValues.clear();
  analysisTypes.clear();
  modelTypes.clear();
  elementTypes.clear();
  materialTypes.clear();

  // now set the start time
  time_t timeT;
  if (time(&timeT) != 0) {
#ifdef _WIN32
    const char *sTime = ctime(&timeT);
    strcpy(startTime, sTime);
#else
    ctime_r(&timeT, &startTime[0]);
#endif

  }

  strcpy(endTime," ");
  numInputFiles = 0;
  
  return 0;
}

int
SimulationInformation::end(void)
{
  time_t timeT;
  if (time(&timeT) != 0) {
#ifdef _WIN32
    const char *eTime = ctime(&timeT);
    strcpy(endTime, eTime);
#else
    ctime_r(&timeT, &endTime[0]);
#endif
  }
  
  return 0;
}

int 
SimulationInformation::setTitle(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (title != 0)
    delete [] title;

  // create new space & copy string
  title = new char[strlen(name)+1];
  if (title == 0)
    return -1;
  strcpy(title, name);

  return 0;
}

int 
SimulationInformation::setDescription(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (description != 0)
    delete [] description;

  // create new space & copy string
  description = new char[strlen(name)+1];
  if (description == 0)
    return -1;
  strcpy(description, name);

  return 0;
}

int 
SimulationInformation::setContact(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (contactName != 0)
    delete [] contactName;

  // create new space & copy string
  contactName = new char[strlen(name)+1];
  if (contactName == 0)
    return -1;
  strcpy(contactName, name);

  return 0;
}

int 
SimulationInformation::setLengthUnit(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (lengthUnit != 0)
    delete [] lengthUnit;

  // create new space & copy string
  lengthUnit = new char[strlen(name)+1];
  if (lengthUnit == 0)
    return -1;
  strcpy(lengthUnit, name);

  return 0;
}

int 
SimulationInformation::setForceUnit(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (forceUnit != 0)
    delete [] forceUnit;

  // create new space & copy string
  forceUnit = new char[strlen(name)+1];
  if (forceUnit == 0)
    return -1;
  strcpy(forceUnit, name);

  return 0;
}

int 
SimulationInformation::setTimeUnit(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (timeUnit != 0)
    delete [] timeUnit;

  // create new space & copy string
  timeUnit = new char[strlen(name)+1];
  if (timeUnit == 0)
    return -1;
  strcpy(timeUnit, name);

  return 0;
}


void 
PrintFiles(OPS_Stream &s, File *theFile) 
{
  if (theFile == 0)
    return;
  
  const char *fileName = theFile->getName();

  if (theFile->isDir() == true) {
    s << fileName << endln;
    FileIter theDirFiles = theFile->getFiles();
    File *theDirFile;
    while ((theDirFile = theDirFiles()) != 0)
      PrintFiles(s, theDirFile);
  } else
    s << "  " << fileName << endln;

}

int 
SimulationInformation::addInputFile(const char *fileName, const char *path)
{
  if (numInputFiles == 0)
    theFiles->addFile(fileName, path, "Main Input File");
  else
    theFiles->addFile(fileName, path, "Input File");


  //  PrintFiles(opserr, theFiles);
  numInputFiles++;

  return 0;
}



int 
SimulationInformation::addOutputFile(const char *fileName, const char *path)
{
  theFiles->addFile(fileName, path, "Output File");
  //  PrintFiles(opserr, theFiles);
  return 0;
}




int 
SimulationInformation::addParameter(const char *name, const char *value)
{
  // check for valid name & value
  if (name == 0 || value == 0) 
    return -1;
  
  // create new array to hold pointers and copy pointers there
  paramNames.addString(name);
  paramValues.addString(value);
  
  return 0;
}

int 
SimulationInformation::addModelType(const char *theType)
{
  return modelTypes.addString(theType);
}

int 
SimulationInformation::addAnalysisType(const char *theType)
{
  return analysisTypes.addString(theType);
}

int 
SimulationInformation::addLoadingType(const char *theType)
{
  return loadingTypes.addString(theType);
}

int 
SimulationInformation::addElementType(const char *theType)
{
  return elementTypes.addString(theType);
}

int 
SimulationInformation::addMaterialType(const char *theType)
{
  return materialTypes.addString(theType);
}


// TclSimulationInformation_defaultUnits()
// to define basic units. the following is based on code provided by S. Mazzoni




void 
SimulationInformation::Print(OPS_Stream &s) const
{
  //  char version[10];
  //  strcpy(version,OPS_VERSION);

  s.tag("Central");
  s.tag("SimulationRun");

  if (title != 0)
    s.tag("title",title);

  if (description != 0)
    s.tag("description",description);

  if (contactName != 0)
    s.tag("contact",contactName);

  if (lengthUnit != 0)
    s.tag("lengthUnit", lengthUnit);

  if (forceUnit != 0)
    s.tag("lengthUnit", forceUnit);

  if (timeUnit != 0)
    s.tag("timeUnit", forceUnit);


  // need anotherTime to get rid of /n
  char *c = (char *)strchr(startTime,'\n');
  if (c != 0)
    strcpy(c,"");
  s.tag("startDate",startTime);

  c = (char *)strchr(endTime,'\n');
  if (c != 0)
    strcpy(c,"");
  s.tag("endDate",endTime);

  int numStrings;

  numStrings = modelTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationModelType",modelTypes.getString(i));    

  numStrings = analysisTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationAnalysisType",analysisTypes.getString(i));    

  numStrings = loadingTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationLoadingType",loadingTypes.getString(i));    

  numStrings = elementTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationElementType",elementTypes.getString(i));    

  numStrings = materialTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationMaterialType",materialTypes.getString(i));    

  s.tag("Software");
  s.tag("program","OpenSees");
  //  s.tag("version",version);
  s.endTag(); // Software


  s.tag("ComputerResource");
  s.tag("OS","Linux");
  s.tag("machine","local");
  s.endTag(); // Computer

  FileIter theFileIter = theFiles->getFiles();
  File *theFile;
  while ((theFile = theFileIter()) != 0) {
    s.tag("DataFile");
    s.attr("name",theFile->getName());
    s.attr("description",theFile->getDescription());
    s.endTag();
  }

  s.endTag(); // SimulationRun
  s.endTag(); // Central

}    

#ifdef _HTTTPS

extern int neesLogin(const char *user,
		     const char *pass,
		     char **cookieRes);


extern int neesSENDTrial(const char *cookie,
			 int projID,
			 int expID,
			 const char *name,
			 const char *title,
			 const char *objective,
			 const char *description);

extern int neesADD_TrialAnalysisFile(const char *cookie,
				     int projID,
				     int expID,
				     int trialID,
				     const char *path,
				     const char *name,
				     const char *description);

extern int neesADD_TrialAnalysisDir(const char *cookie,
				    int projID,
				    int expID,
				    int trialID,
				    const char *path,
				    const char *dirName);

#endif

int 
SimulationInformation::neesUpload(const char *username, 
				  const char *passwd, 
				  int projID, 
				  int expID)
{
#ifdef _HTTTPS

  // call endTime if not already called
  if (strcmp(endTime, " ") == 0)
    this->end();

  int res;
  int trialID;
  char *cookie = 0;

  opserr << "SimulationInformation::neesUpload\n";
  opserr << *this;
  return 0;

  res = neesLogin(username, passwd, &cookie);
  if (res != 0) {
    if (cookie != 0)
      free(cookie);

    return -1;
  }

  trialID = neesSENDTrial(cookie, projID, expID, title, title, " ", description);

  if (trialID < 0) {
    if (cookie != 0)
      opserr << "SimulationInformation::neesUpload() - failed to send Trial information\n";
      free(cookie);

    return -2;
  }

  FileIter theFileIter = theFiles->getFiles();
  File *theFile;
  while ((res == 0) && ((theFile = theFileIter()) != 0)) {

    res =  neesADD_TrialAnalysisFile(cookie,
				     projID,
				     expID,
				     trialID,
				     "/",
				     theFile->getName(),
				     theFile->getDescription());    
    if (res < 0) {
      opserr << "SimulationInformation::neesUpload() - failed to send file: " << theFile->getName() << endln;
    }
  }

  
  free(cookie);
  return res;
#else
  
  opserr << "ERROR: SimulationInformation::neesUpload() - not available in this build\n";
  return -1;

#endif
}


OPS_Stream &operator<<(OPS_Stream &s, const SimulationInformation &E)
{
  E.Print(s);
  return s;
}


