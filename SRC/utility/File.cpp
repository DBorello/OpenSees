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
                                                                        
// $Revision: 1.1 $
// $Date: 2007-09-29 01:59:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/File.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 09/07
//
// Description: This file contains the class definition for File used in SimulationINformation.
//
// What: "@(#) SimulationInformation.h, revA"


#include <map>
#include <string>
using namespace std;

#include <File.h>
#include <FileIter.h>

File::File(const char *theName, const char *theDescription, bool isDir)
  :isDirectory(isDir), name(theName), description(theDescription)
{
  
}

File::~File()
{

}

int 
File::addFile(File *theFile)
{
  if (isDirectory == false)
    return -1;
  
  if (dirFiles.find(theFile->name) == dirFiles.end())
    dirFiles[theFile->name] = theFile;

  return 0;
}

const char *
File::getName(void)
{
  return name.c_str();
}

const char *
File::getDescription(void)
{
  return description.c_str();
}

bool 
File::isDir(void)
{
  return isDirectory;
}

FileIter 
File::getFiles(void)
{
  return FileIter(*this);
}

