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
// $Date: 2004-11-13 00:54:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputFileHandler.cpp,v $
                                                                        
// Written: fmk 
// Date: 10/04
//
// Description: This file contains the class implementation for
// DataOutputFileHandler. 
//
// What: "@(#) DataOutputFileHandler.C, revA"

#include "DataOutputFileHandler.h"
#include <Vector.h>

DataOutputFileHandler::DataOutputFileHandler(const char *theFileName, 
					     echoMode theEMode, 
					     openMode theOMode)
  :fileName(0), theEchoMode(theEMode), theOpenMode(theOMode), numColumns(-1)
{
  if (theFileName != 0) {
    fileName = new char [strlen(theFileName)+1];
    if (fileName == 0) {
    opserr << "DataOutputFileHandler::DataOutputFileHandler() - out of memory\n";      
    } else 
      strcpy(fileName, theFileName);
  }

  if (fileName != 0 && outputFile.setFile(fileName, theOpenMode) < 0) {
    opserr << "DataOutputFileHandler::DataOutputFileHandler() - setFile() failed\n";
    if (fileName != 0) {
      delete [] fileName;
      fileName = 0;
    }
  }
}

DataOutputFileHandler::~DataOutputFileHandler()
{
  if (fileName != 0)
    delete [] fileName;
}

int 
DataOutputFileHandler::open(char **dataDescription, int numData)
{
  if (fileName == 0) {
    opserr << "DataOutputFileHandler::open() - no filename passed in constructor\n";
    return -1;
  }

  if (dataDescription == 0) {
    opserr << "DataOutputFileHandler::open() - no column description data passed\n";
    return -1;
  }

  if (numData < 0) {
    opserr << "DataOutputFileHandler::open() - numColumns (" << numData << ") < 0\n";
    return -1;
  } else
    numColumns = numData;


  if (theEchoMode == DATA_FILE) {
    for (int i=0; i<numData; i++)
      outputFile << dataDescription[i] << " ";
    outputFile << endln;

  } else if (theEchoMode == XML_FILE) {
    
    // create a copy of the  file name
    int res = 0;
    const char *name = outputFile.getFileName();
    int fileNameLength = strlen(name);
    char *xmlFileName = new char[fileNameLength + 5];
    
    if (xmlFileName == 0) {
      opserr << "DataOutputFileHandler::open - out of memory creating copy of string " << xmlFileName << endln;
      return -1;
  }
    
    strcpy(xmlFileName, name);
    if (fileNameLength > 4 && strcmp(".out", &name[fileNameLength-4]) == 0) {
      xmlFileName[fileNameLength-4] = '\0';
    }
    
    strcat(xmlFileName,".xml");
    
    FileStream xmlFile;
    if (xmlFile.setFile(xmlFileName, OVERWRITE) == 0) {
      
      // write the xml data
      xmlFile << "<?xml version=\"1.0\"?>\n";
      xmlFile << "<NumericalFileDataDescription>\n";
      xmlFile << "\t<DataFile>\n";
      xmlFile << "\t\t<DataFileName> " << name << "</DataFileName>\n";
      xmlFile << "\t\t<NumberDataColumns> " << numData << "</NumberDataColumns>\n";
      xmlFile << "\t</DataFile>\n";
      for (int i=0; i<numData; i++) {
	xmlFile << "\t<DataColumnDescription>\n";
	xmlFile << "\t\t<ColumnLocation> " << i+1 << "</ColumnLocation>\n";      
	xmlFile << "\t\t<Description> " << dataDescription[i] << "</Description>\n";
	xmlFile << "\t</DataColumnDescription>\n";
      }
      xmlFile << "</NumericalFileDataDescription>\n";
      xmlFile.close();
      
    } else {
      opserr << "DataOutputFileHandler::open - failed to open cml file: " << xmlFileName << endln;
      delete [] xmlFileName;
      res = -1;
    }
    
    // no longer need xmlFileName
    delete [] xmlFileName;
  }

  return 0;
}

int 
DataOutputFileHandler::write(Vector &data) 
{
  if (fileName == 0 || numColumns < 0) {
    opserr << "DataOutputFileHandler::write() - no filename or data description has been set\n";
    return -1;
  }

  if (data.Size() == numColumns)
    outputFile << data;
  else {
    opserr << fileName;
    opserr << "DataOutputStreamHandler::write() - Vector not of correct size\n";
    return -1;
  }

  return 0;
}
 
