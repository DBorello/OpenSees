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
                                                                        
// $Revision: 1.13 $
// $Date: 2005-12-20 23:37:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/FileDatastore.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 10/98
//
// Description: This file contains the class implementation for FileDatastore.
// FileDatastore is a concrete subclas of FE_Datastore. A FileDatastore 
// object is used in the program to store/restore the geometry and state 
// information in a domain at a particular instance in the analysis. The
// information is stored in text files.
//
// What: "@(#) FileDatastore.C, revA"

#include "FileDatastore.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <bool.h>
#include <iostream>
#include <iomanip>
using std::ios;
using std::setiosflags;


#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>


FileDatastore::FileDatastore(const char *dataBaseName,
			     Domain &theDomain, 
			     FEM_ObjectBroker &theObjBroker) 
  :FE_Datastore(theDomain, theObjBroker), 
   data(0), sizeData(0), currentMaxInt(0), currentMaxDouble(0)
{
  dataBase = new char [strlen(dataBaseName)+1];
  strcpy(dataBase, dataBaseName);

  this->resizeDouble(1024);
}

FileDatastore::~FileDatastore() 
{
  if (data != 0)
    delete [] data;

  //  while (theIDFilesIter != theIDFiles.end()) {
  //    theIDFilesIter++;
  // }
  for (theIDFilesIter = theIDFiles.begin(); theIDFilesIter != theIDFiles.end(); theIDFilesIter++) {
    FileDatastoreOutputFile *theFileStruct = theIDFilesIter->second;
    fstream *theFile = theFileStruct->theFile;
    theFile->seekp(0, ios::beg);
    theIntData->dbTag = theFileStruct->maxDbTag;
    theFile->write(data, sizeof(int));    
    theFile->close();
    delete theFile;
    delete theFileStruct;
  }
  theIDFiles.clear();


  for (theMatFilesIter = theMatFiles.begin(); theMatFilesIter != theMatFiles.end(); theMatFilesIter++) {
    FileDatastoreOutputFile *theFileStruct = theMatFilesIter->second;
    fstream *theFile = theFileStruct->theFile;
    theFile->seekp(0, ios::beg);
    theIntData->dbTag = theFileStruct->maxDbTag;
    theFile->write(data, sizeof(int));    
    theFile->close();
    delete theFile;
    delete theFileStruct;
  }
  theMatFiles.clear();


  for (theVectFilesIter = theVectFiles.begin(); theVectFilesIter != theVectFiles.end(); theVectFilesIter++) {
    FileDatastoreOutputFile *theFileStruct = theVectFilesIter->second;
    fstream *theFile = theFileStruct->theFile;
    theFile->seekp(0, ios::beg);
    theIntData->dbTag = theFileStruct->maxDbTag;
    theFile->write(data, sizeof(int));    
    theFile->close();
    delete theFile;
    delete theFileStruct;
  }
  theVectFiles.clear();
}


int 
FileDatastore::commitState(int commitTag)
{
  int result = FE_Datastore::commitState(commitTag);
  if (result == commitTag) 
    resetFilePointers();
    
  return result;
}


void
FileDatastore::resetFilePointers(void) {

  for (theIDFilesIter = theIDFiles.begin(); theIDFilesIter != theIDFiles.end(); theIDFilesIter++) {
    int tag = theIDFilesIter->first;
    FileDatastoreOutputFile *theFileStruct = theIDFilesIter->second;
    fstream *theFile = theFileStruct->theFile;
    if (theFile != 0) {
      theFile->seekp(0, ios::beg);
      theIntData->dbTag = theFileStruct->maxDbTag;
      theFile->write(data, sizeof(int));    
      if (theFile != 0)
	theFile->close();
      theFileStruct->theFile = 0;
    }
  }

  for (theMatFilesIter = theMatFiles.begin(); theMatFilesIter != theMatFiles.end(); theMatFilesIter++) {
    int tag = theMatFilesIter->first;
    FileDatastoreOutputFile *theFileStruct = theMatFilesIter->second;
    fstream *theFile = theFileStruct->theFile;
    if (theFile != 0) {
      theFile->seekp(0, ios::beg);
      theIntData->dbTag = theFileStruct->maxDbTag;
      theFile->write(data, sizeof(int));    
      if (theFile != 0)
	theFile->close();
      theFileStruct->theFile = 0;
    }
  }

  for (theVectFilesIter = theVectFiles.begin(); theVectFilesIter != theVectFiles.end(); theVectFilesIter++) {
    int tag = theVectFilesIter->first;
    FileDatastoreOutputFile *theFileStruct = theVectFilesIter->second;
    fstream *theFile = theFileStruct->theFile;
    if (theFile != 0) {
      theFile->seekp(0, ios::beg);
      theIntData->dbTag = theFileStruct->maxDbTag;
      theFile->write(data, sizeof(int));    
      if (theFile != 0)
	theFile->close();
      theFileStruct->theFile = 0;
    }
  }
  currentCommitTag = -1;
}


int 
FileDatastore::sendMsg(int dataTag, int commitTag, 
		       const Message &, 
		       ChannelAddress *theAddress)
{
  opserr << "FileDatastore::sendMsg() - not yet implemented\n";
  return -1;
}		       

int 
FileDatastore::recvMsg(int dataTag, int commitTag, 
		       Message &, 
		       ChannelAddress *theAddress)
{
  opserr << "FileDatastore::recvMsg() - not yet implemented\n";
  return -1;
}		       



int 
FileDatastore::sendID(int dataTag, int commitTag, 
		      const ID &theID, 
		      ChannelAddress *theAddress)
{
  if (currentCommitTag != commitTag)
    this->resetFilePointers();

  currentCommitTag = commitTag;
  
  FileDatastoreOutputFile *theFileStruct;
  fstream *theFile = 0;
  int maxDataTag = 0;
  int currentPos = 0;

  
  //
  // next we see if we already have this file; 
  //  if not we need to create data structure & open it
  //  if we have data structure, need to check file is opened (we close in a commit)
  //

  int idSize = theID.Size();
  int stepSize = (1 + idSize)*sizeof(int);

  theIDFilesIter = theIDFiles.find(idSize); 
  if (theIDFilesIter == theIDFiles.end()) {

    // we first check if we need to resize send buffer
    if (idSize > currentMaxInt) {
      if (this->resizeInt(idSize) < 0) {
	opserr << "FileDatastore::sendID() - failed in resizeInt()\n";
	return -1;
      }
    }

    char *fileName = new char[strlen(dataBase)+21];
    theFileStruct = new FileDatastoreOutputFile;

    if (fileName == 0 || theFileStruct == 0) {
      opserr << "FileDatastore::sendID() - out of memory\n";
      return -1;
    }

    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",idSize,commitTag);
    strcat(fileName,".IDs.");
    strcat(fileName,intName);
    
    if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
      opserr << "FileDatastore::sendID() - could not open file\n";
      return -1;
    } else
      theIDFiles.insert(MAP_FILES_TYPE(idSize, theFileStruct));
    
  } else {

    theFileStruct = theIDFilesIter->second;

    // make sure not close from a last commit
    if (theFileStruct->theFile == 0) {

      if (idSize > currentMaxInt) {
	if (this->resizeInt(idSize) < 0) {
	  opserr << "FileDatastore::sendID() - failed in resizeInt()\n";
	  return -1;
	}
      }
      char *fileName = new char[strlen(dataBase)+21];
      char intName[20];
      strcpy(fileName, dataBase);
      sprintf(intName,"%d.%d",idSize,commitTag);
      strcat(fileName,".IDs.");
      strcat(fileName,intName);
      
      if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
	opserr << "FileDatastore::sendID() - could not open file\n";
	return -1;
      }
    }
  }

  fstream *theStream = theFileStruct->theFile;
  long int fileEnd = theFileStruct->fileEnd;

  
  //
  // find location in file to place the data
  //
  
  bool found = false;  
  int pos;
  
  // we first check if the data can go at the end of the file
  // true if commitTag larger than any we have encountered so far
  if (theFileStruct->maxDbTag < dataTag)  {

    pos = fileEnd;
    found = true;
    theFileStruct->maxDbTag = dataTag;

    // we have to search from the beginning of the file

  } else { 
    pos = sizeof(int);
    theStream->seekg(pos, ios::beg);
    while ((pos < fileEnd) && (found == false)) {
      theStream->read(data, stepSize);
      if (theIntData->dbTag == dataTag) 
	found = true;
      else 
	pos += stepSize;
    }
  }


  //
  // we now place the data to be sent into our buffer
  //

  theIntData->dbTag = dataTag;
  for (int i=0; i<idSize; i++)
    theIntData->data[i] = theID(i);

  //
  // we now write the data
  //

  theStream->seekp(pos, ios::beg); // reset so can go write at the end
  theStream->write(data, stepSize);
  
  // update the size of file if we have added to eof
  if (theFileStruct->fileEnd <= pos)
    theFileStruct->fileEnd += stepSize;  

  return 0;
}		       

int 
FileDatastore::recvID(int dataTag, int commitTag, 
		      ID &theID, 
		      ChannelAddress *theAddress)    
{
  if (currentCommitTag != commitTag)
    this->resetFilePointers();

  currentCommitTag = commitTag;
  
  FileDatastoreOutputFile *theFileStruct;
  fstream *theFile = 0;
  long int loc = 0;
  int maxDataTag = 0;
  int currentPos = 0;
  
  //
  // next we see if we already have this file; 
  //  if not we need to create data structure & open it
  //  if we have data structure, need to check file is opened (we close in a commit)
  //

  int idSize = theID.Size();
  int stepSize = (1 + idSize)*sizeof(int);

  theIDFilesIter = theIDFiles.find(idSize); 
  if (theIDFilesIter == theIDFiles.end()) {

    // we first check if we need to resize recv buffer
    if (idSize > currentMaxInt) {
      if (this->resizeInt(idSize) < 0) {
	opserr << "FileDatastore::recvID() - failed in resizeInt()\n";
	return -1;
      }
    }

    char *fileName = new char[strlen(dataBase)+21];
    theFileStruct = new FileDatastoreOutputFile;

    if (fileName == 0 || theFileStruct == 0) {
      opserr << "FileDatastore::recvID() - out of memory\n";
      return -1;
    }

    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",idSize,commitTag);
    strcat(fileName,".IDs.");
    strcat(fileName,intName);
    
    if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
      opserr << "FileDatastore::recvID() - could not open file\n";
      return -1;
    } else
      theIDFiles.insert(MAP_FILES_TYPE(idSize, theFileStruct));
    
  } else {

    theFileStruct = theIDFilesIter->second;

    // make sure not close from a last commit
    if (theFileStruct->theFile == 0) {

      if (idSize > currentMaxInt) {
	if (this->resizeInt(idSize) < 0) {
	  opserr << "FileDatastore::recvID() - failed in resizeInt()\n";
	  return -1;
	}
      }
      char *fileName = new char[strlen(dataBase)+21];
      char intName[20];
      strcpy(fileName, dataBase);
      sprintf(intName,"%d.%d",idSize,commitTag);
      strcat(fileName,".IDs.");
      strcat(fileName,intName);
      
      if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
	opserr << "FileDatastore::recvID() - could not open file\n";
	return -1;
      }
    } 
  }

  fstream *theStream = theFileStruct->theFile;
  long int fileEnd = theFileStruct->fileEnd;
  long int pos = theStream->tellg();
  
  //
  // find location in file to place the data
  //
  
  bool found = false;  

  // we try the current file position first

  if (pos < fileEnd) {
    theStream->read(data, stepSize);
    if ((theIntData->dbTag == dataTag)) {
      found = true;
      pos += stepSize;
    } 
  }

  // we must search from the beginning of the file
  if (found == false) {
    pos = sizeof(int);
    theStream->seekg(pos, ios::beg); 
    while ((pos < fileEnd) && (found == false)) {
      theStream->read(data, stepSize);
      if (theIntData->dbTag == dataTag) 
	found = true;
      else 
	pos += stepSize;
    }
  }


  if (found == false) {
    opserr << "FileDatastore::recvID() - failed\n";
    return -1;
  }

  // we now place the received data into the ID 
  for (int i=0; i<idSize; i++)
    theID(i) = theIntData->data[i];

  return 0;
}		       




int 
FileDatastore::sendMatrix(int dataTag, int commitTag, 
		      const Matrix &theMatrix, 
		      ChannelAddress *theAddress)
{

  if (currentCommitTag != commitTag)
    this->resetFilePointers();

  currentCommitTag = commitTag;
  
  FileDatastoreOutputFile *theFileStruct;
  fstream *theFile = 0;
  int maxDataTag = 0;
  int currentPos = 0;
  
  //
  // next we see if we already have this file; 
  //  if not we need to create data structure & open it
  //  if we have data structure, need to check file is opened (we close in a commit)
  //

  // we first ensure that the Matrix is not too big
  int noMatCols= theMatrix.noCols();
  int noMatRows = theMatrix.noRows();
  int matSize = noMatRows * noMatCols;;
  int stepSize = sizeof(int) + matSize*sizeof(double);

  theMatFilesIter = theMatFiles.find(matSize); 
  if (theMatFilesIter == theMatFiles.end()) {

    // we first check if we need to resize send buffer
    if (matSize > currentMaxDouble) {
      if (this->resizeDouble(matSize) < 0) {
	opserr << "FileDatastore::sendMatrix() - failed in resizeInt()\n";
	return -1;
      }
    }

    char *fileName = new char[strlen(dataBase)+21];
    theFileStruct = new FileDatastoreOutputFile;

    if (fileName == 0 || theFileStruct == 0) {
      opserr << "FileDatastore::sendMatrix() - out of memory\n";
      return -1;
    }

    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",matSize,commitTag);
    strcat(fileName,".MATs.");
    strcat(fileName,intName);
    
    if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
      opserr << "FileDatastore::sendMatrix() - could not open file\n";
      return -1;
    } else
      theMatFiles.insert(MAP_FILES_TYPE(matSize, theFileStruct));
    
  } else {

    theFileStruct = theMatFilesIter->second;

    // make sure not close from a last commit
    if (theFileStruct->theFile == 0) {

      if (matSize > currentMaxDouble) {
	if (this->resizeDouble(matSize) < 0) {
	  opserr << "FileDatastore::sendMatrix() - failed in resizeInt()\n";
	  return -1;
	}
      }
      char *fileName = new char[strlen(dataBase)+21];
      char intName[20];
      strcpy(fileName, dataBase);
      sprintf(intName,"%d.%d",matSize,commitTag);
      strcat(fileName,".MATs.");
      strcat(fileName,intName);
      
      if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
	opserr << "FileDatastore::sendMatrix() - could not open file\n";
	return -1;
      }
    }
  }

  fstream *theStream = theFileStruct->theFile;
  long int fileEnd = theFileStruct->fileEnd;
  
  //
  // find location in file to place the data
  //
  
  bool found = false;  
  int pos;
  
  // we first check if the data can go at the end of the file
  // true if commitTag larger than any we have encountered so far
  if (theFileStruct->maxDbTag < dataTag)  {

    pos = fileEnd;
    found = true;
    theFileStruct->maxDbTag = dataTag;

    // we have to search from the beginning of the file

  } else { 
    pos = sizeof(int);
    theStream->seekg(pos, ios::beg);
    while ((pos < fileEnd) && (found == false)) {
      theStream->read(data, stepSize);
      if (theDoubleData->dbTag == dataTag) 
	found = true;
      else 
	pos += stepSize;
    }
  }


  //
  // we now place the data to be sent into our buffer
  //

  theDoubleData->dbTag = dataTag;
  int loc=0;
  for (int j=0; j<noMatCols; j++)
    for (int k=0; k < noMatRows; k++) {
      theDoubleData->data[loc] = theMatrix(k,j);
      loc++;
    }  

  //
  // we now write the data
  //

  theStream->seekp(pos, ios::beg); // reset so can go write at the end
  theStream->write(data, stepSize);
  
  // update the size of file if we have added to eof
  if (theFileStruct->fileEnd <= pos)
    theFileStruct->fileEnd += stepSize;  

  return 0;
}		       




int 
FileDatastore::recvMatrix(int dataTag, int commitTag, 
		      Matrix &theMatrix, 
		      ChannelAddress *theAddress)    
{

  if (currentCommitTag != commitTag)
    this->resetFilePointers();

  currentCommitTag = commitTag;

  FileDatastoreOutputFile *theFileStruct;
  fstream *theFile = 0;
  int maxDataTag = 0;
  int currentPos = 0;
  
  //
  // next we see if we already have this file; 
  //  if not we need to create data structure & open it
  //  if we have data structure, need to check file is opened (we close in a commit)
  //

  // we first check Matrix not too big
  int noMatCols= theMatrix.noCols();
  int noMatRows = theMatrix.noRows();
  int matSize = noMatRows * noMatCols;;
  int stepSize = sizeof(int) + matSize*sizeof(double);

  theMatFilesIter = theMatFiles.find(matSize); 
  if (theMatFilesIter == theMatFiles.end()) {

    // we first check if we need to resize recv buffer
    if (matSize > currentMaxDouble) {
      if (this->resizeDouble(matSize) < 0) {
	opserr << "FileDatastore::recvMatrix() - failed in resizeDouble()\n";
	return -1;
      }
    }

    char *fileName = new char[strlen(dataBase)+21];
    theFileStruct = new FileDatastoreOutputFile;

    if (fileName == 0 || theFileStruct == 0) {
      opserr << "FileDatastore::recvMatrix() - out of memory\n";
      return -1;
    }

    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",matSize,commitTag);
    strcat(fileName,".MATs.");
    strcat(fileName,intName);
    
    if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
      opserr << "FileDatastore::recvMatrix() - could not open file\n";
      return -1;
    } else
      theMatFiles.insert(MAP_FILES_TYPE(matSize, theFileStruct));
    
  } else {

    theFileStruct = theMatFilesIter->second;

    // make sure not close from a last commit
    if (theFileStruct->theFile == 0) {

      if (matSize > currentMaxDouble) {
	if (this->resizeDouble(matSize) < 0) {
	  opserr << "FileDatastore::recvMatrix() - failed in resizeInt()\n";
	  return -1;
	}
      }
      char *fileName = new char[strlen(dataBase)+21];
      char intName[20];
      strcpy(fileName, dataBase);
      sprintf(intName,"%d.%d",matSize,commitTag);
      strcat(fileName,".MATs.");
      strcat(fileName,intName);
      
      if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
	opserr << "FileDatastore::recvMatrix() - could not open file\n";
	return -1;
      }
    } 
  }

  fstream *theStream = theFileStruct->theFile;
  long int fileEnd = theFileStruct->fileEnd;
  long int pos = theStream->tellg();
  
  //
  // find location in file to place the data
  //
  
  bool found = false;  

  // we try the current file position first

  if (pos < fileEnd) {
    theStream->read(data, stepSize);
    if ((theDoubleData->dbTag == dataTag)) {
      found = true;
      pos += stepSize;
    } 
  }

  // we must search from the beginning of the file
  if (found == false) {
    pos = sizeof(int);
    theStream->seekg(pos, ios::beg); 
    while ((pos < fileEnd) && (found == false)) {
      theStream->read(data, stepSize);
      if (theDoubleData->dbTag == dataTag) 
	found = true;
      else 
	pos += stepSize;
    }
  }

  if (found == false) {
    opserr << "FileDatastore::recvMatrix() - failed\n";
    return -1;
  }

  int loc=0;
  for (int j=0; j<noMatCols; j++)
    for (int k=0; k < noMatRows; k++) {
      theMatrix(k,j) = theDoubleData->data[loc];
      loc++;
    }

  return 0;
}		       



int 
FileDatastore::sendVector(int dataTag, int commitTag, 
		      const Vector &theVector, 
		      ChannelAddress *theAddress)
{

  if (currentCommitTag != commitTag)
    this->resetFilePointers();

  currentCommitTag = commitTag;
  
  FileDatastoreOutputFile *theFileStruct;
  fstream *theFile = 0;
  int maxDataTag = 0;
  int currentPos = 0;
  
  //
  // next we see if we already have this file; 
  //  if not we need to create data structure & open it
  //  if we have data structure, need to check file is opened (we close in a commit)
  //

  // we first ensure that the Matrix is not too big
  int vectSize = theVector.Size();
  int stepSize = sizeof(int) + vectSize*sizeof(double);

  theVectFilesIter = theVectFiles.find(vectSize); 
  if (theVectFilesIter == theVectFiles.end()) {

    // we first check if we need to resize send buffer
    if (vectSize > currentMaxDouble) {
      if (this->resizeDouble(vectSize) < 0) {
	opserr << "FileDatastore::sendVector() - failed in resizeInt()\n";
	return -1;
      }
    }

    char *fileName = new char[strlen(dataBase)+21];
    theFileStruct = new FileDatastoreOutputFile;

    if (fileName == 0 || theFileStruct == 0) {
      opserr << "FileDatastore::sendVector() - out of memory\n";
      return -1;
    }

    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",vectSize,commitTag);
    strcat(fileName,".VECs.");
    strcat(fileName,intName);
    
    if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
      opserr << "FileDatastore::sendVector() - could not open file\n";
      return -1;
    } else
      theVectFiles.insert(MAP_FILES_TYPE(vectSize, theFileStruct));
    
  } else {

    theFileStruct = theVectFilesIter->second;

    // make sure not close from a last commit
    if (theFileStruct->theFile == 0) {

      if (vectSize > currentMaxDouble) {
	if (this->resizeDouble(vectSize) < 0) {
	  opserr << "FileDatastore::sendVector() - failed in resizeInt()\n";
	  return -1;
	}
      }
      char *fileName = new char[strlen(dataBase)+21];
      char intName[20];
      strcpy(fileName, dataBase);
      sprintf(intName,"%d.%d",vectSize,commitTag);
      strcat(fileName,".VECs.");
      strcat(fileName,intName);
      
      if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
	opserr << "FileDatastore::sendVector() - could not open file\n";
	return -1;
      }
    }
  }

  fstream *theStream = theFileStruct->theFile;
  long int fileEnd = theFileStruct->fileEnd;
  
  //
  // find location in file to place the data
  //
  
  bool found = false;  
  int pos;
  
  // we first check if the data can go at the end of the file
  // true if commitTag larger than any we have encountered so far
  if (theFileStruct->maxDbTag < dataTag)  {

    pos = fileEnd;
    found = true;
    theFileStruct->maxDbTag = dataTag;

    // we have to search from the beginning of the file

  } else { 
    pos = sizeof(int);
    theStream->seekg(pos, ios::beg);
    while ((pos < fileEnd) && (found == false)) {
      theStream->read(data, stepSize);
      if (theDoubleData->dbTag == dataTag) 
	found = true;
      else 
	pos += stepSize;
    }
  }


  //
  // we now place the data to be sent into our buffer
  //

  theDoubleData->dbTag = dataTag;
  for (int i=0; i<vectSize; i++)
    theDoubleData->data[i] = theVector(i);

  //
  // we now write the data
  //

  theStream->seekp(pos, ios::beg); // reset so can go write at the end
  theStream->write(data, stepSize);
  
  // update the size of file if we have added to eof
  if (theFileStruct->fileEnd <= pos)
    theFileStruct->fileEnd += stepSize;  

  return 0;
}		       

int 
FileDatastore::recvVector(int dataTag, int commitTag, 
		      Vector &theVector, 
		      ChannelAddress *theAddress)    
{

  if (currentCommitTag != commitTag)
    this->resetFilePointers();

  currentCommitTag = commitTag;

  FileDatastoreOutputFile *theFileStruct;
  fstream *theFile = 0;
  int maxDataTag = 0;
  int currentPos = 0;
  
  //
  // next we see if we already have this file; 
  //  if not we need to create data structure & open it
  //  if we have data structure, need to check file is opened (we close in a commit)
  //

  // we first check Vector not too big
  int vectSize = theVector.Size();
  int stepSize = sizeof(int) + vectSize*sizeof(double);

  theVectFilesIter = theVectFiles.find(vectSize); 
  if (theVectFilesIter == theVectFiles.end()) {

    // we first check if we need to resize recv buffer
    if (vectSize > currentMaxDouble) {
      if (this->resizeDouble(vectSize) < 0) {
	opserr << "FileDatastore::recvVectrix() - failed in resizeDouble()\n";
	return -1;
      }
    }

    char *fileName = new char[strlen(dataBase)+21];
    theFileStruct = new FileDatastoreOutputFile;

    if (fileName == 0 || theFileStruct == 0) {
      opserr << "FileDatastore::recvVectrix() - out of memory\n";
      return -1;
    }

    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",vectSize,commitTag);
    strcat(fileName,".VECs.");
    strcat(fileName,intName);
    
    if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
      opserr << "FileDatastore::recvVectrix() - could not open file\n";
      return -1;
    } else
      theVectFiles.insert(MAP_FILES_TYPE(vectSize, theFileStruct));
    
  } else {

    theFileStruct = theVectFilesIter->second;

    // make sure not close from a last commit
    if (theFileStruct->theFile == 0) {

      if (vectSize > currentMaxDouble) {
	if (this->resizeDouble(vectSize) < 0) {
	  opserr << "FileDatastore::recvVectrix() - failed in resizeInt()\n";
	  return -1;
	}
      }
      char *fileName = new char[strlen(dataBase)+21];
      char intName[20];
      strcpy(fileName, dataBase);
      sprintf(intName,"%d.%d",vectSize,commitTag);
      strcat(fileName,".VECs.");
      strcat(fileName,intName);
      
      if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
	opserr << "FileDatastore::recvVectrix() - could not open file\n";
	return -1;
      }
    } 
  }

  fstream *theStream = theFileStruct->theFile;
  long int fileEnd = theFileStruct->fileEnd;
  long int pos = theStream->tellg();
  
  //
  // find location in file to place the data
  //
  
  bool found = false;  

  // we try the current file position first

  if (pos < fileEnd) {
    theStream->read(data, stepSize);
    if ((theDoubleData->dbTag == dataTag)) {
      found = true;
      pos += stepSize;
    } 
  }

  // we must search from the beginning of the file
  if (found == false) {
    pos = sizeof(int);
    theStream->seekg(pos, ios::beg); 
    while ((pos < fileEnd) && (found == false)) {
      theStream->read(data, stepSize);
      if (theDoubleData->dbTag == dataTag) 
	found = true;
      else 
	pos += stepSize;
    }
  }

  if (found == false) {
    opserr << "FileDatastore::recvVector() - failed\n";
    return -1;
  }

  for (int i=0; i<vectSize; i++)
    theVector(i) = theDoubleData->data[i];

  return 0;
}		       


#include <fstream>
using std::ofstream;
using std::cerr;


int 
FileDatastore::createTable(const char *tableName, int numColumns, char *columns[])
{
  // open the file
  int res = 0;
  char *fileName = new char[strlen(tableName) + strlen(dataBase) + 1];
  if (fileName == 0) {
    opserr << "FileDatastore::insertData - out of memory; failed to open file: " << fileName << endln;
    return -1;
  }

  strcpy(fileName, dataBase);    
  strcat(fileName,".");
  strcat(fileName, tableName);

  ofstream table;
  table.open(fileName, ios::out | ios::trunc); 

  if (table.is_open()) {
    // write the data
    for (int i=0; i<numColumns; i++) {
      table << columns[i] << "\t";
    }
    table << "\n";
    table.close();
    
  } else {
    opserr << "FileDatastore::insertData - failed to open file: " << fileName << endln;
    res = -1;
  }
  
  delete [] fileName;
  return res;
}

int 
FileDatastore::insertData(const char *tableName, char *columns[], 
			  int commitTag, const Vector &data)
{
  // open the file
  char *fileName = new char[strlen(tableName) + strlen(dataBase) + 1];
  if (fileName == 0) {
    opserr << "FileDatastore::insertData - out of memory; failed to open file: " << fileName << endln;
    return -1;
  }

  strcpy(fileName, dataBase);    
  strcat(fileName,".");
  strcat(fileName, tableName);

  ofstream table;
  table.open(fileName, ios::app); 

  table << setiosflags(ios::scientific);
  table << std::setprecision(16);

  if (table.is_open()) {
    // write the data
    for (int i=0; i<data.Size(); i++) {
      table << data(i) << "\t";
    }
    
    table << "\n";
    table.close();

  } else {
    opserr << "FileDatastore::insertData - failed to open file: " << fileName << endln;
    return -1;
  }

  delete [] fileName;
  return 0;
}


int 
FileDatastore::getData(const char *tableName, char *columns[], int commitTag, Vector &data)
{
  return 0;
}



/*******************************************************************
 *              MISC METHODS & FUNCTONS FOR OPENING THE FILE       *
 *******************************************************************/

int
FileDatastore::openFile(char *fileName, FileDatastoreOutputFile *theFileStruct, int dataSize)
{
  fstream *res = new fstream();
  if (res == 0) {
    opserr << "FileDatastore::openFile - out of memory; failed to open file: " << fileName << endln;
    return 0;
  }

  res->open(fileName, ios::in | ios::out | ios::binary); 

  // if file did not exist, need to pass trunc flag to open it
  if (res->bad() == true || res->is_open() == false) {
    // delete & new again for unix gcc compiler to work!
    delete res;
    res = new fstream();
    if (res == 0) {
      opserr << "FileDatastore::openFile - out of memory; failed to open file: " << fileName << endln;
      theFileStruct->theFile = res;
      return -1;
    }
    res->open(fileName, ios::in | ios::out | ios::trunc | ios::binary);   
  }

  if (res->bad() == true || res->is_open() == false) {
    opserr << "FATAL - FileDatastore::openFile() - could not open file " << fileName << endln;
    delete res;
    theFileStruct->theFile = 0;
    return -1;
  }

  // set the position for writing to eof
  res->seekp(0,ios::end);  
  long int fileEnd = res->tellp();
  int maxDataTag = 0;
  
  if (fileEnd == -1 || fileEnd == 0) {
    theIntData->dbTag = maxDataTag;
    res->write(data, sizeof(int));    
    fileEnd = sizeof(int);
    maxDataTag = -1;
  } else {
    long int pos = 0;
    int stepSize = dataSize;
    res->seekg(0, ios::beg);  
    res->read(data, sizeof(int));
    maxDataTag = theIntData->dbTag;
  }

  // fill in the structure data
  theFileStruct->theFile = res;
  theFileStruct->fileEnd = fileEnd;

  theFileStruct->maxDbTag = maxDataTag; 	        

  return 0;
}

int 
FileDatastore::resizeInt(int newSize) {
  int sizeOfChar = sizeof(char);
  int sizeOfInt = sizeof(int);
  int sizeOfDouble = sizeof(double);  
  newSize = (newSize+1)*sizeOfInt/sizeOfChar;

  if (newSize < sizeData)
    return 0; // already big enough

  if (newSize <= 0) {
    opserr << "FileDatastore::resizeInt(int newSize) - invalidSize " << newSize << endln;
    return -1; // invalid size
  }

  if (data != 0)
    delete [] data;

  data = new char[newSize];
  if (data == 0) {
    opserr << "FileDatastore::resizeInt(int newSize) - out of memory for size: " << newSize << endln;
    return -1;
  }

  sizeData = newSize;

  currentMaxInt = (sizeOfChar*sizeData-sizeOfInt)/sizeOfInt;
  currentMaxDouble = (sizeOfChar*sizeData-sizeOfInt)/sizeOfDouble;  
  theIntData = (IntData *)data;
  theDoubleData = (DoubleData *)data;
  return 0;
}

int
FileDatastore::resizeDouble(int newSize) {
  int sizeOfChar = sizeof(char);
  int sizeOfInt = sizeof(int);
  int sizeOfDouble = sizeof(double);  
  newSize = (newSize*sizeOfDouble+sizeOfInt)/sizeOfChar;

  if (newSize < sizeData)
    return 0; // already big enough

  if (newSize <= 0) {
    opserr << "FileDatastore::resizeInt(int newSize) - invalidSize " << newSize << endln;
    return -1; // invalid size
  }

  if (data != 0)
    delete [] data;

  data = new char[newSize];
  if (data == 0) {
    opserr << "FileDatastore::resizeInt(int newSize) - out of memory for size: " << newSize << endln;
    return -1;
  }

  sizeData = newSize;

  currentMaxInt = (sizeOfChar*sizeData-sizeOfInt)/sizeOfInt;
  currentMaxDouble = (sizeOfChar*sizeData-sizeOfInt)/sizeOfDouble;  
  theIntData = (IntData *)data;
  theDoubleData = (DoubleData *)data;
  return 0;
}
