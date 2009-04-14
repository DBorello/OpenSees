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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-04-14 21:12:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataFileStream.cpp,v $


#include <DataFileStream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <ID.h>
#include <Channel.h>
#include <Message.h>
#include <Matrix.h>

using std::cerr;
using std::ios;
using std::setiosflags;
using std::ifstream;
using std::string;
using std::getline;

DataFileStream::DataFileStream(int indent)
  :OPS_Stream(OPS_STREAM_TAGS_DataFileStream), 
   fileOpen(0), fileName(0), indentSize(indent), sendSelfCount(0), theChannels(0)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+5];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");
}


DataFileStream::DataFileStream(const char *file, openMode mode, int indent)
  :OPS_Stream(OPS_STREAM_TAGS_DataFileStream), 
   fileOpen(0), fileName(0), indentSize(indent), sendSelfCount(0), theChannels(0)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");

  this->setFile(file, mode);
}


DataFileStream::~DataFileStream()
{
  if (fileOpen == 1)
    theFile.close();
  
#ifdef _NO_PARALLEL_FILESYSTEM

  if (sendSelfCount < 0) {

    ifstream theFile0;

    theFile0.open(fileName, ios::in);
    char c=theFile0.peek();
    int ordering = -1;
    if (c == 'O') {
      ordering = 0;
    }

    if (ordering == 0) {
      char Ordering[11];
      theFile0.get(Ordering, 11);
      int numColumns =0;
      theFile0 >> numColumns;

      static ID numColumnID(1);
      numColumnID(0) = numColumns;
      theChannels[0]->sendID(0,0,numColumnID);
      
      if (numColumns != 0) {
	
	ID *theColumns = new ID(numColumns);
	double *theData = new double [numColumns];
	Vector data(theData, numColumns);
	
	for (int j=0; j<numColumns; j++)
	  theFile0 >> (*theColumns)[j];

	theChannels[0]->sendID(0, 0, *theColumns);
	delete theColumns;

	bool done = false;

	while (done == false) {
	  
	  for (int j=0; j<numColumns; j++) {
	    theFile0 >> theData[j];
	  }
	  
	  if (theFile0.eof()) {
	    done = true;
	    theFile0.close();
	  } 
	  
	  if (done == false) {
	    theChannels[0]->sendVector(0, 0, data);	  
	  }
	} // while (done == false)
	
	delete [] theData;

      } // if (numColumns != 0)
      else
      theFile0.close();

    } // ordereing == 0

  } else  if (sendSelfCount > 0) {      

      ifstream theFile0;
      theFile0.open(fileName, ios::in);
      int ordering = -1;
      char c=theFile0.peek();
      if (c == 'O') {
	ordering = 0;
      }
      string s;

      // see if ORDERING: first few char of file
      if (ordering == 0) {

	int fileNameLength = strlen(fileName);
	sprintf(&fileName[fileNameLength-2],"");
	
	theFile.open(fileName, ios::out);

	ID sizeColumns(sendSelfCount+1);
	ID **theColumns = new ID *[sendSelfCount+1];
	double **theData = new double *[sendSelfCount+1];
	Vector **theRemoteData = new Vector *[sendSelfCount+1];

	char Ordering[11];
	theFile0.get(Ordering, 11);
	int numColumns =0;
	theFile0 >> numColumns;

	sizeColumns(0) = numColumns;
	if (numColumns != 0) {
	  theColumns[0] = new ID(numColumns);
	  theData[0] = new double [numColumns];
	} else {
	  theColumns[0] = 0;
	  theData[0] = 0;
	}      
	theRemoteData[0] = 0;

	for (int j=0; j<numColumns; j++)
	  theFile0 >> (*theColumns[0])[j];

	getline(theFile0, s);	// clean up to end of line
	

	int maxCount = 0;
	if (numColumns != 0)
	  maxCount = (*theColumns[0])[numColumns-1];

	for (int i=0; i<sendSelfCount; i++) {

	  static ID numColumnID(1);	  

	  theChannels[i]->recvID(0, 0, numColumnID);

	  int numColumns = numColumnID(0);

	  sizeColumns(i+1) = numColumns;
	  if (numColumns != 0) {
	    theColumns[i+1] = new ID(numColumns);
	    theChannels[i]->recvID(0, 0, *theColumns[i+1]);

	    if (numColumns != 0 && (*theColumns[i+1])[numColumns-1] > maxCount)
	      maxCount = (*theColumns[i+1])[numColumns-1];

	    theData[i+1] = new double [numColumns];
	    theRemoteData[i+1] = new Vector(theData[i+1], numColumns);
	  } else {
	    theColumns[i+1] = 0;
	    theData[i+1] = 0;
	    theRemoteData[i+1] = 0;
	  }
	}

	bool done = false;
	ID currentLoc(sendSelfCount+1);
	ID currentCount(sendSelfCount+1);
	
	Matrix printMapping(3,maxCount+1);
	
	for (int i=0; i<=sendSelfCount; i++) {
	  currentLoc(i) = 0;
	  if (theColumns[i] != 0)
	    currentCount(i) = (*theColumns[i])[0];
	  else
	    currentCount(i) = -1;
	}
	
	int count =0;
	while (count <= maxCount) {
	  for (int i=0; i<=sendSelfCount; i++) {
	    if (currentCount(i) == count) {
	      printMapping(0,count) = i;
	      
	      int maxLoc = theColumns[i]->Size();
	      int loc = currentLoc(i);
	      int columnCounter = 0;
	      
	      printMapping(1,count) = loc;
	      
	      while (loc < maxLoc && (*theColumns[i])(loc) == count) {
		loc++;
		columnCounter++;
	      }
	      
	      printMapping(2,count) = columnCounter;
	      
	      currentLoc(i) = loc;
	      
	      if (loc < maxLoc)
		currentCount(i) = (*theColumns[i])(loc);		
	      else
		currentCount(i) = -1; 		
	    }
	  }
	  
	  count++;
	}
	
            
	// go through each file, reading a line & sending to the output file
	done = false;


	while (done == false) {
	  for (int i=0; i<=sendSelfCount; i++) {
	    
	    int numColumns = sizeColumns(i);
	    double *data = theData[i];
	    
	    if (i == 0) {

	      if (numColumns != 0) {
		for (int j=0; j<numColumns; j++) {
		  theFile0 >> data[j];
		}
	      } else {
		getline(theFile0, s);	
	      }

	      if (theFile0.eof()) {
		done = true;
		theFile0.close();
	      }

	    } else { // i == 0

	      if (done == false) {
		if (numColumns != 0) {
		  theChannels[i-1]->recvVector(0, 0, *(theRemoteData[i]));
		}
	      }
	    }
	  }
	  
	  
	  if (done == false) {
	    for (int i=0; i<maxCount+1; i++) {
	      int fileID = printMapping(0,i);
	      int startLoc = printMapping(1,i);
	      int numData = printMapping(2,i);
	      double *data = theData[fileID];
	      for (int j=0; j<numData; j++)
		theFile << data[startLoc++] << " ";
	    }
	  }
	    
	  if (done == false)
	    theFile << "\n";
	}

	for (int i=0; i<=sendSelfCount; i++) {
	  if (theColumns[i] != 0) delete theColumns[i];
	  if (theData[i] != 0) delete [] theData[i];
	  if (theRemoteData[i] != 0) delete theRemoteData[i];
	}

	delete [] theColumns;
	delete [] theData;
	delete [] theRemoteData;

	theFile.close();
	
      } else {
	
	
      }
  }

#else
  if (sendSelfCount > 0) {

    int fileNameLength = strlen(fileName);
    sprintf(&fileName[fileNameLength-2],"");
    
    theFile.open(fileName, ios::out);

    ifstream **theFiles = new ifstream *[sendSelfCount+1];
    
    int ordering = -1;
    // open up the files
    for (int i=0; i<=sendSelfCount; i++) {
      theFiles[i] = new ifstream;
      sprintf(&fileName[fileNameLength-2],".%d",i);
      theFiles[i]->open(fileName, ios::in);
      char c=theFiles[i]->peek();
      if (c == 'O') {
	ordering = 0;
      }
    }

    // see if ORDERING: first few char of file
    if (ordering == 0) {
      ID sizeColumns(sendSelfCount+1);
      ID **theColumns = new ID *[sendSelfCount+1];
      double **theData = new double *[sendSelfCount+1];
      int maxCount = 0;

      for (int i=0; i<=sendSelfCount; i++) {
	char Ordering[11];

	theFiles[i]->get(Ordering, 11);
	int numColumns;
	(*theFiles[i]) >> numColumns;

	sizeColumns(i) = numColumns;
	if (numColumns != 0) {
	  theColumns[i] = new ID(numColumns);
	  theData[i] = new double [numColumns];
	} else {
	  theColumns[i] = 0;
	  theData[i] = 0;
	}

	for (int j=0; j<numColumns; j++)
	  (*theFiles[i]) >> (*theColumns[i])[j];

	if (numColumns != 0 && (*theColumns[i])[numColumns-1] > maxCount)
	  maxCount = (*theColumns[i])[numColumns-1];

      }
      
      bool done = false;
      ID currentLoc(sendSelfCount+1);
      ID currentCount(sendSelfCount+1);

      Matrix printMapping(3,maxCount+1);

      for (int i=0; i<=sendSelfCount; i++) {
	currentLoc(i) = 0;
	if (theColumns[i] != 0)
	  currentCount(i) = (*theColumns[i])[0];
	else
	  currentCount(i) = -1;
      }

      int count =0;
      while (count <= maxCount) {
	for (int i=0; i<=sendSelfCount; i++) {
	  if (currentCount(i) == count) {
	    printMapping(0,count) = i;
	    
	    int maxLoc = theColumns[i]->Size();
	    int loc = currentLoc(i);
	    int columnCounter = 0;

	    printMapping(1,count) = loc;

	    while (loc < maxLoc && (*theColumns[i])(loc) == count) {
	      loc++;
	      columnCounter++;
	    }
	    
	    printMapping(2,count) = columnCounter;

	    currentLoc(i) = loc;
	    
	    if (loc < maxLoc)
	      currentCount(i) = (*theColumns[i])(loc);		
	    else
	      currentCount(i) = -1; 		
	  }
	}

	count++;
      }

      // go through each file, reading a line & sending to the output file
      done = false;
      string s;
      
      while (done == false) {
	for (int i=0; i<=sendSelfCount; i++) {

	  if (done == false) {
	    int numColumns = sizeColumns(i);
	    double *data = theData[i];
	    for (int j=0; j<numColumns; j++) {
	      (*theFiles[i]) >> data[j];
	    }
	  }

	  if (theFiles[i]->eof()) {
	    done = true;
	    theFiles[i]->close();
	    delete theFiles[i];
	  }
	}

	if (done == false) {
	  for (int i=0; i<maxCount+1; i++) {
	    int fileID = printMapping(0,i);
	    int startLoc = printMapping(1,i);
	    int numData = printMapping(2,i);
	    double *data = theData[fileID];
	    for (int j=0; j<numData; j++)
	      theFile << data[startLoc++] << " ";
	  }
	}

	if (done == false)
	  theFile << "\n";
      }

      for (int i=0; i<=sendSelfCount; i++) {
	if (theColumns[i] != 0) delete theColumns[i];
	if (theData[i] != 0) delete [] theData[i];
      }

      delete [] theColumns;
      delete [] theData;
      delete [] theFiles;

      theFile.close();

    } else {

      // go through each file, reading a line & sending to the output file
      bool done = false;
      string s;
      
      while (done == false) {
	for (int i=0; i<=sendSelfCount; i++) {
	  bool eoline = false;
	  getline(*(theFiles[i]), s);	
	  theFile << s;
	  
	  if (theFiles[i]->eof()) {
	    done = true;
	    theFiles[i]->close();
	    delete theFiles[i];
	  }
	}
	theFile << "\n";
      }
      
      delete [] theFiles;
      theFile.close();
    }

  }
#endif

  if (theChannels != 0)
    delete [] theChannels;

  if (indentString != 0)
    delete [] indentString;

  if (fileName != 0)
    delete [] fileName;
}

int 
DataFileStream::setFile(const char *name, openMode mode)
{
  if (name == 0) {
    std::cerr << "DataFileStream::setFile() - no name passed\n";
    return -1;
  }

  // first create a copy of the file name
  if (fileName != 0) {
    if (strcmp(fileName, name) != 0)
      delete [] fileName;
    fileName = 0;
  }
  if (fileName == 0) {
    fileName = new char[strlen(name)+5];
    if (fileName == 0) {
      std::cerr << "DataFileStream::setFile() - out of memory copying name: " << name << std::endl;
      return -1;
    }
    
    // copy the strings
    strcpy(fileName, name);
  }

  // if file already open, close it
  if (fileOpen == 1) {
    theFile.close();
    fileOpen = 0;
  }

  if (mode == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

  return 0;
}

int 
DataFileStream::open(void)
{
  // check setFile has been called
  if (fileName == 0) {
    std::cerr << "DataFileStream::open(void) - no file name has been set\n";
    return -1;
  }

  // if file already open, return
  if (fileOpen == 1) {
    return 0;
  }

  if (sendSelfCount > 0) {
    strcat(fileName, ".0");
  }

  if (theOpenMode == OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  theOpenMode = APPEND;

  if (theFile.bad()) {
    std::cerr << "WARNING - DataFileStream::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;
    fileOpen = 0;
    return -1;
  } else
    fileOpen = 1;

  return 0;
}

int 
DataFileStream::close(void)
{
  if (fileOpen != 0)
    theFile.close();
  fileOpen = 0;

  return 0;
}


int 
DataFileStream::setPrecision(int prec)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << std::setprecision(prec);

  return 0;
}

int 
DataFileStream::setFloatField(floatField field)
{
  if (fileOpen == 0)
    this->open();

  if (field == FIXEDD) {
    if (fileOpen != 0)
      theFile << setiosflags(ios::fixed);
  }
  else if (field == SCIENTIFIC) {
    if (fileOpen != 0)
      theFile << setiosflags(ios::scientific);
  }

  return 0;
}


int 
DataFileStream::tag(const char *tagName)
{
  return 0;
}

int 
DataFileStream::tag(const char *tagName, const char *value)
{
  return 0;
}


int 
DataFileStream::endTag()
{
  return 0;
}

int 
DataFileStream::attr(const char *name, int value)
{
  return 0;
}

int 
DataFileStream::attr(const char *name, double value)
{
  return 0;
}

int 
DataFileStream::attr(const char *name, const char *value)
{
  return 0;
}

int 
DataFileStream::write(Vector &data)
{
  if (fileOpen == 0)
    this->open();

  (*this) << data;  

  return 0;
}



OPS_Stream& 
DataFileStream::write(const char *s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write(s, n);

  return *this;
}

OPS_Stream& 
DataFileStream::write(const unsigned char*s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
DataFileStream::write(const signed char*s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
DataFileStream::write(const void *s, int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}

OPS_Stream& 
DataFileStream::write(const double *s, int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0) {
    int nm1 = n-1;
    for (int i=0; i<n; i++)
      theFile << s[i] << " ";
    //    theFile << s[nm1] << endln;
    theFile << endln;
  }
  return *this;
}


OPS_Stream& 
DataFileStream::operator<<(char c)
{  
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(unsigned char c)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(signed char c)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(const char *s)
{
  if (fileOpen == 0)
    this->open();

  // note that we do the flush so that a "/n" before
  // a crash will cause a flush() - similar to what 
  if (fileOpen != 0) {
    theFile << s;
    theFile.flush();
  }

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(const unsigned char *s)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(const signed char *s)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(const void *p)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << p;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(unsigned int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(long n)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(unsigned long n)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(short n)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(unsigned short n)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(bool b)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << b;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(double n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(float n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << n;

  return *this;
}


int 
DataFileStream::sendSelf(int commitTag, Channel &theChannel)
{
  sendSelfCount++;

  Channel **theNextChannels = new Channel *[sendSelfCount];
  for (int i=0; i<sendSelfCount-1; i++)
    theNextChannels[i] = theChannels[i];
  theNextChannels[sendSelfCount-1] = &theChannel;
  if (theChannels != 0)
    delete [] theChannels;
  theChannels = theNextChannels;

  static ID idData(3);
  int fileNameLength = 0;
  if (fileName != 0)
    fileNameLength = strlen(fileName);

  idData(0) = fileNameLength;

  if (theOpenMode == OVERWRITE)
    idData(1) = 0;
  else
    idData(1) = 1;

  idData(2) = sendSelfCount;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "DataFileStream::sendSelf() - failed to send id data\n";
    return -1;
  }

  if (fileNameLength != 0) {
    Message theMessage(fileName, fileNameLength);
    if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
      opserr << "DataFileStream::sendSelf() - failed to send message\n";
      return -1;
    }
  }

  
  return 0;
}

int 
DataFileStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(3);

  sendSelfCount = -1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "DataFileStream::recvSelf() - failed to recv id data\n";
    return -1;
  }

  int fileNameLength = idData(0);
  if (idData(1) == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

  if (fileNameLength != 0) {
    if (fileName != 0)
      delete [] fileName;
    fileName = new char[fileNameLength+5];
    if (fileName == 0) {
      opserr << "DataFileStream::recvSelf() - out of memory\n";
      return -1;
    }

    Message theMessage(fileName, fileNameLength);
    if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
      opserr << "DataFileStream::recvSelf() - failed to recv message\n";
      return -1;
    }

    int tag = idData(2);

    sprintf(&fileName[fileNameLength],".%d",tag);

    if (this->setFile(fileName, theOpenMode) < 0) {
      opserr << "DataFileStream::DataFileStream() - setFile() failed\n";
      if (fileName != 0) {
	delete [] fileName;
	fileName = 0;
      }
    }
  }
  
  return 0;
}


void
DataFileStream::indent(void)
{
  if (fileOpen != 0)
    for (int i=0; i<numIndent; i++)
      theFile << indentString;
}
