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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputFileHandler.h,v $

#ifndef _DataOutputFileHandler
#define _DataOutputFileHandler

#include <DataOutputHandler.h>
#include <FileStream.h>

enum echoMode  {NONE, DATA_FILE, XML_FILE};

class DataOutputFileHandler : public DataOutputHandler
{
 public:
  DataOutputFileHandler(const char *fileName =0, 
			echoMode = NONE, 
			openMode mode = OVERWRITE);
  virtual ~DataOutputFileHandler();

  virtual int open(char **dataDescription, int numData);
  virtual int write(Vector &data);

 private:
  FileStream outputFile;
  char *fileName;

  echoMode theEchoMode;
  openMode theOpenMode;
  int numColumns;
};

#endif
