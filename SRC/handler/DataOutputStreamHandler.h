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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputStreamHandler.h,v $

#ifndef _DataOutputStreamHandler
#define _DataOutputStreamHandler

#include <DataOutputHandler.h>
#include <StandardStream.h>

class DataOutputStreamHandler : public DataOutputHandler
{
 public:
  DataOutputStreamHandler(bool echoDescription = false);
  virtual ~DataOutputStreamHandler();

  virtual int open(char **dataDescription, int numData);
  virtual int write(Vector &data);

 private:
  StandardStream outputStream;
  bool echoDescription;
  int numColumns;
};

#endif
