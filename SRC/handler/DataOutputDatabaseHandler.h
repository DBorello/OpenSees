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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputDatabaseHandler.h,v $

#ifndef _DataOutputDatabaseHandler
#define _DataOutputDatabaseHandler

#include <DataOutputHandler.h>
class FE_Datastore;

class DataOutputDatabaseHandler : public DataOutputHandler
{
 public:
  DataOutputDatabaseHandler(FE_Datastore *database, const char *tableName);
  virtual ~DataOutputDatabaseHandler();

  virtual int open(char **dataDescription, int numData);
  virtual int write(Vector &data);

  int setDatabase(FE_Datastore &theDatabase, const char *tableName);

 private:
  FE_Datastore *theDatabase;
  char *tableName;

  int numColumns;
  char **columns;
  int commitTag;
};

#endif
