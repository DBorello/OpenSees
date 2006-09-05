/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
// $Date: 2006-09-05 20:21:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/Parameter.h,v $

#ifndef Parameter_h
#define Parameter_h

#include <Information.h>
#include <TaggedObject.h>

class DomainComponent;
class MovableObject;

class Parameter : public TaggedObject
{
 public:
  Parameter(int tag, DomainComponent *theObject,
	    const char **argv, int argc);
  Parameter(const Parameter &param);
  ~Parameter();

  int addObject(DomainComponent *theObject, const char **argv, int argc);
  
  void Print(OPS_Stream &s, int flag =0);
  
  int update(int newValue); 
  int update(double newValue); 
  int activate(bool active);
  
  int addObject(int parameterID, MovableObject *object);

 protected:
  
 private:
  Information theInfo;
  enum {initialSize = 100};
  enum {expandSize = 100};
  MovableObject **theObjects;
  int *parameterID;
  int numObjects;
  int maxNumObjects;
};

#endif
