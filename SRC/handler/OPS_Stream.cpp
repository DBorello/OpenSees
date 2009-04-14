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
                                                                        
// $Revision: 1.5 $
// $Date: 2009-04-14 21:11:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/OPS_Stream.cpp,v $

#include <OPS_Stream.h>

OPS_Stream::OPS_Stream(int classTag)
  :MovableObject(classTag)
{
  
}
OPS_Stream::~OPS_Stream()
{

}


 OPS_Stream& 
 OPS_Stream::write(const char *s, int n) {return *this;}


 OPS_Stream& 
 OPS_Stream::write(const unsigned char *s, int n) {return *this;}

OPS_Stream& 
OPS_Stream::write(const signed char *s, int n) {return *this;}

OPS_Stream& 
OPS_Stream::write(const void *s, int n) {return *this;}

OPS_Stream& 
OPS_Stream::write(const double *d, int n) 
{
  int nm1 = n-1;

  for (int i=0; i<nm1; i++)
    (*this) << d[i] << " ";
  (*this) << d[nm1] << "\n";

  return *this;
}

OPS_Stream& 
OPS_Stream::operator<<(char c) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(unsigned char c) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(signed char c) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(const char *s) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(const unsigned char *s) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(const signed char *s) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(const void *p) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(int n) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(unsigned int n) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(long n) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(unsigned long n) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(short n) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(unsigned short n) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(bool b) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(double n) {return *this;}

OPS_Stream& 
OPS_Stream::operator<<(float n) {return *this;}
