// $Revision: 1.3 $
// $Date: 2001-09-13 19:11:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/MultiYieldSurface.h,v $
                                                                        
// Written: ZHY
// Created: August 2000
//
// MultiYieldSurface.h
// -------------------
//

#ifndef _MultiYieldSurface_H_
#define _MultiYieldSurface_H_

#include <T2Vector.h>


#define LOCK_VALUE  1.0e+30

// global function to find the roots of a second order equation
double secondOrderEqn(double A, double B, double C, int i);

// define yield surface in stress space
class MultiYieldSurface
{
 
public:
  //constructors
  MultiYieldSurface();
  MultiYieldSurface(const Vector & center_init, double size_init, 
                    double plas_modul); 
  ~MultiYieldSurface();
	const Vector & center() const {return theCenter; }
	double size() const {return theSize; }
	double modulus() const {return plastShearModulus; }
  void  setCenter(const Vector & newCenter);
  friend ostream & operator<< (ostream & os, const MultiYieldSurface & );  
	friend istream & operator>> (istream & is, MultiYieldSurface & );

protected:

private:
  double theSize;
  Vector theCenter;  
  double plastShearModulus;

};

#endif
