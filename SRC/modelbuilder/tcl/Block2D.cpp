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
// $Date: 2001-07-12 23:08:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/Block2D.cpp,v $
                                                                        
// Written: Ed Love
// Created: 07/99
//
// Description: This file contains the class definition for Block2D.

//
// What: "@(#) Block2D.cpp, revA"


#include <Block2D.h>


//constructor
Block2D::Block2D(int numx, int numy, 
		 const ID& nodeID, 
		 const Matrix& coorArray ) 
: 
coor(3), 
element(4) 
{
  this->nx = numx;
  this->ny = numy;

  this->setUpXl( nodeID, coorArray );

}


//destructor
Block2D::~Block2D( )
{ }


//set up xl array
void  Block2D::setUpXl( const ID &nodeID, const Matrix &coorArray ) 
{

  int i, j;

  for ( i=0; i<4; i++ ){
    if ( nodeID(i) == -1 ) {
      cerr << "Warning : in Block2D, block node " 
	   << i 
	   << " is not defined.  No Generation will take place."
	   << endl;
      break; 
    }//end if
  }//end for i

  //local storage xl = transpose(coorArray)
  for ( i=0; i<3; i++ ) {
    for ( j=0; j<9; j++ )
      xl[i][j] = coorArray(j,i);
  }//end for i


  if ( nodeID(4) == -1 ) {
    for ( i=0; i<3; i++ )
      xl[i][4] = 0.5*( xl[i][0] + xl[i][1] );
  }//endif

  if ( nodeID(5) == -1 ) {
    for ( i=0; i<3; i++ )
      xl[i][5] = 0.5*( xl[i][1] + xl[i][2] );
  }//endif

  if ( nodeID(6) == -1 ) {
    for ( i=0; i<3; i++ )
      xl[i][6] = 0.5*( xl[i][2] + xl[i][3] );
  }//endif

  if ( nodeID(7) == -1 ) {
    for ( i=0; i<3; i++ )
      xl[i][7] = 0.5*( xl[i][3] + xl[i][0] );
  }//endif

  if ( nodeID(8) == -1 ) {
    for ( i=0; i<3; i++ )
      xl[i][8] = 0.25*( xl[i][0] + xl[i][1] + xl[i][2] + xl[i][3] );
  }//endif

  return;
}


//generate node
const Vector&
Block2D::getNodalCoords( int i, int j )
{

  double hx = 2.0 / nx;

  double hy = 2.0 / ny;

  double x = -1.0 + (i*hx);

  double y = -1.0 + (j*hy);

  coor(0) = x;
  coor(1) = y;
  coor(2) = 0.0;

  this->transformNodalCoordinates( );
  
  return coor;
}


//generate element
const ID&
Block2D::getElementNodes( int i, int j )  
{
  
  int nenx = nx + 1;

  int neny = ny + 1;

  int node1, node2, node3, node4;

  node1 = i + j*nenx;
  node2 = node1 + 1;

  node3 = node2 + nenx;
  node4 = node1 + nenx;

  element(0) = node1;
  element(1) = node2;
  element(2) = node3;
  element(3) = node4;

  return element;
}



//transform to real coordiantes
void  Block2D::transformNodalCoordinates( )
{

  static double shape[9];
  
  static double natCoor[2];

  natCoor[0] = coor(0);
  natCoor[1] = coor(1);

  coor.Zero( );

  this->shape2d( natCoor[0], natCoor[1], shape );

  for ( int j=0; j<9; j++ ) {
      
    for ( int dim=0; dim<3; dim++ )
      coor(dim) += shape[j]*xl[dim][j];

  } //end for j

  return;

}


//shape functions
void  Block2D::shape2d( double x, double y, 
	                double shape[9]     ) 
{
  static double Nx[3];
  static double Ny[3];

  Nx[0] = 0.5 * x * ( x - 1.0 );
  Nx[1] = 1.0 - (x*x);
  Nx[2] = 0.5 * x * ( x + 1.0 );

  Ny[0] = 0.5 * y * ( y - 1.0 );
  Ny[1] = 1.0 - (y*y);
  Ny[2] = 0.5 * y * ( y + 1.0 );

  shape[0] = Nx[0]*Ny[0];
  shape[1] = Nx[2]*Ny[0];
  shape[2] = Nx[2]*Ny[2];
  shape[3] = Nx[0]*Ny[2];

  shape[4] = Nx[1]*Ny[0];
  shape[5] = Nx[2]*Ny[1];
  shape[6] = Nx[1]*Ny[2];
  shape[7] = Nx[0]*Ny[1];

  shape[8] = Nx[1]*Ny[1];

  return;
}
