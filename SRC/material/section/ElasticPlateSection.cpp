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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-05-19 06:20:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticPlateSection.cpp,v $

// Ed "C++" Love
//
//  Elastic Plate Section
//


#include <ElasticPlateSection.h>
#include <Matrix.h>
#include <Vector.h>

//parameters
const double ElasticPlateSection::five6 = 5.0/6.0 ; //shear correction

//static vector and matrices
Vector ElasticPlateSection::stress(5) ;
Matrix ElasticPlateSection::tangent(5,5) ;
ID     ElasticPlateSection::array(5) ;


//null constructor
ElasticPlateSection::ElasticPlateSection( ) : 
SectionForceDeformation( 0, SEC_TAG_ElasticPlateSection ), 
strain(5) 
{ }



//full constructor
ElasticPlateSection::ElasticPlateSection(  int    tag, 
                                           double young,
                                           double poisson,
                                           double thickness ) :
SectionForceDeformation( tag, SEC_TAG_ElasticPlateSection ),
strain(5)
{
  this->E  = young ;
  this->nu = poisson ;
  this->h  = thickness ;
}



//destructor
ElasticPlateSection::~ElasticPlateSection( ) 
{ } 



//make a clone of this material
SectionForceDeformation* ElasticPlateSection::getCopy( ) 
{
  ElasticPlateSection *clone ;   

  clone = new ElasticPlateSection( ) ; //new instance of this class

  *clone = *this ; //assignment to make copy

  return clone ;
}



//send back order of strain in vector form
int ElasticPlateSection::getOrder( ) const
{
  return 5 ;
}


//send back order of strain in vector form
const ID& ElasticPlateSection::getType( ) const 
{
  return array ;
}



//swap history variables
int ElasticPlateSection::commitState( ) 
{
  return 0 ;
}



//revert to last saved state
int ElasticPlateSection::revertToLastCommit( )
{
  return 0 ;
}

//revert to start
int ElasticPlateSection::revertToStart( )
{
  return 0 ;
}


//get the strain 
int ElasticPlateSection ::
setTrialSectionDeformation( const Vector &strain_from_element)
{
  this->strain = strain_from_element ;

  return 0 ;
}


//send back the strain
const Vector& ElasticPlateSection::getSectionDeformation( )
{
  return this->strain ;
}


//send back the stress 
const Vector& ElasticPlateSection::getStressResultant( )
{
  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ; //bending modulus

  double G  =  0.5 * E / ( 1.0 + nu ) ; //shear modulus
 
  G *= five6 ;
  G *= h ;


  stress(0) = -( D*strain(0) + nu*D*strain(1) ) ;
 
  stress(1) = -( nu*D*strain(0) + D*strain(1) ) ;

  stress(2) = -0.5*D*( 1.0 - nu )*strain(2) ;

  stress(3) = G*strain(3) ;

  stress(4) = G*strain(4) ;

 
  return this->stress ;
}


//send back the tangent 
const Matrix& ElasticPlateSection::getSectionTangent( )
{

  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;

  double G  =  0.5 * E / ( 1.0 + nu ) ;


  tangent.Zero() ;

  tangent(0,0) = -D ;
  tangent(1,1) = -D ;

  tangent(0,1) = -nu*D ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(3,3) = five6*G*h ;

  tangent(4,4) = tangent(3,3) ;


  return this->tangent ;
}

//print out data
void  ElasticPlateSection::Print( ostream &s, int flag )
{
  s << "ElasticPlateSection: \n " ;
  s <<  "  Young's Modulus E  = "  <<  E  <<  '\n' ;
  s <<  "  Poisson's Ratio nu = " <<  nu <<  '\n' ;
  s <<  "  Thickness h = "        <<  h  <<  '\n' ;

  return ;
}


int 
ElasticPlateSection::sendSelf(int commitTag, Channel &theChannel) 
{
  return -1;
}


int 
ElasticPlateSection::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}
