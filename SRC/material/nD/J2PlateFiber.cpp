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
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2001-07-16 22:19:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2PlateFiber.cpp,v $

// Written: Ed "C++" Love
//
// J2PlateFiber isotropic hardening material class
// 
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi) 
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_infty + (sigma_0 - sigma_infty)*exp(-delta*xi) + H*xi 
//
//  Flow Rules
//  \dot{epsilon_p} =  gamma * d_phi/d_sigma
//  \dot{xi}        = -gamma * d_phi/d_q 
//
//  Linear Viscosity 
//  gamma = phi / eta  ( if phi > 0 ) 
//
//  Backward Euler Integration Routine 
//  Yield condition enforced at time n+1 
//
//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//                    2 eps_01   }   <--- note the 2
// 
//  set eta := 0 for rate independent case
//

#include <J2PlateFiber.h>

Vector J2PlateFiber :: strain_vec(5) ;
Vector J2PlateFiber :: stress_vec(5) ;
Matrix J2PlateFiber :: tangent_matrix(5,5) ;

//null constructor
J2PlateFiber ::  J2PlateFiber( ) : 
J2Plasticity( ) 
{  }


//full constructor
J2PlateFiber :: 
J2PlateFiber(   int    tag, 
                 double K,
                 double G,
                 double yield0,
                 double yield_infty,
                 double d,
                 double H,
                 double viscosity ) : 
J2Plasticity(tag, ND_TAG_J2PlateFiber, 
             K, G, yield0, yield_infty, d, H, viscosity )
{ }


//elastic constructor
J2PlateFiber :: 
J2PlateFiber(   int    tag, 
                 double K, 
                 double G ) :
J2Plasticity(tag, ND_TAG_J2PlateFiber, K, G )
{ }


//destructor
J2PlateFiber :: ~J2PlateFiber( ) 
{  } 


//make a clone of this material
NDMaterial* J2PlateFiber :: getCopy( ) 
{ 
  J2PlateFiber  *clone;
  clone = new J2PlateFiber( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* J2PlateFiber :: getType( ) const 
{
  return "PlateFiber" ;
}


//send back order of strain in vector form
int J2PlateFiber :: getOrder( ) const 
{ 
  return 5 ; 
} 

//get the strain and integrate plasticity equations
int J2PlateFiber :: setTrialStrain( const Vector &strain_from_element ) 
{
  const double tolerance = 1e-8 ;

  const int max_iterations = 25 ;
  int iteration_counter  = 0 ;

  int i, j, k, l ;
  int ii, jj ;

  double eps22  =  strain(2,2) ;
  strain.Zero( ) ;

  strain(0,0) =        strain_from_element(0) ;
  strain(1,1) =        strain_from_element(1) ;

  strain(0,1) = 0.50 * strain_from_element(2) ;
  strain(1,0) =        strain(0,1) ;

  strain(1,2) = 0.50 * strain_from_element(3) ;
  strain(2,1) =        strain(1,2) ;
  
  strain(2,0) = 0.50 * strain_from_element(4) ;
  strain(0,2) =        strain(2,0) ;

  strain(2,2) =        eps22 ; 

  //enforce the plane stress condition sigma_22 = 0 
  //solve for epsilon_22 
  iteration_counter = 0 ;  
  do {

     this->plastic_integrator( ) ;
    
     strain(2,2) -= stress(2,2) / tangent[2][2][2][2] ;

     //cerr << stress(2,2) << endl ;

     iteration_counter++ ;
     if ( iteration_counter > max_iterations ) {
       cerr << "More than " << max_iterations ;
       cerr << " iterations in setTrialStrain of J2PlateFiber \n" ;
       break ;
     }// end if 

  } while ( fabs(stress(2,2)) > tolerance ) ;

  //modify tangent for plane stress 
  for ( ii = 0; ii < 5; ii++ ) {
    for ( jj = 0; jj < 5; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          tangent[i][j][k][l] -=   tangent[i][j][2][2] 
                                 * tangent[2][2][k][l] 
                                 / tangent[2][2][2][2] ;

          //minor symmetries 
          tangent [j][i][k][l] = tangent[i][j][k][l] ;
          tangent [i][j][l][k] = tangent[i][j][k][l] ;
          tangent [j][i][l][k] = tangent[i][j][k][l] ;

    } // end for jj
  } // end for ii 

  return 0 ;
}


//unused trial strain functions
int J2PlateFiber :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int J2PlateFiber :: setTrialStrainIncr( const Vector &v ) 
{
    return -1 ;
}

int J2PlateFiber :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
    return -1 ;
}



//send back the strain
const Vector& J2PlateFiber :: getStrain( ) 
{

  strain_vec(0) =       strain(0,0) ;
  strain_vec(1) =       strain(1,1) ;

  strain_vec(2) = 2.0 * strain(0,1) ;

  strain_vec(3) = 2.0 * strain(1,2) ;

  strain_vec(4) = 2.0 * strain(2,0) ;

  return strain_vec ;
} 


//send back the stress 
const Vector& J2PlateFiber :: getStress( ) 
{
 
  stress_vec(0) = stress(0,0) ;
  stress_vec(1) = stress(1,1) ;

  stress_vec(2) = stress(0,1) ;

  stress_vec(3) = stress(1,2) ;
  
  stress_vec(4) = stress(2,0) ;

  return stress_vec ;
}

//send back the tangent 
const Matrix& J2PlateFiber :: getTangent( ) 
{

  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 )
  //   3           1 2  ( or 2 1 )
  //   4           2 0  ( or 0 2 ) 
    
  int ii, jj ;
  int i, j, k, l ;

  for ( ii = 0; ii < 5; ii++ ) {
    for ( jj = 0; jj < 5; jj++ ) {

      index_map( ii, i, j ) ;
      index_map( jj, k, l ) ;

      tangent_matrix(ii,jj) = tangent[i][j][k][l] ;

    } //end for j
  } //end for i
       

  return tangent_matrix ;
} 

//this is mike's problem
int J2PlateFiber :: setTrialStrain(const Tensor &v) 
{
  return -1 ;
}

int J2PlateFiber :: setTrialStrain(const Tensor &v, const Tensor &r)     
{
  return -1 ;
}

int J2PlateFiber :: setTrialStrainIncr(const Tensor &v) 
{
  return -1 ;
}

int J2PlateFiber :: setTrialStrainIncr(const Tensor &v, const Tensor &r) 
{
  return -1 ;
}

const Tensor& J2PlateFiber :: getTangentTensor( ) 
{
  return rank4 ;
}

//jeremic@ucdavis.edu 22jan2001const Tensor& J2PlateFiber :: getStressTensor( ) 
//jeremic@ucdavis.edu 22jan2001{
//jeremic@ucdavis.edu 22jan2001  return rank2 ;
//jeremic@ucdavis.edu 22jan2001}
//jeremic@ucdavis.edu 22jan2001
//jeremic@ucdavis.edu 22jan2001const Tensor& J2PlateFiber :: getStrainTensor( ) 
//jeremic@ucdavis.edu 22jan2001{
//jeremic@ucdavis.edu 22jan2001  return rank2 ;
//jeremic@ucdavis.edu 22jan2001}


//this is frank's problem
int J2PlateFiber :: sendSelf(int commitTag, Channel &theChannel)
{
  return -1 ;
}

int J2PlateFiber :: recvSelf(int commitTag, Channel &theChannel, 
	                      FEM_ObjectBroker &theBroker)
{
  return -1 ;
}


//matrix_index ---> tensor indices i,j
// plane stress different because of condensation on tangent
// case 3 switched to 1-2 and case 4 to 3-3 
void J2PlateFiber :: index_map( int matrix_index, int &i, int &j )
{
  switch ( matrix_index+1 ) { //add 1 for standard tensor indices

    case 1 :
      i = 1 ; 
      j = 1 ;
      break ;
 
    case 2 :
      i = 2 ;
      j = 2 ; 
      break ;

    case 3 :
      i = 1 ;
      j = 2 ;
      break ;

    case 4 :
      i = 2 ;
      j = 3 ;
      break ;

    case 5 :
      i = 3 ;
      j = 1 ;
      break ;

    case 6 :
      i = 3 ;
      j = 3 ;
      break ;


    default :
      i = 1 ;
      j = 1 ;
      break ;

  } //end switch

i-- ; //subtract 1 for C-indexing
j-- ;

return ; 
}

