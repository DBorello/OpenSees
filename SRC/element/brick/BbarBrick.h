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
// $Date: 2001-07-11 21:54:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/BbarBrick.h,v $

// Ed "C++" Love
//
// Eight node BbarBrick element 
//

#include <iostream.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class BbarBrick : public Element {

  public :
    
    //null constructor
    BbarBrick( ) ;
  
    //full constructor
    BbarBrick( int tag, 
			int node1,
			int node2,
		        int node3,
			int node4,
			int node5,
			int node6,
			int node7,
			int node8,
			NDMaterial &theMaterial ) ;

    //destructor 
    virtual ~BbarBrick( ) ;

    //set domain because frank is a dumb ass 
    void setDomain( Domain *theDomain ) ;

    //get the number of external nodes
    int getNumExternalNodes( ) const ;
 
    //return connected external nodes
    const ID &getExternalNodes( ) ;

    //return number of dofs
    int getNumDOF( ) ;

    //commit state
    int commitState( ) ;
    
    //revert to last commit 
    int revertToLastCommit( ) ;
    
    //revert to start 
    int revertToStart( ) ;

    //print out element data
    void Print( ostream &s, int flag ) ;
	
    //return stiffness matrix 
    const Matrix &getTangentStiff( ) ;
    
    //return secant matrix 
    const Matrix &getSecantStiff( ) ;
    
    //return damping matrix because frank is a dumb ass 
    const Matrix &getDamp( ) ;
    
    //return mass matrix
    const Matrix &getMass( ) ;

    //zero the load -- what load?
    void zeroLoad( ) ;
 
    //add load -- what load?
    int addLoad( const Vector &addP ) ;

    //get residual
    const Vector &getResistingForce( ) ;
    
    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
      
    //plotting 
    int displaySelf(Renderer &theViewer, int displayMode, float fact);

  private : 

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

    //quadrature data
    static const double root3 ;
    static const double one_over_root3 ;    
    static const double sg[2] ;
    static const double wg[8] ;

  
    //node information
    ID connectedExternalNodes ;  //four node numbers
    Node *nodePointers[8] ;      //pointers to four nodes


    //material information
    NDMaterial *materialPointers[8] ; //pointers to eight materials
					  
    //local nodal coordinates, three coordinates for each of four nodes
    //    static double xl[3][8] ; 
    static double xl[][8] ; 

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //compute coordinate system
    void computeBasis( ) ;

    //compute Bbar matrix
    Matrix computeBbar( int node, 
			const double shp[4][8], 
			const double shpBar[4][8] ) ;
  
    //Matrix transpose
    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

} ; 







