// Ed "C++" Love
// Do not ask Prashant about this code.  He has no clue. 
//
//  Elastic Plate Section
//

#include <iostream.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <SectionForceDeformation.h>

#define SEC_TAG_ElasticPlateSection 4005 


class ElasticPlateSection : public SectionForceDeformation {

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    ElasticPlateSection( ) ;

    //full constructor
    ElasticPlateSection(   int    tag, 
                           double E,
                           double nu,
                           double h = 1.0 ) ;



    //destructor
    ~ElasticPlateSection( ) ;

    //make a clone of this material
    SectionForceDeformation *getCopy( ) ;


    //send back order of strain in vector form
    int getOrder( ) const ;

    //send back order of strain in vector form
    const ID& getType( ) const ;

    //swap history variables
    int commitState( ) ; 

    //revert to last saved state
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    //get the strain and integrate plasticity equations
    int setTrialSectionDeformation( const Vector &strain_from_element ) ;

    //send back the strain
    const Vector& getSectionDeformation( ) ;

    //send back the stress 
    const Vector& getStressResultant( ) ;

    //send back the tangent 
    const Matrix& getSectionTangent( ) ;

    //print out data
    void Print( ostream &s, int flag ) ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  private :

    double E ;  // elastic modulus
    double nu ; // poisson ratio
    double h  ; // plate thickness

    static const double five6 ; // =5/6 = shear correction factor

    Vector strain ;

    static Vector stress ;

    static Matrix tangent ;

    static ID array ;  

} ; //end of ElasticPlateSection declarations


