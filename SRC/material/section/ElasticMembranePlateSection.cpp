// Ed "C++" Love
// Do not ask Prashant about this code.  He has no clue. 
//
//  Elastic Plate Section with membrane
//


#include <ElasticMembranePlateSection.h>


//parameters
const double ElasticMembranePlateSection::five6 = 5.0/6.0 ; //shear correction

//static vector and matrices
Vector  ElasticMembranePlateSection::stress(8) ;
Matrix  ElasticMembranePlateSection::tangent(8,8) ;
ID      ElasticMembranePlateSection::array(8) ;


//null constructor
ElasticMembranePlateSection::ElasticMembranePlateSection( ) : 
SectionForceDeformation( 0, SEC_TAG_ElasticMembranePlateSection ), 
strain(8) 
{ }



//full constructor
ElasticMembranePlateSection::ElasticMembranePlateSection(  
                                           int    tag, 
                                           double young,
                                           double poisson,
                                           double thickness ) :
SectionForceDeformation( tag, SEC_TAG_ElasticMembranePlateSection ),
strain(8)
{
  this->E  = young ;
  this->nu = poisson ;
  this->h  = thickness ;
}



//destructor
ElasticMembranePlateSection::~ElasticMembranePlateSection( ) 
{ } 



//make a clone of this material
SectionForceDeformation*  ElasticMembranePlateSection::getCopy( ) 
{
  ElasticMembranePlateSection *clone ;   

  clone = new ElasticMembranePlateSection( ) ; //new instance of this class

  *clone = *this ; //assignment to make copy

  return clone ;
}



//send back order of strain in vector form
int ElasticMembranePlateSection::getOrder( ) const
{
  return 8 ;
}


//send back order of strain in vector form
const ID& ElasticMembranePlateSection::getType( ) const 
{
  return array ;
}



//swap history variables
int ElasticMembranePlateSection::commitState( ) 
{
  return 0 ;
}



//revert to last saved state
int ElasticMembranePlateSection::revertToLastCommit( )
{
  return 0 ;
}

//revert to start
int ElasticMembranePlateSection::revertToStart( )
{
  return 0 ;
}


//get the strain 
int ElasticMembranePlateSection ::
setTrialSectionDeformation( const Vector &strain_from_element)
{
  this->strain = strain_from_element ;

  return 0 ;
}


//send back the strain
const Vector& ElasticMembranePlateSection::getSectionDeformation( )
{
  return this->strain ;
}


//send back the stress 
const Vector&  ElasticMembranePlateSection::getStressResultant( )
{

  double M  = E / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * E / ( 1.0 + nu ) ; //shear modulus
 
  G *= h ;  //multiply by thickness
  M *= h ;

  //membrane resultants

  stress(0) =  M*strain(0) + (nu*M)*strain(1)  ;
 
  stress(1) =  (nu*M)*strain(0) +  M*strain(1)  ;

  stress(2) =  G*strain(2) ;

 

  G *= five6 ;  //multiply by shear correction factor

  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  //bending resultants

  stress(3) = -( D*strain(3) + nu*D*strain(4) ) ;
 
  stress(4) = -( nu*D*strain(3) + D*strain(4) ) ;

  stress(5) = -0.5*D*( 1.0 - nu )*strain(5) ;

  stress(6) = G*strain(6) ;

  stress(7) = G*strain(7) ;

 
  return this->stress ;
}


//send back the tangent 
const Matrix&  ElasticMembranePlateSection::getSectionTangent( )
{

  double M  = E / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * E / ( 1.0 + nu ) ; //shear modulus

  G *= h ;  //multiply by thickness
  M *= h ;

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = M ;
  tangent(1,1) = M ;

  tangent(0,1) = nu*M ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = G ;



  G *= five6 ;  //multiply by shear correction factor

  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  //bending tangent terms

  tangent(3,3) = -D ;
  tangent(4,4) = -D ;

  tangent(3,4) = -nu*D ;
  tangent(4,3) = tangent(3,4) ;

  tangent(5,5) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(6,6) = G ;

  tangent(7,7) = tangent(6,6) ;


  return this->tangent ;
}

//print out data
void  ElasticMembranePlateSection::Print( ostream &s, int flag )
{
  s << "ElasticMembranePlateSection: \n " ;
  s <<  "  Young's Modulus E = "  <<  E  <<  '\n' ;
  s <<  "  Poisson's Ratio nu = " <<  nu <<  '\n' ;
  s <<  "  Thickness h = "        <<  h  <<  '\n' ;

  return ;
}

int 
ElasticMembranePlateSection::sendSelf(int commitTag, Channel &theChannel) 
{
  return -1;
}


int 
ElasticMembranePlateSection::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}
 
