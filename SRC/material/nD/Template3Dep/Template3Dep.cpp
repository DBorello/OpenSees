/*
################################################################################
# COPYRIGHT (C):     :-))                                                      #
# PROJECT:           Object Oriented Finite Element Program                    #
# PURPOSE:           General platform for elaso-plastic constitutive model     #
#                    implementation                                            #
# CLASS:             Template3Dep (the base class for all material point)     #
#                                                                              #
# VERSION:                                                                     #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
# TARGET OS:         DOS || UNIX || . . .                                      #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                               #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
#                                                                              #
#                                                                              #
# DATE:              08-03-2000                                                #
# UPDATE HISTORY:    09-12-2000                                                #
#                                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
# SHORT EXPLANATION: This file contains the class implementation for           #
#                    Template3Dep.                                             #
#                                                                              #
################################################################################
*/

#ifndef Template3Dep_CPP
#define Template3Dep_CPP

//#include <string.h>

#include "Template3Dep.h"
#include <Channel.h>
#include <G3Globals.h>

#include <DP_YS.h>
#include <DP_PS.h>
#include <EPState.h>

//================================================================================
// Constructor
//================================================================================

Template3Dep::Template3Dep( int tag                 ,
                            YieldSurface     *YS_   ,         
                            PotentialSurface *PS_   , 
		  	    EPState          *EPS_  ,
	       	     	    EvolutionLaw_S   *ELS1_ , 
	       	     	    EvolutionLaw_S   *ELS2_ , 
	       	     	    EvolutionLaw_S   *ELS3_ , 
	       	     	    EvolutionLaw_S   *ELS4_ , 
	       	     	    EvolutionLaw_T   *ELT1_ ,
	       	     	    EvolutionLaw_T   *ELT2_ ,
	       	     	    EvolutionLaw_T   *ELT3_ ,
	       	     	    EvolutionLaw_T   *ELT4_ )
:NDMaterial(tag, ND_TAG_Template3Dep)
{      	      
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
    
    // Evolution laws
    if ( ELS1_ ) 
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;

    if ( ELS2_ ) 
       ELS2 = ELS2_->newObj();
    else
       ELS2 = 0;

    if ( ELS3_ ) 
       ELS3 = ELS3_->newObj();
    else
       ELS3 = 0;

    if ( ELS4_ ) 
       ELS4 = ELS4_->newObj();
    else
       ELS4 = 0;

    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;

    if ( ELT2_ ) 
       ELT2 = ELT2_->newObj();
    else
       ELT2 = 0;

    if ( ELT3_ ) 
       ELT3 = ELT3_->newObj();
    else
       ELT3 = 0;

    if ( ELT4_ ) 
       ELT4 = ELT4_->newObj();
    else
       ELT4 = 0;

}

//================================================================================
// Constructor 0
//================================================================================
Template3Dep::Template3Dep( int tag               ,
                            YieldSurface     *YS_ ,        
                            PotentialSurface *PS_ ,
              	            EPState          *EPS_)
:NDMaterial(tag, ND_TAG_Template3Dep)
{
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }

    // Evolution laws
    ELS1 = 0;
    ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;
    ELT1 = 0;
    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}


//================================================================================
// Constructor 1
//================================================================================

Template3Dep::Template3Dep(   int tag               ,
                              YieldSurface     *YS_ ,        
                              PotentialSurface *PS_ ,
              	              EPState          *EPS_,
	       	              EvolutionLaw_S   *ELS1_ )
:NDMaterial(tag, ND_TAG_Template3Dep)
{

    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
    
    // Evolution laws
    if ( ELS1_ )
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;

    ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;
    ELT1 = 0;
    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}


//================================================================================
// Constructor 2
//================================================================================

Template3Dep::Template3Dep(   int tag               ,
                              YieldSurface     *YS_ ,
                              PotentialSurface *PS_ ,
              	              EPState          *EPS_,
	       	              EvolutionLaw_T   *ELT1_ ) 
:NDMaterial(tag, ND_TAG_Template3Dep)
{

    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
    
    // Evolution laws
    ELS1 = 0;
    ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;
    
    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;

    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}

//================================================================================
// Constructor 3
//================================================================================

Template3Dep::Template3Dep(   int tag               ,
                              YieldSurface     *YS_ ,        
                              PotentialSurface *PS_ ,
              	    	      EPState          *EPS_,
	       	     	      EvolutionLaw_S   *ELS1_, 
	       	     	      EvolutionLaw_T   *ELT1_ )
:NDMaterial(tag, ND_TAG_Template3Dep)
{
    
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
    
    // Evolution laws
    if ( ELS1_ )
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;
    ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;

    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;
    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}

//================================================================================
// Constructor 4
// Two scalar evolution laws and one tensorial evolution law are provided!
//================================================================================

Template3Dep::Template3Dep(   int tag               ,
                              YieldSurface     *YS_ ,
                              PotentialSurface *PS_ ,
              	    	      EPState          *EPS_,
	       	     	      EvolutionLaw_S   *ELS1_, 
     	       	              EvolutionLaw_S   *ELS2_, 
	       	     	      EvolutionLaw_T   *ELT1_ )
:NDMaterial(tag, ND_TAG_Template3Dep)
{
    
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
    
    if ( ELS1_ )
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;

    if ( ELS2_ )
       ELS2 = ELS2_->newObj();
    else
       ELS2 = 0;

    ELS3 = 0;
    ELS4 = 0;

    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;
    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}

//================================================================================
// Constructor 5
// Two scalar evolution laws and two tensorial evolution law are provided!
//================================================================================

Template3Dep::Template3Dep(   int tag               ,
                              YieldSurface     *YS_ ,        
                              PotentialSurface *PS_ ,
              	    	      EPState          *EPS_,
	       	     	      EvolutionLaw_S   *ELS1_, 
     	       	              EvolutionLaw_S   *ELS2_, 
	       	     	      EvolutionLaw_T   *ELT1_,
			      EvolutionLaw_T   *ELT2_)
:NDMaterial(tag, ND_TAG_Template3Dep)
{
    
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
       g3ErrorHandler->fatal("Template3Dep:: Template3Dep failed to construct the template3Dep");
       exit(-1);
    }
    
    if ( ELS1_ )
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;

    if ( ELS2_ )
       ELS2 = ELS2_->newObj();
    else
       ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;

    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;

    if ( ELT2_ ) 
       ELT2 = ELT2_->newObj();
    else
       ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}

//================================================================================
// Constructor 6
// NO parameter is provided
//================================================================================

Template3Dep::Template3Dep()
:NDMaterial(0, ND_TAG_Template3Dep),ELS1(0), ELS2(0),ELS3(0), ELS4(0),
 ELT1(0), ELT2(0), ELT3(0), ELT4(0) 
{
    YS = new DPYieldSurface();       
    PS = new DPPotentialSurface();
    EPS = new EPState();
}


//================================================================================
// Destructor 
//================================================================================

Template3Dep::~Template3Dep()
{

}

////================================================================================
////copy constructor
////================================================================================
//Template3Dep::Template3Dep(const  Template3Dep & rval) {   
//
//    YS = rval.YS->newObj();
//    PS = rval.PS->newObj();
//    EPS = rval.EPS->newObj();
//    
//    // Scalar Evolution Laws
//    if ( rval.getELS1() ) 
//       ELS1  = rval.getELS1()->newObj();
//    else
//       ELS1 = 0;
//
//    if ( !rval.getELS2() ) 
//       ELS2 = 0;
//    else
//       ELS2  = rval.getELS2()->newObj();
//    
//    if ( !rval.getELS3() ) 
//       ELS3 = 0;
//    else
//       ELS3  = rval.getELS3()->newObj();
//    
//    if ( !rval.getELS4() ) 
//       ELS4 = 0;
//    else
//       ELS4  = rval.getELS4()->newObj();
//    
//    // Tensorial Evolution Laws
//    if ( rval.getELT1() ) 
//       ELT1  = rval.getELT1()->newObj();
//    else
//       ELT1 = 0;
//
//    if ( !rval.getELT2() ) 
//       ELT2 = 0;
//    else
//       ELT2  = rval.getELT2()->newObj();
//
//    if ( !rval.getELT3() ) 
//       ELT3 = 0;
//    else
//       ELT3  = rval.getELT3()->newObj();
//
//    if ( !rval.getELT4() ) 
//       ELT4 = 0;
//    else
//       ELT4  = rval.getELT4()->newObj();
//
//}    	   
//    	   



//================================================================================
// Routine used to generate elastic compliance tensor D for this material point
//================================================================================
tensor Template3Dep::ElasticComplianceTensor(void) const
  {
    tensor ret( 4, def_dim_4, 0.0 );

    double Ey = this->EPS->getE();
    if (Ey == 0) {
      cout << endln << "Ey = 0! Can't give you D!!" << endln;
      exit(1);
    }
    double nuP =this->EPS->getnu();
    
    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);

    tensor I_ijkl = I2("ij")*I2("kl");
    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    
    // Building compliance tensor
    ret = (I_ijkl * (-nuP/Ey)) + (I4s * ((1.0+nuP)/Ey));

    return ret;

}

//================================================================================
// Routine used to generate elastic stiffness tensor D for this material point
//================================================================================
tensor Template3Dep::ElasticStiffnessTensor(void) const
  {
    tensor ret( 4, def_dim_4, 0.0 );

    double E = this->EPS->getE();
    double nu =this->EPS->getnu();
    
    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);

    tensor I_ijkl = I2("ij")*I2("kl");


    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;

    //double x = I4s.trace();
    //cout << "xxxxxx " << x << endln;

    //I4s.null_indices();

    // Building elasticity tensor
    ret = I_ijkl*( E*nu / ( (1.0+nu)*(1.0 - 2.0*nu) ) ) + I4s*( E / (1.0 + nu) );
    //ret.null_indices();

    return ret;


}

//================================================================================
int Template3Dep::setTrialStrain(const Vector &v)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int Template3Dep::setTrialStrain(const Vector &v, const Vector &r)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int Template3Dep::setTrialStrainIncr(const Vector &v)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int Template3Dep::setTrialStrainIncr(const Vector &v, const Vector &r)
{
    // Not yet implemented
    return 0;
}

//================================================================================
const Matrix& Template3Dep::getTangent(void)
{
    // Not yet implemented
    Matrix *M = new Matrix();
    return *M;
}

//================================================================================
const Vector& Template3Dep::getStress(void)
{
    // Not yet implemented
    Vector *V = new Vector();
    return *V;
}

//================================================================================
const Vector& Template3Dep::getStrain(void)
{
    // Not yet implemented
    Vector *V = new Vector();
    return *V;
}

//================================================================================
int Template3Dep::setTrialStrain(const Tensor &v)
{
    EPS->setStrain(v);
    return 0;
}


//================================================================================
int Template3Dep::setTrialStrain(const Tensor &v, const Tensor &r)
{
    EPS->setStrain(v);
    return 0;
}

//================================================================================
//??????????? what is the trial strain? Initial strain?

int Template3Dep::setTrialStrainIncr(const Tensor &v)
{
    EPS->setStrain(v + EPS->getStrain() );
    return 0;
}

//================================================================================
//??????????? what is the TrialStrainIncr? strain increment?
int Template3Dep::setTrialStrainIncr(const Tensor &v, const Tensor &r)
{
    EPS->setStrain(v + EPS->getStrain() );
    return 0;
}

//================================================================================
const tensor& Template3Dep::getTangentTensor(void)
{
    return EPS->getEep();
}

//================================================================================
const Tensor& Template3Dep::getStressTensor(void)
{
    return EPS->getStress();
}

//================================================================================
const Tensor& Template3Dep::getStrainTensor(void)
{
    return EPS->getStrain();
}

//================================================================================
int Template3Dep::commitState(void)
{
	// Nothing to commit
	return 0;
}

//================================================================================
int Template3Dep::revertToLastCommit(void)
{
	// Nothing to commit
	return 0;
}
//================================================================================
int Template3Dep::revertToStart(void)
{
	// Nothing to commit
	return 0;
}

//================================================================================
NDMaterial *Template3Dep::getCopy(void)
{
    Template3Dep * tmp = 
            new Template3Dep( this->getTag()  ,
			      this->getYS()   ,
			      this->getPS()   ,
			      this->getEPS()  ,
			      this->getELS1() ,
			      this->getELS2() ,
			      this->getELS3() ,
			      this->getELS4() ,
			      this->getELT1() ,
			      this->getELT2() ,
			      this->getELT3() ,
			      this->getELT4() );
			      
    return tmp;

}


//================================================================================
NDMaterial *Template3Dep::getCopy(const char *code)
{
    if (strcmp(code,"Template3Dep") == 0)
    {
       Template3Dep * tmp = 
            new Template3Dep( this->getTag()  ,
			      this->getYS()   ,
			      this->getPS()   ,
			      this->getEPS()  ,
			      this->getELS1() ,
			      this->getELS2() ,
			      this->getELS3() ,
			      this->getELS4() ,
			      this->getELT1() ,
			      this->getELT2() ,
			      this->getELT3() ,
			      this->getELT4() );

       tmp->EPS->setZero();			      
       return tmp;
    }
    else
    {
	g3ErrorHandler->fatal("Template3Dep::getCopy failed to get model %s", code);
	return 0;
    }

   //Might need this to fool some compilor
   //return tmp;

}

//================================================================================
const char *Template3Dep::getType(void) const
{
    return "Template3Dep";
}

//================================================================================
//??What is the Order?

int Template3Dep::getOrder(void) const
{
    return 9;
}

//================================================================================
int Template3Dep::sendSelf(int commitTag, Channel &theChannel)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int Template3Dep::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // Not yet implemented
    return 0;
}

//================================================================================
void
Template3Dep::Print(ostream &s, int flag)
{
     s << (*this);
}



//private utilities

//================================================================================
// Set new EPState
//================================================================================
void Template3Dep::setEPS( EPState & rval) {   

     EPS = rval.newObj();
}    	   

//================================================================================
// Get the Yield Surface
//================================================================================
YieldSurface * Template3Dep::getYS() const 
{
    return YS;
}


//================================================================================
// Get the Potential Surface
//================================================================================
PotentialSurface * Template3Dep::getPS() const 
{
    return PS; 
}

//================================================================================
// Get the EPState
//================================================================================
EPState * Template3Dep::getEPS() const 
{ 
    return EPS; 
}

//================================================================================
// Get the 1st Salar evolution law
//================================================================================

EvolutionLaw_S * Template3Dep::getELS1() const 
{ 
    return ELS1; 
}

//================================================================================
// Get the 2ndst Salar evolution law
//================================================================================
EvolutionLaw_S * Template3Dep::getELS2() const 
{ 
    return ELS2; 
}

//================================================================================
// Get the 2ndst Salar evolution law
//================================================================================
EvolutionLaw_S * Template3Dep::getELS3() const 
{ 
    return ELS3; 
}
//================================================================================
// Get the 2ndst Salar evolution law
//================================================================================
EvolutionLaw_S * Template3Dep::getELS4() const 
{ 
    return ELS4; 
}


//================================================================================
// Get the 1st tensorial evolution law
//================================================================================

EvolutionLaw_T * Template3Dep::getELT1() const 
{ 
    return ELT1; 
}

//================================================================================
// Get the 2nd tensorial evolution law
//================================================================================

EvolutionLaw_T * Template3Dep::getELT2() const 
{ 
    return ELT2; 
}
//================================================================================
// Get the 3rd tensorial evolution law
//================================================================================

EvolutionLaw_T * Template3Dep::getELT3() const 
{ 
    return ELT3; 
}
//================================================================================
// Get the 4th tensorial evolution law
//================================================================================

EvolutionLaw_T * Template3Dep::getELT4() const 
{ 
    return ELT4; 
}


//================================================================================
// Overloaded Insertion Operator
// prints an EPState's contents 
//================================================================================
ostream& operator<< (ostream& os, const Template3Dep & MP)
{
    os << endln << "Template3Dep: " << endln;
    os << "\ttag: " << MP.getTag() << endln;
    os << "=================================================================" << endln;
    MP.getYS()->print();
    MP.getPS()->print();
    MP.getEPS()->print();

    cout << endln << "Scalar Evolution Laws: " << endln; 
    if ( MP.ELS1 ){
       cout << "\nFor 1st scalar var:";
       MP.ELS1->print();
    }
    
    if ( MP.ELS2 ){
       cout << "\nFor 2nd scalar var:";
       MP.ELS2->print();
    }
    
    if ( MP.ELS3 ){
       cout << "\nFor 3rd scalar var:";
       MP.ELS3->print();
    }
    
    if ( MP.ELS4 ){
       cout << "\nFor 4th scalar var:";
       MP.ELS4->print();
    }
    

    cout << endln << "Tensorial Evolution Laws: " << endln; 
    if ( MP.ELT1 ){
       cout << "\nFor 1st tensorial var:";
       MP.ELT1->print();
    }
    if ( MP.ELT2 ){
       cout << "\nFor 2nd tensorial var:";
       MP.ELT2->print();
    }
    if ( MP.ELT3 ){
       cout << "\nFor 3rd tensorial var:";
       MP.ELT3->print();
    }
    if ( MP.ELT4 ){
       cout << "\nFor 4th tensorial var:";
       MP.ELT4->print();
    }

    os << endln;           
    return os;
}  


#endif

