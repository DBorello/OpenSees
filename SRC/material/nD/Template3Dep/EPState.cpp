/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
# CLASS:             DPEPState (the base class for all Elasto-plastic state)     #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              08-03-2000                                                  #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: The goal is to create a platform for efficient and easy     #
#                    implemetation of any elasto-plastic constitutive model!     #
#                                                                                #
//================================================================================
*/

#ifndef EPState_CPP
#define EPState_CPP

#include "EPState.h"

//================================================================================
//Normal Constructor
//================================================================================

EPState::EPState(double             Eod,
                 double             Ed,
                 double             nu,
                 double             rho,
                 const stresstensor &stressp,       
                 const straintensor &strainp, 
                 const straintensor &Estrainp,
                 const straintensor &Pstrainp,
                 const straintensor &dEstrainp,
                 const straintensor &dPstrainp,
	         int                NScalarp,
	         const double     * Scalarp,
	         int                NTensorp,
	         const tensor     * Tensorp ) {

      Eo               = Eod;	        
      E_Young          = Ed;	        
      nu_Poisson       = nu;	     
      rho_mass_density = rho; 
      CurrentStress    = stressp;
      //cout << "stressp " << stressp;
      //CurrentStress.null_indices();
      CurrentStrain    =  strainp;
      ElasticStrain    =  Estrainp;
      PlasticStrain    =  Pstrainp;
      dElasticStrain   =  dEstrainp;
      dPlasticStrain   =  dPstrainp;
      Eep = tensor( 4, def_dim_4, 0.0 ); // need to be initialized as 4th order tensor
      //cout << "strainp " << strainp;
      //CurrentStrain.null_indices();

      NScalarVar  =  NScalarp;
      ScalarVar = new double[ NScalarVar ]; 
      if ( !ScalarVar ) {
         ::printf("\n\nInsufficient memory for Scalar hardening vars\n");
         ::exit(1);  
      }
      for (int i = 0; i < NScalarp; i++) {
         //cout << Scalarp[i] << endln; 
	 ScalarVar[i] = Scalarp[i];
      }

      NTensorVar = NTensorp;
      TensorVar = new stresstensor[ NTensorVar ];
      if ( !TensorVar ) {
         ::printf("\n\nInsufficient memory for Tensor hardening vars\n");
         ::exit(1);  
      }
      for (int i = 0; i < NTensorVar; i++) {
	 //cout << Tensorp[i];
	 TensorVar[i] = Tensorp[i];
	 //cout << TensorVar[i];
 	 TensorVar[i].null_indices();
      }
}


//================================================================================
//Normal Constructor --no parameters
//================================================================================

EPState::EPState( ) {

      Eo               = 80000.0;
      E_Young          = 20000.0;
      nu_Poisson       = 0.3;	     
      rho_mass_density = 0.0; 
      Eep = tensor( 4, def_dim_4, 0.0 );

      NScalarVar = 1;
      ScalarVar = new double(0.0); 

      NTensorVar = 1;
      TensorVar = new stresstensor(1);
      //TensorVar = &temp;
}

//================================================================================
//create a clone of itself
//================================================================================

EPState* EPState::newObj() {

      EPState * eps = new  EPState(this->getEo(),	       
      				   this->getE(),
      				   this->getnu(),
      				   this->getrho(),
                                   this->getStress(),
      				   this->getStrain(), 
      				   this->getElasticStrain(), 
      				   this->getPlasticStrain(), 
      				   this->getdElasticStrain(), 
      				   this->getdPlasticStrain(), 
      				   this->getNScalarVar(), 
      				   this->getScalarVar(),
      				   this->getNTensorVar(),
      				   this->getTensorVar() );
      return eps;
}      				 


//================================================================================
// Copy constructor
//================================================================================
EPState::EPState( const EPState &rhs ) {         

      Eo               = rhs.getEo();	        
      E_Young          = rhs.getE();	        
      nu_Poisson       = rhs.getnu();	     
      rho_mass_density = rhs.getrho(); 
      CurrentStress    = rhs.getStress();
      CurrentStrain    = rhs.getStrain();
      ElasticStrain    = rhs.getElasticStrain();
      PlasticStrain    = rhs.getPlasticStrain();
      dElasticStrain   = rhs.getdElasticStrain();
      dPlasticStrain   = rhs.getdPlasticStrain();
      
      //cout << Eep.rank() << " ";
      //Eep.printshort("before copy con ");
      Eep = rhs.getEep();
      //cout << Eep.rank() << " ";
      //Eep.printshort("after copy con ");

      NScalarVar  =  rhs.getNScalarVar();
      ScalarVar = new double[ NScalarVar ]; 
      if ( !ScalarVar ) {
         ::printf("\n\nInsufficient memory for Scalar hardening vars\n");
         ::exit(1);  
      }
      for (int i = 0; i < NScalarVar; i++) 
	 ScalarVar[i] = rhs.ScalarVar[ i ];

      NTensorVar = rhs.getNTensorVar();
      TensorVar = new stresstensor[ NTensorVar ];
      if ( !TensorVar ) {
         ::printf("\n\nInsufficient memory for Tensor hardening vars\n");
         ::exit(1);  
      }
      for (int i = 0; i < NTensorVar; i++) {
	 TensorVar[i] = rhs.TensorVar[ i ];
	 //cout << TensorVar[i];
 	 TensorVar[i].null_indices();
      }
     
}      				 


//================================================================================
//Overloading the assignment sign
//================================================================================
const EPState & EPState::operator=(const EPState &rhs ) { 
        
      if ( this != &rhs ) {
         Eo               = rhs.getEo();
         E_Young          = rhs.getE();	        
         nu_Poisson       = rhs.getnu();	     
         rho_mass_density = rhs.getrho(); 
         CurrentStress    = rhs.getStress();
         //cout << "Current stress " << CurrentStress;
         CurrentStrain    = rhs.getStrain();
         //cout << "strainp " << strainp;
         ElasticStrain    = rhs.getElasticStrain();
         PlasticStrain    = rhs.getPlasticStrain();
         dElasticStrain   = rhs.getdElasticStrain();
         dPlasticStrain   = rhs.getdPlasticStrain();
         Eep              = rhs.getEep();
         
         NScalarVar  =  rhs.getNScalarVar();
         ScalarVar = new double[ NScalarVar ]; 
         if ( !ScalarVar ) {
            ::printf("\n\nInsufficient memory for Scalar hardening vars\n");
            ::exit(1);  
         }
         for (int i = 0; i < NScalarVar; i++) {
            //cout << Scalarp[i] << endln; 
         	 ScalarVar[i] = rhs.getScalarVar(i+1);
         }
         
         NTensorVar = rhs.getNTensorVar();
         TensorVar = new stresstensor[ NTensorVar ];
         if ( !TensorVar ) {
            ::printf("\n\nInsufficient memory for Tensor hardening vars\n");
            ::exit(1);  
         }
         for (int i = 0; i < NTensorVar; i++) {
             TensorVar[i] = rhs.getTensorVar(i+1);
             TensorVar[i].null_indices();
         }
      }
      
      return *this;
      
}      				 

//================================================================================
double EPState::getE() const {
      return E_Young; 
}

//================================================================================
double EPState::getEo() const {
      return Eo; 
}

//================================================================================
double EPState::getnu() const {
      return nu_Poisson; 
}

//================================================================================
double EPState::getrho() const {
      return rho_mass_density; 
};

//================================================================================
int    EPState::getNScalarVar() const {
      return NScalarVar; 
}

//================================================================================
int    EPState::getNTensorVar() const {
      return NTensorVar; 
}

//================================================================================
stresstensor EPState::getStress() const {

     return CurrentStress;

}

//================================================================================
straintensor EPState::getStrain() const {

     return CurrentStrain;

}


//================================================================================
straintensor EPState::getElasticStrain() const {

     return ElasticStrain;

}

//================================================================================
straintensor EPState::getPlasticStrain() const {

     return PlasticStrain;

}
//================================================================================
straintensor EPState::getdElasticStrain() const {

     return dElasticStrain;

}


//================================================================================
straintensor EPState::getdPlasticStrain() const {

     return dPlasticStrain;

}


//================================================================================
tensor EPState::getEep() const {

     return Eep;

}


//================================================================================
void EPState::setE( double Ey ) { 
      E_Young = Ey; 
}


//================================================================================
void EPState::setStress(const stresstensor &newstress ) { 
      CurrentStress = newstress; 
}


//================================================================================
void EPState::setStrain(const straintensor &newstrain ) {

      CurrentStrain = newstrain; 

}

//================================================================================
void EPState::setElasticStrain(const straintensor &newEstrain ) {

      ElasticStrain = newEstrain; 

}

//================================================================================
void EPState::setPlasticStrain(const straintensor &newPstrain ) {

      PlasticStrain = newPstrain; 

}


//================================================================================
void EPState::setdElasticStrain(const straintensor &newdEstrain ) {

      dElasticStrain = newdEstrain; 

}

//================================================================================
void EPState::setdPlasticStrain(const straintensor &newdPstrain ) {

      dPlasticStrain = newdPstrain; 

}

//================================================================================
void EPState::setEep(const tensor &newEep )  {

     Eep = newEep;

}


//================================================================================
// Retrun the nth Scalar Var.... Starting from 1!!
//================================================================================
double EPState::getScalarVar( int WhichOne) const { 

      if ( WhichOne <= getNScalarVar() )  
         return ScalarVar[ WhichOne - 1 ]; 
      else 
      {
         cout << " Out of ScalarVar's range!";	  
	 exit(1);
      }

}


//================================================================================
// Retrun the nth Tensor Var.... Starting from 1!!
//================================================================================
stresstensor EPState::getTensorVar(int WhichOne) const { 

      if ( WhichOne <= getNTensorVar() )  
         return TensorVar[ WhichOne - 1 ]; 
      else 
      {
         cout << " Out of TensorVar's range!";	  
	 exit(1);
      }


}

//================================================================================
// Retrun Scalar pointer 
//================================================================================
double * EPState::getScalarVar() const {
      
      return ScalarVar; 
      
}

//================================================================================
// Retrun Tensor pointer 
//================================================================================
tensor * EPState::getTensorVar() const { 
    
      return TensorVar; 

}


//================================================================================
// set nth Scalar Var.... Starting from 1!!
//================================================================================
void EPState::setScalarVar(int WhichOne, double rval) { 

      if ( WhichOne <= getNScalarVar() )  
         ScalarVar[ WhichOne - 1 ] = rval; 
      else 
      {
         cout << " Out of ScalarVar's range!";	  
	 exit(1);
      }
      
      
}

//================================================================================
// set nth Tensor Var.... Starting from 1!!
//================================================================================
void EPState::setTensorVar(int WhichOne, const stresstensor &rval) { 

      if ( WhichOne <= getNTensorVar() )  
         TensorVar[ WhichOne - 1 ] = rval;
      else 
      {
         cout << " Out of Tensor Var's range!";	  
	 exit(1);
      }

}

//================================================================================
void EPState::print() { 
      cout << *this;

}


//================================================================================
// Set all state variables to initials

void EPState::setZero() { 
      
      Eep = tensor( 4, def_dim_4, 0.0 );

      // not finished yet! set internal vars to initial value or zero?
      stresstensor Stress0;
      straintensor Strain0;

      CurrentStress = Stress0;
      CurrentStrain = Strain0;
      ElasticStrain = Strain0;
      PlasticStrain = Strain0;
      dElasticStrain = Strain0;
      dPlasticStrain = Strain0;
			       
}


#endif

