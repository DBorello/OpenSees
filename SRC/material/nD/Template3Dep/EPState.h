/*
################################################################################
# COPYRIGHT (C):     :-))                                                      #
# PROJECT:           Object Oriented Finite Element Program                    #
# PURPOSE:           General platform for elaso-plastic constitutive model     #
#                    implementation                                            #
# CLASS:             DPEPState (the base class for all Elasto-plastic state)   #
#                                                                              #
# VERSION:                                                                     #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
# TARGET OS:         DOS || UNIX || . . .                                      #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                               #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
#                                                                              #
#                                                                              #
# DATE:              08-03-2000                                                #
# UPDATE HISTORY:                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
# SHORT EXPLANATION: The goal is to create a platform for efficient and easy   #
#                    implemetation of any elasto-plastic constitutive model!   #
#                                                                              #
################################################################################
*/

#ifndef EPState_H
#define EPState_H

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>

#include <iostream.h>
#include <iomanip.h>
#define endln "\n"

class EPState
{
  // Elastic parameters
  public:
    double Eo;	                  // Young's modulus when p = p_atmosphere  
    double E_Young;	          // Young's modulus  			    
    double nu_Poisson;	          // Poisson's ratio  			    
    double rho_mass_density;      // Mass density     			    
    
    stresstensor CurrentStress;   // Current stress  --total                
    straintensor CurrentStrain;	  // Current strain  --total 	            
    straintensor ElasticStrain;	  // Current cumulative elastic strain      
    straintensor PlasticStrain;	  // Current cumulative plastic strain      
    straintensor dElasticStrain;  // Current elastic strain increment       
    straintensor dPlasticStrain;  // Current plastic strain increment       
    tensor       Eep;             // Current Elastic plastic stifness tensor
									    
    int NScalarVar;	 //Number of scalar hardening vars		    
    int NTensorVar;	 //Number of tensor hardening vars		    
    									    
    double *ScalarVar;   // scalar variable array for scalar hardening vars 
    tensor *TensorVar;	 // tensor variable array for tensor hardening vars 
    
  public:
    EPState();                         //Normal Constructor--no parameters
    EPState(double              Eod,    //Normal Constructor
            double              Ed,
            double              nu,
	    double              rho,
            const stresstensor  &stressp,       
            const straintensor  &strainp,
            const straintensor  &Estrainp, 
            const straintensor  &Pstrainp,  
            const straintensor  &dEstrainp, 
            const straintensor  &dPstrainp,  
	    int                 NScalarp,
	    const double       *Scalarp,
	    int                 NTensorp,
	    const tensor       *Tensorp);

    EPState *newObj();                 //create a clone of itself
    EPState( const EPState &rhs );     // Copy constructor    
    const EPState & EPState::operator=(const EPState &rhs );  //Overloading assignment

    double getE() const;
    double getEo() const;
    double getnu() const;
    double getrho() const;
    int getNScalarVar() const;
    int getNTensorVar() const;

    stresstensor getStress() const;
    straintensor getStrain() const;
    straintensor getElasticStrain() const;
    straintensor getPlasticStrain() const;
    straintensor getdElasticStrain() const;
    straintensor getdPlasticStrain() const;
    tensor getEep() const;
 
    void setE( double Ey );
    void setStress( const stresstensor &newstress );
    void setStrain( const straintensor &newstrain );
    void setElasticStrain( const straintensor &newstrain );
    void setPlasticStrain( const straintensor &newstrain );
    void setdElasticStrain( const straintensor &newstrain );
    void setdPlasticStrain( const straintensor &newstrain );
    void setEep( const tensor &);

    // WhichOne starts from 1 and ends up to  NScalarVar
    double getScalarVar( int WhichOne) const;
    stresstensor getTensorVar(int WhichOne) const;

    // return the pointers
    double *getScalarVar() const;
    tensor *getTensorVar() const;

    // WhichOne starts from 1 and ends up to NTensorVar
    void setScalarVar(int WhichOne, double rval);
    void setTensorVar(int WhichOne, const stresstensor &rval);
    void setZero();

    void print();

    //================================================================================
    // Overloaded Insertion Operator
    // prints an EPState's contents 
    //================================================================================
    // this is here because g++
    //
    // cml00-newfem4md 25>g++ -v
    // Reading specs from /usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66/specs
    // gcc version egcs-2.91.66 19990314/Linux (egcs-1.1.2 release)
    //
    // is complaining if it is in .cc ???
    
    friend ostream & operator<< (ostream& os, const EPState & EPS)
    {
        os.setf( ios::showpos | ios::scientific);
        os.precision(4);
        os.width(10);       
        os << endln << "Elastic plastic state parameters: "  << endln;

        //os.width(10);       
        os << "Eo = " << EPS.getEo() << ";";
        os << " E_Young = " << EPS.getE() << ";";
        //os.width(10);       
	os << " nu_Poisson = " << EPS.getnu() << ";";
	//os.width(10);       
	os << " rho = " << EPS.getrho() << endln;

        //os.width(10);       
        os << endln << "Current Stress:" << EPS.getStress() << endln;

	os << "Current Strain:" << EPS.getStrain() << endln;
	os << "ElasticStrain :" << EPS.getElasticStrain() << endln;
	os << "PlasticStrain :" << EPS.getPlasticStrain() << endln;
	os << "dElasticStrain:" << EPS.getdElasticStrain() << endln;
	os << "dPlasticStrain:" << EPS.getdPlasticStrain() << endln;
	os << "Eep.rank():" << EPS.getEep().rank() << endln;

        os.unsetf( ios::showpos );
	int NS = EPS.getNScalarVar();
	int NT = EPS.getNTensorVar();
	
	os << endln << "NScalarVar = " << NS << endln; 
    
        for (int i = 0; i < NS; i++) {
            os << "No." << i+1 << " " << EPS.ScalarVar[i] << "; ";
	}
        os << endln << endln;
    
        os << "NTensorVar = " << NT;
        for (int i = 0; i < NT; i++) {
           os.unsetf( ios::showpos);
           os << endln << "No." << i+1 << " tensorial var:";
           os.setf( ios::showpos);
           cout << EPS.TensorVar[i];
        }

        os << endln;           
        return os;
    }  


};


#endif

