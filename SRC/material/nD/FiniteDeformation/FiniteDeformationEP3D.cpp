//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              July 2004
//# UPDATE HISTORY:
//#
//===============================================================================

#ifndef FiniteDeformationEP3D_CPP
#define FiniteDefornationEP3D_CPP

#include "FiniteDeformationEP3D.h"

const int    Max_Iter  = 40;
const double tolerance = 1.0e-8;

tensor tensorZ2(2, def_dim_2, 0.0);
tensor tensorI2("I", 2, def_dim_2);
tensor tensorZ4(4, def_dim_4, 0.0);

//Constructor 00----------------------------------------------------------------------------------
FiniteDeformationEP3D::FiniteDeformationEP3D()
:NDMaterial(0, 0)
{
	fde3d = 0;
	fdy = 0;
	fdf = 0;
	fdEvolutionS = 0;
	fdEvolutionT = 0;
	fdeps = 0;

        int err;  
	err = this->revertToStart();	
}

// Constructor 01---------------------------------------------------------------------------------
FiniteDeformationEP3D::FiniteDeformationEP3D(int tag,
	                                     NDMaterial *fde3d_in,
	                                     fdYield *fdy_in,
	                                     fdFlow *fdf_in,
	                                     fdEvolution_S *fdEvolutionS_in,
					     fdEvolution_T *fdEvolutionT_in)
:NDMaterial(tag, ND_TAG_FiniteDeformationEP3D)
{
	if (fde3d_in)
	  fde3d = fde3d_in->getCopy();
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdElastic3D\n";
	  exit (-1);
	}
	
	if (fdy_in)
	  fdy = fdy_in->newObj(); 
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdYield\n";
	  exit (-1);
	}	  	
	
	if (fdf_in)
	  fdf = fdf_in->newObj(); 
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdFlow\n";
	  exit (-1);
	}
	  		
	if (fdEvolutionS_in)
	  fdEvolutionS = fdEvolutionS_in->newObj(); 
	else
	  fdEvolutionS = 0;

	if (fdEvolutionT_in)
	  fdEvolutionT = fdEvolutionT_in->newObj(); 
	else
	  fdEvolutionT = 0;
	  
	//if (fdeps == 0)
	fdeps = new FDEPState(); 
	  			
}	  	  

// Constructor 02---------------------------------------------------------------------------------
FiniteDeformationEP3D::FiniteDeformationEP3D(int tag,
	                                     NDMaterial *fde3d_in,
	                                     fdYield *fdy_in,
	                                     fdFlow *fdf_in,
	                                     fdEvolution_S *fdEvolutionS_in)
:NDMaterial(tag, ND_TAG_FiniteDeformationEP3D)
{
	if (fde3d_in)
	  fde3d = fde3d_in->getCopy();
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdElastic3D\n";
	  exit (-1);
	}
	
	if (fdy_in)
	  fdy = fdy_in->newObj(); 
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdYield\n";
	  exit (-1);
	}	  	
	
	if (fdf_in)
	  fdf = fdf_in->newObj(); 
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdFlow\n";
	  exit (-1);
	}
	  		
	if (fdEvolutionS_in)
	  fdEvolutionS = fdEvolutionS_in->newObj(); 
	else
	  fdEvolutionS = 0;

	fdEvolutionT = 0;
	  
	//if (fdeps == 0)
	fdeps = new FDEPState(); 
	  			
}	  

// Constructor 03---------------------------------------------------------------------------------
FiniteDeformationEP3D::FiniteDeformationEP3D(int tag,
	                                     NDMaterial *fde3d_in,
	                                     fdYield *fdy_in,
	                                     fdFlow *fdf_in,
					     fdEvolution_T *fdEvolutionT_in)
:NDMaterial(tag, ND_TAG_FiniteDeformationEP3D)
{
	if (fde3d_in)
	  fde3d = fde3d_in->getCopy();
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdElastic3D\n";
	  exit (-1);
	}
	
	if (fdy_in)
	  fdy = fdy_in->newObj(); 
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdYield\n";
	  exit (-1);
	}	  	
	
	if (fdf_in)
	  fdf = fdf_in->newObj(); 
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdFlow\n";
	  exit (-1);
	}
	  		
	fdEvolutionS = 0;

	if (fdEvolutionT_in)
	  fdEvolutionT = fdEvolutionT_in->newObj(); 
	else
	  fdEvolutionT = 0;
	  
	//if (fdeps == 0)
	fdeps = new FDEPState(); 
	  			
}	

// Constructor 04---------------------------------------------------------------------------------
FiniteDeformationEP3D::FiniteDeformationEP3D(int tag,
	                                     NDMaterial *fde3d_in,
	                                     fdYield *fdy_in,
	                                     fdFlow *fdf_in)
:NDMaterial(tag, ND_TAG_FiniteDeformationEP3D)
{
	if (fde3d_in)
	  fde3d = fde3d_in->getCopy();
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdElastic3D\n";
	  exit (-1);
	}
	
	if (fdy_in)
	  fdy = fdy_in->newObj(); 
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdYield\n";
	  exit (-1);
	}	  	
	
	if (fdf_in)
	  fdf = fdf_in->newObj(); 
	else {
	  opserr << "FiniteDeformationEP3D:: FiniteDeformationEP3D failed to construct the fdFlow\n";
	  exit (-1);
	}
	  		
	fdEvolutionS = 0;

	fdEvolutionT = 0;
	  
	//if (fdeps == 0)
	fdeps = new FDEPState(); 
	  			
}	

// Destructor-------------------------------------------------------------------------------------
FiniteDeformationEP3D::~FiniteDeformationEP3D()
{
	if (fde3d)
	  delete fde3d;
	
	if (fdy)
	  delete fdy;	
	
	if (fdf)
	  delete fdf;	
	
	if (fdEvolutionS)
	  delete fdEvolutionS;

	if (fdEvolutionT)
	  delete fdEvolutionT;

	if (fdeps)
	  delete fdeps;	  	   	  
}

//----------------------------------------------------------------------
double FiniteDeformationEP3D::getRho(void)
{
	return fde3d->getRho();
}

//----------------------------------------------------------------------
int FiniteDeformationEP3D::setTrialF(const straintensor &f)
{
        F = f;

	return ImplicitAlgorithm();
}

//----------------------------------------------------------------------
int FiniteDeformationEP3D::setTrialFIncr(const straintensor &df)
{
        return this->setTrialF(this->getF() + df);
}

//----------------------------------------------------------------------
const tensor& FiniteDeformationEP3D::getTangentTensor(void)
{
	return iniTangent;
}

//----------------------------------------------------------------------
const straintensor FiniteDeformationEP3D::getStrainTensor(void)
{
	return iniGreen;
}

//----------------------------------------------------------------------
const stresstensor FiniteDeformationEP3D::getStressTensor(void)
{
	return iniPK2;
}

//----------------------------------------------------------------------
const straintensor FiniteDeformationEP3D::getF(void)
{
	return F;
}

//----------------------------------------------------------------------
const straintensor FiniteDeformationEP3D::getFp(void)
{
	return fdeps->getFpInVar();
}

//----------------------------------------------------------------------
int FiniteDeformationEP3D::commitState(void)
{
	return fdeps->commitState();
}

//----------------------------------------------------------------------
int FiniteDeformationEP3D::revertToLastCommit(void)
{	
	return fdeps->revertToLastCommit();
}

//----------------------------------------------------------------------
int FiniteDeformationEP3D::revertToStart(void)
{	
	return fdeps->revertToStart();
}

//----------------------------------------------------------------------
NDMaterial* FiniteDeformationEP3D::getCopy (void)
{
	NDMaterial* tmp = new FiniteDeformationEP3D(
	  this->getTag(),
	  this->getFDE3D(),
	  this->getFDY(),
	  this->getFDF(),
	  this->getFDEvolutionS(),
	  this->getFDEvolutionT() );
	
	return tmp;                                                     
}

//----------------------------------------------------------------------
NDMaterial* FiniteDeformationEP3D::getCopy (const char *code)
{
	if ( strcmp(code,"FiniteDeformationEP3D") == 0
	     || strcmp(code,"FDEP3D") == 0 )
	{
	  NDMaterial* tmp = new FiniteDeformationEP3D (
	  this->getTag(),
	  this->getFDE3D(),
	  this->getFDY(),
	  this->getFDF(),
	  this->getFDEvolutionS(),
	  this->getFDEvolutionT() );
	
	  return tmp;
	}
	else
	{
	  opserr << "FiniteDeformationEP3D::getCopy fainled:" << code << "\n";
	  exit(-1); 
	}                                                   
}

//----------------------------------------------------------------------  
const char* FiniteDeformationEP3D::getType (void) const
{
	return "ThreeDimentionalFD";
}
  

//----------------------------------------------------------------------
int FiniteDeformationEP3D::getOrder (void) const
{
	return 6;	
}

//----------------------------------------------------------------------
int FiniteDeformationEP3D::sendSelf(int commitTag, Channel &theChannel)
{
	// Not yet implemented
	return 0;	
}

//----------------------------------------------------------------------
int FiniteDeformationEP3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// Not yet implemented
	return 0;	
}

//----------------------------------------------------------------------
void FiniteDeformationEP3D::Print(OPS_Stream &s, int flag)
{
	s << (*this);
}

//----------------------------------------------------------------------
const stresstensor FiniteDeformationEP3D::getCauchyStressTensor(void)
{
	return cauchystress;
}
	
//----------------------------------------------------------------------
NDMaterial * FiniteDeformationEP3D::getFDE3D() const
{
	return fde3d;
}

//----------------------------------------------------------------------
fdYield * FiniteDeformationEP3D::getFDY() const
{
	return fdy;
}

//----------------------------------------------------------------------
fdFlow * FiniteDeformationEP3D::getFDF() const
{
	return fdf;
}

//----------------------------------------------------------------------
fdEvolution_S * FiniteDeformationEP3D::getFDEvolutionS() const
{
	return fdEvolutionS;
}

//----------------------------------------------------------------------
fdEvolution_T * FiniteDeformationEP3D::getFDEvolutionT() const
{
	return fdEvolutionT;
}

//----------------------------------------------------------------------
FDEPState * FiniteDeformationEP3D::getFDEPState() const
{
	return fdeps;
}


//----------------------------------------------------------------------
//
//----------------------------------------------------------------------
//Mandel Version
int FiniteDeformationEP3D::ImplicitAlgorithm()
{
    // *** This is the key function ***


    // Initializing
    double yieldfun = 0.0;              // trial value of yield function
    double D_gamma  = 0.0;              // consistency parameter /Delta_{gamma}
    double d_gamma  = 0.0;              // increment of consistency parameter /delta_{gamma}
    int	 iter_counter = 0;

    straintensor res_Ee = tensorZ2;     // residual of intermediate Ee
    straintensor res_eta = tensorZ2;    // norm of residual of eta

    double  res_xi = 0.0;	        // residual of strain like internal variable
    double  res_norm_eta = 0.0;         // norm of residual of eta    
    double  res_norm_Ee  = 0.0;         // norm of residual of intermediate Ee
    double  residual = 0.0;	        // residual of Ee and xi and eta

    straintensor Fp;	                // Plastic deformation gradient
    double xi = 0.0;		        // strain like internal isotropic variable
    double q = 0.0;		        // stress like internal isotropic variable
    straintensor eta;		        // strain like internal kinematic variable
    stresstensor a;		        // stress like internal kinematic variable

    straintensor Fpinv = tensorI2;      // inverse of Fp
    straintensor Ce = tensorI2;	        // intermediate C
    straintensor Ee = tensorZ2;	        // intermediate Ee
    
    straintensor Ee_n = tensorZ2;       // Ee at the incremental step n, calculated from Fp_n
    double xi_n;			// xi at the incremental step n, known
    straintensor eta_n = tensorZ2;	// eta at the incremental step n, known
        
    stresstensor Mtensor = tensorZ2;    // --> dFl/dT
    stresstensor symMtensor = tensorZ2; // --> dFl/dS
    tensor Ltensor = tensorZ4 ;         // Tangent tensor in the intermediate configuration
    tensor LEPtensor = tensorZ4;	// Consistent tangent tensor in the reference configuration
        
    double nscalar = 0.0;		// --> dFl/d(xi)
    double Kscalar = 0.0;		// Isotropic hardening modoulus
    straintensor ntensor;		// --> dFl/d(eta)
    tensor Ktensor = tensorZ4;	        // Kinematic hardening modoulus
    
    stresstensor dfods = tensorZ2 ;     // --> dY/d(stress)
    double dfodq = 0.0;			// --> dY/d(xi)
    stresstensor dfoda = tensorZ2 ;     // --> dY/d(eta)

    tensor dMods = tensorZ4;		// --> dM/d(stress), from d2Fl/d(stress)d(stress)    
    tensor dsymMods = tensorZ4;		// --> dMs/d(stress)

    double dnsodq = 0.0;	        // --> d2Fl/dqdq
    tensor dntoda = tensorZ4;           // --> d2Fl/dada
    
    straintensor D_Ee = tensorZ2;
    double D_xi  = 0.0;
    straintensor D_eta = tensorZ2;

    tensor A11 = tensorZ4;  tensor A12 = tensorZ2;  tensor A13 = tensorZ4; 		
    tensor A21 = tensorZ2;  double a22 = 0.0;	    tensor A23 = tensorZ2;
    tensor A31 = tensorZ4;  tensor A32 = tensorZ2;  tensor A33 = tensorZ4;   

    BJmatrix C99(9, 9, 0.0);
    BJmatrix CCC(19, 19, 0.0);

    tensor tensorTemp1;
    tensor tensorTemp2;
    tensor tensorTemp3;
    tensor tensorTemp4;
    tensor tensorTemp5;
    tensor tensorTemp6;
    tensor tensorTemp7;
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    double temp4 = 0.0;
    double temp5 = 0.0;

    tensor LM = tensorZ4;	         // For Mandel Tangent Stiffness
    stresstensor  B_Mandel = tensorZ2;   // For Mandel stress

    tensorTemp3 = tensorI2("ij")*tensorI2("kl");
      tensorTemp3.null_indices();
    tensor tensorI4 = tensorTemp3.transpose0110();
    
    // Read the previous incremental step history variables
    Fp = fdeps->getCommitedFpInVar();
    Fpinv = Fp.inverse();
    Fe = F("ij")*Fpinv("jk");   Fe.null_indices();
    Ce = Fe("ki")*Fe("kj");     Ce.null_indices();
    Ee = (Ce - tensorI2) * 0.5;
    Ee_n = Ee;

    xi = fdeps->getCommitedStrainLikeInVar();
    xi_n = xi;
    q = fdeps->getCommitedStressLikeInVar();

    if ( fdEvolutionT ) {
      eta = fdeps->getCommitedStrainLikeKiVar();
      eta_n = eta;
      a = fdeps->getCommitedStressLikeKiVar();
    }
    
    // Return stress from finite deformation elastic model
    fde3d->setTrialC(Ce);         // Note: It is C, not F!!!
    B_PK2 = fde3d->getStressTensor();

    B_Mandel = Ce("ik")*B_PK2("kj");	// Mandel Stress
      B_Mandel.null_indices();

    // Evaluate the value of yield function
    yieldfun = fdy->Yd(B_Mandel, *fdeps);
    //printf("\nY0= %e\n", yieldfun);

    if ( yieldfun > (fdy->getTolerance()) ) { // start of plastic part
      D_gamma = 0.0;
      d_gamma = 0.0;
      iter_counter = 0;
        
      do {             // beginning of do - while

	// Return Symmetric tensor M and n_scalar and n_tensor
        Mtensor = fdf->dFods(B_Mandel, *fdeps);
        Mtensor = Ce("ik")*Mtensor("kj");
	  Mtensor.null_indices();
        symMtensor = (Mtensor + Mtensor.transpose11()) * 0.5;
    
        // Return tangent variables
	Ltensor = fde3d->getTangentTensor();
	
	tensorTemp1 = tensorI2("ij")*B_PK2("mn");
	  tensorTemp1.null_indices();
	tensorTemp2 = Ce("ik")*Ltensor("kjmn");
	  tensorTemp2.null_indices();
        LM = tensorTemp2 + tensorTemp1.transpose0110() + tensorTemp1.transpose0111();  

	if ( fdEvolutionS) {
	  nscalar = fdf->dFodq(B_Mandel, *fdeps);
	  Kscalar = fdEvolutionS->HModulus(B_Mandel, *fdeps);
	  dfodq = fdy->dYodq(B_Mandel, *fdeps);
	  dnsodq = fdf->d2Fodqdq(B_Mandel, *fdeps);
	}

	if ( fdEvolutionT) {
	  ntensor = fdf->dFoda(B_Mandel, *fdeps);
	  Ktensor = fdEvolutionT->HModulus(B_Mandel, *fdeps);
	  dfoda = fdy->dYoda(B_Mandel, *fdeps);
	  dntoda = fdf->d2Fodada(B_Mandel, *fdeps);
	}

        // Return ...
        dfods = fdy->dYods(B_Mandel, *fdeps);        
        dMods = fdf->d2Fodsds(B_Mandel, *fdeps);
        dsymMods = Ce("ik")*dMods("kjmn");
	  dsymMods.null_indices();
        dsymMods = (dsymMods + dsymMods.transpose1100()) * 0.5;
        

        // Return tensor Axx
        A11 = (dsymMods *D_gamma) + tensorI4;

	a22 = (dnsodq *D_gamma) + 1.0;
	
	A12 = fdf->d2Fodsdq(B_Mandel, *fdeps);
        A12 = Ce("ik")*A12("kj");
	  A12.null_indices();
        A12 = (A12 + A12.transpose11()) * (0.5*D_gamma);	  

	if ( fdEvolutionT ) {
    	  A33 = (dntoda *D_gamma) + tensorI4;

          A13 = fdf->d2Fodsda(B_Mandel, *fdeps);
          A13 = Ce("ik")*A13("kjmn");
	    A13.null_indices();
          A13 = (A13 + A13.transpose1100()) * (0.5*D_gamma);

	  A23 = fdf->d2Fodqda(B_Mandel, *fdeps) *D_gamma;
	}

	// CCC: Tensor -> Matrix
	C99 = A11.BJtensor2BJmatrix_2();
	for (int i =1; i <=9; i++) {
	  for (int j =1; j <=9; j++) {
	    CCC.val(i,j) = C99.cval(i,j);
	    CCC.val(10+i,10+j) = 1.0; // For non-singularity
	  }
	}

	CCC.val(1,10)=CCC.val(10,1)=A12.cval(1,1);  
	CCC.val(2,10)=CCC.val(10,2)=A12.cval(1,2);	
	CCC.val(3,10)=CCC.val(10,3)=A12.cval(1,3);
	CCC.val(4,10)=CCC.val(10,4)=A12.cval(2,1);  
	CCC.val(5,10)=CCC.val(10,5)=A12.cval(2,2);	
	CCC.val(6,10)=CCC.val(10,6)=A12.cval(2,3);
	CCC.val(7,10)=CCC.val(10,7)=A12.cval(3,1);  
	CCC.val(8,10)=CCC.val(10,8)=A12.cval(3,2);	
	CCC.val(9,10)=CCC.val(10,9)=A12.cval(3,3);
	
	CCC.val(10,10) = a22;

	if ( fdEvolutionT ) {
	  C99 = A13.BJtensor2BJmatrix_2();
	  for (int i =1; i <=9; i++) {
	    for (int j =1; j <=9; j++) {
	      CCC.val(10+i,j) = C99.cval(i,j);
	      CCC.val(j,10+i) = C99.cval(i,j);
	    }
	  }

	  CCC.val(11,10)=CCC.val(10,11)=A23.cval(1,1);  
	  CCC.val(12,10)=CCC.val(10,12)=A23.cval(1,2);  
	  CCC.val(13,10)=CCC.val(10,13)=A23.cval(1,3);
	  CCC.val(14,10)=CCC.val(10,14)=A23.cval(2,1);  
	  CCC.val(15,10)=CCC.val(10,15)=A23.cval(2,2);  
	  CCC.val(16,10)=CCC.val(10,16)=A23.cval(2,3);
	  CCC.val(17,10)=CCC.val(10,17)=A23.cval(3,1);  
	  CCC.val(18,10)=CCC.val(10,18)=A23.cval(3,2);  
	  CCC.val(19,10)=CCC.val(10,19)=A23.cval(3,3);

	  C99 = A33.BJtensor2BJmatrix_2();
	  for (int i =1; i <=9; i++) {
	    for (int j =1; j <=9; j++) {
	      CCC.val(10+i,10+j) = C99.cval(i,j);
	    }
	  }
	}

	// Inverse of CCC
	//CCC.print("C","\n");
	CCC = CCC.inverse();
	//CCC.print();

	// CCC: Matrix -> Tensor
	for (int i =1; i <=9; i++) {
	  for (int j =1; j <=9; j++) {
	    C99.val(i,j) = CCC.cval(i,j);
	  }
	}
	A11 = C99.BJmatrix2BJtensor_2();

	A12.val(1,1)=CCC.cval(1,10); 
	A12.val(1,2)=CCC.cval(2,10);  
	A12.val(1,3)=CCC.cval(3,10);
	A12.val(2,1)=CCC.cval(4,10); 
	A12.val(2,2)=CCC.cval(5,10);  
	A12.val(2,3)=CCC.cval(6,10);
	A12.val(3,1)=CCC.cval(7,10); 
	A12.val(3,2)=CCC.cval(8,10);  
	A12.val(3,3)=CCC.cval(9,10);

	A21.val(1,1)=CCC.cval(10,1); 
	A21.val(1,2)=CCC.cval(10,2);  
	A21.val(1,3)=CCC.cval(10,3);
	A21.val(2,1)=CCC.cval(10,4); 
	A21.val(2,2)=CCC.cval(10,5);  
	A21.val(2,3)=CCC.cval(10,6);
	A21.val(3,1)=CCC.cval(10,7); 
	A21.val(3,2)=CCC.cval(10,8);  
	A21.val(3,3)=CCC.cval(10,9);

	a22  = CCC.cval(10,10);

	if ( fdEvolutionT ) {
	  for (int i =1; i <=9; i++) {
	    for (int j =1; j <=9; j++) {
	      C99.val(i,j) = CCC.cval(i,10+j);
	    }
	  }
	  A13 = C99.BJmatrix2BJtensor_2();

	  for (int i =1; i <=9; i++) {
	    for (int j =1; j <=9; j++) {
	      C99.val(j,i) = CCC.cval(10+i,j);
	    }
	  }
	  A31 = C99.BJmatrix2BJtensor_2();

	  A32.val(1,1)=CCC.cval(11,10); 
	  A32.val(1,2)=CCC.cval(12,10);  
	  A32.val(1,3)=CCC.cval(13,10);
	  A32.val(2,1)=CCC.cval(14,10); 
	  A32.val(2,2)=CCC.cval(15,10);  
	  A32.val(2,3)=CCC.cval(16,10);
	  A32.val(3,1)=CCC.cval(17,10); 
	  A32.val(3,2)=CCC.cval(18,10);  
	  A32.val(3,3)=CCC.cval(19,10);

	  A23.val(1,1)=CCC.cval(10,11); 
	  A23.val(1,2)=CCC.cval(10,12); 
	  A23.val(1,3)=CCC.cval(10,13);
	  A23.val(2,1)=CCC.cval(10,14); 
	  A23.val(2,2)=CCC.cval(10,15); 
	  A23.val(2,3)=CCC.cval(10,16);
	  A23.val(3,1)=CCC.cval(10,17); 
	  A23.val(3,2)=CCC.cval(10,18); 
	  A23.val(3,3)=CCC.cval(10,19);

	  for (int i =1; i <=9; i++) {
	    for (int j =1; j <=9; j++) {
	      C99.val(i,j) = CCC.cval(10+i,10+j);
	    }
	  }
	  A33 = C99.BJmatrix2BJtensor_2();
	}

        // Return d_gamma
	tensorTemp6 = dfods("ij")*LM("ijkl");
	  tensorTemp6.null_indices();
	if ( fdEvolutionT ) {
	  tensorTemp7 = dfoda("ij")*Ktensor("ijkl");
	    tensorTemp7.null_indices();
	}
	// !!!tensorTemp6&7 will be used later

        tensorTemp1 = tensorTemp6("kl")*A11("klmn"); 
	  tensorTemp1.null_indices();
        tensorTemp2  = A21 *(dfodq *Kscalar);
	tensorTemp4 = tensorTemp1 + tensorTemp2;
	if ( fdEvolutionT ) {
          tensorTemp3 = tensorTemp7("kl")*A31("klmn"); 
	    tensorTemp3.null_indices();
	  tensorTemp4 += tensorTemp3;
	}

        tensorTemp1 = tensorTemp6("kl")*A12("kl"); 
	  tensorTemp1.null_indices();
	temp1 = tensorTemp1.trace();
        temp2 = a22 *(dfodq *Kscalar);
	temp2 += temp1;
	if ( fdEvolutionT ) {
          tensorTemp3 = tensorTemp7("kl")*A32("kl");
	    tensorTemp3.null_indices();
	  temp3 = tensorTemp3.trace();
	  temp2 += temp3;
	}

        tensorTemp1 = tensorTemp6("kl")*A13("klmn");
	  tensorTemp1.null_indices();
        tensorTemp2  = A23 *(dfodq *Kscalar);
	tensorTemp5 = tensorTemp1 + tensorTemp2;
	if ( fdEvolutionT ) {
          tensorTemp3 = tensorTemp7("kl")*A33("klmn");
	    tensorTemp3.null_indices();
	  tensorTemp5 += tensorTemp3;
	}

	tensorTemp1 = tensorTemp4("ij") *res_Ee("ij");
	  tensorTemp1.null_indices();
	temp1 = tensorTemp1.trace();
	temp4 = temp1 + temp2*res_xi;
	if ( fdEvolutionT ) {
	  tensorTemp3 = tensorTemp5("ij") *res_eta("ij");
	    tensorTemp3.null_indices();
	  temp3 = tensorTemp3.trace();
	  temp4 += temp3;
	}

	tensorTemp1 = tensorTemp4("ij") *symMtensor("ij");
	  tensorTemp1.null_indices();
	temp1 = tensorTemp1.trace();
	temp5 = temp1 + temp2*nscalar;
	if ( fdEvolutionT ) {
	  tensorTemp3 = tensorTemp5("ij") *ntensor("ij");
	    tensorTemp3.null_indices();
	  temp3 = tensorTemp3.trace();
	  temp5 += temp3;
	}
	// !!! temp5 will be used later
	     	
	d_gamma = (yieldfun - temp4) / temp5;  // Here is d_gamma!!!
	//if(d_gamma < 0.0) d_gamma = 0.0;
	//printf("d_gamma= %e\n ", d_gamma);
	
	// Calculate incremental variables
	tensorTemp4 = symMtensor*(-d_gamma) - res_Ee;
	temp2 = nscalar*(-d_gamma) - res_xi;
	if ( fdEvolutionT ) {
	  tensorTemp5 = ntensor*(-d_gamma) - res_eta;
	}

	tensorTemp1 = A11("ijkl") *tensorTemp4("kl");
	  tensorTemp1.null_indices();
	D_Ee = tensorTemp1 + A12*temp2;
	if ( fdEvolutionT ) {
	  tensorTemp3 = A13("ijkl") *tensorTemp5("kl");
	    tensorTemp3.null_indices();
	  D_Ee += tensorTemp3;
	}

	tensorTemp1 = A21("kl") *tensorTemp4("kl");
	  tensorTemp1.null_indices();
	temp1 = tensorTemp1.trace();
	D_xi = temp1 + a22*temp2 ;
	if ( fdEvolutionT ) {
	  tensorTemp3 = A23("kl") *tensorTemp5("kl");
	    tensorTemp3.null_indices();
	  temp3 = tensorTemp3.trace();
	  D_xi += temp3;
	}

	if ( fdEvolutionT ) {
	  tensorTemp1 = A31("ijkl") *tensorTemp4("kl");
	    tensorTemp1.null_indices();
	  tensorTemp3 = A33("ijkl") *tensorTemp5("kl");
	    tensorTemp3.null_indices();
	  D_eta = tensorTemp1 + A32*temp2 + tensorTemp3;
	}

        // Update Variables	
        D_gamma += d_gamma;                  // updated D_gamma

	Ee += D_Ee;
	Ce = Ee*2.0 + tensorI2;
	//Ce.print();
	fde3d->setTrialC(Ce);  // Note: It is C, not F!!!
	B_PK2 = fde3d->getStressTensor();    // Updated B_PK2
        B_Mandel = Ce("ik")*B_PK2("kj");     // Update Mandel Stress
          B_Mandel.null_indices();

	xi += D_xi;
	q += Kscalar*D_xi;
	fdeps->setStrainLikeInVar(xi);
	fdeps->setStressLikeInVar(q);

	if ( fdEvolutionT ) {
	  eta += D_eta;
	  tensorTemp2 = Ktensor("ijkl")*D_eta("kl");
	    tensorTemp2.null_indices();
	  a += tensorTemp2;
	  fdeps->setStrainLikeKiVar(eta);
	  fdeps->setStressLikeKiVar(a);
	}

	// Calculate residuals
	res_Ee = (symMtensor*D_gamma) + Ee - Ee_n;		  

	res_xi = (nscalar*D_gamma) + xi - xi_n;
	tensorTemp1 = res_Ee("ij")*res_Ee("ij"); 
	  tensorTemp1.null_indices();
	res_norm_Ee = tensorTemp1.trace();
        residual = sqrt(( res_norm_Ee + res_xi*res_xi)/2.0);
	
	if ( fdEvolutionT ) {
	  res_eta = (ntensor*D_gamma) + eta - eta_n;
	  tensorTemp3 = res_eta("ij")*res_eta("ij"); 
	    tensorTemp3.null_indices();
	  res_norm_eta = tensorTemp3.trace();
	  residual = sqrt(( res_norm_Ee + res_xi*res_xi + res_norm_eta)/3.0);
	} 
	
	yieldfun = fdy->Yd(B_Mandel, *fdeps);        // Updated yieldfun
	//printf("Y= %e\n ", yieldfun);

	iter_counter++;
	
	if ( iter_counter > Max_Iter ) {
	  opserr << "Stop: Iteration More than " << Max_Iter;
	  opserr << " in return mapping algorithm of FD EP model" << "\n";
	  exit (-1);
	}

      } while ( yieldfun > fdy->getTolerance() || residual > tolerance ); // end of do - while

      //printf("Iter# = %d\n", iter_counter);
      // For Numerical stability
      D_gamma *= (1.0 - 0.5*tolerance);  
      if ( D_gamma < 0.0 ) {D_gamma = 0.0;}
      
      // Elastic-plastic tangent tensor
      tensorTemp1 = A11("klmn")*tensorTemp6("mn"); 
	tensorTemp1.null_indices();
      tensorTemp2  = A21 *(dfodq *Kscalar);
      tensorTemp4 = tensorTemp1 + tensorTemp2;
      if ( fdEvolutionT ) {
        tensorTemp3 = A31("klmn")*tensorTemp7("mn"); 
	  tensorTemp3.null_indices();
        tensorTemp4 += tensorTemp3;
      }

      tensorTemp1 = symMtensor("kl")*A11("klmn"); 
	tensorTemp1.null_indices();
      tensorTemp2  = A21 *(nscalar *Kscalar);
      tensorTemp5 = tensorTemp1 + tensorTemp2;
      if ( fdEvolutionT ) {
        tensorTemp3 = ntensor("kl")*A31("klmn"); 
	  tensorTemp3.null_indices();
        tensorTemp5 += tensorTemp3;
      }

      tensorTemp1 = tensorTemp4("ij")*tensorTemp5("kl");
        tensorTemp1.null_indices();
      
      tensorTemp2 = A11 - tensorTemp1*(1.0/temp5);
      LEPtensor = Ltensor("ijkl")*tensorTemp2("klmn");
        LEPtensor.null_indices();        
      
      // Update Fp
      tensorTemp2 = Mtensor("ij")*Fp("jk");  
        tensorTemp2.null_indices(); 
      Fp += tensorTemp2 *D_gamma;
      fdeps->setFpInVar(Fp); 
       
    } // end of plastic part
    
    else { // start of elastic part

      D_gamma = 0.0;
      d_gamma = 0.0;
      iter_counter = 0;
      
    } // end of elastic part 

    // Return iniTangent and iniPK2
    Fpinv = (fdeps->getFpInVar()).inverse();  //Using the iterative FP
    //Fpinv =  fdeps->getCommitedFpInVar().inverse(); //Using the initial Fp
    Fe = F("ij")*Fpinv("jk");   Fe.null_indices();

    Ce = Fe("ki")*Fe("kj");     Ce.null_indices();
    fde3d->setTrialC(Ce);       // Note: It is C, not F!!!
    B_PK2 = fde3d->getStressTensor();
    
    // if elastic, return only elastic tangent tensor
    if ( iter_counter == 0 )
      LEPtensor = fde3d->getTangentTensor(); 

    iniGreen = F("ki")*F("kj");
      iniGreen.null_indices();
    iniGreen = (iniGreen - tensorI2) * 0.5;

    // Pull backward of intermediate stress
    iniPK2 = Fpinv("ip")*Fpinv("jq")*B_PK2("pq");

    // Pull backward of intermediate EP tangent stiffness
    // iniTangent = Fpinv("im")*Fpinv("jn")*Fpinv("kr")*Fpinv("ls")*LEPtensor("mnrs");
    tensorTemp1 = Fpinv("im")*Fpinv("jn");
      tensorTemp1.null_indices();
    tensorTemp2 = Fpinv("kr")*Fpinv("ls");
      tensorTemp2.null_indices();
    tensorTemp3 = tensorTemp1("imjn")*LEPtensor("mnrs");
      tensorTemp3.null_indices();
    iniTangent = tensorTemp3("ijrs")*tensorTemp2("krls");
      iniTangent.null_indices();

    // Push forward of intermediate stress
    // cauchystress = Fe("ip")*B_PK2("pq")*Fe("jq")/F.determinant();
    cauchystress = Fe("ip")*B_PK2("pq");
      cauchystress.null_indices();
    cauchystress = cauchystress("iq")*Fe("jq");
      cauchystress.null_indices();
    cauchystress = cauchystress*(1.0/F.determinant());
    //cauchystress.print("X","\n");

    return 0;    
}  // */

/*//PK2 Version
int FiniteDeformationEP3D::ImplicitAlgorithm()
{
    // *** This is the key function ***

    // Initializing
    double yieldfun = 0.0;              // trial value of yield function
    double D_gamma  = 0.0;              // consistency parameter /Delta_{gamma}
    double d_gamma  = 0.0;              // increment of consistency parameter /delta_{gamma}
    int	 iter_counter = 0;

    straintensor res_Ee = tensorZ2;     // residual of intermediate Ee
    double  res_norm_Ee = 0.0;          // norm of residual of intermediate Ee
    double  res_xi = 0.0;	        // residual of strain like internal variable
    //double  residual = 0.0;	        // residual of Ee and xi

    straintensor Fp;	                // Plastic deformation gradient
    double xi = 0.0;		        // strain like internal variable
    double q = 0.0;		        // stress like internal variable
    straintensor Fpinv = tensorI2;      // inverse of Fp
    straintensor Ce = tensorI2;	        // intermediate C
    straintensor Ee = tensorZ2;	        // intermediate Ee
    
    straintensor Ee_n = tensorZ2;       // Ee at the incremental step n, calculated from Fp_n
    double xi_n;
            
    stresstensor Mtensor = tensorZ2;
    stresstensor symMtensor = tensorZ2;
    tensor Ltensor = tensorZ4 ;
    tensor LAinv = tensorZ4;
    tensor LEPtensor = tensorZ4;		
        
    double nscalar = 0.0;
    double Kscalar = 0.0;
    
    stresstensor dfods = tensorZ2 ;
    double dfodq = 0.0;
    tensor dMods = tensorZ4;
    tensor dsymMods = tensorZ4;
    double dnodq = 0.0;
    tensor At = tensorZ4;
    tensor Atinv = tensorZ4;
    double bs = 0.0;
    double bsinv = 0.0;
    tensor tensorTemp1;
    tensor tensorTemp2;
    tensor tensorTemp3;
    tensor tensorTemp4;
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    double temp4 = 0.0;
    double temp5 = 0.0;
    double D_xi  = 0.0;
    straintensor D_Ee = tensorZ2;

    tensorTemp3 = tensorI2("ij")*tensorI2("kl");
      tensorTemp3.null_indices();
    tensor tensorI4 = tensorTemp3.transpose0110();
    
    // Read the previous incremental step history variables
    Fp = fdeps->getCommitedFpInVar();
    xi = fdeps->getCommitedStrainLikeInVar();
    q = fdeps->getCommitedStressLikeInVar();

    Fpinv = Fp.inverse();
    Fe = F("ij")*Fpinv("jk");   Fe.null_indices();
    Ce = Fe("ki")*Fe("kj");     Ce.null_indices();
    Ee = (Ce - tensorI2) * 0.5;
    
    Ee_n = Ee;
    xi_n = xi;
    
    // Return stress from finite deformation elastic model
    fde3d->setTrialC(Ce);         // Note: It is C, not F!!!
    B_PK2 = fde3d->getStressTensor();

    // Evaluate the value of yield function
    yieldfun = fdy->Yd(B_PK2, *fdeps);
    ////printf("\nY0= %e\n", yieldfun);

    if ( yieldfun > (fdy->getTolerance()) ) { // start of plastic part
      D_gamma = 0.0;
      d_gamma = 0.0;
      iter_counter = 0;
        
      do {       // beginning of do - while	 

	// Return Symmetric tensor M and scalar n
        Mtensor = fdf->dFods(B_PK2, *fdeps);
        Mtensor = Ce("ij")*Mtensor("jk");
	  Mtensor.null_indices();
        symMtensor = (Mtensor + Mtensor.transpose11()) * 0.5;
    
        nscalar = fdf->dFodq(B_PK2, *fdeps);

        // Return tangent variables
	Ltensor = fde3d->getTangentTensor();
	
	Kscalar = fdEvolutionS->HModulus(B_PK2, *fdeps);

        // Return df/ds and df/dq
        dfods = fdy->dYods(B_PK2, *fdeps);
        
	dfodq = fdy->dYodq(B_PK2, *fdeps);

        // Return dsymM/ds and dn/dq
        dMods = fdf->d2Fodsds(B_PK2, *fdeps);
        dsymMods = Ce("ik")*dMods("kjmn");
	  dsymMods.null_indices();
        dsymMods = (dsymMods + dsymMods.transpose1100()) * 0.5;
        
	dnodq = fdf->d2Fodqdq(B_PK2, *fdeps);

        // Return tensor A and scalar b
        At = (dsymMods*D_gamma) + tensorI4;    
	
	bs = (dnodq*D_gamma) + 1.0;
    
        // Return inverse of A and B
        Atinv = At.inverse();
	
	bsinv = 1.0/bs;

        // Note: Only applied for d2fodsdq == tempZ2

        // Return d_gamma
        tensorTemp1 = dfods("ij")*Ltensor("ijkl")*Atinv("klmn");
	  tensorTemp1.null_indices();
        tensorTemp2 = tensorTemp1("mn")*res_Ee("mn");          
	  tensorTemp2.null_indices();
        temp2 = tensorTemp2.trace();
        tensorTemp3 = tensorTemp1("mn")*symMtensor("mn");       
	  tensorTemp3.null_indices();
        temp3 = tensorTemp3.trace();

        temp1 = dfodq*Kscalar*bsinv;
        temp4 = temp1*res_xi;
        temp5 = temp1*nscalar;
    	
	d_gamma = (yieldfun - temp2 - temp4) / (temp3 + temp5);  // Here is d_gamma!!!
	////printf("d_gamma= %e\n ", d_gamma);
	
	// Calculate incremental variables
        D_xi = bsinv * ( -d_gamma*nscalar - res_xi);
        
	tensorTemp1 =  symMtensor*(-d_gamma) - res_Ee;
        D_Ee = Atinv("ijkl") * tensorTemp1("kl");             
	  D_Ee.null_indices();

        // Update Variables	
        D_gamma += d_gamma;                  // updated D_gamma

	Ee += D_Ee;
	Ce = Ee*2.0 + tensorI2;
	fde3d->setTrialC(Ce);  // Note: It is C, not F!!!
	B_PK2 = fde3d->getStressTensor();    // Updated B_PK2

	xi += D_xi;
	q += Kscalar*D_xi;
	fdeps->setStrainLikeInVar(xi);
	fdeps->setStressLikeInVar(q);
	
	// Calculate residuals
	res_Ee = (symMtensor*D_gamma) + Ee - Ee_n;		  
	res_xi = (nscalar*D_gamma) + xi - xi_n;
	
	tensorTemp1 = res_Ee("ij")*res_Ee("ij"); 
	  tensorTemp1.null_indices();
	res_norm_Ee = tensorTemp1.trace();
	residual = sqrt( res_norm_Ee + res_xi*res_xi );  // combined norm of res_Ee & res_xi
	
	yieldfun = fdy->Yd(B_PK2, *fdeps);        // Updated yieldfun
	////printf("Y= %e\n ", yieldfun);

	iter_counter++;
	
	if ( iter_counter > Max_Iter ) {
	  opserr << "Stop: Iteration More than " << Max_Iter;
	  opserr << " in return mapping algorithm of FD EP model" << "\n";
	  exit (-1);
	}

      } while ( yieldfun > fdy->getTolerance() || residual > tolerance ); // end of do - while
      //} while ( yieldfun > fdy->getTolerance() || res_norm_Ee > tolerance || res_xi > tolerance); // end of do - while

      D_gamma *= (1.0 - tolerance);  // For Numerical stability
      if ( D_gamma < 0.0 ) D_gamma = 0.0;
      
      // Elastic-plastic tangent tensor
      LAinv = Ltensor("ijkl")*Atinv("klmn");
        LAinv.null_indices();      

      tensorTemp1 = LAinv("ijmn")*symMtensor("mn");
        tensorTemp1.null_indices();      
      tensorTemp2 = dfods("ij")*LAinv("ijmn");
	tensorTemp2.null_indices();
      tensorTemp3 = tensorTemp1("ij")*tensorTemp2("mn");
        tensorTemp3.null_indices();
      tensorTemp4 = tensorTemp2("mn")*symMtensor("mn");
        tensorTemp4.null_indices();
      
      temp1 = tensorTemp4.trace();
      temp2 = dfodq*Kscalar*bsinv*nscalar;
      temp3 = 1.0 / (temp1 + temp2);

      LEPtensor = LAinv - tensorTemp3*temp3;

      // Update Fp
      tensorTemp1 = Mtensor("ij")*Fp("jk");  
        tensorTemp1.null_indices(); 
      Fp += tensorTemp1*D_gamma;
      fdeps->setFpInVar(Fp); 
       
    } // end of plastic part
    
    else { // start of elastic part

      D_gamma = 0.0;
      d_gamma = 0.0;
      iter_counter = 0;
      
    } // end of elastic part 

    // Return iniTangent and iniPK2
    Fpinv = (fdeps->getFpInVar()).inverse();  //Using the iterative FP
    //Fpinv =  fdeps->getCommitedFpInVar().inverse(); //Using the initial Fp
    Fe = F("ij")*Fpinv("jk");   Fe.null_indices();

    Ce = Fe("ki")*Fe("kj");     Ce.null_indices();
    fde3d->setTrialC(Ce);       // Note: It is C, not F!!!
    B_PK2 = fde3d->getStressTensor();
    
    // if elastic, return only elastic tangent tensor
    if ( iter_counter == 0 )
      LEPtensor = fde3d->getTangentTensor(); 

    iniGreen = F("ki")*F("kj");
      iniGreen.null_indices();
    iniGreen = (iniGreen - tensorI2) * 0.5;

    // Pull backward of intermediate stress
    iniPK2 = Fpinv("ip")*Fpinv("jq")*B_PK2("pq");

    // Pull backward of intermediate EP tangent stiffness
    // iniTangent = Fpinv("im")*Fpinv("jn")*Fpinv("kr")*Fpinv("ls")*LEPtensor("mnrs");
    tensorTemp1 = Fpinv("im")*Fpinv("jn");
      tensorTemp1.null_indices();
    tensorTemp2 = Fpinv("kr")*Fpinv("ls");
      tensorTemp2.null_indices();
    tensorTemp3 = tensorTemp1("imjn")*LEPtensor("mnrs");
      tensorTemp3.null_indices();
    iniTangent = tensorTemp3("ijrs")*tensorTemp2("krls");
      iniTangent.null_indices();

    // Push forward of intermediate stress
    // cauchystress = Fe("ip")*B_PK2("pq")*Fe("jq")/F.determinant();
    cauchystress = Fe("ip")*B_PK2("pq");
      cauchystress.null_indices();
    cauchystress = cauchystress("iq")*Fe("jq");
      cauchystress.null_indices();
    cauchystress = cauchystress*(1.0/F.determinant());

    return 0;    
}  //*/ 

#endif
