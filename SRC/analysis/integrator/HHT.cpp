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
                                                                        
// $Revision: 1.7 $
// $Date: 2002-12-05 22:33:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/HHT.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/HHT.C
// 
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the implementation of the HHT class.
//
// What: "@(#) HHT.C, revA"

#include <HHT.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//-----------------------------------------
//set alpha = 0.5 
//    beta  = 0.5 
//    gamma = 1.0 
//for symplectic midpoint
//-----------------------------------------


HHT::HHT()
:TransientIntegrator(INTEGRATOR_TAGS_HHT),
 alpha(0.5), gamma(1.0), beta(0), 
 alphaM(0.0), betaK(0.0), betaKi(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{
    
}

HHT::HHT(double _alpha)
:TransientIntegrator(INTEGRATOR_TAGS_HHT),
 alpha(_alpha), gamma(1.5-_alpha), beta((2-_alpha)*(2-_alpha)*0.25),
 alphaM(0.0), betaK(0.0), betaKi(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{
    
}

HHT::HHT(double _alpha, double alpham, double betak, double betaki)
:TransientIntegrator(INTEGRATOR_TAGS_HHT),
 alpha(_alpha), gamma(1.5-_alpha), beta((2-_alpha)*(2-_alpha)*0.25),  
 alphaM(alpham), betaK(betak), betaKi(betaki),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{

}


//generalized HHT alpha method
HHT::HHT(double _alpha, double _beta, double _gamma)
:TransientIntegrator(INTEGRATOR_TAGS_HHT),
 alpha(_alpha), gamma(_gamma), beta(_beta),
 alphaM(0.0), betaK(0.0), betaKi(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{


}

HHT::HHT(double _alpha, double _beta, double _gamma,
	 double alpham, double betak, double betaki)
:TransientIntegrator(INTEGRATOR_TAGS_HHT),
 alpha(_alpha), gamma(_gamma), beta(_beta),
 alphaM(alpham), betaK(betak), betaKi(betaki),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{

}

HHT::~HHT()
{
  // clean up the memory created
  if (Ut != 0)
    delete Ut;
  if (Utdot != 0)
    delete Utdot;
  if (Utdotdot != 0)
    delete Utdotdot;
  if (U != 0)
    delete U;
  if (Udot != 0)
    delete Udot;
  if (Udotdot != 0)
    delete Udotdot;
  if (Ualpha != 0)
    delete Ualpha;
  if (Udotalpha != 0)
    delete Udotalpha;
}


int
HHT::initialize(void)
{
  return 0;
}



int
HHT::newStep(double deltaT)
{

  if (beta == 0 || gamma == 0 ) {
    cerr << "HHT::newStep() - error in variable\n";
    cerr << "gamma = " << gamma << " beta= " << beta << endl;
    return -1;
  }
    
  if (deltaT <= 0.0) {
    cerr << "HHT::newStep() - error in variable\n";
    cerr << "dT = " << deltaT << endl;
    return -2;	
  }
  
  c1 = alpha;
  c2 = gamma/(beta*deltaT);
  c3 = 1.0/(beta*deltaT*deltaT);

  AnalysisModel *theModel = this->getAnalysisModelPtr();

  if (U == 0) {
    cerr << "HHT::newStep() - domainChange() failed or hasn't been called\n";
    return -3;	
  }

  // set response at t to be that at t+delta t of previous step
  (*Ut) = *U;        
  (*Utdot) = *Udot;  
  (*Utdotdot) = *Udotdot;  
    
  // set new velocity and accelerations at t + delta t
  double a1 = (1.0 - gamma/beta); 

  double a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
  // (*Udot) *= a1;  
  Udot->addVector(a1,*Utdotdot,a2);

  double a3 = -1.0/(beta*deltaT);
  double a4 = 1 - 0.5/beta;
  //  (*Udotdot) *= a4;  
   Udotdot->addVector(a4,*Utdot,a3);

  (*Ualpha) = *Ut;
  (*Udotalpha) = *Utdot;
  //  (*Udotalpha) *= (1 - alpha);
  Udotalpha->addVector((1-alpha),*Udot, alpha);

  // set the new trial response quantities
  theModel->setResponse(*Ualpha,*Udotalpha,*Udotdot);        

  // increment the time and apply the load
  double time = theModel->getCurrentDomainTime();
  time +=deltaT;
  if (theModel->updateDomain(time, deltaT) < 0) {
    cerr << "HHT::newStep() - failed to update the domain\n";
    return -4;
  }

  return 0;
}


int
HHT::revertToLastStep()
{
  // set response at t+delta t to be that at t .. for next newStep
  if (U != 0) {
    (*U) = *Ut;        
    (*Udot) = *Utdot;  
    (*Udotdot) = *Utdotdot;  
  }
  return 0;
}



int
HHT::formEleTangent(FE_Element *theEle)
{
  theEle->zeroTangent();
  if (statusFlag == CURRENT_TANGENT) {
    theEle->addKtToTang(c1) ;
    theEle->addCtoTang(c2);
    theEle->addMtoTang(c3);
  } else if (statusFlag == INITIAL_TANGENT) {
    theEle->addKiToTang(c1);
    theEle->addCtoTang(c2);
    theEle->addMtoTang(c3);
  }    

  return 0;
}    



int
HHT::formNodTangent(DOF_Group *theDof)
{
  theDof->zeroTangent();
  theDof->addMtoTang(c3);
  theDof->addCtoTang(c2);
  
  return(0);
}    



int 
HHT::domainChanged()
{
  AnalysisModel *myModel = this->getAnalysisModelPtr();
  LinearSOE *theLinSOE = this->getLinearSOEPtr();
  const Vector &x = theLinSOE->getX();
  int size = x.Size();

  // if damping factors exist set them in the ele & node of the domain
  if (alphaM != 0.0 || betaK != 0.0 || betaKi != 0.0)
    myModel->setRayleighDampingFactors(alphaM, betaK, betaKi);

  // create the new Vector objects
  if (Ut == 0 || Ut->Size() != size) {

    // delete the old
    if (Ut != 0)
      delete Ut;
    if (Utdot != 0)
      delete Utdot;
    if (Utdotdot != 0)
      delete Utdotdot;
    if (U != 0)
      delete U;
    if (Udot != 0)
      delete Udot;
    if (Udotdot != 0)
      delete Udotdot;
    if (Ualpha != 0)
      delete Ualpha;
    if (Udotalpha != 0)
      delete Udotalpha;
    
    // create the new
    Ut = new Vector(size);
    Utdot = new Vector(size);
    Utdotdot = new Vector(size);
    U = new Vector(size);
    Udot = new Vector(size);
    Udotdot = new Vector(size);
    Ualpha = new Vector(size);
    Udotalpha = new Vector(size);

    // check we obtained the new
    if (Ut == 0 || Ut->Size() != size ||
	Utdot == 0 || Utdot->Size() != size ||
	Utdotdot == 0 || Utdotdot->Size() != size ||
	U == 0 || U->Size() != size ||
	Udot == 0 || Udot->Size() != size ||
	Udotdot == 0 || Udotdot->Size() != size ||
	Ualpha == 0 || Ualpha->Size() != size ||
	Udotalpha == 0 || Udotalpha->Size() != size) {
  
      cerr << "HHT::domainChanged - ran out of memory\n";

      // delete the old
      if (Ut != 0)
	delete Ut;
      if (Utdot != 0)
	delete Utdot;
      if (Utdotdot != 0)
	delete Utdotdot;
      if (U != 0)
	delete U;
      if (Udot != 0)
	delete Udot;
      if (Udotdot != 0)
	delete Udotdot;
    if (Ualpha != 0)
      delete Ualpha;
    if (Udotalpha != 0)
      delete Udotalpha;

      Ut = 0; Utdot = 0; Utdotdot = 0;
      U = 0; Udot = 0; Udotdot = 0; Udotalpha=0; Ualpha =0;
      return -1;
    }
  }        
    
  // now go through and populate U, Udot and Udotdot by iterating through
  // the DOF_Groups and getting the last committed velocity and accel

  DOF_GrpIter &theDOFs = myModel->getDOFs();
  DOF_Group *dofPtr;

  while ((dofPtr = theDOFs()) != 0) {
    const ID &id = dofPtr->getID();
    int idSize = id.Size();


	int i;
    const Vector &disp = dofPtr->getCommittedDisp();	
    for (i=0; i < idSize; i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*U)(loc) = disp(i);		
      }
    }

    const Vector &vel = dofPtr->getCommittedVel();
    for (i=0; i < idSize; i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*Udot)(loc) = vel(i);
      }
    }

    const Vector &accel = dofPtr->getCommittedAccel();	
    for (i=0; i < idSize; i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*Udotdot)(loc) = accel(i);
      }
    }
    /** NOTE WE CAN't DO TOGETHER BECAUSE DOF_GROUPS USING SINGLE VECTOR ******
    for (int i=0; i < id.Size(); i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*U)(loc) = disp(i);		
 	(*Udot)(loc) = vel(i);
 	(*Udotdot)(loc) = accel(i);
      }
    }
    *******************************************************************************/

  }    

  return 0;
}


int
HHT::update(const Vector &deltaU)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING HHT::update() - no AnalysisModel set\n";
    return -1;
  }	

  // check domainChanged() has been called, i.e. Ut will not be zero
  if (Ut == 0) {
    cerr << "WARNING HHT::update() - domainChange() failed or not called\n";
    return -2;
  }	

  // check deltaU is of correct size
  if (deltaU.Size() != U->Size()) {
    cerr << "WARNING HHT::update() - Vectors of incompatable size ";
    cerr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endl;
    return -3;
  }
    
  //  determine the response at t+delta t
  (*U) += deltaU;
  Udot->addVector(1.0,deltaU,c2);

  Udotdot->addVector(1.0,deltaU,c3);
  Ualpha->addVector(1.0,deltaU, c1 ) ;   //c1 = alpha
  Udotalpha->addVector(1.0,deltaU, c2*alpha);
  
  // update the responses at the DOFs
  theModel->setResponse(*Ualpha,*Udotalpha,*Udotdot);        
  if (theModel->updateDomain() < 0) {
    cerr << "HHT::update() - failed to update the domain\n";
    return -4;
  }
		    
  return 0;
}    


int
HHT::commit(void)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING HHT::commit() - no AnalysisModel set\n";
    return -1;
  }	  

  // update the responses at the DOFs
  theModel->setResponse(*U,*Udot,*Udotdot);        
  theModel->updateDomain();

  return theModel->commitDomain();
}

int
HHT::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(7);
    data(0) = alpha;
    data(1) = beta;
    data(2) = gamma;
    data(3) = 1.0;	
    data(4) = alphaM;
    data(5) = betaK;
    data(6) = betaKi;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
	cerr << "WARNING HHT::sendSelf() - could not send data\n";
	return -1;
    }	
    return 0;
}

int
HHT::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(7);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
	cerr << "WARNING HHT::recvSelf() - could not receive data\n";
	return -1;
    }

    alpha = data(0);
    beta = data(1);
    gamma = data(2);
    
    alphaM = data(4);
    betaK = data(5);
    betaKi = data(6);
    
    return 0;
}

void
HHT::Print(ostream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t HHT - currentTime: " << currentTime << endl ;
        s << "\t alpha = " << alpha << endl ;
	s << "\t beta  = " << beta  << endl ;
        s << "\t gamma = " << gamma << endl ;
	s << "  Rayleigh Damping - alphaM: " << alphaM;
	s << "  betaK: " << betaK << "   betaKi: " << betaKi << endl;	    
    } else 
	s << "\t HHT - no associated AnalysisModel\n";
}

