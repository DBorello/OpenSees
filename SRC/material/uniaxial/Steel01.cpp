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
                                                                        
// $Revision: 1.9 $
// $Date: 2003-02-25 23:33:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel01.cpp,v $
                                                                        
                                                                        
// File: ~/material/Steel01.C
//
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class implementation for 
// Steel01. 
//
// What: "@(#) Steel01.C, revA"


#include <Steel01.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>

Steel01::Steel01
(int tag, double FY, double E, double B,
 double A1, double A2, double A3, double A4,
 double min, double max) :
   UniaxialMaterial(tag,MAT_TAG_Steel01),
   fy(FY), E0(E), b(B), a1(A1), a2(A2), a3(A3), a4(A4),
   epsmin(min), epsmax(max)
{
   // Calculated material parameters
   epsy = fy/E0;         // Yield strain
   Esh = b*E0;           // Hardening modulus

   // Sets all history and state variables to initial values
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;
   Cfailed = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;
   Tfailed = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	gradientIdentifier = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Steel01::Steel01():UniaxialMaterial(0,MAT_TAG_Steel01),
 fy(0.0), E0(0.0), b(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0),
 epsmin(0.0), epsmax(0.0)
{

// AddingSensitivity:BEGIN /////////////////////////////////////
	gradientIdentifier = 0;
// AddingSensitivity:END //////////////////////////////////////

}

Steel01::~Steel01 ()
{
   // Does nothing
}

int Steel01::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tfailed = Cfailed;

   // Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   //   if (fabs(dStrain) > DBL_EPSILON)   
      determineTrialState (dStrain);

   return 0;
}

int Steel01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tfailed = Cfailed;

   // Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain;
   dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   // if (fabs(dStrain) > DBL_EPSILON) 
      determineTrialState (dStrain);

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void Steel01::determineTrialState (double dStrain)
{

   if (Tstrain <= epsmin || Tstrain >= epsmax)
      Tfailed = 1;

   if (Tfailed) {
      Tstress = 0.0;
      Ttangent = 0.0;
      return;
   } else {

      double c, c1, c2, c3, fyOneMinusB, c1c3, c1c2;

      fyOneMinusB = fy * (1.0 - b);

      c1 = Esh*Tstrain;

      c2 = TshiftN*fyOneMinusB;

      c3 = TshiftP*fyOneMinusB;

      c1c3 = c1+c3;

      c = Cstress + E0*dStrain;

      if (c1c3 < c) 
	Tstress = c1c3;
      else
	Tstress = c;

      c1c2=c1-c2;
      if (c1c2 > Tstress) 
         Tstress = c1c2;

      if (fabs(Tstress-c) < DBL_EPSILON)
	  Ttangent = E0;
      else
	  Ttangent = Esh;      
      
      //
      // Determine if a load reversal has occurred due to the trial strain
      //

      // Determine initial loading condition
      if (Tloading == 0 && dStrain != 0.0) {
	  if (dStrain > 0.0)
	    Tloading = 1;
	  else
	    Tloading = -1;
      }

      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment
      if (Tloading == 1 && dStrain < 0.0) {
	  Tloading = -1;
	  if (Cstrain > TmaxStrain)
	    TmaxStrain = Cstrain;
	  TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Tloading == -1 && dStrain > 0.0) {
	  Tloading = 1;
	  if (Cstrain < TminStrain)
	    TminStrain = Cstrain;
	  TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
      }
   }
}

double Steel01::getStrain ()
{
   return Tstrain;
}

double Steel01::getStress ()
{
   return Tstress;
}

double Steel01::getTangent ()
{
   return Ttangent;
}


int Steel01::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;
   Cfailed = Tfailed;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int Steel01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tfailed = Cfailed;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int Steel01::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;
   Cfailed = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;
   Tfailed = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;

   return 0;
}

UniaxialMaterial* Steel01::getCopy ()
{
   Steel01* theCopy = new Steel01(this->getTag(), fy, E0, b,
				  a1, a2, a3, a4, epsmin, epsmax);

   // Calculated material properties
   theCopy->epsy = epsy;
   theCopy->Esh = Esh;

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CshiftP = CshiftP;
   theCopy->CshiftN = CshiftN;
   theCopy->Cloading = Cloading;
   theCopy->Cfailed = Cfailed;

   // Trial history variables
   theCopy->TminStrain = CminStrain;
   theCopy->TmaxStrain = CmaxStrain;
   theCopy->TshiftP = CshiftP;
   theCopy->TshiftN = CshiftN;
   theCopy->Tloading = Cloading;
   theCopy->Tfailed = Cfailed;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   // Trial state variables
   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   return theCopy;
}

int Steel01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(19);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = a1;
   data(5) = a2;
   data(6) = a3;
   data(7) = a4;
   data(8) = epsmin;
   data(9) = epsmax;

   // History variables from last converged state
   data(10) = CminStrain;
   data(11) = CmaxStrain;
   data(12) = CshiftP;
   data(13) = CshiftN;
   data(14) = Cloading;
   data(15) = Cfailed;

   // State variables from last converged state
   data(16) = Cstrain;
   data(17) = Cstress;
   data(18) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Steel01::sendSelf() - failed to send data\n";

   return res;
}

int Steel01::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(19);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "Steel01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      b = data(3);
      a1 = data(4);
      a2 = data(5);
      a3 = data(6);
      a4 = data(7);
      epsmin = data(8);
      epsmax = data(9);

      epsy = fy/E0;
      Esh = b*E0;

      // History variables from last converged state
      CminStrain = data(10);
      CmaxStrain = data(11);
      CshiftP = data(12);
      CshiftN = data(13);
      Cloading = int(data(14));
      Cfailed = int(data(15));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;
      Tfailed = Cfailed;

      // State variables from last converged state
      Cstrain = data(16);
      Cstress = data(17);
      Ctangent = data(18);      

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}

void Steel01::Print (OPS_Stream& s, int flag)
{
   s << "Steel01 tag: " << this->getTag() << endln;
   s << "  fy: " << fy << " ";
   s << "  E0: " << E0 << " ";
   s << "  b:  " << b << " ";
   s << "  a1: " << a1 << " ";
   s << "  a2: " << a2 << " ";
   s << "  a3: " << a3 << " ";
   s << "  a4: " << a4 << " ";
   if (epsmin != NEG_INF_STRAIN)
     s << "  epsmin: " << epsmin << " ";
   if (epsmax != POS_INF_STRAIN)
     s << "  epsmax: " << epsmax << endln;
}

int
Steel01::setParameter(const char **argv, int argc, Information &info)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"fy") == 0) {
		info.theType = DoubleType;
		return 1;
	}
	if (strcmp(argv[0],"E0") == 0) {
		info.theType = DoubleType;
		return 2;
	}
	if (strcmp(argv[0],"b") == 0) {
		info.theType = DoubleType;
		return 3;
	}
	if (strcmp(argv[0],"a1") == 0) {
		info.theType = DoubleType;
		return 4;
	}
	if (strcmp(argv[0],"a2") == 0) {
		info.theType = DoubleType;
		return 5;
	}
	if (strcmp(argv[0],"a3") == 0) {
		info.theType = DoubleType;
		return 6;
	}
	if (strcmp(argv[0],"a4") == 0) {
		info.theType = DoubleType;
		return 7;
	}
	if (strcmp(argv[0],"epsmin") == 0) {
		info.theType = DoubleType;
		return 8;
	}
	if (strcmp(argv[0],"epsmax") == 0) {
		info.theType = DoubleType;
		return 9;
	}
	else
		opserr << "WARNING: Could not set parameter in Steel01. " << endln;
                
	return -1;
}

int
Steel01::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->fy = info.theDouble;
		break;
	case 2:
		this->E0 = info.theDouble;
		break;
	case 3:
		this->b = info.theDouble;
		break;
	case 4:
		this->a1 = info.theDouble;
		break;
	case 5:
		this->a2 = info.theDouble;
		break;
	case 6:
		this->a3 = info.theDouble;
		break;
	case 7:
		this->a4 = info.theDouble;
		break;
	case 8:
		this->epsmin = info.theDouble;
		break;
	case 9:
		this->epsmax = info.theDouble;
		break;
	default:
		return -1;
	}

	epsy = fy/E0;           // Yield strain
	Esh = b*E0;             // Hardening modulus
	Ttangent = E0;          // Initial stiffness

	return 0;
}


// AddingSensitivity:BEGIN ///////////////////////////////////
int
Steel01::gradient(bool compute, int identifier, double & gradient)
{
/*	The gradient method can be called with four different purposes:
	1) To clear the sensitivity flag so that the object does not contribute:
			gradient(false, 0)
	2) To set the sensitivity flag so that the object contributes
	   (the sensitivity flag is stored as the value of parameterID):
			gradient(false, parameterID)
	3) To obtain the gradient vector from the object (like for the residual):
			gradient(true, 0)
	4) To commit unconditional sensitivities for path-dependent problems:
			gradient(true, gradNumber)
*/

// COMPUTE GRADIENTS
	if (compute) {

// IF "PHASE 1" IN THE GRADIENT COMPUTATIONS (RETURN GRADIENT VECTOR)
		if (identifier == 0) {

			// For now, simply return the derivatives of these equations:
			// sigma = E0 * epsilon                 for abs(epsilon) < fy/E0
			// sigma = fy + (epsilon-fy/E0)*b*E0    for abs(epsilon) > fy/E0
			// (That is, this method only works for path-independent problems for now.)

			if ( gradientIdentifier == 0 ) {
				gradient = 0.0;
			}
			else if ( fabs(Tstrain) < fy/E0 ) {				// below the yield limit
				if ( gradientIdentifier == 1 ) {		// d{sigma}d{fy}
					gradient = 0.0;
				}
				else if ( gradientIdentifier == 2  ) {	// d{sigma}d{E0}
					gradient = Tstrain;
				}
				else if ( gradientIdentifier == 3  )	// d{sigma}d{b}
					gradient = 0.0;
				else
					gradient = 0.0;
			}
			else {										// past the yield limit
				if ( gradientIdentifier == 1 ) {		// d{sigma}d{fy}
					gradient = 1.0 - b;
				}
				else if ( gradientIdentifier == 2  ) {	// d{sigma}d{E0}
					gradient = Tstrain*b;
				}
				else if ( gradientIdentifier == 3  ) {	// d{sigma}d{b}
					gradient = Tstrain*E0 - fy;
				}
				else {
					gradient = 0.0;
				}
			}

			return 0;
		}

// IF "PHASE 2" IN THE GRADIENT COMPUTATIONS (COMMIT UNCONDITIONAL GRADIENT)
		else {
			// Not treated yet. 
		}
	}
	
// DO NOT COMPUTE GRADIENTS, JUST SET FLAG
	else {
		gradientIdentifier = identifier;
	}

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////
