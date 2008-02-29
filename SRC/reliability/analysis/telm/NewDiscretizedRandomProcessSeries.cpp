/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewDiscretizedRandomProcessSeries.cpp,v $

#include <NewDiscretizedRandomProcessSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <ModulatingFunction.h>
#include <Filter.h>
#include <classTags.h>


NewDiscretizedRandomProcessSeries::NewDiscretizedRandomProcessSeries(int num, 
								     ModulatingFunction **theModFuncs,
								     double p_mean,
								     double p_maxStdv)
:TimeSeries(TSERIES_TAG_DiscretizedRandomProcessSeries)
{
  numRandVar=0;
  randomVariables = 0;
  kickInTimes = 0;
  theModulatingFunctions = theModFuncs;
  numModFuncs = num;
  mean = p_mean;
  maxStdv = p_maxStdv;
  arrayID = 0;
  MaxRanVarID = 0;
  parameterID = -1;
  c = maxStdv;
  active = 0;
  
  output.open("excitation.txt", ios::out);
  
}


NewDiscretizedRandomProcessSeries::~NewDiscretizedRandomProcessSeries()
{
  if (randomVariables != 0)  {
    delete randomVariables;
    randomVariables = 0;
  }
  
  if (kickInTimes != 0) {
    delete kickInTimes;
    kickInTimes = 0;
  }
  if (arrayID !=0 ) {
    delete [] arrayID;
    arrayID = 0;
  }
  if (active !=0 ) {
    delete [] active;
    active = 0;
  }
  if ( theModulatingFunctions != 0 ){
    delete [] theModulatingFunctions;
    theModulatingFunctions = 0;
  }
}


int
NewDiscretizedRandomProcessSeries::setParameter (const char **argv, int argc, Information &info)
{

    if (argc < 1)
        return -1;

	int* arrayID_tmp=0;
	int i;
	Vector* temp=0;

	int rvNumber = info.theInt; // ID# of the random variable

	// The second argument tells when the random variable "kicks in".
	// Store this in a table...
	// In case the vector doesn't exist
	if (numRandVar == 0 ) {
		numRandVar++;
		if( kickInTimes != 0 ){
			delete kickInTimes;
			kickInTimes=0;
		}
		kickInTimes = new Vector(numRandVar);
		double ktime=(double)atof(argv[0])+1.0e-8;
		(*kickInTimes)(numRandVar-1)=(double)atof(argv[0])+1.0e-8;

		if( active != 0 ) {
			delete [] active;
			active = 0;
		}
		active = new bool[numRandVar];
		active[numRandVar-1]=false;

		if( randomVariables != 0 ){
			delete randomVariables;
			randomVariables = 0;
		}
		randomVariables = new Vector(numRandVar);
		for( i=0; i<numRandVar; i++) (*randomVariables)(i)=0.0; 

		if( arrayID != 0 ){
			delete [] arrayID;
			arrayID = 0;
		}
		arrayID = new int[rvNumber+1];
		for ( i=0; i<rvNumber; i++ ) arrayID[i]=-1; 
		arrayID[rvNumber]=numRandVar-1;
		MaxRanVarID = rvNumber;
	}
	else {
		if( rvNumber <= MaxRanVarID ){
			if(arrayID[rvNumber]==0){
				numRandVar++;
				if(numRandVar==750){
					double tmpp=0.0;
				}
				arrayID[rvNumber]=numRandVar-1;

				if( active != 0 ) {
					delete [] active;
					active = 0;
				}
				active = new bool[numRandVar];
				for ( i=0; i< numRandVar; i++ )active[i]=false;

				if(temp != 0){
					delete temp;
					temp = 0;
				}
				temp = new Vector(kickInTimes->Size());
				(*temp) = (*kickInTimes);
				if( kickInTimes != 0 ) {
					delete kickInTimes;
					kickInTimes = 0;
				}
				kickInTimes = new Vector(numRandVar);
				for ( i=0; i< numRandVar-1; i++ ) (*kickInTimes)(i) = (*temp)(i);
			    double ktime=(double)atof(argv[0])+1.0e-8;
				(*kickInTimes)(numRandVar-1)=(double)atof(argv[0])+1.0e-8;

				if( randomVariables != 0 ) {
					delete randomVariables;
					randomVariables = 0;
				}
				randomVariables = new Vector(numRandVar);
				for( i=0; i<numRandVar; i++) (*randomVariables)(i)=0.0; 

			}else{
				opserr<< " error in discretized time sereise \n";
				opserr<< " random variable "<<rvNumber<< " is defined already \n";
				exit(-1);
			}
		}else{
			numRandVar++;
			if(numRandVar==750){
				double tmpp=0.0;
			}
//			opserr<< "numRandVar=" <<numRandVar <<"\n";
//			opserr<< "arrayID_tmp bef delete =" <<arrayID_tmp<<"\n";
			if( arrayID_tmp != 0 ){
				delete [] arrayID_tmp;
				arrayID_tmp = 0;
			}
//			opserr<< "arrayID_tmp aft delete =" <<arrayID_tmp<<"\n";
			arrayID_tmp = new int[rvNumber+1];
//			opserr<< "arrayID_tmp aft new =" <<arrayID_tmp<<"\n";
			for ( i=0; i<=MaxRanVarID; i++ ) arrayID_tmp[i]=arrayID[i];
//			opserr << " check arrayID 1 \n";
			if( rvNumber > MaxRanVarID+1) 
			for ( i=MaxRanVarID+1; i<rvNumber; i++ ) arrayID[i]=-1;
//			opserr << " check arrayID 2 \n";
			arrayID_tmp[rvNumber]=numRandVar-1;
//			opserr << " check arrayID 3 \n";
//			opserr<< "arrayID bef delete =" <<arrayID<<"\n";
			if( arrayID != 0 ) {
				delete [] arrayID;
				arrayID = 0;
			}
//			opserr<< "arrayID aft delete =" <<arrayID<<"\n";
			arrayID = new int[rvNumber+1];
//			opserr<< "arrayID aft new =" <<arrayID<<"\n";
			for ( i=0; i<=rvNumber; i++ ) arrayID[i]=arrayID_tmp[i];
//			opserr << " check arrayID 4 \n";
			MaxRanVarID = rvNumber;

//			opserr<< "active bef delete =" <<active<<"\n";
			if ( active != 0 ) { 
				delete [] active;
				active = 0;
			}
//			opserr<< "active aft delete =" <<active<<"\n";
			active = new bool[numRandVar];
//			opserr<< "active aft new =" <<active<<"\n";
			for ( i=0; i< numRandVar; i++ )active[i]=false;
//			opserr << " check active 4 \n";

//			opserr<< "temp bef delete =" <<temp<<"\n";
			if( temp !=0) {
				delete temp;
				temp=0;
			}
//			opserr<< "temp aft delete =" <<temp<<"\n";
			temp = new Vector(kickInTimes->Size());
//			opserr<< "temp aft new =" <<temp<<"\n";
			(*temp) = (*kickInTimes);
//			opserr << " check temp1 \n";
//			opserr << "kickInTimes bef delete =" <<kickInTimes<<"\n";
			if ( kickInTimes !=0 ) {
				delete kickInTimes;
				kickInTimes=0;
			}
//			opserr << "kickInTimes aft delete =" <<kickInTimes<<"\n";
			kickInTimes = new Vector(numRandVar);
//			opserr << "kickInTimes aft new =" <<kickInTimes<<"\n";
			for ( i=0; i< numRandVar-1; i++ )(*kickInTimes)(i) = (*temp)(i);
//			opserr << " check kickintime1 \n";
		    double ktime=(double)atof(argv[0])+1.0e-8;
			(*kickInTimes)(numRandVar-1)=(double)atof(argv[0])+1.0e-8;
//			opserr << " check kickintime2 \n";

//			opserr << "randomVariables bef delete =" <<randomVariables<<"\n";
			if( randomVariables != 0) {
				delete randomVariables;
				randomVariables=0;
			}
//			opserr << "randomVariables aft delete =" <<randomVariables<<"\n";
			randomVariables = new Vector(numRandVar);
//			opserr << "randomVariables aft new =" <<randomVariables<<"\n";
			for( i=0; i<numRandVar; i++) (*randomVariables)(i)=0.0; 
//			opserr << " check randomVariables \n";
		}
	}

	if(arrayID_tmp !=0) {
		delete [] arrayID_tmp;
		arrayID_tmp=0;
	}
	if(temp!=0) {
		delete temp;
		temp=0;
	}

	// The random variable number is returned as a parameter ID
	return rvNumber;
}
double
NewDiscretizedRandomProcessSeries::getFactor(double time)
{
	/* fmk
	double r1=(*randomVariables)(0);
	double r2=(*randomVariables)(1);
	double r3=(*randomVariables)(2);
	double r4=(*randomVariables)(3);
	double r5=(*randomVariables)(4);
	double r6=(*randomVariables)(5);
	double r7=(*randomVariables)(6);
	double r8=(*randomVariables)(7);
	double r9=(*randomVariables)(8);
	double r10=(*randomVariables)(9);
	*/

	if (time == 0.0) {
		return 0.0;
	}
	else if (randomVariables == 0 || kickInTimes == 0) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactor(): " << endln
			<< " random variables or kick-in times vector(s) do not exist. " << endln;
		return 0.0;
	}
	else if (kickInTimes->Size() != randomVariables->Size() ) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactor(): " << endln
			<< " number of random variables is not the same as kick-in times. " << endln;
		return 0.0;
	}
	else {
		double sum1;
		double sum2;
//		int nrv = 0;
		double modFuncAmplitude, filterAmplitude;
		Filter *theFilter;

//		nrv = randomVariables->Size();
		for (int i=0; i<numRandVar; i++) {
			active[i]=false;
			if(time>(*kickInTimes)(i)-1.0e-7){
				active[i]=true;
			}
		}

		// Loop over all modulating functions
		sum1 = 0.0;
		double dtime;
		for (int k=0; k<numModFuncs; k++) {

			// Get value of modulating function number k at time t
			modFuncAmplitude = theModulatingFunctions[k]->getAmplitude(time);
			theFilter = theModulatingFunctions[k]->getFilter();

			// Number of discretizing random variables
//			nrv = randomVariables->Size();

			// Loop over all active rv's 
			sum2 = 0.0;
			for (int i=0; i<numRandVar; i++) {

				// Get value of filter for argument (t-ti)
				dtime=time-(*kickInTimes)(i);
				theFilter->setKickTime((*kickInTimes)(i));
				filterAmplitude = theFilter->getAmplitude(dtime);
				
				// Add contribution 'ui * hi'
				sum2 += (*randomVariables)(i) * filterAmplitude;

				// Break when we get to inactive rv's
				if (dtime <= -1.0e-7) {
					break;
				}
			}

			sum1 += sum2*modFuncAmplitude;
		}

		double result = mean + c*sum1;

		output << "time... ," << time << ","<< "value.. ," << sum1 << "\n"; 
		output.flush();

		return result;
	}
}


double
NewDiscretizedRandomProcessSeries::getFactorSensitivity(double time)
{
	// The parameterID has been set to the number of 
	// the random variable in question

	// So, do the same thing as above, just set x(i-1) equal to 1.0
	// for i==parameterID

	if (time == 0.0 || parameterID<0) {
		return 0.0;
	}
	else if (randomVariables == 0 || kickInTimes == 0) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " random variables or kick-in times vector(s) do not exist. " << endln;
		return 0.0;
	}
	else if (kickInTimes->Size() != randomVariables->Size() ) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " number of random variables is not the same as kick-in times. " << endln;
		return 0.0;
	}
	else {

		double sum1;
		double sum2;
//		int nrv = 0;
		double modFuncAmplitude;
		Filter *theFilter;

		// Loop over all modulating functions
		double dtime;
		sum1 = 0.0;
		for (int k=0; k<numModFuncs; k++) {

			// Get value of modulating function number k at time t
			modFuncAmplitude = theModulatingFunctions[k]->getAmplitude(time);
			theFilter = theModulatingFunctions[k]->getFilter();

			// Number of discretizing random variables
			//nrv = randomVariables->Size();

			// Loop over all rv's (even though some may be zero at this time)
			dtime=time-(*kickInTimes)(parameterID);
//			if(fabs(dtime)<=1.0e-7) dtime=0.0;
			sum2 = theFilter->getAmplitude(dtime);
			sum1 += sum2*modFuncAmplitude;
		}

//		double result=0.0;
//		if(fabs(time-(*kickInTimes)(parameterID-1))<= 1.0e-8) result=1.0;
		double result = mean + c*sum1;
		return result;
	}
}


double
NewDiscretizedRandomProcessSeries::getFactorSensitivity(double time, double ktime)
{
	// The parameterID has been set to the number of 
	// the random variable in question

	// So, do the same thing as above, just set x(i-1) equal to 1.0
	// for i==parameterID

	if (time == 0.0) {
		return 0.0;
	}
	else if (randomVariables == 0 || kickInTimes == 0) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " random variables or kick-in times vector(s) do not exist. " << endln;
		return 0.0;
	}
	else if (kickInTimes->Size() != randomVariables->Size() ) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " number of random variables is not the same as kick-in times. " << endln;
		return 0.0;
	}
	else {

		double sum1;
		double sum2;
//		int nrv = 0;
		double modFuncAmplitude;
		Filter *theFilter;

		// Loop over all modulating functions
		double dtime;
		sum1 = 0.0;
		for (int k=0; k<numModFuncs; k++) {

			// Get value of modulating function number k at time t
			modFuncAmplitude = theModulatingFunctions[k]->getAmplitude(time);
			theFilter = theModulatingFunctions[k]->getFilter();

			// Number of discretizing random variables
			//nrv = randomVariables->Size();

			// Loop over all rv's (even though some may be zero at this time)
			dtime=time-ktime;
			sum2 = theFilter->getAmplitude(dtime);
			sum1 += sum2*modFuncAmplitude;
		}

//		double result=0.0;
//		if(fabs(time-(*kickInTimes)(parameterID-1))<= 1.0e-8) result=1.0;
		double result = mean + c*sum1;
		return result;
	}
}

int
NewDiscretizedRandomProcessSeries::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}


int 
NewDiscretizedRandomProcessSeries::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return 0;    
}


void
NewDiscretizedRandomProcessSeries::Print(OPS_Stream &s, int flag)
{
}


int
NewDiscretizedRandomProcessSeries::updateParameter (int parameterID, Information &info)
{
	// In case the vector doesn't exist
	if (randomVariables == 0) {
		opserr << "DiscretizedRandomProcessSeries::updateParameter \n";
		opserr << "randomVariables ==0 require check \n";
		exit(-1);
//		randomVariables = new Vector(parameterID);
//		(*randomVariables)(parameterID-1) = info.theDouble;
	}
	// In case the vector isn't big enough
	else if (MaxRanVarID < parameterID) {
//	else if (randomVariables->Size() < parameterID) {
		// Store old values in a temporary vector
//		Vector temp = (*randomVariables);

		// Create a large enough vector
//		delete randomVariables;
//		randomVariables = new Vector(parameterID);

//		// Put in old values
//		for (int i=0; i<temp.Size(); i++) {
//			(*randomVariables)(i) = temp(i);
//		}

		// Put in new value
//		(*randomVariables)(parameterID-1) = info.theDouble;
		opserr << "DiscretizedRandomProcessSeries::updateParameter \n";
		opserr << "MaxRanVarID < parameterID require check \n";
		opserr << "MaxRanVarID="<< MaxRanVarID
			   << "parameterID="<< parameterID<<"\n";
		exit(-1);
//
	}
	else {
		int ii=arrayID[parameterID];
//		(*randomVariables)(parameterID-1) = info.theDouble;
		(*randomVariables)(ii) = info.theDouble;
	}

	return 0;
}

int
NewDiscretizedRandomProcessSeries::activateParameter(int passedParameterID)
{
	if( passedParameterID == 0 ){
//      ------ inactivate ---------
		parameterID=-1;
		for( int i=0; i<numRandVar; i++) active[i]=false; 
		return 0;
	}else if( passedParameterID < 0 ){
//      ------ inactivate ---------
		parameterID=-1;
		int itemp=arrayID[-passedParameterID];
		active[itemp]=false; 
		return 0;
	}else{
//      ------ activate ---------
		parameterID=arrayID[passedParameterID];
		if(active[parameterID]) return 1;
		else return 0;
	}

}
double NewDiscretizedRandomProcessSeries::getkickInTimes(int rvNum) const
{ 
	if( rvNum < 0 ){
		return (*kickInTimes)(rvNum);
	}else{
		if( rvNum >MaxRanVarID) return -1.0;  
		int ii=arrayID[rvNum];
		return (*kickInTimes)(ii);
	}
}
int NewDiscretizedRandomProcessSeries::getPulseSequentialID( int rvNum ) const
{
	if( rvNum >MaxRanVarID) return -1;  
	return arrayID[rvNum];
}


int
NewDiscretizedRandomProcessSeries::updateRV (int nrv, double value)
{
	// In case the vector doesn't exist
	if (randomVariables == 0) {
		opserr << "DiscretizedRandomProcessSeries::updateParameter \n";
		opserr << "randomVariables ==0 require check \n";
		exit(-1);
	}
	else if (MaxRanVarID < nrv) {
		opserr << "DiscretizedRandomProcessSeries::updateParameter \n";
		opserr << "MaxRanVarID < parameterID require check \n";
		opserr << "MaxRanVarID="<< MaxRanVarID
			   << "parameterID="<< parameterID<<"\n";
		exit(-1);
	}
	else {
		int ii=arrayID[nrv];
		(*randomVariables)(ii) = value;
	}

	return 0;
}
/*void 
DiscretizedRandomProcessSeries::sortPulse (void)
{
	if(numRandVar>1){
		for( int i=0; i<numRandVar-1; i++){
			for (int j=i+1; j<numRandVar; j++){
				double ki=(*kickInTimes)(i);
				double kj=(*kickInTimes)(j);
				if(ki==kj){
					opserr << " !! Warning in DiscretizedRandomProcessSeries::setParameter /n";
					opserr << " two random pulses are defined on a identical time instance /n";
					opserr << " time ="<<ki;
				}else if(ki>kj){
					double rvi = 
					(*kickInTimes)(i)=kj;
					(*kickInTimes)(j)=ki;
					int kki;
					for(int k=1;k<=MaxRanVarID;k++){
						if( arrayID[k] == i ){
							kki=k;
							break;
						}
					}
					int kkj;
					for(k=1;k<=MaxRanVarID;k++){
						if( arrayID[k] == j ){
							kkj=k;
							break;
						}
					}
					arrayID[kki]=j;
					arrayID[kkj]=i;
				}
			}
		}
	}
}
*/


