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
// $Date: 2001-06-13 05:06:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunction.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <LimitStateFunction.h>
#include <Vector.h>
#include <string.h>
#include <fstream.h>
#include <iomanip.h>


LimitStateFunction::LimitStateFunction(	int passedTag, 
									    char *passedExpression)
:ReliabilityDomainComponent(passedTag, 14727)
{
	tag = passedTag;
	expression = new char[100];
	strcpy(expression,passedExpression);
	FORMAnalysisPerformed = false;
	SimulationAnalysisPerformed = false;
	EvaluateLimitStateFunctionAtStartPointPerformed = false;
	PointFittingSORMAnalysisPerformed = false;
	CurvatureFittingSORMAnalysisPerformed = false;
	CurvaturesFromSearchAlgorithmSORMAnalysisPerformed = false;
}


LimitStateFunction::~LimitStateFunction()
{
	delete expression;
}


void
LimitStateFunction::Print(ostream &s, int flag)  
{
}



char *
LimitStateFunction::getExpression()
{
	return expression;
}



void
LimitStateFunction::printSummaryOfResults(ofstream &outputFile)
{

	if (FORMAnalysisPerformed || SimulationAnalysisPerformed || CurvaturesFromSearchAlgorithmSORMAnalysisPerformed) {

		outputFile << "#  LIMIT-STATE FUNCTION NUMBER "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<tag <<"                                   #" << endl;
		if (FORMAnalysisPerformed) {
		outputFile << "#   FORM Analysis:                    " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(19)<<FORMReliabilityIndexBeta 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<FORMProbabilityOfFailure_pf1 
			<< " #" << endl;
		}
		if (SimulationAnalysisPerformed) {
		outputFile << "#   Simulation Analysis:              " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(19)<<SimulationReliabilityIndexBeta 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<SimulationProbabilityOfFailure_pfsim 
			<< " #" << endl;
		}
		if (CurvaturesFromSearchAlgorithmSORMAnalysisPerformed) {
		outputFile << "#   SORM(1) (Breitung)                "
			<<setiosflags(ios::left)<<setprecision(5)<<setw(19)<<SORMUsingSearchBetaBreitung 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<SORMUsingSearchPf2Breitung 
			<< " #" << endl;
		}
		outputFile << "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #" << endl;
	}





/*  ideas for future output options

#######################################################################
#  SORM - POINT FITTING - ANALYSIS RESULTS                            #
#                                                                     #
#  Reliability index beta: ............................... 2.37544    #
#  Estimated probability of failure pf2: ................. 0.00876409 #
#                                                                     #
#######################################################################


#######################################################################
#  SORM - CURVATURE FITTING - ANALYSIS RESULTS                        #
#                                                                     #
#  Reliability index beta: ............................... 2.37544    #
#  Estimated probability of failure pf2: ................. 0.00876409 #
#                                                                     #
#######################################################################


#######################################################################
#  SENSITIVITY ANALYSIS RESULTS                                       #
#                                                                     #
#  d(beta)/d(parameter)                                               #
#  rv#    mean      stdv      par1     par2    par3     par4          #
#  1   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  2   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  3   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  4   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  5   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  6   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  7   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#                                                                     #
#  d(pf1)/d(parameter)                                                #
#  rv#    mean      stdv      par1     par2    par3     par4          #
#  1   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  2   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  3   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  4   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  5   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  6   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#  7   109019    -1.38242  -0.581826 -0.581826 -0.581826 -0.581826    #
#                                                                     #
#######################################################################


#######################################################################
#  INPUT DETAILS                                                      #
#                                                                     #
# rv#     mean        stdv      startPt                               #
#  1     200000.0    80000.0   200000.0                               #
#  2     200000.0    80000.0   200000.0                               #
#  3     200000.0    80000.0   200000.0                               #
#  4     200000.0    80000.0   200000.0                               #
#  5     200000.0    80000.0   200000.0                               #
#  6     200000.0    80000.0   200000.0                               #
#  7     200000.0    80000.0   200000.0                               #
#  8     200000.0    80000.0   200000.0                               #
#  9     200000.0    80000.0   200000.0                               #
#                                                                     #
#  correlation?  rvPositioners?  limitStateFunctions?                 #
#  analysisSetUp?                                                     #
#######################################################################



*/


}

void
LimitStateFunction::printResults(ofstream &outputFile)
{


	// PRINT EvaluateLimitStateFunctionAtStartPoint RESULTS
	if (EvaluateLimitStateFunctionAtStartPointPerformed) {
		outputFile << "#######################################################################" << endl;
		outputFile << "#  VALUE OF LIMIT-STATE FUNCTION "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<tag <<" AT START POINT                  #" << endl;
		outputFile << "#                                                                     #" << endl;
		outputFile << "#  Limit-state function value at start point: ......... " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<GFunValueAtStartPt 
			<< "  #" << endl;
		outputFile << "#                                                                     #" << endl;
		outputFile << "#######################################################################" << endl << endl << endl;
	}


	// PRINT FORM RESULTS
	if (FORMAnalysisPerformed) {
		outputFile << "#######################################################################" << endl;
		outputFile << "#  FORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<tag <<"            #" << endl;
		outputFile << "#                                                                     #" << endl;
		outputFile << "#  Limit-state function value at start point: ......... " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<GFunValueAtStartPt 
			<< "  #" << endl;
		outputFile << "#  Number of iterations: .............................. " 
			<<setiosflags(ios::left)<<setw(12)<<numberOfIterationsToFindDesignPoint 
			<< "  #" << endl;
		outputFile << "#  Reliability index beta: ............................ " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<FORMReliabilityIndexBeta 
			<< "  #" << endl;
		outputFile << "#  Estimated probability of failure pf1: .............. " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<FORMProbabilityOfFailure_pf1 
			<< "  #" << endl;
		outputFile << "#                                                                     #" << endl;
		outputFile << "# rv#     x*          u*        alpha     gamma     delta     eta     #" << endl;
		outputFile.setf( ios::scientific, ios::floatfield );
		for (int i=0;  i<designPoint_x_inOriginalSpace.Size(); i++) {
		outputFile << "#  " <<setw(3)<<(i+1)<<" ";
		outputFile.setf(ios::scientific, ios::floatfield);
		if (designPoint_x_inOriginalSpace(i)<0.0) { outputFile << "-"; }
		else { outputFile << " "; }
		outputFile <<setprecision(3)<<setw(11)<<fabs(designPoint_x_inOriginalSpace(i));
		if (designPoint_u_inStdNormalSpace(i)<0.0) { outputFile << "-"; }
		else { outputFile << " "; }
		outputFile <<setprecision(3)<<setw(11)<<fabs(designPoint_u_inStdNormalSpace(i));
		outputFile.unsetf( ios::scientific );
		outputFile.setf(ios::fixed, ios::floatfield);

		if (normalizedNegativeGradientVectorAlpha(i)<0.0) { outputFile << "-"; }
		else { outputFile << " "; }
		outputFile<<setprecision(5)<<setw(8)<<fabs(normalizedNegativeGradientVectorAlpha(i));

		if (importanceVectorGamma(i)<0.0) { outputFile << "-"; }
		else { outputFile << " "; }
		outputFile<<setprecision(5)<<setw(8)<<fabs(importanceVectorGamma(i));		
		
		outputFile<<"   -     ";
		outputFile<<"   -     ";
		outputFile<<"   #" << endl;
		}
		outputFile << "#                                                                     #" << endl;
		outputFile << "#######################################################################" << endl << endl << endl;
	}


	// PRINT SIMULATION RESULTS
	if (SimulationAnalysisPerformed) {
		outputFile << "#######################################################################" << endl;
		outputFile << "#  SIMULATION ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<tag <<"      #" << endl;
		outputFile << "#                                                                     #" << endl;
		outputFile << "#  Reliability index beta: ............................ " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<SimulationReliabilityIndexBeta 
			<< "  #" << endl;
		outputFile << "#  Estimated probability of failure pf_sim: ........... " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<SimulationProbabilityOfFailure_pfsim 
			<< "  #" << endl;
		outputFile << "#  Number of simulations: ............................. " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<NumberOfSimulations 
			<< "  #" << endl;
		outputFile << "#  Coefficient of variation (of pf): .................. " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<CoefficientOfVariationOfPfFromSimulation 
			<< "  #" << endl;
		outputFile << "#                                                                     #" << endl;
		outputFile << "#######################################################################" << endl << endl << endl;

	}

	// PRINT CurvaturesFromSearchAlgorithmSORMAnalysis RESULTS
	if (CurvaturesFromSearchAlgorithmSORMAnalysisPerformed) {
		outputFile << "#######################################################################" << endl;
		outputFile << "#  SORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<tag <<"            #" << endl;
		outputFile << "#  (Curvatures found from search algorithm.)                          #" << endl;
		outputFile << "#                                                                     #" << endl;
		outputFile << "#  Number of principal curvatures used: ............... " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<numberOfCurvatauresUsed
			<< "  #" << endl;
		outputFile << "#  Reliability index beta (impr. Breitung's formula):.. " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<SORMUsingSearchBetaBreitung 
			<< "  #" << endl;
		outputFile << "#  Corresponding estimated probability of failure pf2:.." 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<SORMUsingSearchPf2Breitung 
			<< "  #" << endl;
		outputFile << "#                                                                     #" << endl;
		outputFile << "#######################################################################" << endl << endl << endl;

	}




}
