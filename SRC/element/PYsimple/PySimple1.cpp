/* *********************************************************************
**    Module:	PySimple1.cpp 
**
**    Purpose:	Provide a simple p-y spring for OpenSees.
**
**    Developed by Ross W. Boulanger and Boris Jeremic
**    (C) Copyright 2001, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2001/12/15
// $Source: /OpenSees/SRC/material/uniaxial/PySimple1.cpp

// Written: RWB
// Created: Dec 2001
// Revision: A
//
// Description: This file contains the class implementation for PySimple1

#include <PySimple1.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

/////////////////////////////////////////////////////////////////////
//	Constructor with data

PySimple1::PySimple1(int tag, int soil, double p_ult, double y_50,
				 double dragratio, double dash_pot)
:UniaxialMaterial(tag,MAT_TAG_PySimple1),
 soilType(soil), pult(p_ult), y50(y_50), drag(dragratio), dashpot(dash_pot)
{
	// Initialize PySimple variables and history variables
	//
    this->revertToStart();
}

/////////////////////////////////////////////////////////////////////
//	Default constructor

PySimple1::PySimple1()
:UniaxialMaterial(0,MAT_TAG_PySimple1),
 soilType(0), pult(0.0), y50(0.0), drag(0.0), dashpot(0.0)
{
	// Initialize variables
	//
	this->revertToStart();
}

/////////////////////////////////////////////////////////////////////
//	Default destructor
PySimple1::~PySimple1()
{
    // Does nothing
}

/////////////////////////////////////////////////////////////////////
void PySimple1::getGap(double ylast, double dy, double dy_old)
{
	// For stability in Closure spring, may limit "dy" step size to avoid
	// overshooting on the closing of this gap.
	//
	TGap_y = ylast + dy;
	if(TGap_y > TClose_yright) {dy = 0.75*(TClose_yright - ylast);}
	if(TGap_y < TClose_yleft)  {dy = 0.75*(TClose_yleft  - ylast);}

	// Limit "dy" step size if it is oscillating in sign and not shrinking
	//
	if(dy*dy_old < 0.0 && fabs(dy/dy_old) > 0.5) dy = -dy_old/2.0;
	
	// Combine the Drag and Closure elements in parallel, starting by
	// resetting TGap_y in case the step size was limited.
	//
	TGap_y   = ylast + dy;
	getClosure(ylast,dy);
	getDrag(ylast,dy);
	TGap_p = TDrag_p + TClose_p;
	TGap_tang = TDrag_tang + TClose_tang;

	return;
}

/////////////////////////////////////////////////////////////////////
void PySimple1::getFarField(double y)
{
	TFar_y   = y;
	TFar_tang= TFar_tang;
	TFar_p   = TFar_tang * TFar_y;

	return;
}

/////////////////////////////////////////////////////////////////////
void PySimple1::getClosure(double ylast, double dy)
{
	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TClose_yleft != CClose_yleft)  TClose_yleft = CClose_yleft;
	if(TClose_yright!= CClose_yright) TClose_yright= CClose_yright;

	// Check if plastic deformation in Near Field should cause gap expansion
	//
	TClose_y = ylast + dy;
	double yrebound=1.5*y50;
	if(TNF_y+TClose_y > -TClose_yleft + yrebound)
		TClose_yleft=-(TNF_y+TClose_y) + yrebound;
	if(TNF_y+TClose_y < -TClose_yright - yrebound)
		TClose_yright=-(TNF_y+TClose_y) - yrebound;

	// Spring force and tangent stiffness
	//
	TClose_p=1.8*pult*(y50/50.0)*(pow(y50/50.0 + TClose_yright - TClose_y,-1.0)
		-pow(y50/50.0 + TClose_y - TClose_yleft,-1.0));
	TClose_tang=1.8*pult*(y50/50.0)*(pow(y50/50.0+ TClose_yright - TClose_y,-2.0)
		+pow(y50/50.0 + TClose_y - TClose_yleft,-2.0));

	// Ensure that tangent not zero or negative.
	//	
	if(TClose_tang <= 1.0e-2*pult/y50) {TClose_tang = 1.0e-2*pult/y50;}

	return;
}

/////////////////////////////////////////////////////////////////////
void PySimple1::getDrag(double ylast, double dy)
{
	TDrag_y = ylast + dy;
	double pmax=drag*pult;
	double dyTotal=TDrag_y - CDrag_y;

	// Treat as elastic if dyTotal is below tolerance
	//
	if(fabs(dyTotal*TDrag_tang/pult) < 10.0*tolerance) 
	{
		TDrag_p = TDrag_p + dy*TDrag_tang;
		if(fabs(TDrag_p) >=pmax) TDrag_p =(TDrag_p/fabs(TDrag_p))*(1.0-1.0e-8)*pmax;
		return;
	}
	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TDrag_pin != CDrag_pin)
	{
		TDrag_pin = CDrag_pin;
		TDrag_yin = CDrag_yin;
	}

	// Change from positive to negative direction
	//
	if(CDrag_y > CDrag_yin && dyTotal < 0.0)
	{
		TDrag_pin = CDrag_p;
		TDrag_yin = CDrag_y;
	}
	// Change from negative to positive direction
	//
	if(CDrag_y < CDrag_yin && dyTotal > 0.0)
	{
		TDrag_pin = CDrag_p;
		TDrag_yin = CDrag_y;
	}
	
	// Positive loading
	//
	if(dyTotal >= 0.0)
	{
		TDrag_p=pmax-(pmax-TDrag_pin)*pow(y50/2.0,nd)
					*pow(y50/2.0 + TDrag_y - TDrag_yin,-nd);
		TDrag_tang=nd*(pmax-TDrag_pin)*pow(y50/2.0,nd)
					*pow(y50/2.0 + TDrag_y - TDrag_yin,-nd-1.0);
	}
	// Negative loading
	//
	if(dyTotal < 0.0)
	{
		TDrag_p=-pmax+(pmax+TDrag_pin)*pow(y50/2.0,nd)
					*pow(y50/2.0 - TDrag_y + TDrag_yin,-nd);
		TDrag_tang=nd*(pmax-TDrag_pin)*pow(y50/2.0,nd)
					*pow(y50/2.0 - TDrag_y + TDrag_yin,-nd-1.0);
	}
	// Ensure that |p|<pmax and tangent not zero or negative.
	//
	if(fabs(TDrag_p) >=pmax) {
		TDrag_p =(TDrag_p/fabs(TDrag_p))*(1.0-tolerance)*pmax;}
	if(TDrag_tang <=1.0e-2*pult/y50) TDrag_tang = 1.0e-2*pult/y50;

	return;
}

/////////////////////////////////////////////////////////////////////
void PySimple1::getNearField(double ylast, double dy, double dy_old)
{
	// Limit "dy" step size if it is oscillating in sign and not shrinking
	//
	if(dy*dy_old < 0.0 && fabs(dy/dy_old) > 0.5) dy = -dy_old/2.0;

	// Set "dy" so "y" is at middle of elastic zone if oscillation is large.
	// Note that this criteria is based on the min step size in setTrialStrain.
	//
	if(dy*dy_old < -y50*y50) dy = (TNFyinr + TNFyinl)/2.0 - ylast;
	
	// Establish trial "y" and direction of loading (with NFdy) for entire step
	//
	TNF_y = ylast + dy;
	double NFdy = TNF_y - CNF_y;

	// Treat as elastic if NFdy is below tolerance
	//
	if(fabs(NFdy*TNF_tang/pult) < 10.0*tolerance) 
	{
		TNF_p = TNF_p + dy*TNF_tang;
		if(fabs(TNF_p) >=pult) TNF_p=(TNF_p/fabs(TNF_p))*(1.0-tolerance)*pult;
		return;
	}

	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TNFpinr != CNFpinr || TNFpinl != CNFpinl)
	{
		TNFpinr = CNFpinr;
		TNFpinl = CNFpinl;
		TNFyinr = CNFyinr;
		TNFyinl = CNFyinl;
	}

	// For stability, may have to limit "dy" step size if direction changed.
	//
	bool changeDirection = false;
	
	// Direction change from a yield point triggers new Elastic range
	//
	double minE = 0.25;		// The min Elastic range on +/- side of p=0
	if(CNF_p > CNFpinr && NFdy <0.0){				// from pos to neg
		changeDirection = true;
		TNFpinr = CNF_p;
		TNFpinl = TNFpinr - 2.0*pult*Elast;
		if (TNFpinl > -minE*pult) {TNFpinl = -minE*pult;}
		TNFyinr = CNF_y;
		TNFyinl = TNFyinr - (TNFpinr-TNFpinl)/NFkrig; 
	}
	if(CNF_p < CNFpinl && NFdy > 0.0){				// from neg to pos
		changeDirection = true;
		TNFpinl = CNF_p;
		TNFpinr = TNFpinl + 2.0*pult*Elast;
		if (TNFpinr < minE*pult) {TNFpinr = minE*pult;}
		TNFyinl = CNF_y;
		TNFyinr = TNFyinl + (TNFpinr-TNFpinl)/NFkrig; 
	}
	// Now if there was a change in direction, limit the step size "dy"
	//
	if(changeDirection == true) {
		double maxdy = 0.25*pult/NFkrig;
		if(fabs(dy) > maxdy) dy = (dy/fabs(dy))*maxdy;
	}

	// Now, establish the trial value of "y" for use in this function call.
	//
	TNF_y = ylast + dy;

	// Postive loading
	//
	if(NFdy >= 0.0){
		// Check if elastic using y < yinr
		if(TNF_y <= TNFyinr){							// stays elastic
			TNF_tang = NFkrig;
			TNF_p = TNFpinl + (TNF_y - TNFyinl)*NFkrig;
		}
		else {
			TNF_tang = np * (pult-TNFpinr) * pow(yref,np) 
				* pow(yref - TNFyinr + TNF_y, -np-1.0);
			TNF_p = pult - (pult-TNFpinr)* pow(yref/(yref-TNFyinr+TNF_y),np);
		}
	}

	// Negative loading
	//
	if(NFdy < 0.0){
		// Check if elastic using y < yinl
		if(TNF_y >= TNFyinl){							// stays elastic
			TNF_tang = NFkrig;
			TNF_p = TNFpinr + (TNF_y - TNFyinr)*NFkrig;
		}
		else {
			TNF_tang = np * (pult+TNFpinl) * pow(yref,np) 
				* pow(yref + TNFyinl - TNF_y, -np-1.0);
			TNF_p = -pult + (pult+TNFpinl)* pow(yref/(yref+TNFyinl-TNF_y),np);
		}
	}

	// Ensure that |p|<pult and tangent not zero or negative.
	//
	if(fabs(TNF_p) >=pult) TNF_p=(TNF_p/fabs(TNF_p))*(1.0-tolerance)*pult;
	if(TNF_tang <= 1.0e-2*pult/y50) TNF_tang = 1.0e-2*pult/y50;

    return;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple1::setTrialStrain (double newy, double yRate)
{
	// Set trial values for displacement and load in the material
	// based on the last Tangent modulus.
	//
	double dy = newy - Ty;
	double dp = Ttangent * dy;
	TyRate    = yRate;

	// Limit the size of step (dy or dp) that can be imposed. Prevents
	// numerical difficulties upon load reversal at high loads
	// where a soft loading modulus becomes a stiff unloading modulus.
	//
	int numSteps = 1;
	double stepSize = 1.0;
	if(fabs(dp/pult) > 0.1) numSteps = 1 + int(fabs(dp/(0.1*pult)));
	if(fabs(dy/y50)  > 0.25 ) numSteps = 1 + int(fabs(dy/(0.25*y50)));
	stepSize = 1.0/float(numSteps);
	if(numSteps > 100) numSteps = 100;

	dy = stepSize * dy;

	// Main loop over the required number of substeps
	//
	for(int istep=1; istep <= numSteps; istep++)
	{
		Ty = Ty + dy;
		dp = Ttangent * dy;
		
	// May substep within Gap or NearField element if oscillating, which can happen
	// when they jump from soft to stiff.
	//
		double dy_gap_old = ((Tp + dp) - TGap_p)/TGap_tang;
		double dy_nf_old  = ((Tp + dp) - TNF_p) /TNF_tang;

	// Iterate to distribute displacement among the series components.
	// Use the incremental iterative strain & iterate at this strain.
	//
	for (int j=1; j < maxIterations; j++)
	{
		Tp = Tp + dp;
		if(fabs(Tp) >(1.0-tolerance)*pult) Tp=(1.0-tolerance)*pult*(Tp/fabs(Tp));

		// Stress & strain update in Near Field element
		double dy_nf = (Tp - TNF_p)/TNF_tang;
		getNearField(TNF_y,dy_nf,dy_nf_old);
		
		// Residuals in Near Field element
		double p_unbalance = Tp - TNF_p;
		double yres_nf = (Tp - TNF_p)/TNF_tang;
		dy_nf_old = dy_nf;

		// Stress & strain update in Gap element
		double dy_gap = (Tp - TGap_p)/TGap_tang;
		getGap(TGap_y,dy_gap,dy_gap_old);

		// Residuals in Gap element
		double p_unbalance2 = Tp - TGap_p;
		double yres_gap = (Tp - TGap_p)/TGap_tang;
		dy_gap_old = dy_gap;

		// Stress & strain update in Far Field element
		double dy_far = (Tp - TFar_p)/TFar_tang;
		TFar_y = TFar_y + dy_far;
		getFarField(TFar_y);

		// Residuals in Far Field element
		double p_unbalance3 = Tp - TFar_p;
		double yres_far = (Tp - TFar_p)/TFar_tang;

		// Update the combined tangent modulus
		Ttangent = pow(1.0/TGap_tang + 1.0/TNF_tang + 1.0/TFar_tang, -1.0);

		// Residual deformation across combined element
		double dv = Ty - (TGap_y + yres_gap)
			- (TNF_y + yres_nf) - (TFar_y + yres_far);

		// Residual "p" increment 
		dp = Ttangent * dv;

		// Test for convergence
		double psum = fabs(p_unbalance) + fabs(p_unbalance2) + fabs(p_unbalance3);
		if(psum/pult < tolerance) break;
	}
	}

	return 0;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple1::getStress(void)
{
	// Dashpot force is only due to velocity in the far field.
	// If converged, proportion by Tangents.
	// If not converged, proportion by ratio of displacements in components.
	//
	double ratio_disp =(1.0/TFar_tang)/(1.0/TFar_tang + 1.0/TNF_tang + 1.0/TGap_tang);
	if(Ty != Cy) {
		ratio_disp = (TFar_y - CFar_y)/(Ty - Cy);
		if(ratio_disp > 1.0) ratio_disp = 1.0;
		if(ratio_disp < 0.0) ratio_disp = 0.0;
	}
	double dashForce = dashpot * TyRate * ratio_disp;

	// Limit the combined force to pult.
	//
	if(fabs(Tp + dashForce) >= (1.0-tolerance)*pult)
		return (1.0-tolerance)*pult*(Tp+dashForce)/fabs(Tp+dashForce);
	else return Tp + dashForce;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple1::getTangent(void)
{
    return this->Ttangent;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple1::getDampTangent(void)
{
	// Damping tangent is produced only by the far field component.
	// If converged, proportion by Tangents.
	// If not converged, proportion by ratio of displacements in components.
	//
	double ratio_disp =(1.0/TFar_tang)/(1.0/TFar_tang + 1.0/TNF_tang + 1.0/TGap_tang);
	if(Ty != Cy) {
		ratio_disp = (TFar_y - CFar_y)/(Ty - Cy);
		if(ratio_disp > 1.0) ratio_disp = 1.0;
		if(ratio_disp < 0.0) ratio_disp = 0.0;
	}

	double DampTangent = dashpot * ratio_disp;

	// Minimum damping tangent referenced against Farfield spring
	//
	if(DampTangent < TFar_tang * 1.0e-12) DampTangent = TFar_tang * 1.0e-12;

	return DampTangent;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple1::getStrain(void)
{
    return this->Ty;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple1::getStrainRate(void)
{
    return this->TyRate;
}
/////////////////////////////////////////////////////////////////////
int 
PySimple1::commitState(void)
{
	// Commit trial history variable -- Combined element
    Cy       = Ty;
    Cp       = Tp;
    Ctangent = Ttangent;
    
	// Commit trial history variables for Near Field component
	CNFpinr   = TNFpinr;
	CNFpinl   = TNFpinl; 
	CNFyinr   = TNFyinr;
	CNFyinl   = TNFyinl;	
	CNF_p     = TNF_p;
	CNF_y     = TNF_y;
	CNF_tang  = TNF_tang;

	// Commit trial history variables for Drag component
	CDrag_pin = TDrag_pin;
	CDrag_yin = TDrag_yin;
	CDrag_p   = TDrag_p;
	CDrag_y   = TDrag_y;
	CDrag_tang= TDrag_tang;

	// Commit trial history variables for Closure component
	CClose_yleft  = TClose_yleft;
	CClose_yright = TClose_yright;
	CClose_p      = TClose_p;
	CClose_y      = TClose_y;
	CClose_tang   = TClose_tang;

	// Commit trial history variables for the Gap
	CGap_y    = TGap_y;
	CGap_p    = TGap_p;
	CGap_tang = TGap_tang;
    
	// Commit trial history variables for the Far Field
	CFar_y    = TFar_y;
	CFar_p    = TFar_p;
	CFar_tang = TFar_tang;
    
    return 0;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple1::revertToLastCommit(void)
{
	// Nothing to do here
    return 0;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple1::revertToStart(void)
{
	maxIterations = 40;
	tolerance     = 1.0e-12;

	// Reset gap "drag" if zero (or negative).
	//
	if(drag <= tolerance) drag = tolerance;

	// Only allow zero or positive dashpot values
	//
	if(dashpot < 0.0) dashpot = 0.0;

	// Do not allow zero or negative values for y50 or pult.
	//
	if(pult <= 0.0 || y50 <= 0.0)
		g3ErrorHandler->fatal("%s -- only accepts positive nonzero pult and y50",
			"PySimple1::PySimple1");

	// Initialize variables for Near Field rigid-plastic spring
	//
	if(this->soilType ==1) {
		this->yref  = 10.0*this->y50;
		this->np    = 5.0;
		this->Elast = 0.35;
		this->nd    = 1.0;
		TFar_tang   = pult/(8.0*pow(Elast,2.0)*y50);
	}
	else if (this->soilType == 2){
		this->yref  = 0.5*this->y50;
		this->np    = 2.0;
		this->Elast = 0.2;
		this->nd    = 1.0;
		// This TFar_tang assumes Elast=0.2, but changes very little for other
		// reasonable Elast values. i.e., API curves are quite linear initially.
		TFar_tang   = 0.542*pult/y50;
	}
	else{
		g3ErrorHandler->fatal("%s -- only accepts soilType of 1 or 2",
			"PySimple1::PySimple1");
	}

	// Far Field components: TFar_tang was set under "soil type" statements.
	//
	TFar_p  = 0.0;
	TFar_y  = 0.0;

	// Near Field components
	//
	this->NFkrig  = 100.0 * (0.5 * this->pult) / this->y50;
    this->TNFpinr = this->Elast*this->pult;
	this->TNFpinl = -this->TNFpinr;
	this->TNFyinr = this->TNFpinr / this->NFkrig;
	this->TNFyinl = -this->TNFyinr;
	this->TNF_p   = 0.0;
	this->TNF_y   = 0.0;
	this->TNF_tang= this->NFkrig;

	// Drag components
	//
	TDrag_pin = 0.0;
	TDrag_yin = 0.0;
	TDrag_p   = 0.0;
	TDrag_y   = 0.0;
	TDrag_tang= nd*(pult*drag-TDrag_p)*pow(y50/2.0,nd)
					*pow(y50/2.0 - TDrag_y + TDrag_yin,-nd-1.0);

	// Closure components
	//
	TClose_yleft = -y50/100.0;
	TClose_yright=  y50/100.0;
	TClose_p     = 0.0; 
	TClose_y     = 0.0;
	TClose_tang  = 1.8*pult*(y50/50.0)*(pow(y50/50.0+ TClose_yright - TClose_y,-2.0)
		+pow(y50/50.0 + TClose_y - TClose_yleft,-2.0));

	// Gap (Drag + Closure in parallel)
	//
	TGap_y   = 0.0;
	TGap_p   = 0.0;
	TGap_tang= TClose_tang + TDrag_tang;

	// Entire element (Far field + Near field + Gap in series)
	//
	Ty       = 0.0;
	Tp       = 0.0;
	Ttangent = pow(1.0/TGap_tang + 1.0/TNF_tang + 1.0/TFar_tang, -1.0);
	TyRate   = 0.0;

	// Now get all the committed variables initiated
	//
	this->commitState();

    return 0;
}

/////////////////////////////////////////////////////////////////////
UniaxialMaterial *
PySimple1::getCopy(void)
{
    PySimple1 *theCopy =
	new PySimple1(this->getTag(),soilType,pult,y50,drag,dashpot);

	// Copy parameters
	theCopy->yref    = yref;
	theCopy->np      = np;
	theCopy->Elast   = Elast;
	theCopy->nd      = nd;

	// Copy internal parameters or constants
	theCopy->NFkrig  = NFkrig;
    
	// Copy committed history variables for Near Field
    theCopy->CNFpinr = CNFpinr;
    theCopy->CNFpinl = CNFpinl;
    theCopy->CNFyinr = CNFyinr;
    theCopy->CNFyinl = CNFyinl;
    theCopy->CNF_p   = CNF_p;
    theCopy->CNF_y   = CNF_y;
    theCopy->CNF_tang= CNF_tang;

	// Copy trial history variables for Near Field
    theCopy->TNFpinr = TNFpinr;
    theCopy->TNFpinl = TNFpinl;
    theCopy->TNFyinr = TNFyinr;
    theCopy->TNFyinl = TNFyinl;
    theCopy->TNF_p   = TNF_p;
    theCopy->TNF_y   = TNF_y;
    theCopy->TNF_tang= TNF_tang;

	// Copy committed history variables for Drag component
    theCopy->CDrag_pin = CDrag_pin;
    theCopy->CDrag_yin = CDrag_yin;
	theCopy->CDrag_p   = CDrag_p;
	theCopy->CDrag_y   = CDrag_y;
	theCopy->CDrag_tang= CDrag_tang;

	// Copy trial history variables for Drag component
    theCopy->TDrag_pin = TDrag_pin;
    theCopy->TDrag_yin = TDrag_yin;
	theCopy->TDrag_p   = TDrag_p;
	theCopy->TDrag_y   = TDrag_y;
	theCopy->TDrag_tang= TDrag_tang;

	// Copy committed history variables for Closure component
	theCopy->CClose_yleft = CClose_yleft;
	theCopy->CClose_yright= CClose_yright;
	theCopy->CClose_p     = CClose_p; 
	theCopy->CClose_y     = CClose_y;
	theCopy->CClose_tang  = CClose_tang;
	
	// Copy trail history variables for Closure component	
	theCopy->TClose_yleft = TClose_yleft;
	theCopy->TClose_yright= TClose_yright;
	theCopy->TClose_p     = TClose_p; 
	theCopy->TClose_y     = TClose_y;
	theCopy->TClose_tang  = TClose_tang;

	// Copy committed history variables for Gap component
	theCopy->CGap_y   = CGap_y;
	theCopy->CGap_p   = CGap_p;
	theCopy->CGap_tang= CGap_tang;	

	// Copy trial history variables for Gap component
	theCopy->TGap_y   = TGap_y;
	theCopy->TGap_p   = TGap_p;
	theCopy->TGap_tang= TGap_tang;	

	// Copy committed history variables for Far Field component
	theCopy->CFar_y   = CFar_y;
	theCopy->CFar_p   = CFar_p;
	theCopy->CFar_tang= CFar_tang;	

	// Copy trial history variables for Far Field component
	theCopy->TFar_y   = TFar_y;
	theCopy->TFar_p   = TFar_p;
	theCopy->TFar_tang= TFar_tang;	

	// Copy committed history variables for Entire Material
	theCopy->Cy       = Cy;
	theCopy->Cp       = Cp;
	theCopy->Ctangent = Ctangent;

	// Copy trial history variables for Entire Material
	theCopy->Ty       = Ty;
	theCopy->Tp       = Tp;
	theCopy->Ttangent = Ttangent;
	theCopy->TyRate   = TyRate;

    return theCopy;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple1::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(38);
  
  data(0) = this->getTag();
  data(1) = soilType;
  data(2) = pult;
  data(3) = y50;
  data(4) = drag;
  data(5) = dashpot;
  data(6) = yref;
  data(7) = np;
  data(8) = Elast;
  data(9) = nd;
  data(10)= NFkrig;

  data(11) = CNFpinr;
  data(12) = CNFpinl;
  data(13) = CNFyinr;
  data(14) = CNFyinl;
  data(15) = CNF_p;
  data(16) = CNF_y;
  data(17) = CNF_tang;

  data(18) = CDrag_pin;
  data(19) = CDrag_yin;
  data(20) = CDrag_p;
  data(21) = CDrag_y;
  data(22) = CDrag_tang;

  data(23) = CClose_yleft;
  data(24) = CClose_yright;
  data(25) = CClose_p;
  data(26) = CClose_y;
  data(27) = CClose_tang;

  data(28) = CGap_y;
  data(29) = CGap_p;
  data(30) = CGap_tang;

  data(31) = CFar_y;
  data(32) = CFar_p;
  data(33) = CFar_tang;

  data(34) = Cy;
  data(35) = Cp;
  data(36) = Ctangent;
  data(37) = TyRate;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    cerr << "PySimple1::sendSelf() - failed to send data\n";

  return res;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple1::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(38);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      cerr << "PySimple1::recvSelf() - failed to receive data\n";
      CNF_tang = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
	soilType = data(1);
	pult     = data(2);
	y50      = data(3);
	drag     = data(4);
	dashpot  = data(5);
	yref     = data(6);
	np       = data(7);
	Elast    = data(8);
	nd       = data(9);
	NFkrig   = data(10);

	CNFpinr  = data(11);
	CNFpinl  = data(12);
	CNFyinr  = data(13);
	CNFyinl  = data(14);
	CNF_p    = data(15);
	CNF_y    = data(16);
	CNF_tang = data(17);

	CDrag_pin = data(18);
    CDrag_yin = data(19);
	CDrag_p   = data(20);
	CDrag_y   = data(21);
	CDrag_tang= data(22);

	CClose_yleft = data(23);
	CClose_yright= data(24);
	CClose_p     = data(25);
	CClose_y     = data(26);
	CClose_tang  = data(27);

	CGap_y    = data(28);
	CGap_p    = data(29);
	CGap_tang = data(30);

	CFar_y    = data(31);
	CFar_p    = data(32);
	CFar_tang = data(33);

	Cy        = data(34);
	Cp        = data(35);
	Ctangent  = data(36);
	TyRate    = data(37);
  }
    
  return res;
}

/////////////////////////////////////////////////////////////////////
void 
PySimple1::Print(ostream &s, int flag)
{
    s << "PySimple1, tag: " << this->getTag() << endl;
    s << "  soilType: " << soilType << endl;
    s << "  pult: " << pult << endl;
    s << "  y50: " << y50 << endl;
    s << "  drag: " << drag << endl;
	s << "  dashpot: " << dashpot << endl;
}

/////////////////////////////////////////////////////////////////////

