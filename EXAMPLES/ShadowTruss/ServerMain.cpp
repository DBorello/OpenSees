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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-08-29 09:24:04 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/ShadowTruss/ServerMain.cpp,v $


// Written: fmk 08/99
//
// Purpose: this file contains a C++ main procedure to perform the analysis
// of example1 (found in most documents). One of the elements does it's 
// processing remotely. This is the Client program.

#include <stdlib.h>

#include <OPS_Globals.h>
#include <StandardStream.h>

#include <TCP_Socket.h>
#include <FEM_ObjectBroker.h>
#include "ActorTruss.h"

// init the global variabled defined in OPS_Globals.h
StandardStream sserr;
OPS_Stream &opserr = sserr;

double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;



// main routine
int main(int argc, char **argv)
{
    TCP_Socket *theChannel = new TCP_Socket(atoi(argv[1]));
    FEM_ObjectBroker *theBroker = new FEM_ObjectBroker();

    // create the actor & run it
    
    Actor *truss3 = new ActorTruss(*theChannel, *theBroker);        
    truss3->run();
    delete truss3;

    exit(0);
}	
	
