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
//# DATE:              19AUg2003
//# UPDATE HISTORY:    Sept 2003
//#
//#
//===============================================================================
#include <tcl.h>
#include <OPS_Globals.h>

//#include <stdlib.h>
//#include <string.h>

#include <Vector.h>

#include <Domain.h>

#include <ErrorHandler.h>
#include <TclModelBuilder.h>

#include <FiniteDeformationElastic3D.h>
#include <W.h>
#include <LogWEnergy.h>
#include <MooneyRivlinWEnergy.h>
#include <NeoHookeanWEnergy.h>
#include <OgdenWEnergy.h>
#include <SimoPisterWEnergy.h>
#include <OgdenSimoWEnergy.h>
#include <MooneyRivlinSimoWEnergy.h>


// the functions to create the component objects (defined at eof)
WEnergy *EvaluateWEnergyArgs(ClientData, Tcl_Interp *, TCL_Char *tclString);

// little function to free memory after invoke Tcl_SplitList
// note Tcl_Split list stores the array of pointers and the strings in
// one array, which is why Tcl_Free needs only be called on the array.

static void cleanup(TCL_Char **argv)
{
    Tcl_Free((char *) argv);
}

FiniteDeformationElastic3D *
TclModelBuilder_addFiniteDeformationElastic3D(ClientData clientData, Tcl_Interp *interp,  int argc,
          TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // create some empty pointers which we fill in as parse the command line
  int tag = 0;
  double rho_in = 0.0;
  WEnergy  *wenergy =0;

  int loc = eleArgStart;

  if (argc < 6) {
    opserr << "WARNING FiniteDeformation3DElastic -insufficient number of arguments\n";
    return 0;
  }

  if (Tcl_GetInt(interp, argv[loc+1], &tag) != TCL_OK)
  {
    opserr << "nDMaterial FiniteDeformationElastic3D - invalid tag " << argv[loc+1] << endln;
    return 0;
  }

  if ( (strcmp(argv[loc+2],"-wenergy") == 0) || (strcmp(argv[loc+2],"-WEnergy") == 0)
          || (strcmp(argv[loc+2],"-w") == 0) || (strcmp(argv[loc+2],"-W") == 0) )
  {
      wenergy = EvaluateWEnergyArgs(clientData, interp, argv[loc+3]);
      if (wenergy == 0)
      {
        opserr << "nDMaterial FiniteDeformationElastic3D - could not create a WEnergy from" << argv[loc+3] << endln;
        return 0;
      }
  }

  if (Tcl_GetDouble(interp, argv[loc+4], &rho_in) != TCL_OK)
  {
    opserr << "nDMaterial FiniteDeformationElastic3D - invalid rho " << argv[loc+4] << endln;
    return 0;
  }

  FiniteDeformationElastic3D *theMaterial = 0;
  if ( wenergy != 0)
    theMaterial = new FiniteDeformationElastic3D(tag, wenergy, rho_in);

  else
    opserr << "invalid number of args used to create a FiniteDeformationElastic3D material\n";

  return theMaterial;
}



WEnergy *EvaluateWEnergyArgs(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK)
    return 0;

  if (argc == 0)
    return 0;

  // now parse the list & construct the required object
  WEnergy *wenergy = 0;

  // #1 Logarithmic
  if ((strcmp(argv[0],"-Log") == 0) || (strcmp(argv[0],"-log") == 0) || (strcmp(argv[0],"-Logarithmic") == 0) )
    {
      double E_in = 0.0;
      double nu_in = 0.0;
      if (argc > 2)
      {
        if (Tcl_GetDouble(interp, argv[1], &E_in) != TCL_OK)
        {
          opserr << "Warning: invalid Young's Module for Logarithmic Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[2], &nu_in) != TCL_OK)
        {
          opserr << "Warning: invalid Poisson Ratio for Logarithmic Strain Energy Function.";
          return 0;
        }

      wenergy = new LogWEnergy(E_in, nu_in);
     }
   }

  // #2 Mooney Rivlin
    if ((strcmp(argv[0],"-MR") == 0) || (strcmp(argv[0],"-MooneyRivlin") == 0) )
    {
      double E_in = 0.0;
      double nu_in = 0.0;
      double c1_in = 0.0;
      double c2_in = 0.0;
      if (argc > 4)
      {
        if (Tcl_GetDouble(interp, argv[1], &E_in) != TCL_OK)
        {
          opserr << "Warning: invalid 1st material parameter for Mooney-Rivlin Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[2], &nu_in) != TCL_OK)
        {
          opserr << "Warning: invalid 2nd material parameter for Mooney-Rivlin Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[3], &c1_in) != TCL_OK)
        {
          opserr << "Warning: invalid 1st material parameter for Mooney-Rivlin Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[4], &c2_in) != TCL_OK)
        {
          opserr << "Warning: invalid 2nd material parameter for Mooney-Rivlin Strain Energy Function.";
          return 0;
        }

      wenergy = new MooneyRivlinWEnergy(E_in, nu_in, c1_in, c2_in);
     }
   }

  // #3 Neo Hookean
    if ((strcmp(argv[0],"-NH") == 0) || (strcmp(argv[0],"-NeoHookean") == 0) )
    {
      double E_in = 0.0;
      double nu_in = 0.0;
      if (argc > 2)
      {
        if (Tcl_GetDouble(interp, argv[1], &E_in) != TCL_OK)
        {
          opserr << "Warning: invalid Young's Module for Neo-Hookean Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[2], &nu_in) != TCL_OK)
        {
          opserr << "Warning: invalid Poisson Ratio for Neo-Hookean Strain Energy Function.";
          return 0;
        }

      wenergy = new NeoHookeanWEnergy(E_in, nu_in);
     }
   }

  // #4 Ogden
  if ((strcmp(argv[0],"-Ogden") == 0) || (strcmp(argv[0],"-ogden") == 0) )
    {  
      double E_in = 0.0;
      double nu_in = 0.0; 
      int N_in = 0;

      if (argc > 2*N_in+3)
      {
        if (Tcl_GetDouble(interp, argv[1], &E_in) != TCL_OK)
        {
          opserr << "Warning: invalid vector parameter number for Ogden Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[2], &nu_in) != TCL_OK)
        {
          opserr << "Warning: invalid vector parameter number for Ogden Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetInt(interp, argv[3], &N_in) != TCL_OK)
        {
          opserr << "Warning: invalid vector parameter number for Ogden Strain Energy Function.";
          return 0;
        }
      
      double *cr_in = new double[N_in];
      double *mur_in = new double[N_in];
         
      for (int i=0; i<N_in; i++)
        {
          if (Tcl_GetDouble(interp, argv[4+i], &cr_in[i]) != TCL_OK)
          {
            opserr << "Warning: invalid parameter for Ogden Strain Energy Function.";
            return 0;
          }
          if (Tcl_GetDouble(interp, argv[4+N_in+i], &mur_in[i]) != TCL_OK)
          {
            opserr << "Warning: invalid parameter for Ogden Strain Energy Function.";
            return 0;
          }
        }

      wenergy = new OgdenWEnergy(E_in, nu_in, N_in, cr_in, mur_in);
     }
   }

  // #5 Simo Pister
    if ((strcmp(argv[0],"-SP") == 0) || (strcmp(argv[0],"-SimoPister") == 0) )
    {
      double E_in = 0.0;
      double nu_in = 0.0;
      if (argc > 2)
      {
        if (Tcl_GetDouble(interp, argv[1], &E_in) != TCL_OK)
        {
          opserr << "Warning: invalid Young's Module for Simo-Pister Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[2], &nu_in) != TCL_OK)
        {
          opserr << "Warning: invalid Poisson Ratio for Simo-Pister Strain Energy Function.";
          return 0;
        }

      wenergy = new SimoPisterWEnergy(E_in, nu_in);
     }
   }

  // #6 Ogden-Simo
  if ((strcmp(argv[0],"-OgdenSimo") == 0) || (strcmp(argv[0],"-OS") == 0) )
    {  
      double E_in = 0.0;
      double nu_in = 0.0; 
      int N_in = 0;

      if (argc > 2*N_in+3)
      {
        if (Tcl_GetDouble(interp, argv[1], &E_in) != TCL_OK)
        {
          opserr << "Warning: invalid vector parameter number for Ogden Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[2], &nu_in) != TCL_OK)
        {
          opserr << "Warning: invalid vector parameter number for Ogden Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetInt(interp, argv[3], &N_in) != TCL_OK)
        {
          opserr << "Warning: invalid vector parameter number for Ogden Strain Energy Function.";
          return 0;
        }
      
      double *cr_in = new double[N_in];
      double *mur_in = new double[N_in];
         
      for (int i=0; i<N_in; i++)
        {
          if (Tcl_GetDouble(interp, argv[4+i], &cr_in[i]) != TCL_OK)
          {
            opserr << "Warning: invalid parameter for Ogden Strain Energy Function.";
            return 0;
          }
          if (Tcl_GetDouble(interp, argv[4+N_in+i], &mur_in[i]) != TCL_OK)
          {
            opserr << "Warning: invalid parameter for Ogden Strain Energy Function.";
            return 0;
          }
        }

      wenergy = new OgdenSimoWEnergy(E_in, nu_in, N_in, cr_in, mur_in);
     }
   }

  // #7 Mooney Rivlin Simo
    if ((strcmp(argv[0],"-MRS") == 0) || (strcmp(argv[0],"-MooneyRivlinSimo") == 0) )
    {
      double E_in = 0.0;
      double nu_in = 0.0;
      double c1_in = 0.0;
      double c2_in = 0.0;
      if (argc > 4)
      {
        if (Tcl_GetDouble(interp, argv[1], &E_in) != TCL_OK)
        {
          opserr << "Warning: invalid 1st material parameter for Mooney-Rivlin Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[2], &nu_in) != TCL_OK)
        {
          opserr << "Warning: invalid 2nd material parameter for Mooney-Rivlin Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[3], &c1_in) != TCL_OK)
        {
          opserr << "Warning: invalid 1st material parameter for Mooney-Rivlin Strain Energy Function.";
          return 0;
        }

        if (Tcl_GetDouble(interp, argv[4], &c2_in) != TCL_OK)
        {
          opserr << "Warning: invalid 2nd material parameter for Mooney-Rivlin Strain Energy Function.";
          return 0;
        }

      wenergy = new MooneyRivlinSimoWEnergy(E_in, nu_in, c1_in, c2_in);
     }
   }

  cleanup(argv);
  return wenergy;

}




