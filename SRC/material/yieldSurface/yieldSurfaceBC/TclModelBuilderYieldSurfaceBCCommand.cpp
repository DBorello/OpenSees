#include "YieldSurface_BC.h"
#include <TclModelBuilder.h>
#include <string.h>
#include <Vector.h>

#include "Attalla2D.h"
#include "Orbison2D.h"
#include "Hajjar2D.h"
#include "ElTawil2D.h"
#include "ElTawil2DUnSym.h"


static void printCommand(int argc, char **argv)
{
    cerr << "Input command: ";
    for (int i=0; i<argc; i++)
	cerr << argv[i] << " ";
    cerr << endl;
}

int
TclModelBuilderYieldSurface_BCCommand (ClientData clienData, Tcl_Interp *interp, int argc,
					char **argv, TclModelBuilder *theBuilder)
{
    // Make sure there is a minimum number of arguments
    if (argc < 3) {
	cerr << "WARNING insufficient number of uniaxial material arguments\n";
	cerr << "Want: yieldSurfaceBC type? tag? <specific material args>" << endl;
	return TCL_ERROR;
    }

    // Pointer to a ys that will be added to the model builder
    YieldSurface_BC *theYS = 0;

    // Check argv[2] for uniaxial material type
    if (strcmp(argv[1],"Orbison2D") == 0)
	{
		if (argc < 6)
		{
			cerr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			// Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
			cerr << "Want: yieldSurfaceBC Orbison2D tag? xCap? yCap? ys_model_tag?" << endl;
			return TCL_ERROR;
		}

		int tag;
		double xCap, yCap;
		// int matID1, matID2;
		int modelID;
//		double isoRatio;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC Orbison2D tag" << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &xCap) != TCL_OK)
		{
			cerr << "WARNING invalid xCap\n";
			cerr << "yieldSurfaceBC Orbison2D tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &yCap) != TCL_OK)
		{
			cerr << "WARNING invalid yCap\n";
			cerr << "yieldSurfaceBC Orbison2D tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[5], &modelID) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC Orbison2D matID1" << modelID << endl;
			return TCL_ERROR;
		}

		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			cerr << "WARNING yieldSurfaceBC Orbison2D no ys_model exixts with tag: " << modelID << endl;
			return TCL_ERROR;
		}

		// Parsing was successful, allocate the material
		theYS = new Orbison2D(tag, xCap, yCap, *theModel);
    }

    else if (strcmp(argv[1],"ElTawil2D") == 0)
	{
		if (argc < 7)
		{
			cerr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			// Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
			cerr << "Want: yieldSurfaceBC ElTawil2D tag? xCap? yCap? ys_model_tag?" << endl;
			return TCL_ERROR;
		}

		int tag;
		double xBal, yBal;
		double yPos, yNeg;

		int modelID;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC ElTawil2D tag" << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &xBal) != TCL_OK)
		{
			cerr << "WARNING invalid xBal\n";
			cerr << "yieldSurfaceBC ElTawil2D tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &yBal) != TCL_OK)
		{
			cerr << "WARNING invalid yBal\n";
			cerr << "yieldSurfaceBC ElTawil2D tag: " << tag << endl;
			return TCL_ERROR;
		}
		
		if (Tcl_GetDouble(interp, argv[5], &yPos) != TCL_OK)
		{
			cerr << "WARNING invalid xPos\n";
			cerr << "yieldSurfaceBC ElTawil2D tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[6], &yNeg) != TCL_OK)
		{
			cerr << "WARNING invalid yNeg\n";
			cerr << "yieldSurfaceBC ElTawil2D tag: " << tag << endl;
			return TCL_ERROR;
		}
		

		if (Tcl_GetInt(interp, argv[7], &modelID) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC ElTawil2D matID1" << modelID << endl;
			return TCL_ERROR;
		}

		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			cerr << "WARNING yieldSurfaceBC ElTawil2D no ys_model exixts with tag: " << modelID << endl;
			return TCL_ERROR;
		}

		// Parsing was successful, allocate the material
		theYS = new ElTawil2D(tag, xBal, yBal, yPos, yNeg, *theModel);
    }
	
	    else if (strcmp(argv[1],"ElTawil2DUnSym") == 0)
	{
		if (argc < 9)
		{
			cerr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			// Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
			cerr << "Want: yieldSurfaceBC ElTawil2DUnSym tag? xPosBal? yPosBal? "
			     << "xNegBal? yPos? yNeg? ys_model_tag?" << endl;
			return TCL_ERROR;
		}

		int tag;
		double xPosBal, yPosBal;
		double xNegBal, yNegBal;
		double yPos, yNeg;
		
		int modelID;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC ElTawil2DUnSym tag" << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &xPosBal) != TCL_OK)
		{
			cerr << "WARNING invalid xPosBal\n";
			cerr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &yPosBal) != TCL_OK)
		{
			cerr << "WARNING invalid yPosBal\n";
			cerr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endl;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[5], &xNegBal) != TCL_OK)
		{
			cerr << "WARNING invalid xNegBal\n";
			cerr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[6], &yNegBal) != TCL_OK)
		{
			cerr << "WARNING invalid yNegBal\n";
			cerr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endl;
			return TCL_ERROR;
		}
		
		if (Tcl_GetDouble(interp, argv[7], &yPos) != TCL_OK)
		{
			cerr << "WARNING invalid xPos\n";
			cerr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[8], &yNeg) != TCL_OK)
		{
			cerr << "WARNING invalid yNeg\n";
			cerr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endl;
			return TCL_ERROR;
		}
		

		if (Tcl_GetInt(interp, argv[9], &modelID) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC ElTawil2DUnSym matID1" 
			     << modelID << endl;
			return TCL_ERROR;
		}

		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			cerr << "WARNING yieldSurfaceBC ElTawil2D no ys_model exixts with tag: " << modelID << endl;
			return TCL_ERROR;
		}

		// Parsing was successful, allocate the material
		theYS = new ElTawil2DUnSym(tag, xPosBal, yPosBal, xNegBal, yNegBal, yPos, yNeg, *theModel);
    }

	else if (strcmp(argv[1],"Attalla2D") == 0)
	{
// 	     Attalla2D( int tag, double xmax, double ymax, YS_HardeningModel &model,
// 					double x_offset=0, double y_offset=0,
// 					double a01=0.19,  double a02=0.54, double a03=-1.4,
// 					double a04=-1.64, double a05=2.21, double a06=2.10);

		if (argc < 6 || argc > 14)
		{
			cerr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			cerr << "Want: yieldSurfaceBC Attalla2D tag? xCap? yCap? matXTag? maxYTag? isoRatio? <..>" << endl;
			return TCL_ERROR;
		}

		int tag;
		double xCap, yCap;
		//int matID1, matID2;
		int modelID;
		//double isoRatio;
		//double x_offset = 0, y_offset = 0;
		Vector param(6);

		param[0] = 0.19;
		param[1] = 0.54;
		param[2] =-1.40;
		param[3] =-1.64;
		param[4] = 2.21;
		param[5] = 2.10;


		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC Attalla2D tag" << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &xCap) != TCL_OK)
		{
			cerr << "WARNING invalid xCap\n";
			cerr << "yieldSurfaceBC Attalla2D tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &yCap) != TCL_OK)
		{
			cerr << "WARNING invalid yCap\n";
			cerr << "yieldSurfaceBC Attalla2D tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[5], &modelID) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC Attalla2D modelID" << modelID << endl;
			return TCL_ERROR;
		}
		
		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			cerr << "WARNING yieldSurfaceBC Orbison2D no ys_model exixts with tag: " << modelID << endl;
			return TCL_ERROR;
		}



		if(argc > 6)
		{
			int count = 6;
			double temp;

			for(int i=0; i < 6; i++)
			{
				if (Tcl_GetDouble(interp, argv[count], &temp) != TCL_OK)
				{
					cerr << "WARNING invalid parameter " << i+1 << "\n";
					cerr << "yieldSurfaceBC Attalla2D tag: " << tag << endl;
					return TCL_ERROR;
				}
				param(i) = temp;
				count++;
			}

		}

		// Parsing was successful, allocate the material
		theYS = new Attalla2D(tag, xCap, yCap, *theModel,
							  param(0),param(1),param(2),param(3),param(4),param(5));

	}

	else if (strcmp(argv[1],"Hajjar2D") == 0)
	{
		if (argc < 9)
		{
			cerr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			cerr << "Want: yieldSurfaceBC Hajjar2D tag? ysModelTag? D? b? t? fc? fy?" << endl;
			return TCL_ERROR;
		}

		int tag;
		//int matID1, matID2;
		int modelID;
//		double isoRatio;
		double D, b, t, fc, fy;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC Hajjar2D  tag" << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3], &modelID) != TCL_OK)
		{
			cerr << "WARNING invalid yieldSurfaceBC Hajjar2D  matID1" << modelID << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &D) != TCL_OK)
		{
			cerr << "WARNING invalid D \n";
			cerr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK)
		{
			cerr << "WARNING invalid b \n";
			cerr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[6], &t) != TCL_OK)
		{
			cerr << "WARNING invalid t \n";
			cerr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[7], &fc) != TCL_OK)
		{
			cerr << "WARNING invalid fc \n";
			cerr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[8], &fy) != TCL_OK)
		{
			cerr << "WARNING invalid fy \n";
			cerr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endl;
			return TCL_ERROR;
		}
		
		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			cerr << "WARNING yieldSurfaceBC Orbison2D no ys_model exixts with tag: " << modelID << endl;
			return TCL_ERROR;
		}

		// Parsing was successful, allocate the material
		theYS = new Hajjar2D(tag, *theModel, D, b, t, fc, fy);
    }

	else
	{
		cerr << "Warning - unknown yield surface type \n";
		printCommand(argc,argv);
	}

	///////////////////////////////////////////////////////////////
	// Now add the ys to the modelBuilder
	///////////////////////////////////////////////////////////////

	if (theBuilder->addYieldSurface_BC(*theYS) < 0)
	{
		cerr << "WARNING could not add YieldSurfaceBC to the domain\n";
		cerr << *theYS << endl;
		delete theYS; // invoke the material objects destructor, otherwise mem leak
		return TCL_ERROR;
	}


    return TCL_OK;
}







