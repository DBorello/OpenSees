//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):                 
//                 ``This    source  code is Copyrighted in 
//                 U.S.,   by   the   The  Regents  of  the 
//                 University   of   California,   for   an 
//                 indefinite  period,  and  anybody caught 
//                 using it without our permission, will be 
//                 mighty  good  friends  of ourn, cause we 
//                 don't  give a darn. Hack it. Compile it. 
//                 Debug it. Run it. Yodel it. Enjoy it. We 
//                 wrote  it, that's all we wanted to do.'' 
//     
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Plastic Bowl (aka Domain Reduction) implementation:
//#                    This file contains the class definition 
//#                    for PBowlLoading. 
//#                    PBowlLoading is an subclass of loadPattern.
//# CLASS:             PBowlLoading
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic
//# PROGRAMMER(S):     Jinxou Liao, Zhaohui Yang, Boris Jeremic
//#
//#
//# DATE:              21Oct2002
//# UPDATE HISTORY:
//#
//#
//#
//#
//#
//===============================================================================

// Purpose: This file contains the class definition for PBowlLoading.
// PBowlLoading is an subclass of loadPattern.

#include <PBowlLoading.h>

#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <stdlib.h>
#include <Channel.h>
#include <ErrorHandler.h>

#include <fstream.h>
#include <iostream.h>


PBowlLoading::PBowlLoading()
:LoadPattern(0, PATTERN_TAG_PBowlLoading),
PBowlElements(0),
PBowlNodes(0),
PBowlLoads(0),
U(0),
Udd(0)
{

}

PBowlLoading::PBowlLoading(int tag,
                           const ID &_PBElements,
                           char *DispfName,
                           char *AccefName,
                           double theTimeIncr,
                           double theFactor)
:LoadPattern(tag, 
             PATTERN_TAG_PBowlLoading),
             PBTimeIncr(theTimeIncr), 
             cFactor(theFactor)
  {
    // determine the number of data points .. open file and count num entries
    int timeSteps1, timeSteps2;
    int numDataPoints =0;
    double dataPoint;


    ifstream theFile;

  //--------------------------------
  //Input displacement data
  //--------------------------------
    theFile.open(DispfName);
    if (theFile.bad()) 
      {
        cerr << "WARNING - PBowlLoading::PBowlLoading()";
        cerr << " - could not open file " << DispfName << endl;
      } 
    else 
      {
        //Input the number of time steps first
        theFile >> timeSteps1;
        //Loop to count the number of data points
        while (theFile >> dataPoint)
          numDataPoints++;
      }
    theFile.close();
    UnumDataPoints = numDataPoints;
    thetimeSteps = timeSteps1;
    if ( timeSteps1 == 0) 
      {
        cerr << "WARNING - PBowlLoading::PBowlLoading()";
        cerr << " - Time steps for displacement equal to zero... " << AccefName << endl;
        exit(0);
      }

  // create a vector and read in the data
    if (numDataPoints != 0) {

    // first open the file
    theFile.open(DispfName, ios::in);
    if (theFile.bad()) {
      cerr << "WARNING - PBowlLoading::PBowlLoading()";
      cerr << " - could not open file " << DispfName << endl;
    } else {

      // now create the vector
      if ( numDataPoints - thetimeSteps*(numDataPoints/thetimeSteps) != 0 ) {
         cerr << "WARNING - PBowlLoading::PBowlLoading()";
         cerr << " - Input data not sufficient...Patched with zeros " << DispfName << endl;
      }
      int cols = numDataPoints/thetimeSteps;
      U = new Matrix(cols, thetimeSteps);

      // ensure we did not run out of memory
      if (U == 0 || U->noRows() == 0 || U->noCols() == 0) {
  cerr << "PBowlLoading::PBowlLoading() - ran out of memory constructing";
  cerr << " a Matrix of size (cols*rows): " << cols << " * " << thetimeSteps << endl;

  if (U != 0)
    delete U;
  U = 0;
      }

      // read the data from the file
      else {
  theFile >> timeSteps1;
  for (int t=0; t< thetimeSteps; t++)
           for  (int j=0;j<cols; j++) {
        theFile >> dataPoint;
        (*U)(j, t) = dataPoint;
           }
  }
      }

      // finally close the file
      theFile.close();
  }

  //--------------------------------
  //Input acceleration data
  //--------------------------------
  theFile.open(AccefName);
  numDataPoints = 0;

  if (theFile.bad()) {
    cerr << "WARNING - PBowlLoading::PBowlLoading()";
    cerr << " - could not open file " << AccefName << endl;
  } else {
    //Input the number of time steps first
    theFile >> timeSteps2;
    //Loop to count the number of data points
    while (theFile >> dataPoint)
      numDataPoints++;
  }
  theFile.close();
  UddnumDataPoints = numDataPoints;
  if ( timeSteps1 !=  timeSteps2) {
    cerr << "WARNING - PBowlLoading::PBowlLoading()";
    cerr << " - Time steps for displacement not equal to that for acceleration... " << AccefName << endl;
  }


  // create a vector and read in the data
  if (numDataPoints != 0) {

    // first open the file
    theFile.open(AccefName, ios::in);
    if (theFile.bad()) {
      cerr << "WARNING - PBowlLoading::PBowlLoading()";
      cerr << " - could not open file " << AccefName << endl;
    } else {

      // now create the vector
      if ( numDataPoints - thetimeSteps*(numDataPoints/thetimeSteps) != 0 ) {
         cerr << "WARNING - PBowlLoading::PBowlLoading()";
         cerr << " - Input data not sufficient...Patched with zeros " << AccefName << endl;
      }
      int cols = numDataPoints/thetimeSteps;
      Udd = new Matrix(cols, thetimeSteps);

      // ensure we did not run out of memory
      if (Udd == 0 || Udd->noRows() == 0 || Udd->noCols() == 0) {
  cerr << "PBowlLoading::PBowlLoading() - ran out of memory constructing";
  cerr << " a Matrix of size (cols*rows): " << cols << " * " << thetimeSteps << endl;

  if (Udd != 0)
    delete Udd;
  Udd = 0;
      }

      // read the data from the file
      else {
  theFile >> timeSteps2;
  for (int t=0; t< thetimeSteps; t++)
           for  (int j=0;j<cols; j++) {
        theFile >> dataPoint;
        (*Udd)(j, t) = dataPoint;
           }
  }
      }

      // finally close the file
      theFile.close();
  }

  //Adding the plastic bowl node ID
  this->addPBElements(_PBElements);
  LoadComputed = false;
}



PBowlLoading::~PBowlLoading()
{
  // invoke the destructor on all ground motions supplied

  if ( PBowlNodes != 0)
    delete PBowlNodes;

}


void
PBowlLoading::setDomain(Domain *theDomain)
{
  this->LoadPattern::setDomain(theDomain);

  // // now we go through and set all the node velocities to be vel0
  // if (vel0 != 0.0) {
  //   NodeIter &theNodes = theDomain->getNodes();
  //   Node *theNode;
  //   Vector newVel(1);
  //   int currentSize = 1;
  //   while ((theNode = theNodes()) != 0) {
  //     int numDOF = theNode->getNumberDOF();
  //     if (numDOF != currentSize)
  //   newVel.resize(numDOF);
  //
  //     newVel = theNode->getVel();
  //     newVel(theDof) = vel0;
  //     theNode->setTrialVel(newVel);
  //     theNode->commitState();
  //   }
  // }
}


void
PBowlLoading::applyLoad(double time)
{

  Domain *theDomain = this->getDomain();
  if (theDomain == 0)
    return;

  //Finding the all the nodes in the plastic bowl and sort it ascendingly
  if ( !LoadComputed )
     this->CompPBLoads();

  // see if quick return, i.e. no plastic bowl nodes or domain set
  int numPBnodes = PBowlNodes->Size();
  if (numPBnodes == 0)
    return;

  //Apply loads on each plastic bowl nodes
  Node *theNode;
  for (int i=0; i<numPBnodes; i++) {
    const Vector &load=this->getNodalLoad((*PBowlNodes)[i], time);
    theNode = theDomain->getNode( (*PBowlNodes)[i] );
    theNode->addUnbalancedLoad(load);
  }

}

int
PBowlLoading::sendSelf(int commitTag, Channel &theChannel)
{
  cerr << "PBowlLoading::sendSelf() - not yet implemented\n";
  return 0;
}

int
PBowlLoading::recvSelf(int commitTag, Channel &theChannel,
       FEM_ObjectBroker &theBroker)
{
  cerr << "PBowlLoading::recvSelf() - not yet implemented\n";
  return 0;
}

/* **********************************************************************************************
int
PBowlLoading::sendSelf(int commitTag, Channel &theChannel)
{
  // first send the tag and info about the number of ground motions
  int myDbTag = this->getDbTag();
  ID theData(2);
  theData(0) = this->getTag();
  theData(1) = numMotions;

  if (theChannel.sendID(myDbTag, commitTag, theData) < 0) {
    g3ErrorHandler->warning("PBowlLoading::sendSelf - channel failed to send the initial ID");
    return -1;
  }

  // now for each motion we send it's classsss tag and dbtag
  ID theMotionsData(2*numMotions);
  for (int i=0; i<numMotions; i++) {
    theMotionsData[i] = theMotions[i]->getClassTag();
    int motionsDbTag = theMotions[i]->getDbTag();
    if (motionsDbTag == 0) {
      motionsDbTag = theChannel.getDbTag();
      theMotions[i]->setDbTag(motionsDbTag);
    }
    theMotionsData[i+numMotions] = motionsDbTag;
  }

  if (theChannel.sendID(myDbTag, commitTag, theMotionsData) < 0) {
    g3ErrorHandler->warning("PBowlLoading::sendSelf - channel failed to send the motions ID");
    return -1;
  }

  // now we send each motion
  for (int j=0; j<numMotions; j++)
    if (theMotions[j]->sendSelf(commitTag, theChannel) < 0) {
      g3ErrorHandler->warning("PBowlLoading::sendSelf - motion no: %d failed in sendSelf", j);
      return -1;
    }

  // if get here successfull
  return 0;
}

int
PBowlLoading::recvSelf(int commitTag, Channel &theChannel,
       FEM_ObjectBroker &theBroker)
{
  // first get the tag and info about the number of ground motions from the Channel
  int myDbTag = this->getDbTag();
  ID theData(2);
  if (theChannel.recvID(myDbTag, commitTag, theData) < 0) {
    g3ErrorHandler->warning("PBowlLoading::recvSelf - channel failed to recv the initial ID");
    return -1;
  }

  // set current tag
  this->setTag(theData(0));

  // now get info about each channel
  ID theMotionsData (2*theData(1));
  if (theChannel.recvID(myDbTag, commitTag, theMotionsData) < 0) {
    g3ErrorHandler->warning("PBowlLoading::recvSelf - channel failed to recv the motions ID");
    return -1;
  }


  if (numMotions != theData(1)) {

    //
    // we must delete the old motions and create new ones and then invoke recvSelf on these new ones
    //

    if (numMotions != 0) {
      for (int i=0; i<numMotions; i++)
  delete theMotions[i];
      delete [] theMotions;
    }
    numMotions = theData[1];
    theMotions = new (GroundMotion *)[numMotions];
    if (theMotions == 0) {
      g3ErrorHandler->warning("PBowlLoading::recvSelf - out of memory creating motion array of size %d\n", numMotions);
      numMotions = 0;
      return -1;
    }

    for (int i=0; i<numMotions; i++) {
      theMotions[i] = theBroker.getNewGroundMotion(theMotionsData[i]);
      if (theMotions[i] == 0) {
  g3ErrorHandler->warning("PBowlLoading::recvSelf - out of memory creating motion array of size %d\n", numMotions);
  numMotions = 0;
  return -1;
      }
      theMotions[i]->setDbTag(theMotionsData[i+numMotions]);
      if (theMotions[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
  g3ErrorHandler->warning("PBowlLoading::recvSelf - motion no: %d failed in recvSelf", i);
  numMotions = 0;
  return -1;
      }
    }

  } else {

    //
    // we invoke rrecvSelf on the motions, note if a motion not of correct class
    // we must invoke the destructor on the motion and create a new one of correct type
    //

    for (int i=0; i<numMotions; i++) {
      if (theMotions[i]->getClassTag() == theMotionsData[i]) {
  if (theMotions[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
    g3ErrorHandler->warning("PBowlLoading::recvSelf - motion no: %d failed in recvSelf", i);
    return -1;
  }
      } else {
  // motion not of correct type
  delete theMotions[i];
  theMotions[i] = theBroker.getNewGroundMotion(theMotionsData[i]);
  if (theMotions[i] == 0) {
    g3ErrorHandler->warning("PBowlLoading::recvSelf - out of memory creating motion array of size %d\n", numMotions);
    numMotions = 0;
    return -1;
  }
  theMotions[i]->setDbTag(theMotionsData[i+numMotions]);
  if (theMotions[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
    g3ErrorHandler->warning("PBowlLoading::recvSelf - motion no: %d failed in recvSelf", i);
    numMotions = 0;
    return -1;
  }
      }
    }
  }

  // if get here successfull
  return 0;
}


***************************************************************************************** */

void
PBowlLoading::Print(ostream &s, int flag)
{
  cerr << "PBowlLoading::Print() - not yet implemented\n";
}



// method to obtain a blank copy of the LoadPattern
LoadPattern *
PBowlLoading::getCopy(void)
{

  cerr << "PBowlLoading::getCopy() - not yet implemented\n";
  return 0;
}

/* ****************************************************************************************
// method to obtain a blank copy of the LoadPattern
LoadPattern *
PBowlLoading::getCopy(void)
{
  PBowlLoading *theCopy = new PBowlLoading(0, 0);
  theCopy->setTag(this->getTag);
  theCopy->numMotions = numMotions;
  theCopy->theMotions = new (GroundMotion *)[numMotions];
  for (int i=0; i<numMotions; i++)
    theCopy->theMotions[i] = theMotions[i];

  return 0;
}
***************************************************************************************** */

void
PBowlLoading::addPBElements(const ID &PBEle)
{
  // create a copy of the vector containg plastic bowl elements
  PBowlElements = new ID(PBEle);
  // ensure we did not run out of memory
  if (PBowlElements == 0 || PBowlElements->Size() == 0) {
    cerr << "PBowlLoading::addPBElements() - ran out of memory constructing";
    cerr << " a Vector of size: " <<  PBowlElements->Size() << endl;
    if (PBowlElements != 0)
      delete PBowlElements;
    PBowlElements = 0;
  }

}


void
PBowlLoading::CompPBLoads()
{
  //===========================================================
  // Finding all plastic bowl nodes
  //===========================================================
  Domain *theDomain = this->getDomain();
  Element *theElement = theDomain->getElement( (*PBowlElements)(0) );
  int NIE = theElement->getNumExternalNodes();

  int max_bnode = PBowlElements->Size() * NIE;
  //int *Bowl_node;
  ID *Bowl_node = new ID(max_bnode);
  ID Resort = theElement->getExternalNodes();

  int i, j, k;
  //Inital node list from the first plastic bowl element
  for (i = 0; i<NIE; i++)
     (*Bowl_node)(i) = Resort(i);

  //cout << Resort;

  //... int min;
  //... int min_nb ;
  //... for (int i=0;i<NIE; i++)
  //...   {
  //...     min = max_bnode;
  //...     for (int j=0;j<NIE; j++)
  //...       {
  //...         if (Resort(j) < min)
  //...           {
  //...             min= Resort(j);
  //...             min_nb =j;
  //...           }
  //...       }
  //...     Bowl_node(i)=min;
  //...     Resort(min_nb) = max_bnode;
  //...   }
  //...
  //... //for (int bi=0;bi<20;bi++)
  //... //        cout<<Bowl_node[bi]<<"  ";
  //... //   cout<<endl<<endl;

  //  insert the node number to the array of Bowl_node

  int no_bnode = NIE;
  int Bowl_elem_nb = PBowlElements->Size();
  ID Temp;

  //... for ( i=1; i<Bowl_elem_nb; i++)
  //... {
  //...     theElement = theDomain->getElement( (*PBowlElements)(i) );
  //...     Temp = theElement->getExternalNodes();
  //...     for ( j=0;j<NIE;j++)
  //...     {
  //...         if ( Temp(j) < Bowl_node(0) )
  //...         {
  //...             for (k= no_bnode;k>0 ; k--)
  //...             {
  //...                 Bowl_node(k+1) = Bowl_node(k);
  //...             }
  //...             Bowl_node(0) = Temp(j);
  //...             no_bnode++;
  //...         }
  //...         else  if ( Temp(j) > Bowl_node(no_bnode-1) )
  //...         {
  //...              Bowl_node(no_bnode) = Temp(j);
  //...              no_bnode++;
  //...         }
  //...         else
  //...         {
  //...             int not_set = 1;
  //...             int k = 0;
  //...             int ss;
  //...             while ((not_set==1)  && ( k<no_bnode))
  //...             {
  //...                 if( (Temp(j) > Bowl_node(k)) && (Temp(j) < Bowl_node(k+1)) )
  //...                 {
  //...                     for ( ss=no_bnode;ss> k; ss--)
  //...                     {
  //...                         Bowl_node(ss+1) = Bowl_node(ss);
  //...                     }
  //...                     Bowl_node(k+1) = Temp(j);
  //...                     not_set = 0;
  //...                     no_bnode++;
  //...                 }
  //...                 else k++;
  //...             }
  //...
  //...         }
  //...     }
  //... }

  //--------------------------------------------------
  //Just make a list of all plastic bowl nodes, no need to sort??? Joey Yang Oct. 18, 2002
  //--------------------------------------------------
  for ( i=1; i<Bowl_elem_nb; i++)
  {
      theElement = theDomain->getElement( (*PBowlElements)(i) );
      Temp = theElement->getExternalNodes();
      for ( j=0;j<NIE;j++)
      {
  for (k = 0; k< no_bnode; k++) {
    if ( Temp(j) == (*Bowl_node)(k) )
      //cout << Temp(j) << "  " << Bowl_node(k) << endl;
      break;
    } //endofIF

    if ( k == no_bnode) {
      (*Bowl_node)(no_bnode) = Temp(j);
      no_bnode ++;
        } //end of for (k=0...)

      }//end of for (j=0...)

  }
  //--Joey------------------------------------------------



  //check the bowl_node
  cout << "\nCheck all plastic bowl nodes...\n";
  for (int bi=0;bi<no_bnode;bi++)
       cout<< (*Bowl_node)(bi) <<"  ";
  cout<< endl << "# of pbowl nodes = " << no_bnode<<endl;
  //cout<<"finish inputting  and organizing the Bowl_node array"<<endl;


  //Adding all plastic bowl nodes
  PBowlNodes = new ID(no_bnode);

  if (PBowlNodes == 0 || PBowlNodes->Size() == 0) {
    cerr << "PBowlLoading::PBowlLoading() - ran out of memory constructing";
    cerr << " a Vector of size: " <<  PBowlNodes->Size() << endl;
    if (PBowlNodes != 0)
      delete PBowlNodes;
    PBowlNodes = 0;
  }
  for (i =0; i < no_bnode; i++)
   (*PBowlNodes)(i) = (*Bowl_node)(i);

  //===========================================================
  // Computing the equivalent force for all plastic bowl nodes
  //===========================================================

  int cols = Udd->noRows();
  Matrix *F = new Matrix(cols, thetimeSteps);

  Node *theNode = theDomain->getNode((*PBowlNodes)(0));
  int NDOF = theNode->getNumberDOF();

  Vector *Fm = new Vector(NIE*NDOF);
  Vector *Fk  = new Vector(NIE*NDOF);
  Matrix *Ke= new Matrix(NIE*NDOF,NIE*NDOF);
  Matrix *Me= new Matrix(NIE*NDOF,NIE*NDOF);
  Vector *u_e = new Vector(NIE*NDOF);
  Vector *udd_e = new Vector(NIE*NDOF);


  // intialize the F()
  for ( i=0;i<cols; i++)
      for ( j=0;j<thetimeSteps; j++)
         (*F)(i,j)=0;

  Element *theBowlElements;

  for ( i=0; i<Bowl_elem_nb; i++)
  {
   // get the Brick;
   theBowlElements = theDomain->getElement( (*PBowlElements)(i) );
   ID *nd = new ID(NIE);
   (*nd) = theBowlElements->getExternalNodes();

   (*Ke) = theBowlElements ->getTangentStiff();
   (*Me) = theBowlElements ->getMass();


   //   get the u and u_dotdot for this element
   for ( int t=0;t<thetimeSteps; t++)
   {
     //cout << "element: " << i << "" << " Time step: " << t << endl;
     for (int j=0;j<NIE;j++)  //BJ make it into # of nodes per element (2D or 3D...)
     {
       for (int d=0;d<NDOF;d++)
       {
         (*u_e)(j*NDOF+d)      = (*U)( (*nd)(j)*NDOF-NDOF+d,t);
         (*udd_e)(j*NDOF+d)    = (*Udd)( (*nd)(j)*NDOF-NDOF+d,t);
          //if ( (*nd)(j) == 109 )
   //  cout << "----nd(J) u udd t= " << t << " dim =" << d << " " << (*nd)(j) << " " <<  (*u_e)(j*NDOF+d) << " " << (*udd_e)(j*NDOF+d) << endl;
   //fflush(stdout);
       }
     }

     (*Fm) = (*Me) * (*udd_e);
     //cout<<"Fm = \n"<<Fm<<endl;
     //cout<<Ke;

     (*Fk) = (*Ke) * (*u_e);
     //cout<<"Fk = \n"<<Fk<<endl;

     for (int k=0;k<NIE; k++)
        for (int d=0;d<NDOF;d++)
            (*F)( (*nd)(k)*NDOF-NDOF+d,t) = (*F)( (*nd)(k)*NDOF-NDOF+d,t) + (*Fm)(k*NDOF+d)   + (*Fk)(k*NDOF+d);

   } //end for timestep

  }  // end for bowl element

  PBowlLoads = new Matrix(*F);
  //cout <<  PBowlLoads->noCols() << " " << PBowlLoads->noRows() << endl;

  // ensure we did not run out of memory
  if (PBowlLoads->noRows() == 0 || PBowlLoads->noCols() == 0 ) {
    cerr << "PBowlLoading::PBowlLoads() - ran out of memory";
    cerr << " a Matrix of size: " <<  PBowlLoads->noRows() << " * " << PBowlLoads->noCols() << endl;
  }

  cout<<"\nFinish calculating the forces..." << endl << endl;
  LoadComputed = true;

  delete Fm;
  delete Fk;
  delete Ke;
  delete Me;
  delete u_e;
  delete udd_e;
  delete F;

}


const Vector &
PBowlLoading::getNodalLoad(int nodeTag, double time)
{
  Vector *dummy = new Vector(0);
  //Get the node
  Domain *theDomain = this->getDomain();
  Node *theNode = theDomain->getNode(nodeTag);
  if (theNode == 0) {
     cerr << "PBowlLoading::getNodalLoad() - no nodes associtated to the nodeTag " << nodeTag << "\n";
     return ( *dummy );
  }

  delete dummy;

  //Create the nodal load vector accoding to the DOFs the node has
  int numDOF = theNode->getNumberDOF();
  Vector *nodalLoad = new Vector(numDOF);


  //Get the nodal loads associated to the nodeTag
  // check for a quick return
  if (time < 0.0 || PBowlLoads == 0)
    return (*nodalLoad);

  // determine indexes into the data array whose boundary holds the time
  double incr = time/PBTimeIncr;
  int incr1 = (int) floor(incr)-1;
  int incr2 = incr1 + 1;
  double value1=0, value2=0;

  int i;
  if ( incr2 == thetimeSteps )
    for (i = 0; i < numDOF; i++)
       (*nodalLoad)(i) = (*PBowlLoads)(i, incr1);
  //If beyond time step, return 0 loads
  else if (incr2 > thetimeSteps ) {
//test    if ( nodeTag == 109)
//test       cout << "Time = " << time << " Node # " << nodeTag  << " " << (*nodalLoad)(0) << " "<< (*nodalLoad)(1) << " "<< (*nodalLoad)(2) << endl;
    return (*nodalLoad);
  }

  //If within time step, return interpolated values
  else {
    for (i = 0; i < numDOF; i++){
       value1 = (*PBowlLoads)(numDOF*(nodeTag-1)+i, incr1);
       value2 = (*PBowlLoads)(numDOF*(nodeTag-1)+i, incr2);
       (*nodalLoad)(i) = cFactor*(value1 + (value2-value1)*(time/PBTimeIncr - incr2));
    }
  }

//test  if ( nodeTag == 109)
//test    cout << "Time = " << time << " Node # " << nodeTag  << " " << (*nodalLoad)(0) << " "<< (*nodalLoad)(1) << " "<< (*nodalLoad)(2) << endl;

  return (*nodalLoad);

}
