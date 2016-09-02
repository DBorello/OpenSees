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

// Written: Minjie

// Description: A python wrapper for OpenSees commands
//

#include "PythonWrapper.h"
#include "OpenSeesCommands.h"
#include <OPS_Globals.h>

PythonWrapper* wrapper = 0;

PythonWrapper::PythonWrapper()
    :currentArgv(0), currentArg(0), numberArgs(0),
     methodsOpenSees(), opensees_docstring(""), currentResult(0)
{
    wrapper = this;
}

PythonWrapper::~PythonWrapper()
{
}

void
PythonWrapper::resetCommandLine(int nArgs, int cArg, PyObject* argv)
{
    numberArgs = nArgs;
    currentArg = cArg-1;
    if (currentArg < 0) currentArg = 0;
    currentArgv = argv;
}

void
PythonWrapper::resetCommandLine(int cArg)
{
    if (cArg < 0) {
	currentArg += cArg;
    } else {
	currentArg = cArg-1;
    }
    if (currentArg < 0) currentArg = 0;
}

void
PythonWrapper::addCommand(const char* name, PyCFunction proc)
{
    PyMethodDef method = {name,proc,METH_VARARGS,opensees_docstring};
    methodsOpenSees.push_back(method);
}

PyMethodDef*
PythonWrapper::getMethods()
{
    if (methodsOpenSees.empty()) {
	return 0;
    }
    
    return &methodsOpenSees[0];
}

void
PythonWrapper::setOutputs(int* data, int numArgs)
{
    if (numArgs == 0) return;
    if (numArgs == 1) {
	currentResult = Py_BuildValue("i", data[0]);
	return ;
    }
    currentResult = PyList_New(numArgs);
    for (int i=0; i<numArgs; i++) {
	PyList_SET_ITEM(currentResult, i, Py_BuildValue("i", data[i]));
    }
}

void
PythonWrapper::setOutputs(double* data, int numArgs)
{
    if (numArgs == 0) return;
    if (numArgs == 1) {
	currentResult = Py_BuildValue("d", data[0]);
	return ;
    }
    currentResult = PyList_New(numArgs);
    for (int i=0; i<numArgs; i++) {
	PyList_SET_ITEM(currentResult, i, Py_BuildValue("d", data[i]));
    }
}

void
PythonWrapper::setOutputs(const char* str)
{
    currentResult = Py_BuildValue("s", str);
}

PyObject*
PythonWrapper::getResults()
{
    PyObject* result = currentResult;
    currentResult = 0;

    if (result == 0) {
	Py_INCREF(Py_None);
	result = Py_None;
    }

    return result;
}


//////////////////////////////////////////////
/////// Python wrapper functions  ////////////
/////////////////////////////////////////////
static PyObject *Py_ops_UniaxialMaterial(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_UniaxialMaterial() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_testUniaxialMaterial(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_testUniaxialMaterial() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_setStrain(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setStrain() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getStrain(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getStrain() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getStress(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getStress() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getTangent(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getTangent() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getDampTangent(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getDampTangent() < 0) return NULL;

    return wrapper->getResults();
}


static PyObject *Py_ops_wipe(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_wipe() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_model(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_model() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_node(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Node() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_fix(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_HomogeneousBC() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_element(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Element() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_timeSeries(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_TimeSeries() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_pattern(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Pattern() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodalLoad(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_NodalLoad() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_system(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_System() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_numberer(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Numberer() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_constraints(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_ConstraintHandler() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_integrator(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Integrator() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_algorithm(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Algorithm() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_analysis(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Analysis() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_analyze(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_analyze() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeDisp(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeDisp() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_test(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_CTest() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_section(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Section() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_fiber(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Fiber() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_patch(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Patch() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_layer(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Layer() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_geomTransf(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_CrdTransf() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_beamIntegration(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_BeamIntegration() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_loadConst(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_loadConst() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_eleLoad(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_ElementalLoad() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_reactions(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_calculateNodalReactions() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeReaction(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeReaction() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_eigen(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eigenAnalysis() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nDMaterial(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_NDMaterial() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_block2d(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_doBlock2D() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_block3d(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_doBlock3D() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_rayleigh(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_rayleighDamping() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_wipeAnalysis(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_wipeAnalysis() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_setTime(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setTime() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_remove(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_removeObject() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_mass(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_addNodalMass() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_equalDOF(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_EqualDOF() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeEigenvector(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeEigenvector() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getTime(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getTime() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_eleResponse(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eleResponse() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_SP(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_SP() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_fixX(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_HomogeneousBC_X() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_fixY(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_HomogeneousBC_Y() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_fixZ(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_HomogeneousBC_Z() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_reset(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_resetModel() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_initialize(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_initializeAnalysis() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getLoadFactor(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getLoadFactor() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_build(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_buildModel() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_print(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_printModel() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_printA(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_printA() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_printB(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_printB() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_printGID(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_printModelGID() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getCTestNorms(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getCTestNorms() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getCTestIter(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getCTestIter() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_recorder(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Recorder() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_database(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Database() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_save(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_save() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_restore(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_restore() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_eleForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eleForce() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_eleDynamicalForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eleDynamicalForce() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeUnbalance(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeUnbalance() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeVel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeVel() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_setNodeVel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setNodeVel() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeAccel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeAccel() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeResponse(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeResponse() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeCoord(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeCoord() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_setNodeCoord(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setNodeCoord() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_updateElementDomain(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_updateElementDomain() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_eleNodes(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eleNodes() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeMass(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeMass() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodePressure(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodePressure() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeBounds(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeBounds() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_startTimer(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_startTimer() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_stopTimer(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_stopTimer() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_modalDamping(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_modalDamping() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_modalDampingQ(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_modalDampingQ() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_setElementRayleighDampingFactors(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setElementRayleighDampingFactors() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_region(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_MeshRegion() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_setPrecision(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setPrecision() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_searchPeerNGA(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_peerNGA() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_domainChange(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_domainChange() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_metaData(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_neesMetaData() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_neesUpload(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_neesUpload() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_stripXML(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_stripOpenSeesXML() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_convertTextToBinary(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_convertTextToBinary() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_convertBinaryToText(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_convertBinaryToText() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getEleTags(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getEleTags() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getNodeTags(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getNodeTags() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getParamTags(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getParamTags() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getParamValue(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getParamValue() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionForce() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionDeformation(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionDeformation() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionStiffness(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionStiffness() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionFlexibility(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionFlexibility() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionLocation(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionLocation() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionWeight(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionWeight() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_basicDeformation(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_basicDeformation() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_basicForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_basicForce() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_basicStiffness(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_basicStiffness() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_InitialStateAnalysis(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_InitialStateAnalysis() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_totalCPU(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_totalCPU() < 0) return NULL;

    return wrapper->getResults();
}


static PyObject *Py_ops_solveCPU(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_solveCPU() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_accelCPU(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_accelCPU() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_numFact(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_numFact() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_numIter(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_numIter() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_systemSize(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_systemSize() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_version(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_version() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_setMaxOpenFiles(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_maxOpenFiles() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_background(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_BackgroundMesh() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_limitCurve(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_LimitCurve() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_imposedMotion(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_ImposedMotionSP() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_groundMotion(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_groundMotion() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_equalDOF_Mixed(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_EqualDOF_Mixed() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_rigidLink(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_RigidLink() < 0) return NULL;

    return wrapper->getResults();
}


static PyObject *Py_ops_rigidDiaphragm(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_RigidDiaphragm() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_ShallowFoundationGen(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_ShallowFoundationGen() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_setElementRayleighFactors(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_addElementRayleigh() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_mesh(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_mesh() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_remesh(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_remesh() < 0) return NULL;

    return wrapper->getResults();
}

/////////////////////////////////////////////////
////////////// Add Python commands //////////////
/////////////////////////////////////////////////
void
PythonWrapper::addOpenSeesCommands()
{
    addCommand("uniaxialMaterial", &Py_ops_UniaxialMaterial);
    addCommand("testUniaxialMaterial", &Py_ops_testUniaxialMaterial);
    addCommand("setStrain", &Py_ops_setStrain);
    addCommand("getStrain", &Py_ops_getStrain);
    addCommand("getStress", &Py_ops_getStress);
    addCommand("getTangent", &Py_ops_getTangent);
    addCommand("getDampTangent", &Py_ops_getDampTangent);
    addCommand("wipe", &Py_ops_wipe);
    addCommand("model", &Py_ops_model);
    addCommand("node", &Py_ops_node);
    addCommand("fix", &Py_ops_fix);
    addCommand("element", &Py_ops_element);
    addCommand("timeSeries", &Py_ops_timeSeries);
    addCommand("pattern", &Py_ops_pattern);
    addCommand("load", &Py_ops_nodalLoad);
    addCommand("system", &Py_ops_system);
    addCommand("numberer", &Py_ops_numberer);
    addCommand("constraints", &Py_ops_constraints);
    addCommand("integrator", &Py_ops_integrator);
    addCommand("algorithm", &Py_ops_algorithm);
    addCommand("analysis", &Py_ops_analysis);
    addCommand("analyze", &Py_ops_analyze);
    addCommand("nodeDisp", &Py_ops_nodeDisp);
    addCommand("test", &Py_ops_test);
    addCommand("section", &Py_ops_section);
    addCommand("fiber", &Py_ops_fiber);
    addCommand("patch", &Py_ops_patch);
    addCommand("layer", &Py_ops_layer);
    addCommand("geomTransf", &Py_ops_geomTransf);
    addCommand("beamIntegration", &Py_ops_beamIntegration);
    addCommand("loadConst", &Py_ops_loadConst);
    addCommand("eleLoad", &Py_ops_eleLoad);
    addCommand("reactions", &Py_ops_reactions);
    addCommand("nodeReaction", &Py_ops_nodeReaction);
    addCommand("eigen", &Py_ops_eigen);
    addCommand("nDMaterial", &Py_ops_nDMaterial);
    addCommand("block2D", &Py_ops_block2d);
    addCommand("block3D", &Py_ops_block3d);
    addCommand("rayleigh", &Py_ops_rayleigh);
    addCommand("wipeAnalysis", &Py_ops_wipeAnalysis);
    addCommand("setTime", &Py_ops_setTime);
    addCommand("remove", &Py_ops_remove);
    addCommand("mass", &Py_ops_mass);
    addCommand("equalDOF", &Py_ops_equalDOF);
    addCommand("nodeEigenvector", &Py_ops_nodeEigenvector);
    addCommand("getTime", &Py_ops_getTime);
    addCommand("eleResponse", &Py_ops_eleResponse);
    addCommand("sp", &Py_ops_SP);
    addCommand("fixX", &Py_ops_fixX);
    addCommand("fixY", &Py_ops_fixY);
    addCommand("fixZ", &Py_ops_fixZ);
    addCommand("reset", &Py_ops_reset);
    addCommand("initialize", &Py_ops_initialize);
    addCommand("getLoadFactor", &Py_ops_getLoadFactor);
    addCommand("build", &Py_ops_build);
    addCommand("Print", &Py_ops_print);
    addCommand("printA", &Py_ops_printA);
    addCommand("printB", &Py_ops_printB);
    addCommand("printGID", &Py_ops_printGID);
    addCommand("getCTestNorms", &Py_ops_getCTestNorms);
    addCommand("getCTestIter", &Py_ops_getCTestIter);
    addCommand("recorder", &Py_ops_recorder);
    addCommand("database", &Py_ops_database);
    addCommand("save", &Py_ops_save);
    addCommand("restore", &Py_ops_restore);
    addCommand("eleForce", &Py_ops_eleForce);
    addCommand("eleDynamicalForce", &Py_ops_eleDynamicalForce);
    addCommand("nodeUnbalance", &Py_ops_nodeUnbalance);
    addCommand("nodeVel", &Py_ops_nodeVel);
    addCommand("setNodeVel", &Py_ops_setNodeVel);
    addCommand("nodeAccel", &Py_ops_nodeAccel);
    addCommand("nodeResponse", &Py_ops_nodeResponse);
    addCommand("nodeCoord", &Py_ops_nodeCoord);
    addCommand("setNodeCoord", &Py_ops_setNodeCoord);
    addCommand("updateElementDomain", &Py_ops_updateElementDomain);
    addCommand("eleNodes", &Py_ops_eleNodes);
    addCommand("nodeMass", &Py_ops_nodeMass);
    addCommand("nodePressure", &Py_ops_nodePressure);
    addCommand("nodeBounds", &Py_ops_nodeBounds);
    addCommand("start", &Py_ops_startTimer);
    addCommand("stop", &Py_ops_stopTimer);
    addCommand("modalDamping", &Py_ops_modalDamping);
    addCommand("modalDampingQ", &Py_ops_modalDampingQ);
    addCommand("setElementRayleighDampingFactors", &Py_ops_setElementRayleighDampingFactors);
    addCommand("region", &Py_ops_region);
    addCommand("setPrecision", &Py_ops_setPrecision);
    addCommand("searchPeerNGA", &Py_ops_searchPeerNGA);
    addCommand("domainChange", &Py_ops_domainChange);
    addCommand("metaData", &Py_ops_metaData);
    addCommand("neesUpload", &Py_ops_neesUpload);
    addCommand("stripXML", &Py_ops_stripXML);
    addCommand("convertBinaryToText", &Py_ops_convertBinaryToText);
    addCommand("convertTextToBinary", &Py_ops_convertTextToBinary);
    addCommand("getEleTags", &Py_ops_getEleTags);
    addCommand("getNodeTags", &Py_ops_getNodeTags);
    addCommand("getParamTags", &Py_ops_getParamTags);
    addCommand("getParamValue", &Py_ops_getParamValue);
    addCommand("sectionForce", &Py_ops_sectionForce);
    addCommand("sectionDeformation", &Py_ops_sectionDeformation);
    addCommand("sectionStiffness", &Py_ops_sectionStiffness);
    addCommand("sectionFlexibility", &Py_ops_sectionFlexibility);
    addCommand("sectionLocation", &Py_ops_sectionLocation);
    addCommand("sectionWeight", &Py_ops_sectionWeight);
    addCommand("basicDeformation", &Py_ops_basicDeformation);
    addCommand("basicForce", &Py_ops_basicForce);
    addCommand("basicStiffness", &Py_ops_basicStiffness);
    addCommand("InitialStateAnalysis", &Py_ops_InitialStateAnalysis);
    addCommand("totalCPU", &Py_ops_totalCPU);
    addCommand("solveCPU", &Py_ops_solveCPU);
    addCommand("accelCPU", &Py_ops_accelCPU);
    addCommand("numFact", &Py_ops_numFact);
    addCommand("numIter", &Py_ops_numIter);
    addCommand("systemSize", &Py_ops_systemSize);
    addCommand("version", &Py_ops_version);
    addCommand("setMaxOpenFiles", &Py_ops_setMaxOpenFiles);
    addCommand("background", &Py_ops_background);
    addCommand("limitCurve", &Py_ops_limitCurve);
    addCommand("imposedMotion", &Py_ops_imposedMotion);
    addCommand("imposedSupportMotion", &Py_ops_imposedMotion);
    addCommand("groundMotion", &Py_ops_groundMotion);
    addCommand("equalDOF_Mixed", &Py_ops_equalDOF_Mixed);
    addCommand("rigidLink", &Py_ops_rigidLink);
    addCommand("rigidDiaphragm", &Py_ops_rigidDiaphragm);
    addCommand("ShallowFoundationGen", &Py_ops_ShallowFoundationGen);
    addCommand("setElementRayleighFactors", &Py_ops_setElementRayleighFactors);
    addCommand("mesh", &Py_ops_mesh);
    addCommand("remesh", &Py_ops_remesh);
    
    PyMethodDef method = {NULL,NULL,0,NULL};
    methodsOpenSees.push_back(method);

    
}
