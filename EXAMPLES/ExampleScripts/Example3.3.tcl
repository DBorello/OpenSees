# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Portal Frame Example 3.3
# ------------------------
#  Reinforced concrete one-bay, one-story frame
#  Distributed vertical load on girder
#  Uniform excitation acting at fixed nodes in horizontal direction
#  
# 
# Example Objectives
# -----------------
#  Nonlinear dynamic analyis using Portal Frame Example 1 as staring point
#  Using Tcl Procedures 
#
# 
# Units: kips, in, sec
#
# Written: GLF/MHS/fmk
# Date: January 2001


# ----------------------------------------------------
# Start of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------

# Do operations of Example3.1 by sourcing in the tcl file
source Example3.1.tcl
puts "Gravity load analysis completed"

# Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0

# ----------------------------------------------------
# End of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------


# ----------------------------------------------------
# Start of additional modelling for dynamic loads
# ----------------------------------------------------

# Define nodal mass in terms of axial load on columns
set g 386.4
set m [expr $P/$g];       # expr command to evaluate an expression

#    tag   MX   MY   RZ
mass  3    $m   $m    0
mass  4    $m   $m    0


# Define dynamic loads
# --------------------

# Set some parameters
set outFile ARL360.g3

# Source in TCL proc to read PEER SMD record
source ReadSMDFile.tcl

# Permform the conversion from SMD record to OpenSees record
#              inFile     outFile dt
ReadSMDFile ARL360.at2 $outFile dt

# Set time series to be passed to uniform excitation
set accelSeries "Path -filePath $outFile -dt $dt -factor $g"

# Create UniformExcitation load pattern
#                         tag dir 
pattern UniformExcitation  2   1  -accel $accelSeries


# ----------------------------------------------------
# End of additional modelling for dynamic loads
# ----------------------------------------------------


# ---------------------------------------------------------
# Start of modifications to analysis for transient analysis
# ---------------------------------------------------------

# Delete the old analysis and all it's component objects
wipeAnalysis

# Create the system of equation, a banded general storage scheme
system BandGeneral

# Create the constraint handler, a plain handler as homogeneous boundary
constraints Plain

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
test NormDispIncr 1.0e-12  10 

# Create the solution algorithm, a Newton-Raphson algorithm
algorithm Newton

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM

# Create the integration scheme, the Newmark with alpha =0.5 and beta =.25
integrator Newmark  0.5  0.25 0.0 0.0 0.0 0.000625

# Create the analysis object
analysis Transient

# ---------------------------------------------------------
# End of modifications to analysis for transient analysis
# ---------------------------------------------------------


# ------------------------------
# Start of recorder generation
# ------------------------------

# Create a recorder to monitor nodal displacements
recorder Node node33.out disp -time -node 3 4 -dof 1 2 3

# Create recorders to monitor section forces and deformations
# at the base of the left column
recorder Element 1 -time -file ele1secForce.out section 1 force
recorder Element 1 -time -file ele1secDef.out   section 1 deformation

# --------------------------------
# End of recorder generation
# ---------------------------------


# ------------------------------
# Finally perform the analysis
# ------------------------------

# Perform an eigenvalue analysis
eigen 2

# Perform the transient analysis
#         N   dt
analyze 2000 0.01
puts " "; puts "Transient analysis completed";

# Perform an eigenvalue analysis
eigen 2

# Print state of node 3
print node 3


