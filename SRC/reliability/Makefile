include ../../Makefile.def

RELIABILITY_OBJS = $(FE)/reliability/analysis/analysis/FORMAnalysis.o \
		$(FE)/reliability/analysis/analysis/ReliabilityAnalysis.o \
		$(FE)/reliability/analysis/analysis/SimulationAnalysis.o \
		$(FE)/reliability/analysis/analysis/SORMAnalysis.o \
		$(FE)/reliability/analysis/analysis/SystemAnalysis.o \
		$(FE)/reliability/analysis/curvature/CurvaturesBySearchAlgorithm.o \
		$(FE)/reliability/analysis/curvature/FindCurvatures.o \
		$(FE)/reliability/analysis/curvature/FirstPrincipalCurvature.o \
		$(FE)/reliability/analysis/designPoint/FindDesignPoint.o \
		$(FE)/reliability/analysis/designPoint/SearchWithStepSizeAndStepDirection.o \
		$(FE)/reliability/analysis/direction/HLRFSearchDirection.o \
		$(FE)/reliability/analysis/direction/SearchDirection.o \
		$(FE)/reliability/analysis/gFunction/BasicGFunEvaluator.o \
		$(FE)/reliability/analysis/gFunction/GFunEvaluator.o \
		$(FE)/reliability/analysis/gFunction/OpenSeesGFunEvaluator.o \
		$(FE)/reliability/analysis/misc/MatrixOperations.o \
		$(FE)/reliability/analysis/randomNumber/CStdLibRandGenerator.o \
		$(FE)/reliability/analysis/randomNumber/RandomNumberGenerator.o \
		$(FE)/reliability/analysis/sensitivity/OpenSeesSensitivityEvaluator.o \
		$(FE)/reliability/analysis/sensitivity/SensitivityByFiniteDifference.o \
		$(FE)/reliability/analysis/sensitivity/SensitivityEvaluator.o \
		$(FE)/reliability/analysis/stepSize/ArmijoRule.o \
		$(FE)/reliability/analysis/stepSize/FixedStepSizeRule.o \
		$(FE)/reliability/analysis/stepSize/StepSizeRule.o \
		$(FE)/reliability/analysis/transformation/NatafXuTransformation.o \
		$(FE)/reliability/analysis/transformation/XuTransformation.o \
		$(FE)/reliability/domain/components/CorrelationCoefficient.o \
		$(FE)/reliability/domain/components/LimitStateFunction.o \
		$(FE)/reliability/domain/components/RandomVariable.o \
		$(FE)/reliability/domain/components/RandomVariablePositioner.o \
		$(FE)/reliability/domain/components/ReliabilityDomainComponent.o \
		$(FE)/reliability/domain/components/ReliabilityDomain.o \
		$(FE)/reliability/domain/distributions/BetaRV.o \
		$(FE)/reliability/domain/distributions/GumbelRV.o \
		$(FE)/reliability/domain/distributions/ParetoRV.o \
		$(FE)/reliability/domain/distributions/Type1LargestValueRV.o \
		$(FE)/reliability/domain/distributions/UniformRV.o \
		$(FE)/reliability/domain/distributions/ChiSquareRV.o \
		$(FE)/reliability/domain/distributions/LaplaceRV.o \
		$(FE)/reliability/domain/distributions/RayleighRV.o \
		$(FE)/reliability/domain/distributions/Type1SmallestValueRV.o \
		$(FE)/reliability/domain/distributions/WeibullRV.o \
		$(FE)/reliability/domain/distributions/ExponentialRV.o \
		$(FE)/reliability/domain/distributions/LognormalRV.o \
		$(FE)/reliability/domain/distributions/ShiftedExponentialRV.o \
		$(FE)/reliability/domain/distributions/Type2LargestValueRV.o \
		$(FE)/reliability/domain/distributions/GammaRV.o \
		$(FE)/reliability/domain/distributions/NormalRV.o \
		$(FE)/reliability/domain/distributions/ShiftedRayleighRV.o \
		$(FE)/reliability/domain/distributions/Type3SmallestValueRV.o \
		$(FE)/reliability/FEsensitivity/SensitivityAlgorithm.o \
		$(FE)/reliability/FEsensitivity/SensitivityIntegrator.o \
		$(FE)/reliability/FEsensitivity/StaticSensitivityIntegrator.o \
		$(FE)/reliability/tcl/TclReliabilityBuilder.o

# Compilation control
all:         
	     @$(CD) $(FE)/reliability/analysis; $(MAKE);
	     @$(CD) $(FE)/reliability/domain; $(MAKE);
	     @$(CD) $(FE)/reliability/FEsensitivity; $(MAKE);
	     @$(CD) $(FE)/reliability/tcl; $(MAKE);
	     $(AR) $(ARFLAGS) $(RELIABILITY_LIBRARY) $(RELIABILITY_OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	     @$(CD) $(FE)/reliability/analysis; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/FEsensitivity; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/tcl; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
