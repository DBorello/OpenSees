include ../../../Makefile.def

OBJS       = NDMaterial.o ElasticIsotropic3D.o ElasticIsotropicMaterial.o \
	ElasticIsotropicPlaneStrain2D.o ElasticIsotropicPlaneStress2D.o \
	TclModelBuilderNDMaterialCommand.o J2PlaneStrain.o \
	J2PlaneStress.o  J2Plasticity.o J2AxiSymm.o J2ThreeDimensional.o 

all:	$(OBJS)
	@$(CD) $(FE)/material/nD/Template3Dep; $(MAKE);
	@$(CD) $(FE)/material/nD/ucsd; $(MAKE);

# Miscellaneous

tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/material/nD/Template3Dep; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/ucsd; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
