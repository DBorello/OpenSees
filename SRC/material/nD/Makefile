include ../../../Makefile.def

OBJS       = NDMaterial.o ElasticIsotropic3D.o ElasticIsotropicMaterial.o \
	ElasticIsotropicPlaneStrain2D.o ElasticIsotropicPlaneStress2D.o \
	ElasticIsotropicPlateFiber.o ElasticIsotropicAxiSymm.o\
	TclModelBuilderNDMaterialCommand.o J2PlaneStrain.o \
	J2PlaneStress.o  J2Plasticity.o J2AxiSymm.o J2ThreeDimensional.o \
	BidirectionalMaterial.o J2PlateFiber.o  PressureDependentElastic3D.o


TEMPLATE_OBJS =              \
                MatPoint3D.o         \
                Template3Dep.o						 \
                DP_PS.o													 \
                DP_YS.o													 \
                DP_YS01.o											 \
                EL_LEeq.o											 \
                EL_LEij.o											 \
                EL_NLEeq.o										 \
                EL_NLEij.o										 \
                EL_NLEijMD.o								 \
                EL_S.o														 \
                EL_T.o														 \
                EPState.o											 \
                MD_PS.o													 \
                MD_PS01.o											 \
                MD_YS.o													 \
                VM_PS.o														\
                VM_YS.o														\
                YS.o


all:	$(OBJS) 
	@$(CD) $(FE)/material/nD/soil; $(MAKE);
	@$(CD) $(FE)/material/nD/Template3Dep; $(MAKE);  $(AR) -rv $(MATMOD_LIBRARY) $(TEMPLATE_OBJS)

# Miscellaneous

tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map

spotless: clean
	@$(CD) $(FE)/material/nD/soil; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/Template3Dep; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



