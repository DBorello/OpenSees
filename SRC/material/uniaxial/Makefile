include ../../../Makefile.def

FEAP_OBJS = 	FedeasMaterial.o \
	FedeasBond1Material.o \
	FedeasBond2Material.o \
	FedeasConcr1Material.o \
	FedeasConcr2Material.o \
	FedeasConcr3Material.o \
	FedeasHard1Material.o \
	FedeasHyster1Material.o \
	FedeasHyster2Material.o \
	FedeasSteel1Material.o \
	FedeasSteel2Material.o \
	Bond1.o  Concr1.o  Concr3.o  Hyster1.o  Steel1.o \
	Bond2.o  Concr2.o  Hard1.o   Hyster2.o  Steel2.o

DRAIN_OBJS = DrainMaterial.o \
	Drain00Material.o anal00.o common00.o \
	Drain17Material.o anal17.o common17.o

OBJS       = UniaxialMaterial.o ElasticMaterial.o \
	ElasticPPMaterial.o ParallelMaterial.o \
	HardeningMaterial.o Concrete01.o Steel01.o \
	HystereticMaterial.o EPPGapMaterial.o \
	ViscousMaterial.o SeriesMaterial.o \
	PathIndependentMaterial.o \
	TclModelBuilderUniaxialMaterialCommand.o \
	Steel1.o Steel2.o FedeasSteel2Material.o FedeasMaterial.o \
	FedeasSteel1Material.o Concr1.o FedeasConcr1Material.o Hyster2.o FedeasHyster2Material.o \
	Concr2.o FedeasConcr2Material.o \
	DrainMaterial.o Drain00Material.o anal00.o common00.o \
	Drain17Material.o anal17.o common17.o EPP_Material.o CableMaterial.o


all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
