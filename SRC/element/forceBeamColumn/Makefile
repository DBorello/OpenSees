include ../../../Makefile.def

OBJS       = ForceBeamColumn2d.o \
	TclForceBeamColumnCommand.o \
	BeamIntegration.o \
	LobattoBeamIntegration.o UserDefinedBeamIntegration.o \
	HingeMidpointBeamIntegration2d.o HingeRadauBeamIntegration2d.o \
	UserDefinedHingeIntegration2d.o

# Compilation control

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.









