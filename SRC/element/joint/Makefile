include ../../../Makefile.def

OBJS       = TclJoint2dCommand.o \
	Joint2D.o \
	MP_Joint2D.o \
	BeamColumnJoint2d.o \
	BeamColumnJoint3d.o \
	TclBeamColumnJointCommand.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
