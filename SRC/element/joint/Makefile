include ../../../Makefile.def

OBJS       = InternalSpring.o SimpleJoint2D.o TclJoint2dCommand.o \
		Joint2D.o MP_Joint2D.o MP_SimpleJoint2D.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
