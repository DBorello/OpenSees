include ../../../Makefile.def

OBJS       = InternalSpring.o MP_ForJoint2D.o SimpleJoint2D.o TclSimpleJoint2dCommand.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
