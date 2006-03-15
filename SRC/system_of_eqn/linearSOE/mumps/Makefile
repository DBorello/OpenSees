# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM         = test



OBJS       = MumpsSOE.o \
	MumpsSolver.o 

ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS       = MumpsSOE.o \
	MumpsSolver.o \
	MumpsParallelSOE.o \
	MumpsParallelSolver.o 

endif


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
