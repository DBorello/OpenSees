include ../../Makefile.def

OBJS = ActorTruss.o ShadowTruss.o

all: $(OBJS)

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(PROGRAM)

# DO NOT DELETE THIS LINE -- make depend depends on it.
