include ../../../Makefile.def

OBJS       = 

# Compilation control
all:         
	     @$(CD) $(FE)/reliability/analysis/analysis; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/curvature; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/designPoint; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/direction; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/gFunction; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/misc; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/randomNumber; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/sensitivity; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/stepSize; $(MAKE);
	     @$(CD) $(FE)/reliability/analysis/transformation; $(MAKE);


# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	     @$(CD) $(FE)/reliability/analysis/analysis; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/curvature; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/designPoint; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/direction; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/gFunction; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/misc; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/randomNumber; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/sensitivity; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/stepSize; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/analysis/transformation; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
