include ../../Makefile.def


################### TARGETS ########################
all:
	@$(CD) $(FE)/optimization/tcl; $(MAKE);
	@$(CD) $(FE)/optimization/domain; $(MAKE);
	@$(CD) $(FE)/optimization/SNOPT; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test

clean: tidy

spotless: clean

wipe: spotless
	@$(CD) $(FE)/optimization/tcl; $(MAKE) wipe;
	@$(CD) $(FE)/optimization/domain; $(MAKE) wipe;
	@$(CD) $(FE)/optimization/SNOPT; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
