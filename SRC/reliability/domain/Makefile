include ../../../Makefile.def

OBJS       = 

# Compilation control
all:         
	     @$(CD) $(FE)/reliability/domain/components; $(MAKE);
	     @$(CD) $(FE)/reliability/domain/distributions; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	     @$(CD) $(FE)/reliability/domain/components; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain/distributions; $(MAKE) wipe;


# DO NOT DELETE THIS LINE -- make depend depends on it.
