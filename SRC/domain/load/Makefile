include ../../../Makefile.def

OBJS       = Load.o ElementalLoadIter.o NodalLoadIter.o Beam2dUniformLoad.o \
	Beam2dPointLoad.o BrickSelfWeight.o

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
