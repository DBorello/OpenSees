include ../../../../Makefile.def

OBJS       = FiniteDeformationElastic3D.o \
             W.o \
             LogWEnergy.o \
             MooneyRivlinWEnergy.o \
             NeoHookeanWEnergy.o \
             OgdenWEnergy.o \
             SimoPisterWEnergy.o

all: $(OBJS)

#all: $(OBJS)
#	@$(RM) $(RMFLAGS) $(NDARRAY_LIBRARY);
#	$(AR) $(ARFLAGS) $(NDARRAY_LIBRARY) $(NDARRAY)




fd_test: $(OBJS)   fd_test.o
	$(LINKER) $(LINKFLAGS) fd_test.o $(OBJS)  \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o fd_test

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS)  *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.





