#  Written: fmk 
#  Created: 10/97

include ../../Makefile.def

NDARRAY  = BJtensor.o \
           BJmatrix.o \
	   BJvector.o \
	   nDarray.o  \
	   basics.o   \
	   stresst.o  \
	   straint.o  \


OBJS       = Vector.o ID.o Matrix.o  $(NDARRAY)

################### TARGETS ########################
all: $(OBJS); $(AR) -rv $(NDARRAY_LIBRARY) $(NDARRAY)

libndarray.a: $(NDARRAY)
	$(AR) -rv $(NDARRAY_LIBRARY) $(NDARRAY)


test: $(OBJS) main.o
	$(LINKER) main.o $(OBJS) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(FE_LIBRARY) \
	-o test


ndarray: libndarray.a

math_tst: $(NDARRAY) math_tst.o
	$(LINKER) math_tst.o $(NDARRAY) -o math_tst

elast_tst: $(NDARRAY) elast_tst.o
	$(LINKER) elast_tst.o $(NDARRAY) -o elast_tst

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test
	@$(RM) $(RMFLAGS) math_tst elast_tst 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core test
	@$(RM) $(RMFLAGS) math_tst  
	@$(RM) $(RMFLAGS) elast_tst 
	@$(RM) $(RMFLAGS) test.\$\$\$ 
	@$(RM) $(RMFLAGS) *.ti 
	@$(RM) $(RMFLAGS) tca.map 
	@$(RM) $(RMFLAGS) .inslog*

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
