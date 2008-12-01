include ../../../Makefile.def

LIBRARY         = elasticppf_.so

all: $(LIBRARY)

OBJS = materialtypes.mod materialapi.mod elasticPPf.o 


$(LIBRARY): $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY)   elasticPPf.o -lc

elasticPPC.o: elasticPPC.c
	$(CC++) -fPIC $(INCLUDES)  -c -Wall elasticPPC.c

materialapi.mod: materialAPI.f materialtypee.mod
	$(FC)  -c  materialAPI.f

materialtypes.mod: materialTypes.f
	$(FC)  -c  materialTypes.f

elasticPPf.o: elasticPPC.f materialapi.mod
	$(FC) -shared -fPIC $(INCLUDES)  -c -Wall elasticPPf.f

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
