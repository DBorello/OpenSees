include ../../Makefile.def

CLIENT = client
SERVER = server

OBJS = ActorTruss.o ShadowTruss.o ClientMain.o ServerMain.o

all: $(OBJS)
	@$(LINKER) $(LINKFLAGS) ClientMain.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o $(CLIENT)
	@$(LINKER) $(LINKFLAGS) ServerMain.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o $(SERVER)

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
