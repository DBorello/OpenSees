# Makefile for fe objects

include ../../../../Makefile.def


PROGRAM       = test


#OBJS       = PetscSOE.o PetscSolver.o ShadowPetscSOE.o ActorPetscSOE.o
OBJS       = PetscSOE.o PetscSolver.o PetscSparseSeqSolver.o

all:    $(OBJS)

test:    main.o $(OBJS)
	$(LINKER) $(LINKFLAGS) $(OBJS) main.o  $(MATRIX_LIBS) $(ERROR_LIBS) $(GRAPH_LIBS)\
        $(SysOfEqn_LIBS) $(MACHINE_NUMERICAL_LIBS) $(PETSC_LIB) -o $(PROGRAM)

ERROR_LIBS = $(FE)/handler/StandardStream.o \
	$(FE)/handler/FileStream.o \
	$(FE)/handler/OPS_Stream.o

GRAPH_LIBS   = $(FE)/graph/graph/Graph.o \
	$(FE)/graph/graph/Vertex.o \
	$(FE)/graph/graph/VertexIter.o \
	$(FE)/tagged/TaggedObject.o \
	$(FE)/tagged/storage/ArrayOfTaggedObjects.o \
	$(FE)/tagged/storage/ArrayOfTaggedObjectsIter.o 	

MATRIX_LIBS   = $(FE)/matrix/Matrix.o \
        $(FE)/matrix/Vector.o \
        $(FE)/matrix/ID.o \
	$(FE)/nDarray/basics.o \
	$(FE)/nDarray/nDarray.o \
	$(FE)/nDarray/BJmatrix.o \
	$(FE)/nDarray/BJvector.o \
	$(FE)/nDarray/BJtensor.o \
	$(FE)/nDarray/straint.o  \
	$(FE)/nDarray/stresst.o

SysOfEqn_LIBS = $(FE)/system_of_eqn/SystemOfEqn.o \
		$(FE)/system_of_eqn/linearSOE/LinearSOE.o \
		$(FE)/system_of_eqn/linearSOE/DomainSolver.o \
		$(FE)/system_of_eqn/linearSOE/LinearSOESolver.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSOE.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSolver.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSOE.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSolver.o \
		$(FE)/system_of_eqn/linearSOE/sparseGEN/SuperLU.o \
		$(FE)/system_of_eqn/Solver.o \
		$(FE)/actor/actor/MovableObject.o 

# Miscellaneous

tidy:	
	@rm -f Makefile.bak *~ #*# core

clean: tidy
	@rm -f $(OBJS) *.o

spotless: clean
	@rm -f $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.


