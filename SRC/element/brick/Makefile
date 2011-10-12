include ../../../Makefile.def

OBJS       = Brick.o \
	shp3d.o \
	TclBrickCommand.o \
	BbarBrick.o \
	BbarBrickWithSensitivity.o \
	Twenty_Node_Brick.o \
	TclTwenty_Node_Brick.o


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
