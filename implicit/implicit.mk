SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= af3f.o af3dia.o jacs.o solvr.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

.F.o:
	$(FC) -c $(FFLAGS) $*.F 

#
af3f.o:		af3f.F		$(SRC)/header/common.h  
af3dia.o:	af3dia.F	$(SRC)/header/common.h  
jacs.o:		jacs.F		$(SRC)/header/common.h  
solvr.o:	solvr.F		$(SRC)/header/common.h  
