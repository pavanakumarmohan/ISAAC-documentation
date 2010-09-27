SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= rij.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

.F.o:
	$(FC) -c $(FFLAGS) $*.F 

#
rij.o:       rij.F	$(SRC)/header/common.h  
