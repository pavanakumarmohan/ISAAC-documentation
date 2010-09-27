SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= kepsln.o komega.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

.F.o:
	$(FC) -c $(FFLAGS) $*.F 

#
kepsln.o:       kepsln.F	$(SRC)/header/common.h  
komega.o:       komega.F	$(SRC)/header/common.h  
