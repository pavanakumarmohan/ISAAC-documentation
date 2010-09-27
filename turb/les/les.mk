SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= les.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

.F.o:
	$(FC) -c $(FFLAGS) $*.F 

#
les.o:       les.F	$(SRC)/header/common.h  
