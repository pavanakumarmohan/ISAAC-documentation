SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= rk.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

.F.o:
	$(FC) -c $(FFLAGS) $*.F 

#
rk.o:		rk.F		$(SRC)/header/common.h  
