SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= march.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

.F.o:
	$(FC) -c $(FFLAGS) $*.F 

#
march.o:	march.F		$(SRC)/header/common.h  $(SRC)/header/histry.h  
