SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= bldlmx.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

.F.o:
	$(FC) -c $(FFLAGS) $*.F 

#
bldlmx.o:         bldlmx.F     $(SRC)/header/common.h  
