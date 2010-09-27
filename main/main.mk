SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= main.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

#.F.o:
#	$(FC) -c $(FFLAGS) $*.F 

#
main.o:         main.F     $(SRC)/header/common.h  $(SRC)/header/histry.h
