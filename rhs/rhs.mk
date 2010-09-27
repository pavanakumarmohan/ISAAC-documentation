SRC= $(HOME)/ISAAC/src

FFLAGS= -ansi -fast

OBJS= bc.o fourth.o io.o invscd.o mg.o misc.o resid.o viscos.o

#
all: $(OBJS)

.f,v.f:
	$(CO) $*.f

.F.o:
	$(FC) -c $(FFLAGS) $*.F 

#
bc.o:		bc.F		$(SRC)/header/common.h  
fourth.o:	fourth.F	$(SRC)/header/common.h  
io.o:		io.F		$(SRC)/header/common.h  $(SRC)/header/histry.h  
invscd.o:	invscd.F	$(SRC)/header/common.h  
mg.o:		mg.F		$(SRC)/header/common.h  
misc.o:		misc.F		$(SRC)/header/common.h  
resid.o:	resid.F		$(SRC)/header/common.h  
viscos.o:	viscos.F	$(SRC)/header/common.h  
