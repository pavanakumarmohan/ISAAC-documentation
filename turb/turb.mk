MAKE=make

SRC=  $(HOME)/ISAAC/src
TURB= $(SRC)/turb

FFLAGS= -ansi -fast

all:
#
	cd $(TURB)/algebraic;\
	$(MAKE) -f algebraic.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)";\
	ln *.o $(TURB);
#
	cd $(TURB)/2eqn;\
	$(MAKE) -f 2eqn.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)";\
	ln *.o $(TURB);
#
	cd $(TURB)/rij;\
	$(MAKE) -f rij.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)";\
	ln *.o $(TURB);
#
	cd $(TURB)/les;\
	$(MAKE) -f les.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)";\
	ln *.o $(TURB);

clean:
	/bin/rm *.o
