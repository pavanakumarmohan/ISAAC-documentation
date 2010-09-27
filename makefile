MAKE=make

FC = gfortran
DOXYGEN=doxygen

#SRC=  $(HOME)/CFD/ISAAC/src
SRC=  $(CURDIR)
SRCTARFILE= isaacsrc.tar
MKTARFILE=  isaacmk.tar
RCSTARFILE= isaacrcs.tar

##FFLAGS= -ansi -O2 
##FFLAGS= -ansi -O2 -trapuv -C -Olimit 3000
#FFLAGS= -ansi -trapuv -C 
########FFLAGS= -ansi -trapuv -C -mips2 -O2 -noisam
#normalFFLAGS= -ansi -r8 -trapuv -C -mips2 -O2 -noisam -Nn20000
#too slowFFLAGS= -ansi -r8 -trapuv -C -mips2 -O2 -Olimit 5000 -noisam -Nn20000

#SGI Production flags
#FFLAGS= -ansi -r8 -trapuv -C -mips2 -O2 -noisam -Nn20000
#####FFLAGS= -mips2 -O2 -Olimit 3000 -noisam

#g77 flags
FFLAGS= -O2
CPPFLAGS = -DCPU_TIME

OBJ= main.o \
	bc.o \
	fourth.o \
	io.o \
	invscd.o \
	mg.o \
	misc.o \
	resid.o \
	viscos.o \
	rk.o \
	march.o \
	af3f.o \
	af3dia.o \
	jacs.o \
	solvr.o \
	bldlmx.o \
	kepsln.o \
	komega.o \
	les.o \
	rij.o

all:
#
	cd $(SRC)/main;\
	$(MAKE) -f main.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)" "CPPFLAGS=$(CPPFLAGS)";\
	ln *.o $(SRC)
#
	cd $(SRC)/rhs;\
	$(MAKE) -f rhs.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)" "CPPFLAGS=$(CPPFLAGS)";\
	ln *.o $(SRC)
#
	cd $(SRC)/explicit;\
	$(MAKE) -f explicit.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)" "CPPFLAGS=$(CPPFLAGS)";\
	ln *.o $(SRC)
#
	cd $(SRC)/implicit;\
	$(MAKE) -f implicit.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)" "CPPFLAGS=$(CPPFLAGS)";\
	ln *.o $(SRC)
#
	cd $(SRC)/march;\
	$(MAKE) -f march.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)" "CPPFLAGS=$(CPPFLAGS)";\
	ln *.o $(SRC)
#
	cd $(SRC)/turb;\
	$(MAKE) -f turb.mk "HOME=$(HOME)" "SRC=$(SRC)" "FFLAGS=$(FFLAGS)" \
		"FC=$(FC)" "CPPFLAGS=$(CPPFLAGS)";\
	ln *.o $(SRC);\
        /bin/rm *.o
#
	$(FC) $(FFLAGS) $(OBJ) -o isaac
	/bin/rm *.o

#
# delete all binary files
#
clean:
	/bin/rm *.o

# 
# create a tar file of the source files 
# use binary mode when using ftp to move this to a new machine 
#
srctar:
	tar -cvf $(SRCTARFILE) `find . \
				-name "*.F" -print \
				-o -name "*.h" -print \
				`

# 
# create a tar file of the source files 
# use binary mode when using ftp to move this to a new machine 
#
mktar:
	tar -cvf $(MKTARFILE)  `find . \
				-name "*.mk" -print \
				-o -name "makefile" -print \
				`


# 
# 
# create a tar file of the RCS files 
# use binary mode when using ftp to move this to a new machine 
#
rcstar:
	tar -cvf $(RCSTARFILE) `find . -name "RCS" -a -type d -print

# Create Doxygen documentation	
html:
	$(DOXYGEN) 
