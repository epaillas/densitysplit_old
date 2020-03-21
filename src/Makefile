# Compiler choice
FC = gfortran
MPIFC = mpif90

MLIBS = -lm

#################

EXES = CCF_monopole.exe

all: $(EXES)

CCF_monopole.exe: CCF_monopole.o
	$(FC) CCF_monopole.o -o CCF_monopole.exe $(MLIBS)

CCF_monopole.o: CCF_monopole.f90
	$(FC) -c CCF_monopole.f90

install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
