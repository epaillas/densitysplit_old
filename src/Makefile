# Compiler choice
FC = gfortran
MPIFC = mpif90

MLIBS = -lm

#################

EXES = density_profiles.exe

all: $(EXES)

density_profiles.exe: density_profiles.o
	$(FC) density_profiles.o -o density_profiles.exe $(MLIBS)

density_profiles.o: density_profiles.f90
	$(FC) -c density_profiles.f90

install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
