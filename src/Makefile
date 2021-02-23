# Compiler choice
FC = gfortran
MPIFC = mpif90

MLIBS = -lm

#################

EXES = gaussian_filter.exe tophat_filter.exe tophat_velocity.exe

all: $(EXES)

gaussian_filter.exe: gaussian_filter.o
	$(FC) gaussian_filter.o -o gaussian_filter.exe $(MLIBS)

tophat_filter.exe: tophat_filter.o
	$(FC) tophat_filter.o -o tophat_filter.exe $(MLIBS)

tophat_velocity.exe: tophat_velocity.o
	$(FC) tophat_velocity.o -o tophat_velocity.exe $(MLIBS)


gaussian_filter.o: gaussian_filter.f90
	$(FC) -c gaussian_filter.f90

tophat_filter.o: tophat_filter.f90
	$(FC) -c tophat_filter.f90

tophat_velocity.o: tophat_velocity.f90
	$(FC) -c tophat_velocity.f90

install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
