# Compiler choice
FC = gfortran
MPIFC = mpif90

MLIBS = -lm

#################

EXES = CF_monopole.exe CCF_monopole.exe CCF_rmu.exe CCF_spi.exe\
gaussian_filter.exe tophat_filter.exe

all: $(EXES)

CF_monopole.exe: CF_monopole.o
	$(FC) CF_monopole.o -o CF_monopole.exe $(MLIBS)

CCF_monopole.exe: CCF_monopole.o
	$(FC) CCF_monopole.o -o CCF_monopole.exe $(MLIBS)

CCF_rmu.exe: CCF_rmu.o
	$(FC) CCF_rmu.o -o CCF_rmu.exe $(MLIBS)

CCF_spi.exe: CCF_spi.o
	$(FC) CCF_spi.o -o CCF_spi.exe $(MLIBS)

gaussian_filter.exe: gaussian_filter.o
	$(FC) gaussian_filter.o -o gaussian_filter.exe $(MLIBS)

tophat_filter.exe: tophat_filter.o
	$(FC) tophat_filter.o -o tophat_filter.exe $(MLIBS)


CF_monopole.o: CF_monopole.f90
	$(FC) -c CF_monopole.f90

CCF_monopole.o: CCF_monopole.f90
	$(FC) -c CCF_monopole.f90

CCF_rmu.o: CCF_rmu.f90
	$(FC) -c CCF_rmu.f90

CCF_spi.o: CCF_spi.f90
	$(FC) -c CCF_spi.f90

gaussian_filter.o: gaussian_filter.f90
	$(FC) -c gaussian_filter.f90

tophat_filter.o: tophat_filter.f90
	$(FC) -c tophat_filter.f90

install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
