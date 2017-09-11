FC90 = gfortran
MPI90 = mpif90
FLAGS = -03

rwhet_files = rwhet4.1.f90

all: rwhet_exe

rwhet_exe: $(rwhet_files)
	$(FC90) $(rwhet_files) -o $@

clean: 
	rm *_exe *.mod *.o *.out *.err
