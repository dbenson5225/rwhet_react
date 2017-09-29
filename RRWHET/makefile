FC90 = gfortran
MPI90 = mpif90
FLAGS = -03

rwhet_files = rwhet4.1.f90

all: rwhet_exe

rwhet_exe: $(rwhet_files)
	$(FC90) $(rwhet_files) -o $@ -fcheck=all -Wall > compile.log 2>&1

clean: 
	rm *_exe *.mod *.o *.err compile.log
