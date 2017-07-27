FC90 = gfortran
MPI90 = mpif90
FLAGS = -03

react_files = reaction_test.f90 react.f90 kdtree2.f90

all: react_exe

react_exe: $(react_files)
	$(FC90) kdtree2.f90 -c
	$(FC90) react.f90 -c
	$(FC90) $(react_files) -o $@

clean: 
	rm react_exe *.mod x.txt closeguys.txt *.o *.out
