FC90 = gfortran
MPI90 = mpif90
FLAGS = -03

react_files = reaction_test.f90 react.f90 kdtree2.f90 RM_interface.F90
adre_files = dolo_ADRE.f90 ADRE_mod.f90 RM_interface.F90

all: react_exe adre_exe

react_exe: $(react_files)
	$(FC90) RM_interface.F90 -c
	$(FC90) kdtree2.f90 -c
	$(FC90) react.f90 -c
	$(FC90) $(react_files) -o $@ -L/usr/local/lib -lphreeqcrm
# -I/usr/local/include

adre_exe: $(adre_files)
	$(FC90) RM_interface.F90 -c
	$(FC90) ADRE_mod.f90 -c
	$(FC90) $(adre_files) -o $@ -L/usr/local/lib -lphreeqcrm

clean: 
	rm *_exe *.mod *.o *.out
