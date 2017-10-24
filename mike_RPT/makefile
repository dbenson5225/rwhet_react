FC90 = gfortran
MPI90 = mpif90
FLAGS = -03

RPT_files = RPT_dolo.f90 RPT_mod.f90 RM_interface.F90 kdtree2.f90

all: RPT_exe

RPT_exe: $(RPT_files)
	$(FC90) RM_interface.F90 -c
	$(FC90) kdtree2.f90 -c
	$(FC90) RPT_mod.f90 -c
	$(FC90) $(RPT_files) -o $@ -L/usr/local/lib -lphreeqcrm -g -fcheck=all

clean:
	rm *_exe *.mod *.o *.out *.err dolomite_chem.in myrun*


 # -Wall > compile.log 2>&1
