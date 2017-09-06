FC90 = gfortran
MPI90 = mpif90
FLAGS = -03

adre_files = dolo_ADRE.f90 ADRE_mod.f90 RM_interface.F90

all: adre_exe

adre_exe: $(adre_files)
	$(FC90) RM_interface.F90 -c
	$(FC90) ADRE_mod.f90 -c
	$(FC90) $(adre_files) -o $@ -L/usr/local/lib -lphreeqcrm

clean: 
	rm *_exe *.mod *.o *.out *.err dolomite_chem.in myrun*
