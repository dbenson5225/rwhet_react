FC90 = gfortran
MPI90 = mpif90
# different library locations on different systems
# local library location
# LIB = -L/usr/local/lib -lphreeqcrm
# mio library location
# LIB = -L/u/au/er/mschmidt1/phreeqcrm/lib -lphreeqcrm
# lab machine library location
# LIB = -L/u/au/er/mschmidt1/lib -lphreeqcrm
# optimization level
# FLAGS = -O3
# for easier debugging
# NOTE: -g only works for gfortran... i think
FLAGS = -g -fcheck=all -fbacktrace
WARN = -Wall > compile.log 2>&1

# make on mio couldn't take multiple variables for compiling below, so LIB
# gets concatenated with FLAGS
# FLAGS += $(LIB)

rrwhet_files = blah.f90 kdtree2.f90 rpt_mod.f90

all: rrwhet_exe

rrwhet_exe: $(rrwhet_files)
	$(FC90) kdtree2.f90 -c $(FLAGS)
	-$(FC90) blah.f90 -c $(FLAGS) 2> /dev/null
	# NOTE: the above generates an error since it needs RPT_mod, but the error
		# is piped to /dev/null--change this to see warnings/errors
	$(FC90) rpt_mod.f90 -c $(FLAGS)
	$(FC90) $(rrwhet_files) -o $@ $(FLAGS)

clean:
	rm -r *_exe *.mod *.o *.out *.err dolomite_chem.in myrun* *.dSYM *.inp compile.log

# get number or cores on current machine
# sysctl -n hw.ncpu
