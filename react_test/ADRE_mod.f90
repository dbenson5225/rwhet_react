module ADRE_mod
    implicit none

    type species_list
        character(:), allocatable :: spec
    end type
    type component_list
        character(:), allocatable :: comp
    end type
    type selectout_list
        character(:), allocatable :: head
    end type

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    contains

    ! subroutine to initialize the random number generator seed from clock time
    subroutine init_random_seed()
        integer                            :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))
        call system_clock(count = clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
        deallocate(seed)
    end subroutine init_random_seed

    subroutine phreeqc_input(calcite_in, na_in, mg_in, ca_in, cl_in, co2_in,&
                             dolomite_in, quartz_in)
        double precision, intent(in   ) :: calcite_in, na_in, mg_in, ca_in, cl_in, co2_in,&
                                           dolomite_in, quartz_in
        character*1                     :: tab

        tab = char(9)

        open (unit=11, file='dolomite_chem.in', action='write')
        write (11, *) '# Brine-CO2-Calcite-Quartz system'
        write (11, *) 'SOLUTION 0 Brine'
        write (11, *) tab, 'pH 7.000000 charge'
        write (11, *) tab, 'units mol/L'
        write (11, *) tab, 'temp 60.000000'
        write (11, *) tab, 'pressure 100.000000 atm'
        write (11, '(A, A, f8.6)') tab, 'Na ', Na_in
        write (11, '(A, A, f8.6)') tab, 'Mg ', mg_in
        write (11, '(A, A, f8.6)') tab, 'Ca ', ca_in
        write (11, '(A, A, f8.6)') tab, 'Cl ', cl_in
        write (11, *) 'EQUILIBRIUM_PHASES 0'
        write (11, '(A, A, f8.6)') tab, 'Calcite 0.000000 ', calcite_in
        write (11, '(A, A, f8.6)') tab, 'Dolomite 0.000000 ', dolomite_in
        write (11, '(A, A, f8.6)') tab, 'CO2(g) -3.500000 ', co2_in
        write (11, '(A, A, f8.6)') tab, 'Quartz 0.000000 ', quartz_in
        write (11,*) 'SAVE solution 0'
        write (11,*) 'SOLUTION 1 Domain_brine'
        write (11,*) tab, 'pH 7 charge'
        write (11,*) tab, 'units mol/L'
        write (11,*) tab, 'temp 60.000000'
        write (11,*) tab, 'pressure 100.000000 atm'
        write (11,*) tab, 'Na 0.000000'
        write (11,*) tab, 'Mg 0.000000'
        write (11,*) tab, 'Ca 0.000000'
        write (11,*) tab, 'Cl 0.000000'
        write (11,*) 'EQUILIBRIUM_PHASES 1'
        write (11,*) tab, 'Calcite 0.000000 0.270000'
        write (11,*) tab, 'Dolomite 0.000000 0.000000'
        write (11,*) tab, 'CO2(g) -3.500000 0.000000'
        write (11,*) tab, 'Quartz 0.000000 22.000000'
        write (11,*) 'SELECTED_OUTPUT'
        write (11,*) tab, '-simulation false'
        write (11,*) tab, '-state false'
        write (11,*) tab, '-solution false'
        write (11,*) tab, '-distance false'
        write (11,*) tab, '-time false'
        write (11,*) tab, '-step false'
        write (11,*) tab, '-ph true'
        write (11,*) tab, '-pe false'
        write (11,*) tab, '-equilibrium_phases Calcite Dolomite'
        write (11,*) 'END'
        close (unit=11, status='keep')
    end subroutine

end module ADRE_mod