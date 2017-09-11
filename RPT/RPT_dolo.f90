! reactive partilce tracking simulation for calcite/dolimite system,
! modeling dissolution and precipitation thereof

program RPT_dolo
    use ADRE_mod
    use PhreeqcRM
    implicit none

! ======== Simulation Parameters ========
! double precision, parameter :: maxtime = 60e3 ! 1000 MIN
! integer, parameter          :: maxtime = 15e3 ! 250 MIN
integer, parameter          :: maxtime = 1e3 ! less, for testing
double precision, parameter :: dt = 1e0
double precision, parameter :: dx = 1e-3 ! for grid of solids
integer, parameter          :: npart = 1e4 ! number of mobile particles
integer, parameter          :: nsol = 1e4 ! number of solid particles
integer, parameter          :: nsteps = nint(maxtime/dt)
double precision, parameter :: save_dt = dt * 100.0d0
integer, parameter          :: save_steps = nint(maxtime/(save_dt)) + 1

! ======== Physical Parameters ========
integer, parameter          :: omega = 0.5 ! domain size
double precision, parameter :: darvel = 1.2e-5 ! Darcy velocity [m/s]
double precision, parameter :: init_porosity = 0.5
double precision, parameter :: init_v = darvel/init_porosity
double precision, parameter :: alpha_l = 0.005 ! longitudinal dispersivity
double precision, parameter :: init_D = alpha_l*init_v
! vectors for v and D
double precision            :: v(ntrans - 1) = init_v, D(ntrans - 1) = init_D

! ======== General Variables ========
double precision, allcocatable     :: conc(:, :) ! will be dimension npart x ncomp
double precision, allcocatable     :: solconc(:, :) ! will be dimension nsolids x (1/dt)
double precision                   :: loc(npart) ! 1D locations
double precision                   :: solloc(npart) ! 1D locations
double precision                   :: init_calcite, na_inflow, mg_inflow,&
                                      ca_inflow, cl_inflow, co2_inflow,&
                                      dolomite_inflow, quartz_inflow
double precision, dimension(npart) :: den_0, prs_0, tmp_0, sat_0, por_0, vol_0,&
                                      solvol_post, por_v
double precision                   :: cur_time
integer                            :: i, j, m, id, status, save_on, ngrd, nsolids
integer                            :: ncomp, nchem, so_col
integer                            :: ic1(npart, 7)
integer, dimension(npart)          :: mask
double precision, allocatable      :: bc_conc(:, :), comp_conc(:, :),&
                                      sout(:, :), concs(:, :),&
                                      plot_concs(: , :), plot_times(:)
integer                            :: bc1(1) ! this is for one boundary condition

cur_time = 0
nsolids = 2 ! calcite and dolomite

! initial concentration [mol/L]
init_calcite = .270865 ! = 270.865 mol/m^3

! inflow concentration pCO2 [atm]
! currently assumes molal = mol/kg SOLVENT
! pco2_in = 10.0788
! inflow concentration pCO2 [mol/L]
! currently assumes molal = mol/kg SOLVENT
! co2_inflow = 0.36867
co2_inflow = 10.0 ! this is from Leal's input file

! inflow concentrations [mol/L]
! currently assumes molal = mol/kg SOLVENT
na_inflow = 0.884882
mg_inflow = 0.0491601
ca_inflow = 0.00983202
cl_inflow = 1.00287


end program RPT_dolo
