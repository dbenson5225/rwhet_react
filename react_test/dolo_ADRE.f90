program dolo_ADRE
    use ADRE_mod
    use PhreeqcRM
    implicit none

! ======== Parameters ========
integer, parameter                 :: nthreads = 1 ! number of openmp threads for reaction module
double precision, parameter        :: maxtime = 60e3 ! 1000 MIN
! integer, parameter                 :: maxtime = 15e3 ! 250 MIN
! double precision, parameter        :: dx = 1e-3
double precision, parameter        :: dx = 1e-1 ! ****for faster testing
integer, parameter                 :: ncell = nint(1.0d0/dx) - 1
    ! this could go wrong if dx is a weird number
    ! subtract 1 because won't be calculating chemistry for boundary cell
integer, parameter                 :: ntrans = ncell + 1
    ! number of cells for transport
double precision, parameter        :: dt = 1e0
integer, parameter                 :: nsteps = nint(maxtime/dt) + 1
double precision, parameter        :: save_dt = dt * 100.0d0
integer, parameter                 :: save_steps = nint(maxtime/(save_dt)) + 1
    ! time step for saving concentrations to plot them
double precision, parameter        :: Omega = 0.5
double precision, parameter        :: darvel = 1.2e-5 ! Darcy velocity [m/s]
double precision, parameter        :: init_porosity = 0.5
double precision, parameter        :: init_v = darvel/init_porosity
    ! advective velocity
    ! ****update this on the fly as porosity changes
double precision, parameter        :: alpha_l = 0.005 ! longitudinal dispersivity
double precision, parameter        :: init_D = alpha_l*init_v
double precision, parameter        :: v(ntrans - 1) = init_v, D(ntrans - 1) = init_D
    ! dispersion coefficient in long. direction (transverse is zero)
    ! ****this might be wrong, since units in paper give concs as mol/m^3

! ======== general variables ========
double precision                   :: init_calcite, na_inflow, mg_inflow,&
                                           ca_inflow, cl_inflow, co2_inflow,&
                                           dolomite_inflow, quartz_inflow
double precision, dimension(ncell) :: den_0, prs_0, tmp_0, sat_0, por_0, vol_0
double precision                   :: cur_time
integer                            :: i, j, m, id, status, save_on, ngrd
integer                            :: ncomp, nchem, so_col
integer                            :: ic1(ncell, 7)
type(component_list), allocatable  :: comp_list(:)
type(selectout_list), allocatable  :: head_list(:)
integer, dimension(ncell)          :: mask
character(25)                      :: tempname
double precision, allocatable      :: bc_conc(:, :), comp_conc(:, :),&
                                      sout(:, :), concs(:, :),&
                                      plot_concs(: , :, :), plot_times(:)
integer                            :: bc1(1) ! this is for one boundary condition

double precision :: testmat(2, 2, 2)

cur_time = 0

print *, '======================================='
print *, 'grid_Pe =', dx/alpha_l
print *, 'CFL = ', (init_v*dt)/dx
print *, '======================================='
! pause

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

call phreeqc_input(init_calcite, na_inflow, mg_inflow, ca_inflow, cl_inflow,&
                   co2_inflow, dolomite_inflow, quartz_inflow)

! ======== DEFINE PHYSICAL CONDITIONS ========
! These don't make much difference for this problem
den_0 = 1.0 ! Density
prs_0 = 98.6923 ! Pressure 98.6923 atm = 100 bar
tmp_0 = 60.0 ! Temperature
sat_0 = 1.0 ! Saturation
por_0 = init_porosity ! Porosity
vol_0 = (Omega * 1e3)/dble(ncell) ! Representative Volume of cell (L), default = 1

! print chemistry mask to print only first element
mask = 0
mask(1) = 1

id = RM_Create(ncell, nthreads)
status = RM_LoadDatabase(id, '/usr/local/share/doc/phreeqcrm/database/phreeqc.dat')
if (status < 0) then
    print *, 'Database Load Error'
    call exit(status)
endif

status = RM_OpenFiles(id)
status = RM_SetRepresentativeVolume(id, vol_0)
status = RM_SetSaturation(id, sat_0)
status = RM_SetPorosity(id, por_0)
status = RM_SetTime(id, cur_time)
status = RM_SetTimeStep(id, dt)
status = RM_SetDensity(id, den_0)
status = RM_SetTemperature(id, tmp_0)
status = RM_SetPressure(id, prs_0)
status = RM_SetSelectedOutputOn(id, 1) ! turn on selected output to get dolomite/calcite values
status = RM_SetUnitsSolution(id, 2) ! 2 = mol/L
status = RM_SetUnitsPPassemblage(id, 0) ! 0 = mol/L of representative volume
status = RM_SetPrintChemistryMask(id, mask)
status = RM_SetPrintChemistryOn(id, 0, 0, 0)

status = RM_RunFile(id, 1, 1, 1, 'dolomite_chem.in')

! for some reason, species count needs to be done before the initial condition
ncomp = RM_FindComponents(id)
! get number of grid cells in model (this can change, so it needs to be checked)
ngrd = RM_GetGridCellCount(id)
if (ngrd /= ncell) then
    print *, '****** Cell count differs from particle number ******'
    call exit(status)
endif

nchem = RM_GetChemistryCellCount(id)
if (nchem /= ncell) then
    print *, '****** Chem count differs from particle number ******'
    call exit(status)
endif

! =======================================================================
                ! INITIAL/BOUNDARY CONDITIONS
! =======================================================================
! ===  (1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE, 
! ===  (4) SURFACE, (5) GAS_PHASE, (6) SOLID_SOLUTIONS, and (7) KINETICS
 
ic1 = -1
ic1(:, 1) = 1
ic1(:, 2) = 1

allocate(bc_conc(1, ncomp))
    ! must be 2D array for module--requires singleton dimension in this case
bc1 = 0 ! corresponds to solution zero in input file

status = RM_InitialPhreeqc2Module(id, ic1)
!            >>>>> End of initial condition <<<<<
! =======================================================================

allocate(comp_list(ncomp))
! print *, 'ncomp = ', ncomp
do i = 1, ncomp
    status = RM_GetComponent(id, i, tempname)
    allocate(character(len_trim(tempname)) :: comp_list(i)%comp)
    comp_list(i)%comp = trim(tempname)
    ! print *, comp_list(i)%comp
enddo

status = RM_RunCells(id)

allocate(comp_conc(ncell, ncomp))
status = RM_GetConcentrations(id, comp_conc)

status = RM_InitialPhreeqc2Concentrations(id, bc_conc, 1, bc1)

so_col = RM_GetSelectedOutputcolumncount(id)
allocate(sout(ncell, so_col))
status = RM_GetSelectedOutput(id, sout)

! dimension is ncell x ncomp - 1 because we will not be doing transport
    ! calculations for H2O
! For this specific case:
! (1) H, (2) O, (3) charge, (4) C, (5) Ca, (6) Cl, (7) Mg, (8) Na, (9) Si
! **** plot_concs ****
! (10) pH, (11) calcite, (12) dolomite
allocate(concs(ntrans, ncomp - 1), plot_concs(ntrans, ncomp + 2, save_steps),&
         plot_times(save_steps))
    ! *****probably will want to only save certain species later, but for now,
        ! let's keep them all
concs(1, :) = bc_conc(1, 2 : ncomp)
concs(2 : ntrans, 1 : 9) = comp_conc(:, 2 : ncomp)
concs(2 : ntrans, 10 : 12) = sout
plot_concs(:, :, 1) = concs
plot_times = (/((i - 1) * save_dt, i = 1, save_steps)/)
print *, '188'

j = 2
! time stepping
do m = 1, nsteps
    call advect(concs(:, 1 : 9), v, dx, dt, ntrans)
    call diffuse(concs(:, 1 : 9), D, dx, dt, ntrans)

    if (cur_time >= plot_times(j)) then
        plot_concs(:, :, j) = concs
        j = j + 1
    endif

    cur_time = cur_time + dt
    status = RM_SetTime(id, cur_time)
    print *, '203'

    comp_conc(:, 2 : ncomp) = concs(2 : ntrans, 1 : 9)
    print *, '206'
    status = RM_SetConcentrations(id, comp_conc)
    print *, '208'
    status = RM_RunCells(id)
    print *, '210'
    status = RM_GetConcentrations(id, comp_conc)
    print *, '212'
    concs(2 : ntrans, 1 : 9) = comp_conc(:, 2 : ncomp)
    print *, '214'
    status = RM_GetSelectedOutput(id, sout)
    print *, '216'
    concs(2 : ntrans, 10 : 12) = sout
    print *, '218'

enddo

open (unit=12, file='time_concs.txt', action='write')
write (12, *) shape(plot_concs)
write (12, *) plot_concs
close (unit=12, status='keep')

status = RM_Destroy(id)

end program dolo_ADRE
