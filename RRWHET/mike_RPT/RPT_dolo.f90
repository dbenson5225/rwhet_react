! reactive partilce tracking simulation for calcite/dolimite system,
! modeling dissolution and precipitation thereof

program RPT_dolo
    use RPT_mod
    use PhreeqcRM
    implicit none

! ======== Simulation Parameters ========
integer, parameter          :: nthreads = 0
    ! number of openmp threads for reaction module. <= 0 implies the number of
    ! threads equals the number of processors of the computer
! double precision, parameter :: maxtime = 60e3_dp ! 1000 MIN
! double precision, parameter :: maxtime = 15e3_dp ! 250 MIN
double precision, parameter :: maxtime = 5e3_dp ! less, for testing
double precision, parameter :: lowx = 0.0d0, upx = 0.5d0, omega = upx - lowx
    ! domain size, lower, upper bounds
double precision, parameter :: dt = 1e0_dp
integer, parameter          :: nsteps = nint(maxtime/dt)
    ! avg spacing for immobile particles
integer, parameter          :: nipart = 1e3, ncellv = 5e2
    ! number of immobile particles--keeping these the same should make porosity
    ! calculations easier
    ! also, making this 500 for the 0.5m domain makes each cell a liter and
    ! means you don't have to recalculate calcite input, since EQUILIBRIUM_PHASES
    ! are done in moles, not mol/L
double precision, parameter :: dxv = omega/dble(ncellv) ! for velocity grid
integer, parameter          :: npinject = 1e1
    ! number of particles to be injected each time step at boundary
integer, parameter          :: npbnd = npinject * nsteps, npinit = 1e3,&
                               nptot = npbnd + npinit
    ! total number of boundary, initial, and all particles
double precision, parameter :: save_dt = dt * 1e1_dp
integer, parameter          :: save_steps = nint(maxtime/(save_dt)) + 1

! ======== Physical Parameters ========
double precision, parameter :: darvel = 1.2e-5_dp ! Darcy velocity [m/s]
double precision, parameter :: init_porosity = 0.5d0
double precision, parameter :: init_v = darvel/init_porosity
double precision, parameter :: alpha_l = 0.005d0 ! longitudinal dispersivity
double precision, parameter :: init_D = alpha_l*init_v
double precision            :: v(ncellv) = init_v, D(ncellv) = init_D
    ! vectors for v and D

! ======== General Variables ========
type(mparticle)      :: mparts(nptot) ! total mobile particles
type(iparticle)      :: iparts(nipart) ! total immobile particles
integer              :: indices(nptot), lastinject
integer, allocatable :: alive(:)
    ! arrays for indexing to alive particles
integer              :: nactive ! number of active mobile particles
double precision     :: oactive ! portion of active domain


! ======== PHREEQCRM Variables ========
double precision                      :: init_calcite, na_inflow, mg_inflow,&
                                         ca_inflow, cl_inflow, co2_inflow,&
                                         dolomite_inflow, quartz_inflow, vol_0s,&
                                         vol_0b
double precision, dimension(nipart)   :: den_0, prs_0, tmp_0, sat_0, por_0, vol_0,&
                                         solvol_post, por_v
double precision                      :: cur_time
integer                               :: i, j, m, id, status, save_on, ngrd, nsolids
integer                               :: ncomp, nchem, so_col, grid2chem(nipart)
integer                               :: ic1(nipart, 7)
integer, dimension(nipart)            :: mask
double precision, allocatable         :: bc_conc(:, :), comp_conc(:, :),&
                                         sout(:, :), concs(:, :),&
                                         plot_concs(:, :), plot_times(:)
integer                               :: bc1(1) ! this is for one boundary condition

integer :: start, finish, clock_rate, clock_max, new

! call system_clock (start, clock_rate, clock_max)
! call system_clock (finish, clock_rate, clock_max)
! print *, 'time = ', real(finish - start) / real(clock_rate)
! pause


cur_time = 0.0d0
nsolids = 2 ! calcite and dolomite

! make first npinject boundary particles and all initial particles active
mparts%active = .false.
nactive = npinit
mparts(1 : npinit)%active = .true.
indices = (/(i, i = 1, nptot)/) ! this is for easy array-indexing of mparts

! scatter the initial and immobile particles randomly throughout domain
call init_random_seed()
call random_number(mparts(1 : npinit)%loc)
mparts(1 : npinit)%loc = (upx - lowx) * mparts(1 : npinit)%loc + lowx
call random_number(iparts%loc)
iparts%loc = (upx - lowx) * iparts%loc + lowx
! assign boundary particles loc = 0
mparts (npinit + 1 : nptot)%loc = 0.0d0
! initialize all mobile particle's bin to -999 for error catching
mparts%bin = -999

! initial concentration [mol]
init_calcite = 0.270865d0 ! = 270.865 mol/m^3
! since EQUILIBRIUM_PHASES is total moles per reaction cell, this needs to be
! scaled based on initial particles
init_calcite = init_calcite / ((omega * 1e3_dp)/ dble(nipart))

! inflow concentration pCO2 [atm]
! currently assumes molal = mol/kg SOLVENT
! pco2_in = 10.0788
! inflow concentration pCO2 [mol/L]
! currently assumes molal = mol/kg SOLVENT
! co2_inflow = 0.36867
co2_inflow = 10.0d0 ! this is from Leal's input file
! *****maybe rescale this ahead of time, too*****

! inflow concentrations [mol/L]
! currently assumes molal = mol/kg SOLVENT
na_inflow = 0.884882d0
mg_inflow = 0.0491601d0
ca_inflow = 0.00983202d0
cl_inflow = 1.00287d0

call phreeqc_input(init_calcite, na_inflow, mg_inflow, ca_inflow, cl_inflow,&
                   co2_inflow, dolomite_inflow, quartz_inflow)

! ======== DEFINE PHYSICAL CONDITIONS ========
den_0 = 1.0d0 ! Density
prs_0 = 98.6923d0 ! Pressure 98.6923 atm = 100 bar
tmp_0 = 60.0d0 ! Temperature
sat_0 = 1.0d0 ! Saturation
por_0 = init_porosity ! Porosity
! *****this is going to be tricky with particles*****
! Representative Volume of cell (L), default = 1
    ! ***** maybe change this *****
vol_0 = (omega * 1e3_dp)/dble(nipart)
    ! for initial particles: rep vol = avg inter-partile spacing
    ! factor of 1e3 is conversion from m^3 to L

id = RM_Create(nipart, nthreads)
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

status = RM_SetComponentH2O(id, 0)
    ! 0 implies that H2O is not given as a separate component--just H and O
    ! 1 is default with H2O as separate component
    ! ****be careful, setting to 0 might cause problems****

status = RM_RunFile(id, 1, 1, 1, 'dolomite_chem.in')

! for some reason, species count needs to be done before the initial condition
ncomp = RM_FindComponents(id)
! get number of grid cells in model (this can change, so it needs to be checked)
ngrd = RM_GetGridCellCount(id)
if (ngrd /= nipart) then
    print *, '****** Cell count differs from particle number ******'
    call exit(status)
endif

nchem = RM_GetChemistryCellCount(id)
if (nchem /= nipart) then
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

! this makes a list of the names of components
! allocate(comp_list(ncomp))
! ! print *, 'ncomp = ', ncomp
! do i = 1, ncomp
!     status = RM_GetComponent(id, i, tempname)
!     allocate(character(len_trim(tempname)) :: comp_list(i)%comp)
!     comp_list(i)%comp = trim(tempname)
!     ! print *, comp_list(i)%comp
! enddo

status = RM_RunCells(id)

! so_col = RM_GetSelectedOutputcolumncount(id)
! so_row = RM_GetSelectedOutputrowcount(id)
! allocate(sout(so_row, so_col))
! status = RM_GetSelectedOutput(id, sout)

! this makes a list of the names of selected output elements
! allocate(head_list(so_col))
! do i = 1, so_col
!     status = RM_GetSelectedOutputheading(id, i, tempname)
!     allocate(character(len_trim(tempname)) :: head_list(i)%head)
!     head_list(i)%head = trim(tempname)
!     print *, head_list(i)%head
! enddo

allocate(comp_conc(nipart, ncomp))
status = RM_GetConcentrations(id, comp_conc)
call concs_to_iparts(comp_conc, iparts, nipart)

status = RM_InitialPhreeqc2Concentrations(id, bc_conc, 1, bc1)

so_col = RM_GetSelectedOutputcolumncount(id)
allocate(sout(nipart, so_col))
status = RM_GetSelectedOutput(id, sout)
! **** selected output ****
! (1) pH, (2) calcite, (3) change calcite, (4) dolomite, (5), change dolomite

! For this specific case, concs on particles:
! (1) H, (2) O, (3) charge, (4) C, (5) Ca, (6) Cl, (7) Mg, (8) Na, (9) Si
! **** plot_concs ****
! (1) pH, (2) calcite, (3) dolomite
allocate(plot_concs(nipart, 3), plot_times(save_steps))

! find out which bins the initial particles are in
allocate (alive(nactive))
alive = (/(i, i = 1, nactive)/)
mparts(alive)%bin = floor(mparts(alive)%loc/dxv) + 1 ! bin based on vel grid
! do mass/conc balance for all initial particles
call mass_balance(iparts, mparts, nactive, alive, nipart, D, dt, omega)
deallocate (alive)

plot_concs(: , 1) = sout(:, 1)
plot_concs(: , 2) = sout(:, 2)
plot_concs(: , 3) = sout(:, 4)
plot_times = (/((i - 1) * save_dt, i = 1, save_steps)/)

open (unit=12, file='time_concs.txt', action='write')
write (12, *) shape(plot_concs), save_steps
write (12, *) plot_concs

! activate the first injection of boundary particles
nactive = nactive + npinject
lastinject = nactive ! will use this for adding at boundary later
mparts(npinit + 1 : nactive)%active = .true.
! initialize concs for all boundary particles
oactive = init_v * dt
do i = npinit + 1, nptot
    mparts(i)%concs = bc_conc(1, :)
    ! adjust for support volume being distance traveled in first time step in
    ! domain (oactive) divided by npinject, then adjust from m^3 to L
    mparts(i)%concs = mparts(i)%concs * ((oactive * 1e3_dp) / npinject)
enddo

j = 2
! time stepping
do m = 1, 500
    print *, 'step = ', m, ' of ', nsteps
    ! get indices of alive particles
    nactive = count(mparts%active)
    allocate (alive(nactive)) ! maybe preallocate to avoid repeatedly doing this
    alive = pack(indices, mparts%active)
    mparts(alive)%bin = floor(mparts(alive)%loc/dxv) + 1 ! bin based on vel grid

    ! ***** i suspect that passing the whole particles to the transport subs
        ! may not be as efficient as passing only the locations as arrays
    call advect(mparts, v, dt, alive)
    ! print *, 'mparts(1 : 3)%loc = ', mparts(1 : 3)%loc
    ! call absorbhigh(mparts, upx, alive)
    call diffuse(mparts, nactive, D, dt, alive)
    ! print *, 'mparts(1 : 3)%loc = ', mparts(1 : 3)%loc
    call reflectlow(mparts, lowx, alive)
    call absorbhigh(mparts, upx, alive)
    ! get indices of alive particles again, since some could have exited
    deallocate (alive)
    nactive = count(mparts%active)
    allocate (alive(nactive)) ! maybe preallocate to avoid repeatedly doing this
    alive = pack(indices, mparts%active)
    mparts(alive)%bin = floor(mparts(alive)%loc/dxv) + 1

    call mass_balance(iparts, mparts, nactive, alive, nipart, D, dt, omega)
    call iparts_to_concs(comp_conc, iparts, nipart)

    cur_time = cur_time + dt
    status = RM_SetTime(id, cur_time)
    status = RM_SetConcentrations(id, comp_conc)
    status = RM_RunCells(id)
    status = RM_GetConcentrations(id, comp_conc)

    call concs_to_iparts(comp_conc, iparts, nipart)

    ! new = RM_FindComponents(id)
    ! if (new /= ncomp) then
    !     print *, '****** Component count has changed ******'
    ! endif

    if (cur_time >= plot_times(j)) then
    ! if (m >= 30) then
        status = RM_GetSelectedOutput(id, sout)
        plot_concs(: , 1) = sout(:, 1)
        plot_concs(: , 2) = sout(:, 2)
        plot_concs(: , 3) = sout(:, 4)
        write (12, *) plot_concs
        j = j + 1
    endif

    ! ! recalculate porosity and reset cell saturation and solution volume
    ! ! since, as defined by PHREEQCRM,
    ! ! porosity = (solution volume)/(saturation * representative volume)
    ! ! we want to maintain fully saturated cells (sat := 1, the initial value)
    ! ! this implies that porosity = (solution volume)/(representative volume)
    ! ! status = RM_GetSaturation(id, sat_post)
    ! status = RM_GetSolutionVolume(id, solvol_post)
    ! por_v = solvol_post/vol_0
    ! status = RM_SetSaturation(id, sat_0)
    ! status = RM_SetPorosity(id, por_v)

    ! ! alter velocity and diffusion coefficient vectors, based on change in
    ! ! porosity
    ! v = darvel/por_v; ! advective velocity
    ! D = alpha_l * v; ! dispersion coefficient in long. direction (transverse is zero)

    if (lastinject + npinject <= nptot) then
        mparts(lastinject + 1 : lastinject + npinject)%active = .true.
        lastinject = lastinject + npinject
    endif

    deallocate (alive)

enddo

close (unit=12, status='keep')

status = RM_Destroy(id)
deallocate (bc_conc, comp_conc, sout, plot_concs, plot_times)

end program RPT_dolo
