program reaction_test
    use ADRE_mod
    use PhreeqcRM
    implicit none

integer, parameter                 :: n = 1e1, d = 3, corr = n/2, nthreads = 1
double precision, parameter        :: r2 = 0.01**2
integer, parameter                 :: cgsize = ceiling(n*((4.0d0/3.0d0)*pi*r2))
double precision, allocatable      :: x(:, :)
integer                            :: correltime = 1 ! not quite sure what this does 
integer                            :: i, id, status, save_on, ngrd
double precision                   :: start, finish_tree, finish
double precision                   :: calcite_in, Na_in, mg_in, ca_in, Cl_in, co2_in,&
                                      dolomite_in, quartz_in
double precision, dimension(n)     :: den_0, prs_0, tmp_0, sat_0, por_0, vol_0
integer, dimension(n)              :: mask
double precision                   :: cur_time, dt
integer                            :: ncomp, nspec, nchem, so_col, so_row
integer                            :: ic1(n, 7)
type(species_list), allocatable    :: spec_list(:)
type(component_list), allocatable  :: comp_list(:)
type(selectout_list), allocatable  :: head_list(:)
character(25)                      :: tempname
double precision, allocatable      :: bc_conc(:, :), comp_conc(:, :), spec_conc(:, :), sout(:, :)
integer                            :: bc1(1) ! this is for one boundary condition

cur_time = 0
dt = 0.1

! allocate (closeguys(n))
! allocate (x(d, n))

! call init_random_seed
! call random_number(x)
! call cpu_time(start)
! call maketree(tree, x, d, n)
! call cpu_time(finish_tree)

! this is initial injection brine speciation
calcite_in = 0.0;
Na_in = 0.884882;
mg_in = 0.0491601;
ca_in = 0.00983202;
Cl_in = 1.00287;
co2_in = 0.36867;
dolomite_in = 0.0;
quartz_in = 0.0;

! ======== DEFINE PHYSICAL CONDITIONS ========
! These don't make much difference for this problem
den_0 = 1.0;  ! Density
prs_0 = 98.6923;  ! Pressure 98.6923 atm = 100 bar
tmp_0 = 60.0; ! Temperature
sat_0 = 1.0;  ! Saturation
por_0 = 0.5;  ! Porosity
vol_0 = 1.0;  ! Representative Volume of cell (L), default = 1

! print chemistry mask to print only first element
mask = 0
mask(1) = 1

call phreeqc_input(calcite_in, na_in, mg_in, ca_in, cl_in, co2_in, dolomite_in, quartz_in)

id = RM_Create(n, nthreads)
status = RM_LoadDatabase(id, '/usr/local/share/doc/phreeqcrm/database/phreeqc.dat')
if (status < 0) then
    print *, 'Database Load Error'
    call exit(status)
endif

save_on = RM_SetSpeciesSaveOn(id, 1)
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
nspec = RM_GetSpeciesCount(id)
! get number of grid cells in model (this can change, so it needs to be checked)
ngrd = RM_GetGridCellCount(id)
if (ngrd /= n) then
    print *, '****** Cell count differs from particle number ******'
    call exit(status)
endif

nchem = RM_GetChemistryCellCount(id)

! =======================================================================
                ! INITIAL CONDITION
! =======================================================================
! ===  (1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE, 
! ===  (4) SURFACE, (5) GAS_PHASE, (6) SOLID_SOLUTIONS, and (7) KINETICS
 
ic1 = -1

ic1(:, 1) = 1
ic1(:, 2) = 1

status = RM_InitialPhreeqc2Module(id, ic1)
!            >>>>> End of initial condition <<<<<
! =======================================================================

allocate(spec_list(nspec))
print *, 'nspec = ', nspec
do i = 1, nspec
    status = RM_GetSpeciesName(id, i, tempname)
    allocate(character(len_trim(tempname)) :: spec_list(i)%spec)
    spec_list(i)%spec = trim(tempname)
    print *, spec_list(i)%spec
enddo

allocate(comp_list(ncomp))
print *, 'ncomp = ', ncomp
do i = 1, ncomp
    status = RM_GetComponent(id, i, tempname)
    allocate(character(len_trim(tempname)) :: comp_list(i)%comp)
    comp_list(i)%comp = trim(tempname)
    print *, comp_list(i)%comp
enddo

allocate(bc_conc(1, ncomp))
bc1 = 0 ! corresponds to solution zero in input file--must be 2D array
status = RM_InitialPhreeqc2Concentrations(id, bc_conc, 1, bc1)

allocate(comp_conc(n, ncomp), spec_conc(n, nspec))
! not sure if this initial get is necessary
status = RM_GetConcentrations(id, comp_conc)
status = RM_GetSpeciesConcentrations(id, spec_conc)

status = RM_RunCells(id)
status = RM_GetConcentrations(id, comp_conc)
so_col = RM_GetSelectedOutputcolumncount(id)
so_row = RM_GetSelectedOutputrowcount(id)
allocate(sout(so_row, so_col))
status = RM_GetSelectedOutput(id, sout)

allocate(head_list(so_col))
do i = 1, so_col
    status = RM_GetSelectedOutputheading(id, i, tempname)
    allocate(character(len_trim(tempname)) :: head_list(i)%head)
    head_list(i)%head = trim(tempname)
    print *, head_list(i)%head
enddo

status = RM_Destroy(id)

print *, 'status = ', status

! do i = 1, n
!     call find_neighbs(tree, i, correltime, r2, cgsize, closeguys(i)%indices, closeguys(i)%dists)
! end do

! do i = 1, n
!     write (*,*) 'i = ', i, 'closeguys = ', closeguys(i)%indices
!     write (*,*) 'distances = ', closeguys(i)%dists
! end do

! call cpu_time(finish)

! write (*,*) 'tree build time = ', finish_tree - start
! write (*,*) 'total time = ', finish - start

! call kdtree2_destroy(tree)
! deallocate(x)

end program reaction_test
