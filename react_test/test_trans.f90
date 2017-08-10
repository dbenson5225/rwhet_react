program test_trans
    use ADRE_mod
    use PhreeqcRM
    implicit none

integer, parameter                 :: maxtime = 15e3 ! less, for testing
double precision, parameter        :: Omega = 0.5
double precision, parameter        :: dx = 1e-3
! double precision, parameter        :: dx = 1e-1 ! ****for faster testing
integer, parameter                 :: ncell = nint(Omega/dx) - 1
    ! this could go wrong if dx is a weird number
    ! subtract 1 because won't be calculating chemistry for boundary cell
integer, parameter                 :: ntrans = ncell + 1
    ! number of cells for transport
double precision, parameter        :: dt = 1e0
integer, parameter                 :: nsteps = nint(maxtime/dt) + 1
double precision, parameter        :: save_dt = dt * 100.0d0
integer, parameter                 :: save_steps = nint(maxtime/(save_dt)) + 1
    ! time step for saving concentrations to plot them
double precision, parameter        :: darvel = 1.2e-5 ! Darcy velocity [m/s]
double precision, parameter        :: init_porosity = 0.5
double precision, parameter        :: init_v = darvel/init_porosity
    ! advective velocity
    ! ****update this on the fly as porosity changes
double precision, parameter        :: alpha_l = 0.005 ! longitudinal dispersivity
double precision, parameter        :: init_D = alpha_l*init_v
double precision                   :: v(ntrans - 1) = init_v, D(ntrans - 1) = init_D

double precision, allocatable      :: bc_conc(:, :), comp_conc(:, :),&
                                      sout(:, :), concs(:, :),&
                                      plot_concs(: , :), plot_times(:)
integer                            :: i, j, m
double precision                   :: cur_time

allocate(concs(ntrans, 1), plot_concs(ntrans, 1), plot_times(save_steps))

concs(1, 1) = 1
plot_concs = concs
plot_times = (/((i - 1) * save_dt, i = 1, save_steps)/)

open (unit=12, file='time_concs.txt', action='write')
write (12, *) shape(plot_concs), save_steps + 1
write (12, *) plot_concs

j = 2
! time stepping
do m = 1, nsteps
    call advect(concs(:, :), v, dx, dt, ntrans)
    ! call diffuse(concs(:, :), D, dx, dt, ntrans)

    cur_time = cur_time + dt

    if (cur_time >= plot_times(j)) then
        plot_concs = concs
        write (12, *) plot_concs
        j = j + 1
    endif

enddo
print *, '****done****'

! write (12, *) plot_concs
close (unit=12, status='keep')

end program test_trans