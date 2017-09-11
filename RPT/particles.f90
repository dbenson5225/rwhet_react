program particles
    use ADRE_mod
    use PhreeqcRM
    implicit none

! ======== Parameters ========
integer, parameter          :: npart = 1e4 ! number of particles
integer, parameter          :: omega = 0.5 ! domain size
double precision, parameter :: darvel = 1.2e-5 ! Darcy velocity [m/s]
double precision, parameter :: init_porosity = 0.5
double precision, parameter :: init_v = darvel/init_porosity
double precision, parameter :: alpha_l = 0.005 ! longitudinal dispersivity
double precision, parameter :: init_D = alpha_l*init_v
double precision            :: v(ntrans - 1) = init_v, D(ntrans - 1) = init_D ! vectors for v and D

! ======== general variables ========
double precision                   :: init_calcite, na_inflow, mg_inflow,&
                                           ca_inflow, cl_inflow, co2_inflow,&
                                           dolomite_inflow, quartz_inflow
double precision, dimension(npart) :: den_0, prs_0, tmp_0, sat_0, por_0, vol_0,&
                                      solvol_post, por_v
double precision                   :: cur_time
integer                            :: i, j, m, id, status, save_on, ngrd
integer                            :: ncomp, nchem, so_col
integer                            :: ic1(npart, 7)
integer, dimension(npart)          :: mask
double precision, allocatable      :: bc_conc(:, :), comp_conc(:, :),&
                                      sout(:, :), concs(:, :),&
                                      plot_concs(: , :), plot_times(:)
integer                            :: bc1(1) ! this is for one boundary condition

cur_time = 0



end program particles
