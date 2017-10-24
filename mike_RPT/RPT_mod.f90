module RPT_mod
    use kdtree2_precision_module
    use kdtree2_module
    implicit none

    ! these are char lists for getting the names of components from PHREEQCRM
    ! type component_list
    !     character(:), allocatable :: comp
    ! end type
    ! type selectout_list
    !     character(:), allocatable :: head
    ! end type

    integer, parameter          :: sp = kind(1.0), dp = kind(1.d0)
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    integer, parameter          :: nspec = 8, numsd = 3, num_alloc = 1e4
        ! *****hard-coded nspec--fix to make this general*****
        ! numsd is distances over which particle reactions are considered
        ! num_alloc is the maximum number of nearby particles expected to be
            ! found--may need to adust this up or calculate on the fly in
            ! mass_balance() subroutine, as i don't love having it hard-coded

    ! mobile particle type
    ! NOTE: these are currently for a 1D problem
    type mparticle
        double precision :: loc ! real-valued spatial location
        double precision :: concs(nspec) ! vector of chemical concentrations
        ! *****hard-coded nspec--fix to make this general*****
        logical          :: active
            ! indicates whether particle is active and within the domain
        integer          :: bin
            ! used to indicate which spatial gridpoint a particle resides within
            ! i keep it with the particle to avoid calculating multiple times
            ! in a single time step
    end type
    ! immobile particle type
    type iparticle
        double precision :: loc
        double precision :: concs(nspec)
            ! *****hard-coded nspec--fix to make this general*****
        ! logical        :: active
        ! integer        :: bin
    end type

    ! a couple of derived types for the kd tree search

    ! holds indices of nearby particles
    type index_array
        integer, allocatable :: indices(:)
    end type
    ! holds the distances to the corresponding particle held by index_array
    type dist_array
        double precision, allocatable :: dists(:)
    end type

    contains

    ! subroutine to initialize the random number generator seed from clock time
    subroutine init_random_seed()
        integer              :: i, n, clock
        integer, allocatable :: seed(:)

        call random_seed(size = n)
        allocate (seed(n))
        call system_clock(count = clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
        deallocate(seed)
    end subroutine init_random_seed

    ! this generates the PHREEQC input file for this specific problem
    subroutine phreeqc_input(calcite_in, na_in, mg_in, ca_in, cl_in, co2_in,&
                             dolomite_in, quartz_in)
        double precision, intent(in   ) :: calcite_in, na_in, mg_in, ca_in,&
                                           cl_in, co2_in, dolomite_in, quartz_in
        character*1                     :: tab

        tab = char(9)

        open (unit=11, file='dolomite_chem.in', action='write')
        write (11, *) '# Brine-CO2-Calcite-Quartz system'
        write (11, *) 'SOLUTION 0 Brine'
        write (11, *) tab, 'pH 7.0'
        write (11, *) tab, 'units mol/L'
        write (11, *) tab, 'temp 60.000000'
        write (11, *) tab, 'pressure 98.6923 atm' ! = 100 bar
        write (11, '(A, A, f8.6)') tab, 'Na ', Na_in
        write (11, '(A, A, f8.6)') tab, 'Mg ', mg_in
        write (11, '(A, A, f8.6)') tab, 'Ca ', ca_in
        write (11, '(A, A, f8.6, A)') tab, 'Cl ', cl_in, ' charge'
        write (11, *) 'EQUILIBRIUM_PHASES 0'
        write (11, '(A, A, f9.6)') tab, 'CO2(g) 2.0 ', co2_in
        write (11, *) 'SAVE solution 0'
        write (11, *) 'END'
        write (11, *) 'SOLUTION 1 Domain'
        write (11, *) tab, 'pH 7.0'!
        write (11, *) tab, 'units mol/L'
        write (11, *) tab, 'temp 60.000000'
        write (11, *) tab, 'pressure 98.6923 atm' ! = 100 bar
        write (11, *) tab, 'Cl 0.0100 charge'
        write (11, *) 'EQUILIBRIUM_PHASES 1'
        write (11, '(A, A, f9.6)') tab, 'Calcite 0.000000', calcite_in
        write (11, *) tab, 'Dolomite 0.000000 0.000000'
        ! write (11, *) tab, 'Quartz 0.000000 22.000000'
        write (11, *) 'SAVE solution 1'
        write (11, *) 'SELECTED_OUTPUT'
        write (11, *) tab, '-simulation false'
        write (11, *) tab, '-state false'
        write (11, *) tab, '-solution false'
        write (11, *) tab, '-distance false'
        write (11, *) tab, '-time false'
        write (11, *) tab, '-step false'
        write (11, *) tab, '-ph true' ! get pH output
        write (11, *) tab, '-pe false'
        write (11, *) tab, '-equilibrium_phases Calcite Dolomite'
            ! get calcite and dolomite concentrations
        write (11, *) 'END'
        close (unit=11, status='keep')
    end subroutine

    ! moves active particles via advection
    subroutine advect(p, v, dt, alive)
        type(mparticle),  intent(inout) :: p(:) ! mobile particle array
        double precision, intent(in   ) :: v(:), dt ! veclocity grid and time step
        integer,          intent(in   ) :: alive(:)
            ! array of indices of active particles

        ! use the velocity of the relevant cell, based on bin value
        p(alive)%loc = p(alive)%loc + v(p(alive)%bin) * dt
    end subroutine advect

    ! moves active particles via diffusion
    subroutine diffuse(p, np, D, dt, alive)
        type(mparticle),  intent(inout) :: p(:) ! mobile particle array
        double precision, intent(in   ) :: D(:), dt
            ! diffusion coefficient grid and time step
        integer,          intent(in   ) :: np, alive(:)
        ! number and array of indices of active particles
        double precision                :: normvec(np)
            ! vector which will hold Normal(0, 1) values

        ! call N(0, 1) generator
        call box_mullerp(np, normvec)

        ! use the diffusion coeff of the relevant cell, based on bin value
        p(alive)%loc = p(alive)%loc + sqrt(2.0d0 * D(p(alive)%bin) * dt) * normvec
    end subroutine diffuse

    ! relfective lower boundary
    subroutine reflectlow(p, low, alive)
        type(mparticle),  intent(inout) :: p(:) ! mobile particle array
        double precision, intent(in   ) :: low ! lower spatial boundary
        integer,          intent(in   ) :: alive(:)
            ! array of indices of active particles

        ! if particle has exited lower boundary, flip the negative location
        ! to the same positive value
        where (p(alive)%loc < low) p(alive)%loc = -p(alive)%loc
    end subroutine reflectlow

    ! absorbing upper boundary
    subroutine absorbhigh(p, high, alive)
        type(mparticle),  intent(inout) :: p(:) ! mobile particle array
        double precision, intent(in   ) :: high ! upper spatial boundary
        integer,          intent(in   ) :: alive(:)
            ! array of indices of active particles

        ! if particle has exited upper boundary, make the particle inactive
        ! also, make bin and loc -999 to catch any errors
        where (p(alive)%loc > high) p(alive)%active = .false.
        where (p(alive)%loc > high) p(alive)%bin = -999
        where (p(alive)%loc > high) p(alive)%loc = -999
    end subroutine absorbhigh

    ! since PHREEQCRM can't accept the particle array as input, this subroutine
    ! assigns the values in the 2D concs array to its corresponding immobile
    ! particle
    subroutine concs_to_iparts(c, p, np)
        double precision, intent(in   ) :: c(:, :)
            ! 2D concentration array used by PHREEQCRM
        type(iparticle),  intent(inout) :: p(:) ! immobile particle array
        integer,          intent(in   ) :: np ! number of immobile particles
        integer                         :: i ! iteration variable

        do i = 1, np
            p(i)%concs = c(i, :)
        enddo
    end subroutine concs_to_iparts

    ! this does the opposite of above
    subroutine iparts_to_concs(c, p, np)
        double precision, intent(inout) :: c(:, :)
            ! 2D concentration array used by PHREEQCRM
        type(iparticle),  intent(in   ) :: p(:) ! immobile particle array
        integer,          intent(in   ) :: np ! number of immobile particles
        integer                         :: i ! iteration variable

        do i = 1, np
            c(i, :) = p(i)%concs
        enddo
    end subroutine iparts_to_concs

    ! this subroutine does the probabilistic mass balancing according to the
    ! algorithm set forth in Benson and Bolster, "Arbitrarily complex reactions
    ! on particles," WRR 2016
    subroutine mass_balance(ip, mp, na, alive, ni, D, dt, omega)
        type(iparticle),  intent(inout) :: ip(:) ! immobile particle array
        type(mparticle),  intent(inout) :: mp(:) ! mobile particle array
        integer,          intent(in   ) :: na, alive(:), ni
            ! number and array of indices of active mobile particles and number
            ! of immobile particles
        double precision, intent(in   ) :: D(:), dt, omega
            ! diffusion coefficient array, time step, and domain length
        type(kdtree2), pointer          :: tree ! this is the KD tree
        integer                         :: ntot, dim = 1, aindex, bindex, i, j
            ! total number of particles (active mobile + immobile), number
            ! of spatial dimensions, index of 'B' particle for mass balance loop,
            ! and a couple of loop interators
            !****Note: hard coded one spatial dimension
        real(kdkind)                    :: locs(na + ni), r2
            ! array holding locations of immobile and active immobile particles
            ! and value of squared search radius for KD search
        type(index_array), allocatable  :: closeguys(:)
            ! this holds the indices of nearby particles
        type(dist_array), allocatable   :: close_dists(:)
            ! this holds the distances to the corresponding nearby particle
        double precision                :: ds, const(na), denom(na), v_s,&
                                           dmass(nspec), dmA(nspec), denom1,&
                                           const1
        ! ds, const, and demom are used in v(s) calculation but are precalculated for
        ! efficiency. v_s is encounter probability between particles
        ! (techinically, it is v(s)ds). dmass is the mass change for a particle pair
        ! dmA is the storage value for saving up the mass change in an A particle
        ! if the mass reduction is to be done AFTER the B paticle loop
        ! denom1 and const1 are scalar versions of denom and const

        ! calculate total number of particles to be considered for mass balance
        ntot = na + ni
        ! build locs array--immobile particles will be at the beginning of the
        ! array, and mobile will be at the end
        locs(1 : ni) = real(ip%loc, kdkind)
        locs(ni + 1 : ntot) = real(mp(alive)%loc, kdkind)
        ! calculate interaction distance to be numsd standard deviations of the
        ! Brownian Motion process--r2 is this distance squared
        ! ****NOTE: numsd is a global variable that is hard-coded above
        r2 = (real(numsd, kdkind) * sqrt(2.0_kdkind * maxval(real(D, kdkind)) *&
                                         real(dt, kdkind)))**2

        ! build the KD tree and search it
        ! ****NOTE: num_alloc is a global variable that is hard-coded above
        call maketree(tree, locs, dim, ntot)
        call search(ntot, dim, tree, r2, num_alloc, closeguys, close_dists)
        ! NOTE: this search returns the SQUARED distance between two points
            ! also, the point itself is included in the closeguys list
        call kdtree2_destroy(tree)

        ! do i = 1, ni
        !     print *, 'i = ', i, 'closeguys = ', closeguys(i)%indices
        !     do j = 1, size(closeguys(i)%indices)
        !         print *, 'closeguys(i)%indices(j) = ', closeguys(i)%indices(j)
        !         print *, 'loc A, loc B = ', locs(i), locs(closeguys(i)%indices(j))
        !         print *, 'distance = ', sqrt(close_dists(i)%dists(j))
        !         ! pause
        !     enddo
        ! enddo
        ! pause

        ! calculate ds, the 'support volume' of a particle
        ! ****should immobile particles be a part of this calculation?
        ds = omega / dble(ntot)
            ! ****NOTE: this is average inter-particle spacing of alive particles
        denom = -4.0d0 * D(mp(alive)%bin) * dt
            ! denom is part of the normalizing constant as well as the
            ! the denominator of the exponential argument
        const = ds / sqrt(pi * (-denom))
            ! normalizing constant for v(s), multiplied by ds
            ! this is for mobile/immobile interaction, since immobile aren't
            ! moving via diffusion (note there is only one D above, for the
            ! mobile particle, instead of two)
            ! NOTE: indices of denom and const correspond to indices of alive
            ! array--not the index in the mp array
            ! i.e., particle x which is alive(i) has denom(i) and const(i)

        ! loop over immobile particles
        do i = 1, ni ! immobile particle index--this is the 'A' particle loop
            ! if reducing mass of A particle after B loop set dmA to zero
            ! dmA = 0.0d0
            do j = 1, size(closeguys(i)%indices) ! mobile particle loop
                ! this is the 'B' particle loop
                ! note that the closeguys array is indexed to the loc array,
                ! and thus its true index in the mp array is calculated below

                ! current B particle's index in locs array
                if (closeguys(i)%indices(j) <= ni) cycle
                    ! if B is an immobile particle, skip this calculation
                    ! also prevents calculating with self
                ! make bindex equal to index in alive array so as to correspond
                ! to the const and denom arrays
                ! i.e., alive(bindex) = index in mp array, and const(bindex) =
                ! constant for particle alive(bindex)
                bindex = closeguys(i)%indices(j) - ni
                ! ****get rid of this when confident it's working****
                if (bindex > na .or. bindex < 1) then
                    print *, '****ERROR IN INDEX CALC****'
                    print *, 'bindex = ', bindex
                endif

                ! calculate encounter probability for A and B particle pair
                ! NOTE: distance is not squared since the search already
                ! returned the squared value
                v_s = const(bindex) * exp(close_dists(i)%dists(j)/denom(bindex))
                ! calculate change in mass for A and B particle pair
                dmass = 0.5d0 * (ip(i)%concs - mp(alive(bindex))%concs) * v_s
                ! change mass of A particle
                ! ****should it be changed here or add it up and do it after?****
                ip(i)%concs = ip(i)%concs - dmass
                ! save mass change of A particle until after B loop, so add up
                ! all of the dmasses
                ! dmA = dmA + dmass
                ! change mass of B particle
                mp(alive(bindex))%concs = mp(alive(bindex))%concs + dmass
            enddo
            ! subtract the saved up dmA from the A particle
            ! ip(i)%concs = ip(i)%concs - dmA
        enddo

        ! loop over mobile particles
        do i = ni + 1, ntot ! this is the 'A' particle loop
            ! make aindex equal to index in mp array
            aindex = alive(i - ni)
            ! if reducing mass of A particle after B loop set dmA to zero
            ! dmA = 0.0d0
            do j = 1, size(closeguys(i)%indices) ! this is the 'B' particle loop
                ! note that the closeguys array is indexed to the loc array,
                ! and thus its true index in the mp array is calculated below

                ! current B particle's index in locs array
                if (closeguys(i)%indices(j) <= i) cycle
                    ! want to ignore any indices lower than i to avoid immobile
                    ! particles, as well as double counting interactions

                ! make bindex equal to index in alive array
                ! i.e., alive(bindex) = index in mp array
                ! ****get rid of this when confident it's working****
                bindex = closeguys(i)%indices(j) - ni
                if (bindex > na .or. bindex < 1) then
                    print *, '****ERROR IN INDEX CALC****'
                    print *, 'bindex = ', bindex
                endif
                ! make bindex equal to index in mp array
                bindex = alive(closeguys(i)%indices(j) - ni)

                ! we calculate denom and const on the fly here since A and B
                ! do not necessarily have the same D value
                ! precaculating seems unnecessary, since that would require
                ! conisdering every mobile particle pair, and not every pair
                ! will interact
                denom1 = -4.0d0 * (D(mp(aindex)%bin) +&
                                  D(mp(bindex)%bin)) * dt
                const1 = ds / sqrt(pi * (-denom1))

                ! calculate encounter probability for A and B particle pair
                ! NOTE: distance is not squared since the search already
                ! returned the squared value
                v_s = const1 * exp(close_dists(i)%dists(j)/denom1)
                ! calculate change in mass for A and B particle pair
                dmass = 0.5d0 * (mp(aindex)%concs -&
                        mp(bindex)%concs) * v_s
                ! change mass of A particle immediately--this imposes a sort of
                ! ordering scheme
                ! ****should it be changed here or add it up and do it after?****
                mp(aindex)%concs = mp(aindex)%concs - dmass
                ! save mass change of A particle until after B loop, so add up
                ! all of the dmasses
                ! dmA = dmA + dmass
                ! change mass of B particle
                mp(bindex)%concs = mp(bindex)%concs + dmass
            enddo
            ! subtract the saved up dmA from the A particle
            ! mp(aindex)%concs = mp(aindex)%concs - dmA
        enddo
        deallocate (closeguys, close_dists)
    end subroutine mass_balance

    ! this builds a KD tree
    subroutine maketree(tree2, locs, d, n)
        type(kdtree2), pointer :: tree2 ! this is the KD tree
        real(kdkind)           :: locs(d, n)
            ! location array for particles, with dimension d x n (number of
            ! spatial dimensions x number of particles)
        integer                :: d, n
            ! number of spatial dimensions, number of particles

        ! build the tree
        tree2 => kdtree2_create(locs, sort=.false., rearrange=.false.)
            ! currently don't see a need to sort or rearrange
    end subroutine maketree

    ! this searches an already built KD tree
    subroutine search(n, d, tree, r2, num_alloc, closeguys, close_dists)
        integer,                        intent(in   ) :: n, d, num_alloc
            ! number of particles, number of spatial dimensions, how large to
            ! to preallocate results array within KD tree module
        type(kdtree2), pointer,         intent(in   ) :: tree ! the KD tree
        real(kdkind),                   intent(in   ) :: r2 ! squared search radius
        type(index_array), allocatable, intent(  out) :: closeguys(:)
            ! this holds the indices of nearby particles
        type(dist_array), allocatable,  intent(  out) :: close_dists(:)
            ! this holds the distances to the corresponding nearby particle
        integer                                       :: i, nf
            ! loop iterator and number of particles found by search
        type(kdtree2_result), allocatable             :: results(:)
            ! results array from KD tree module

        allocate (closeguys(n), close_dists(n), results(num_alloc))

        ! loop over all particles
        do i = 1, n
            ! the type of search used here finds all the particles within
            ! squared distance r2 from the i^th particle in the list
            ! the hard-coded 0 is the 'correlation time' of the search
            call kdtree2_r_nearest_around_point(tree, i, 0, r2, nf, num_alloc, results)

            ! allocate these based on how many nearby particles were found
            allocate (closeguys(i)%indices(nf), close_dists(i)%dists(nf))

            closeguys(i)%indices = results(1 : nf)%idx
            close_dists(i)%dists = results(1 : nf)%dis
        end do

        deallocate (results)
    end subroutine search

    ! subroutine bin_concs(p, np, alive, row, col, c)
    !     type(mparticle),   intent(inout) :: p(:)
    !     integer,          intent(in   ) :: np, alive(:), row, col
    !     ! double precision, intent(in   ) :: dx
    !     double precision, intent(  out) :: c(row, col)
    !     integer                         :: i

    !     c = 0.0d0

    !     ! print *, 'np = ', np, size(alive)
    !     ! pause

    !     do i = 1, np
    !         ! p(alive(i))%bin = floor(p(alive(i))%loc/dx) + 1
    !         ! print *, 'i = ', i
    !         ! print *, 'p(alive(i))%bin = ', p(alive(i))%bin
    !         if (p(alive(i))%bin < 1 .or. p(alive(i))%bin > 500) then
    !             print *, 'i = ', i
    !             print *, 'p(alive(i))%bin = ', p(alive(i))%bin
    !             print *, 'p(alive(i))%loc = ', p(alive(i))%loc
    !         endif
    !         c(p(alive(i))%bin, :) = c(p(alive(i))%bin, :) +&
    !                                 p(alive(i))%concs * (dble(np) / dble(row))
    !         ! adjust from particle to cell concs (note that row = ncellsol)
    !         ! **** also, dividing here by total alive particles, not sure if
    !         ! this is smart, given the higher conc near inflow boundary
    !     enddo
    ! end subroutine bin_concs

    ! subroutine unbin_concs(p, np, alive, row, col, c)
    !     type(mparticle),   intent(inout) :: p(:)
    !     integer,          intent(in   ) :: np, alive(:), row, col
    !     double precision, intent(in   ) :: c(row, col)
    !     double precision                :: temp(col)
    !     integer                         :: i, j, npart
    !     integer                         :: idx(np) ! made this overly large to avoid repeated allocation/deallocation
    !     logical                         :: inbin(np)

    !     do i = 1, row
    !         if (p(alive(i))%bin < 1 .or. p(alive(i))%bin > 500) then
    !             print *, 'i = ', i
    !             print *, 'p(alive(i))%bin = ', p(alive(i))%bin
    !             print *, 'p(alive(i))%loc = ', p(alive(i))%loc
    !         endif
    !         inbin = p(alive)%bin == i
    !         npart = count(inbin)
    !             ! ****distributing evenly--not sure if this is a good idea or not
    !         temp = c(i, :)/dble(npart)
    !         temp = temp * (dble(row) / dble(np))
    !             ! adjust from cell to particle concs (note that row = ncellsol)
    !         idx(1 : npart) = pack(alive, inbin)
    !         do j = 1, npart
    !             p(idx(j))%concs = temp
    !         enddo
    !     enddo
    ! end subroutine unbin_concs

    ! these next two subroutines use the Box-Muller transform to generate N(0,1)
    ! random numbers from U(0,1)
    ! https://goo.gl/DQgmMu
    ! Note: this polar formulation seems to be consistently ~20% faster than the
    ! version below that uses trig functions
    ! reference for polar version (and standard version):
    ! https://www.taygeta.com/random/gaussian.html
    subroutine box_mullerp(n, z)
        integer,          intent(in   ) :: n ! size of random vector to be generated
        double precision, intent(  out) :: z(n)
        integer                         :: i, j
        double precision                :: w, x1, x2

        call init_random_seed()
        i = 1

        do j = 1, n/2
            w = 1
            do while (w >= 1)
                x1 = 2.0d0 * rand() - 1.0d0
                x2 = 2.0d0 * rand() - 1.0d0
                w = x1**2 + x2**2
            enddo
            w = sqrt((-2.0d0 * log(w)) / 2.0d0)
            z(i : i + 1) = (/x1 * w, x2 * w/)
            i = i + 2
        enddo
    end subroutine box_mullerp

    ! subroutine box_muller()
    !     integer, parameter :: n = 1e8
    !     integer            :: i, j
    !     double precision   :: x1, x2, y1, y2, z(n)

    !     call init_random_seed()
    !     i = 1

    !     do j = 1, n/2
    !         x1 = rand()
    !         x2 = rand()
    !         y1 = sqrt(-2.0d0 * log(x1)) * cos(2.0d0 * pi * x2)
    !         y2 = sqrt(-2.0d0 * log(x1)) * sin(2.0d0 * pi * x2)
    !         z(i : i + 1) = (/y1, y2/)
    !         i = i + 2
    !     enddo
    ! end subroutine box_muller

end module RPT_mod
