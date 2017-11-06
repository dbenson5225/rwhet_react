!  Modified from Mike Schmidt by David Benson for inclusion into Reactive RWHet
!  For now, do not include interface to Phreeqc.

module RPT_mod
    use global                      ! From rwhet
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
    integer, parameter           :: sp = kind(1.0), dp = kind(1.d0)
    double precision, parameter  :: pi = 4.0d0 * atan(1.0d0)
    integer, parameter          :: numsd = 3
        ! numsd is distances over which particle reactions are considered
 
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

    ! this subroutine does the probabilistic mass balancing according to the
    ! algorithm set forth in Benson and Bolster, "Arbitrarily complex reactions
    ! on particles," WRR 2016
    subroutine mix_particles(imp, pat, cat, ddiff, dt, closeguys, close_dist,Dloc)
        type(particle),  intent(inout) :: pat(:) ! immobile particle array
        type(imparticle),intent(inout) :: imp(:) ! mobile particle array
        type(cell), intent(in)         :: cat(nx,ny,nz)   ! cell info (incl. diff)
        type(kdtree2), pointer         :: tree ! this is the KD tree
        double precision, intent(in   ):: ddiff(:), dt
        double precision               :: ds(:)
        integer                        :: na, alive(:), ni
            ! number and array of indices of active mobile particles and number
            ! of immobile particles
        integer                        :: ntot, dim = 3, aindex, bindex,&
                                           i, j, k, iloop, jloop
            ! total number of particles (active mobile + immobile), number
            ! of spatial dimensions, index of 'B' particle for mass balance loop,
            ! and a couple of loop iterators
            !****Note: hard coded one spatial dimension
        real(kdkind)                   :: locs(na + ni,dim), r2
            ! array holding locations of immobile and active immobile particles
            ! and value of squared search radius for KD search
        type(index_array), intent(inout), allocatable  :: closeguys(:)
            ! this holds the indices of nearby particles
        type(dist_array), intent(inout), allocatable   :: close_dists(:)
            ! this holds the distances to the corresponding nearby particle
!        double precision              :: const(na), denom(na), v_s,&
!                                           dmass(nspec), dmA(nspec), denom1,&
!                                           const1, D(na)
        double precision, allocatable  :: const(:),denom(:),dmass(:),dmA(:),Dloc(:)
        double precision               :: v_s,denom1,const1
 
        logical::im_transfer
        im_transfer=.false.   ! For now, mixing only among mobile particles

        ! ds, const, and demom are used in v(s) calculation but are precalculated for
        ! efficiency. v_s is encounter probability between particles
        ! (technically, it is v(s)ds). dmass is the mass change for a particle pair
        ! dmA is the storage value for saving up the mass change in an A particle
        ! if the mass reduction is to be done AFTER the B particle loop
        ! denom1 and const1 are scalar versions of denom and const

        ! calculate total number of particles to be considered for mass balance
    indices = (/(iloop, iloop = 1, na)/) ! this is for easy array-indexing of pat
    na = count(pat%active)
    ni = count(imp%active)
    allocate (alive(na)) ! maybe preallocate to avoid repeatedly doing this
    allocate (Dloc(na),denom(na),const(na))
    alive = pack(indices, pat%active)
    ntot = na+ni
        ! build locs array--mobile particles will be at the beginning of the
        ! array, and immobile will be at the end
    locs(1 : na) = real(pat(alive)%xyz, kdkind)
    locs(ni + 1 : ntot) = real(imp(alive)%xyz, kdkind)
    ds(1:na)=real(pat(alive)%ds, kdkind)
 
!     if(im_transfer) then
!         ntot=na+ni
!     endif

     allocate(const(ntot),denom(ntot),D(ntot))
     Dloc=0.0D0
!  This seems super slow - I should find a better place to get each particle's Diff?
!  DAB Remember to subtract this Dloc from the random walk part !!!!
!  Dloc should contain Dm and possibly dtran*|v|

   do iloop=1,na
          aindex=alive(iloop)
          i=pat(aindex)%ijk(1); j=pat(aindex)%ijk(2); k=pat(aindex)%ijk(3) 
          Dloc(iloop)=ddiff(cat(i,j,k)%zone)
   enddo
!
        ! calculate interaction distance to be numsd standard deviations of the
        ! Brownian Motion process--r2 is this distance squared
        ! ****NOTE: numsd is a global variable that is hard-coded above

        r2 = (real(numsd, kdkind) * sqrt(2.0_kdkind * maxval(real(Dloc, kdkind)) *&
                                         real(dt, kdkind)))**2

        ! build the KD tree and search it
        ! ****NOTE: num_alloc is a global variable that is hard-coded above

        call maketree(tree, locs, dim, ntot)
        call search(ntot, dim, tree, r2, num_alloc, closeguys, close_dists)

        ! NOTE: this search returns the SQUARED distance between two points
            ! also, the point itself is included in the closeguys list

        call kdtree2_destroy(tree)

        ! calculate ds, the 'support volume' of a particle
        ! ****should immobile particles be a part of this calculation?
!        ds = omega / dble(ntot)
            ! ****NOTE: this is average inter-particle spacing of alive particles
        denom = -4.0d0 * Dloc * dt    ! This should be a vector na big
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

        ! loop over mobile particles
    allocate(dmass(nspec))

        do iloop = 1,na ! this is the 'A' particle loop
            ! make aindex equal to index in mobile pat array
            aindex = alive(iloop)
            ! if reducing mass of A particle after B loop set dmA to zero
            ! dmA = 0.0d0
            do jloop = 1, size(closeguys(iloop)%indices) ! this is the 'B' particle loop
                ! note that the closeguys array is indexed to the loc array,
                ! and thus its true index in the pat array is calculated below

                ! current B particle's index in locs array
                if (closeguys(iloop)%indices(jloop) <= iloop .or. &
                    closeguys(iloop)%indices(jloop) > na) cycle
            ! want to ignore any indices lower than iloop to avoid double 
            ! counting interactions as well as >na, which are immobile particles 
           
           bindex = alive(closeguys(iloop)%indices(jloop))  ! DAB .. right?
                if (bindex > na .or. bindex < 1) then
                    print *, '****ERROR IN INDEX CALC****'
                    print *, 'bindex = ', bindex
                endif
                ! make bindex equal to index in pat array

 
                ! we calculate denom and const on the fly here since A and B
                ! do not necessarily have the same D value
                ! precalculating seems unnecessary, since that would require
                ! considering every mobile particle pair, and not every pair
                ! will interact
                denom1 = -4.0d0 * (D(iloop)+D(closeguys(iloop)%indices(jloop))) * dt  
                const1 = ds(iloop) / sqrt(pi * (-denom1))

                ! calculate encounter probability for A and B particle pair
                ! NOTE: distance is not squared since the search already
                ! returned the squared value
                v_s = const1 * exp(close_dists(iloop)%dists(jloop)/denom1)
                ! calculate change in mass for A and B particle pair
                dmass = 0.5d0 * (pat(aindex)%pmass - pat(bindex)%pmass) * v_s
                ! change mass of A particle immediately--this imposes a sort of
                ! ordering scheme
                ! ****should it be changed here or add it up and do it after?****
                pat(aindex)%pmass = pat(aindex)%pmass - dmass
                ! change mass of B particle
                pat(bindex)%pmass = pat(bindex)%pmass + dmass
            enddo
         enddo

!   Loop over immobile particles

! if(im_transfer) then  ! DAB I have not yet fixed this
! print*,' DAB Immobile mass transfer not yet done correctly!!!'

!        allocate(dmass(inspec))
!        do iloop = ni+1,ntot ! immobile particle index--this is the 'A' particle loop
!            ! if reducing mass of A particle after B loop set dmA to zero
!
!            do jloop = 1, size(closeguys(i)%indices) ! mobile particle loop
!                ! this is the 'B' particle loop
!                ! note that the closeguys array is indexed to the loc array,
!                ! and thus its true index in the mp array is calculated below
!
!                ! current B particle's index in locs array
!                if (closeguys(iloop)%indices(jloop) <= ni) cycle
!                    ! if B is an immobile particle, skip this calculation
!                    ! also prevents calculating with self
!                ! make bindex equal to index in alive array so as to correspond
!                ! to the const and denom arrays
!                ! i.e., alive(bindex) = index in mp array, and const(bindex) =
!                ! constant for particle alive(bindex)
!               bindex = closeguys(i)%indices(j) - ni
!                ! ****get rid of this when confident it's working****
!                if (bindex > na .or. bindex < 1) then
!                    print *, '****ERROR IN INDEX CALC****'
!                    print *, 'bindex = ', bindex
!                endif
!
!                ! calculate encounter probability for A and B particle pair
!                ! NOTE: distance is not squared since the search already
!                ! returned the squared value
!                v_s = const(bindex) * exp(close_dists(i)%dists(j)/denom(bindex))
!                ! calculate change in mass for A and B particle pair
!                dmass = 0.5d0 * (ip(i)%concs - mp(alive(bindex))%concs) * v_s
!                ! change mass of A particle
!                ! ****should it be changed here or add it up and do it after?****
!                ip(i)%concs = ip(i)%concs - dmass
!                ! save mass change of A particle until after B loop, so add up
!                ! all of the dmasses
!                ! dmA = dmA + dmass
!                ! change mass of B particle
!                mp(alive(bindex))%concs = mp(alive(bindex))%concs + dmass
!            enddo
!            ! subtract the saved up dmA from the A particle
!            ! ip(i)%concs = ip(i)%concs - dmA
!        enddo
!   deallocate(dmA)
!endif   ! Mass transfer with immobile particles?

    end subroutine mix_particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine abc_react(imp,pat,closeguys,close_dists,Dloc,dt)

!  This subroutine will have reaction mA (aq) + nB (aq) <--> pC (solid) for debugging
!  The forward rxn has hard-coded rate kf, backwards kr
!  Any solid made will be put on nearby imm particles. If none near, create one
!  If solid dissolves, place on nearby pats.  If none near, create one.

       type(particle),  intent(inout) :: pat(:) ! immobile particle array
       type(imparticle),intent(inout) :: imp(:) ! mobile particle array
       type(cell),      intent(in)    :: cat(nx,ny,nz)   ! cell info (incl. diff)
       type(index_array), intent(in)  :: closeguys(:)
            ! this holds the indices of nearby particles
       type(dist_array), intent(in)   :: close_dists(:)
            ! this holds the distances to the corresponding nearby particle
       double precision, intent(in)   :: Dloc(:)
       integer:: iloop,numjs
       double precision:: kr,kf,stoA,stoB,stoC,conc(nspec),iconc(inspec),ds(:),pmass()
       double precision:: dm, dmA, dmB, dmC(:), dmAtemp, dmBtemp, v_s(:), totwgt

!  Define stoichiometry here
stoA=1.d0
stoB=1.d0
stoC=1.d0
allocate(pmass(inspec))

!  ds=pat%ds
na = count(pat%active)
ni= count(imp%active)
allocate v_s(na+ni)  ! Maximum possible size
v_s=0.d0
pmass=0.d0 
do iloop=1,na   !  These are the mobiles

     conc=(pat(iloop)%pmass)/pat(iloop)%ds   
     dm = kf*dt*(conc(1)/ds)**msto*(conc(2)/ds)**nsto

     dmAtemp=min(stoA*dm,pat(iloop)%pmass(1))
     dmBtemp=min(stoB*dm,pat(iloop)%pmass(2))
     dm=min(dmAtemp/msto,dmBtemp/nsto)
    
     pat(iloop)%pmass(1)=pat(iloop)%pmass(1)-msto*dm
     pat(iloop)%pmass(2)=pat(iloop)%pmass(2)-nsto*dm
     dmC(iloop)=stoC*dm             !  need to put this on a solid particle
!  If I do it here I lose future paralellism, so beware
!  Find out how many immobile particle are nearby
!  indys holds the indexes for all particles (1:na then na+1:not)
!  The index in the imp array is the indy - na

     numjs=size(closeguys(iloop)%indices>na)
      if(numjs.lt.1) then
! make an imp here
         i=pat(iloop)%ijk(1); j=pat(iloop)%ijk(2); k=pat(iloop)%ijk(3) 
         x=pat(iloop)%xyz(1); y=pat(iloop)%xyz(2); z=pat(iloop)%xyz(1)
         pmass(1)=dmC(iloop)
         call addimp(single(curtime),x,y,z,i,j,k,pmass,imp,cat) 
     else
! Place solids on nearby imp particles weighted by prob. of collision
        indys(1:numjs)=pack(closeguys(iloop)%indices,closeguys(iloop)%indices>na)     
        v_s=0.d0

        do jloop = 1, numjs

                denom1 = -4.0d0 * (D(iloop)) * dt  

                ! calculate encounter probability for A and B particle pair
                ! NOTE: distance is not squared since the search already
                ! returned the squared value

                v_s(jloop) = exp(close_dists(iloop)%dists(indys(jloop))/denom1)

        enddo
        totwgt=sum(v_s(1:numjs))
        do jloop=1, numjs 
                weight=v_s(jloop)/totwgt
                bindex = closeguys(iloop)%indices(indys(jloop)) - na 
                imp(bindex)%pmass(1)=weight*stoC*dmC(iloop)
        enddo
enddo  !!!!  end of loop through mobile reactions and precipitation

!  Do the reverse reaction

do iloop=1,ni    !  This will stay parallel

!  Should enforce that kr*dt < 0.1, but fuck it
  if (kr*dt>0.1) print*,' Timesteps probably too big.  kr*dt = ',kr*dt

     dm = min(imp(iloop)%pmass(1),kr*dt*imp(iloop)%pmass(1))

     dmA(iloop)=(stoA/stoC)*dm        ! Store for later delivery to mobile particles
     dmB(iloop)=(stoB/stoC)*dm    
     imp(iloop)%pmass(1)=imp(iloop)%pmass(1)-dm
enddo

do iloop=1,ni

! find close mobile particles to give newly dissolved mass to
!  Remember that the immobile particles site after ni mobiles in the closeguys array

    numjs=size(closeguys(ni+iloop)%indices<=na)
    indys(1:numjs)=pack(closeguys(ni+iloop)%indices,closeguys(ni+iloop)%indices<=na)     
    if(numjs.lt.1) then
! make a new pat here
         i=imp(iloop)%ijk(1); j=imp(iloop)%ijk(2); k=imp(iloop)%ijk(3) 
         x=imp(iloop)%xyz(1); y=imp(iloop)%xyz(2); z=imp(iloop)%xyz(1)
         pmass(1)=dmA(iloop); pmass(2)=dmB(iloop)    ! etc. etc. if needed
         call addp(single(curtime),x,y,z,i,j,k,pmass,pat,cat) 
     else
! Place dissolved solids on nearby pat particles weighted by prob. of collision
        v_s=0.d0

        do jloop = 1, numjs

                denom1 = -4.0d0 * (D(indys(jloop))) * dt  

                ! calculate encounter probability for A and B particle pair
                ! NOTE: distance is not squared since the search already
                ! returned the squared value

                v_s(jloop) = exp(close_dists(ni+iloop)%dists(indys(jloop))/denom1)

        enddo
        totwgt=sum(v_s(1:numjs))
        do jloop=1, numjs 
                weight=v_s(jloop)/totwgt
                bindex = closeguys(ni+iloop)%indices(indys(jloop)) 
                pat(bindex)%pmass=weight*stoC*pmass
        enddo

enddo   !  loop through the immobile for dissolution

deallocate (closeguys, close_dists,pmass,v_s)

end subroutine abc_react
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
