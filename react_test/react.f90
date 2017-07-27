module react
    use kdtree2_precision_module
    use kdtree2_module
    implicit none

type index_array
    integer, allocatable :: indices(:)
end type
type dist_array
    double precision, allocatable :: dists(:)
end type
type result_array
    type(kdtree2_result), allocatable :: results(:)
end type
type treetree
    type(kdtree2), pointer :: tree
end type

double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

contains

subroutine find_neighbs(tree, target, cortime, rad2, nalloc, indices, dists)
    type(kdtree2), pointer,        intent(in   ) :: tree
    integer,                       intent(in   ) :: target, cortime, nalloc
    real(kdkind),                  intent(in   ) :: rad2
    integer, allocatable,          intent(  out) :: indices(:)
    double precision, allocatable, intent(  out) :: dists(:)
    
    type(kdtree2_result), target   :: results(nalloc)
    integer                        :: nfound

    ! allocate (results(nalloc))

    call kdtree2_r_nearest_around_point(tree, target, cortime, rad2,&
                                        nfound, nalloc, results)

    allocate (indices(nfound), dists(nfound))

    indices = results(1 : nfound)%idx
    dists = results(1 : nfound)%dis
end subroutine find_neighbs

subroutine maketree(tree2,data,d,n)
    integer :: n, d
    !   real(kdkind), dimension(:,:), allocatable :: data
    real(kdkind) data(d,n)
    type(kdtree2), pointer     :: tree2

    !    allocate(data(d,n))
    tree2 => kdtree2_create(data,sort=.false.,rearrange=.false.)  ! this is how you create a tree. 
end subroutine maketree

! subroutine to initialize the random number generator seed from clock time
subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count = clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine init_random_seed


end module react