program reaction_test
    use react
    implicit none

integer, parameter              :: n = 1e4, d = 3, corr = n/2
real(kdkind), parameter         :: r2 = 0.0001**2
integer, parameter              :: cgsize = ceiling(n*((4.0d0/3.0d0)*pi*r2))
real(kdkind), allocatable       :: my_array(:, :)
type (index_array), allocatable :: closeguys(:)
type (dist_array), allocatable  :: close_dists(:)
! type(kdtree2), pointer          :: tree
type(treetree) :: bigtree(n)
integer                         :: nn      ! number of neighbors 
integer                         :: i
integer                         :: nf
double precision                :: start, finish_tree, finish, omp_get_wtime,&
                                   ostart, ofinish
! type(result_array) :: resultarr(n)
type(kdtree2_result), allocatable :: results(:)

allocate (closeguys(n), close_dists(n), results(cgsize))
allocate (my_array(d,n))

call init_random_seed
call random_number(my_array)
call cpu_time(start)
!~ ostart = omp_get_wtime()
! call maketree(tree,my_array,d,n)
call cpu_time(finish_tree)

!$omp parallel
!$omp do
do i = 1, n
    call maketree(bigtree(i)%tree,my_array,d,n)
    call find_neighbs(bigtree(i)%tree, i, 1, r2, cgsize, closeguys(i)%indices, close_dists(i)%dists)
    print *, 'i = ', i
    call kdtree2_destroy(bigtree(i)%tree)
end do
!$omp end do
!$omp end parallel

! !$omp parallel
! !$omp do
! do i = 1, n
!     ! allocate (resultarr(i)%results(n))
!     call kdtree2_r_nearest_around_point(tree,i,1,r2,&
!         nf,cgsize,results)

!     allocate (closeguys(i)%indices(nf), close_dists(i)%dists(nf))

!     closeguys(i)%indices = results(1:nf)%idx
!     close_dists(i)%dists = results(1:nf)%dis
! end do
! !$omp end do
! !$omp end parallel

! do i = 1, n
!     write (*,*) 'i = ', i, 'closeguys = ', closeguys(i)%indices
!     write (*,*) 'distances = ', close_dists(i)%dists
! end do

call cpu_time(finish)

write (*,*) 'tree build time = ', finish_tree - start
write (*,*) 'total time = ', finish - start

! call kdtree2_destroy(tree)
deallocate(my_array)

end program reaction_test