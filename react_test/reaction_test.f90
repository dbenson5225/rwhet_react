program reaction_test
    use react
    implicit none

integer, parameter                 :: n = 1e6, d = 3, corr = n/2
real(kdkind), parameter            :: r2 = 0.01**2
integer, parameter                 :: cgsize = ceiling(n*((4.0d0/3.0d0)*pi*r2))
real(kdkind), allocatable          :: my_array(:, :)
type (search_results), allocatable :: closeguys(:)
type(kdtree2), pointer             :: tree
integer                            :: correltime = 1 ! not quite sure what this does 
integer                            :: i
double precision                   :: start, finish_tree, finish

allocate (closeguys(n))
allocate (my_array(d, n))

call init_random_seed
call random_number(my_array)
call cpu_time(start)
call maketree(tree,my_array, d, n)
call cpu_time(finish_tree)

do i = 1, n
    call find_neighbs(tree, i, correltime, r2, cgsize, closeguys(i)%indices, closeguys(i)%dists)
end do

! do i = 1, n
!     write (*,*) 'i = ', i, 'closeguys = ', closeguys(i)%indices
!     write (*,*) 'distances = ', closeguys(i)%dists
! end do

call cpu_time(finish)

write (*,*) 'tree build time = ', finish_tree - start
write (*,*) 'total time = ', finish - start

call kdtree2_destroy(tree)
deallocate(my_array)

end program reaction_test
