module globalblah
!integer::nspec
!open(30,file='nspec.dat')
!read(30,*)nspec
integer, parameter:: nspec=2
integer::iloop

type particle
  real::xyz(3)
  real::pmass(2)
end type particle
type imparticle
  real::xyz(3)
  integer::pmass(3)
  logical::active
end type imparticle


end module globalblah
!-----------------------
program main
use globalblah

implicit none
integer::i, size1, size2
real::blah(nspec)
character(80)::filename
type(particle), allocatable   ::pat(:)
type(imparticle), allocatable ::impat(:,:)
integer :: indys(10),numindys
integer,allocatable:: alive(:)
logical :: blahmask
real :: pmass(3)


indys=(/(i, i = 1, 10)/)

allocate(pat(3))
allocate(impat(1,3))

do i=1,3
pat(i)%xyz=0.0
pat(i)%pmass=[1, 2]
enddo
do i=1,3
impat(1,i)%pmass=[i+1, i+2, i+3]
enddo
impat%active=[.false., .true.,.true.]

numindys=count(impat%active)
allocate(alive(numindys))
print*, count(impat%active)
!print*, impat%active

alive=pack(indys,impat%active)
!print*,alive

print*,impat, indys, alive
!print*, pat
call temp_print(pat)

pmass=impat(1,1)%pmass

stop 
end program main

subroutine temp_print(pat)
use globalblah
implicit none
type(particle), intent(in)  ::pat(3)
print*,pat
return 
end subroutine temp_print