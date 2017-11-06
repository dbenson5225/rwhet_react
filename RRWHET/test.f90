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
type(imparticle), allocatable ::impat(:)
integer :: indys(10),numindys
logical :: blahmask
real :: pmass(3)


indys=(/ 0,0,0,0,0,0,0,0,0,0 /)

allocate(pat(3))
allocate(impat(2))

do i=1,3
pat(i)%xyz=0.0
pat(i)%pmass=[1, 2]
enddo
do i=1,2
impat(i)%pmass=[4, 5, 6]
enddo
pmass=impat(1)%pmass
numindys=count(pmass>4)

indys(1:numindys)=pack(pmass,pmass>4)

!print*,impat, indys
print*, pat
call temp_print(pat)

stop 
end program main

subroutine temp_print(pat)
use globalblah
implicit none
type(particle), intent(in)  ::pat(3)
print*,pat
return 
end subroutine temp_print