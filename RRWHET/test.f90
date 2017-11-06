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
type(particle), allocatable ::pat(:)
type(imparticle), allocatable ::impat(:)
integer :: indys(10),numindys
logical :: blahmask
real :: pmass(3)


indys=(/ 0,0,0,0,0,0,0,0,0,0 /)

allocate(pat(3))
allocate(impat(2))

do i=1,3
pat(i)%pmass=[1, 2]
enddo
do i=1,2
impat(i)%pmass=[4, 5, 6]
enddo
pmass=impat(1)%pmass
numindys=count(pmass>4)

indys(1:numindys)=pack(pmass,pmass>4)

print*,impat, indys

!open(10,file='test.dat')
!read(10,*)i,pat(1)%pmass
!print*, pat

!blah=0.0
!print*,' blah',blah
!do i=1,3
!pat(i)%pmass=[1.0, 2.0]
!enddo

!print*,maxval(max(blah,pat(1)%pmass))

!blah=pat(1)%pmass
!print*,' blah',blah
!do iloop=1,nspec
! write(filename,'(a,i3.3)')'file',iloop
!print*,'1'//trim(filename)//'2'
!  write(*,'(a5,20e10.3)'),pat(iloop)
!enddo
!deallocate(pat)

stop 
end program main

