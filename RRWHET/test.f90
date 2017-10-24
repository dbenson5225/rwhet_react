module globalblah
!integer::nspec
!open(30,file='nspec.dat')
!read(30,*)nspec
integer, parameter:: nspec=2
integer::iloop

type particle
  real::xyz(3)
  real::pmass(nspec)
end type particle
end module globalblah
!-----------------------
program main
use globalblah

implicit none
integer::i
real::blah(nspec)
character(80)::filename

type(particle) pat(3)

open(10,file='test.dat')
read(10,*)i,pat(1)%pmass
print*, pat

blah=0.0
print*,' blah',blah
do i=1,3
pat(i)%pmass=[1.0, 2.0]
enddo

print*,maxval(max(blah,pat(1)%pmass))

blah=pat(1)%pmass
print*,' blah',blah
do iloop=1,nspec
 write(filename,'(a,i3.3)')'file',iloop
print*,'1'//trim(filename)//'2'
!  write(*,'(a5,20e10.3)'),pat(iloop)
enddo
!deallocate(pat)

stop 
end program main

