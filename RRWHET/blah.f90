!   Reactive RWHET.  Original RWHET modified 10/2017 by David Benson to allow chemical reactions.
!
! This code solves arbitrary reactions and transport by particles.
! There are two logical ways to add reactions.  The first, in which each particle is composed
! of one component and changing mass (Bolster et al. WRR, 2016) is best handled by taking the 
! particle array pat(:,:) of species n=1,nspec as species 1 = pat(1,np), species 2 = pat(2,np) ...
! 
! This program follows the algorithm of Benson et al. (WRR 2016), in which each particle carries mass 
! of any or all species.  In this case, it is best handled by changing each particle's attributes,
! specifically the particle mass now has nspec components: particle%pmass is now 
! particle%pmass(nspec).  This requires changing all things that track mass or concentration 
! to arrays of size nspec, e.g., sample%mass(nspec).
!
! RWHET4.1
! AUTHOR: ERIC M. LABOLLE
!
! RWHet was developed for the benefit of the hydrologic sciences. 
! RWHet is officially unsupported by the author, although he welcomes 
! questions and will generally seek to answer them provided that you've 
! read this document and he has the time.  As with any scientific software, 
! the user should thoroughly understand its applicability, limitations and 
! foundational theory before using it. Use of RWHet is at your own risk; 
! the author cannot guarantee RWHet's accuracy for your application. 
! The author assumes no liability for the consequences of using RWHet.  
! RWHet is the sole ownership of the author and RWHet or any of the code 
! therein may not be distributed without this document and the author's permission.
! Email: emlabolle@ucdavis.edu
!
! Version: 2.0 --> 2.1
!
! ploti		- added statement to return if no active boundaries, nbounds=0
!			- previously the code would fail if there were no boundaries
! main		- removed call to output at end of program
!
! Version from 2.1 --> 2.2
!   Handle type 1 constant c boundaries with mass injection rate on cell edge, instead of refresh
!   Minor bug fixes in int output
! Version 2.2 --> 3.0 
!   Add a unique particle number pnumber to pat
!   Add pnumber as an incremental particle number counter
!   Add birth_place to pat
!   Add borehole flow and transport support
!   Add new int output for borehole flow and transport + easy to calculate concentration
! Version 3.2
!   Add new output for sampling at MW
! Version 3.2.1
!   Fix the problem with boundary-type 7 input
! Version 3.2.2
!   Subroutinized the velupdt routine
!   Read recharge from cbc flow file and adjust water table boundary velocities accordingly.
!   Store recharge nodes for recharge BC type.
!   Add recharge BC type 8 and 9
!   Fixed problem with tcntrl that could slow down the solution significantly in some cases
! Version 3.2.3
!   Read CHD and GHB fluxes from cbc flow file.
!   Add recharge BC type 10 and 11 for CHD and GHB
!   Add a line of output identifying the version used.
! Version 3.2.4
!   Fixed recharge and CHD boundaries for the case where 
!   irevz=1
! Version 3.2.5
!   Bug fix: condition dtmin< near still results in particles getting stuck
!   change to displacement < near
! Version 3.2.5
!   Bug fix: condition dtmin< near still results in particles getting stuck
!   change to displacement < near
! Version 3.2.6
!   Switched back to the double precision reflect to avoid stranded particles.
!   Added capability for ASCII TSIM output.
!	Fully include retardation in all aspects of the code.
!	Changed meaning of dtcntrl, now scales linearly with step size; use 1.0 to limit to cell
! Version 4.1
module global
! DAB hardwire number of mobile species (nspec) and immobile species(inspec) right now.
  integer, parameter:: nspec=2,inspec=2
  type sample
    real:: xyz(4),radius,vol                    ! x,y,zbot,ztop,radius of sample location
    real:: conc(nspec),mass(nspec)              ! concentration for all time
    integer*2:: ijk(4),itype,nzone,zone(100)    ! cell location and type (1 = x,y,z, 0 = i,j,k) of sample location
  end type sample            
  type particle
    real:: xyz(3),birth_day,death_day !,birth_place(3)  ! location,birth time
    double precision:: pmass(nspec)             ! mass of each species
    integer*2:: ijk(3)                      ! cell location
    integer:: pnumber                       ! cell location
    logical*2:: active                        ! active particle?
  end type particle

  type imparticle
    real:: xyz(3),birth_day,death_day !,birth_place(3)  ! location,birth time
    double precision:: pmass(inspec)             ! mass of each immobile species
    integer*2:: ijk(3)                      ! cell location
    integer:: pnumber                       ! cell location
    logical*2:: active                        ! active particle?
  end type imparticle


! define a cell
  type cell
    real::             tc                   ! time step control 
    double precision:: cmass(nspec),mass_remove(nspec)    ! mass in cell, mass lost since last report
    double precision:: cimmass(inspec)      ! Immobile mass in cell
    integer::          np_cell,np_remove    ! # of particles per cell, # lost since last report 
    integer::          nimp_cell,nimp_remove   ! # of IMMOBILE particles per cell, # lost since last report 
    integer::          bc_number,zone       ! bc condition #, zone number                
  end type cell            
! boundary  
  type boundary
    integer:: bc_type,ijk(3),kbot,ktop,np,group       ! boundary type,location,group number
    integer:: np_remove                     ! # of particles per cell, # lost since last report 
    integer:: ijkm(3),ijkp(3)               ! ijk loc. of beginning and end of source
    real:: conc(nspec),refresh,nptime       ! aqueous phase concentrations,refresh rate control
    real:: flux,fluxin,fluxout              ! flux in cell, total flux into MAW, total flux out of MAW
    real:: bc_xyzm(3),bc_xyzs(3)            ! xyz loc. of beginning of source type 2, size of source
    double precision:: pmass(nspec)         ! particle mass for all species
    double precision:: massrate(nspec)      ! rate of injection M/T for each species
    double precision:: mass_remove(nspec)   ! mass removed
    double precision:: maw(nspec)           
 !   real,allocatable:: flux_ijk(:)         ! flux for ijk locations bounded by ijkm and ijkp
    real:: tbeg,tend,dt_type1               ! beginning time, ending time
  end type boundary
  type recharge
    real:: flow                             ! recharge flux 
    integer:: k                             ! layer
  end type recharge            
  type constant_head
    real:: flow                             ! constant head flux 
  end type constant_head            

!
! nx,ny,nz		-	problem dimensions
! mx,my,mz		-	1 for active dimension, 0 otherwise
! nxyz		    -	nx*ny*nz
! np			    -   total number of particles
! npack           -   number of packages
! noc             -   number of output control packages (options)

  integer, parameter:: npack=8,nopc=7,nbtype=11
  real, parameter:: large=1.0e+10,small=1.0e-4,bcdxyz=0.01
  character (len=80), parameter:: version=' ------------ R e a c t i v e  R W H e t   V e r s i o n  1.0 ------------------'
!

! xyzl2(3)          - half the length of the domain
! xyzc(3)           - the center of the domain
! xyzmax(3)         - the maximum extant of the domain
! xyzmin(3)         - the minimum extant of the domain
! mass              - total mass in the domain 
! netmassp(3)  - mass passing through maximum extents of domain
! netmassm(3)  - mass passing through minimum extents of domain 
! pnumber      - incremental particle number 
  integer:: nx,ny,nz,nxy,mx,my,mz,nxyz,np,maxnp,pnumber,nline
  integer:: nimp,impnumber
  integer:: nbounds,maxbnd,maxsource,nsource,npntsrc,nchd
  integer:: netxyzp(3),netxyzm(3),netsplit,ixyzp(3),ixyzm(3)
  integer:: ixyzmax(3),ixyzmin(3),npbc(nbtype+1)
  integer:: ibug
  integer:: inbas,iout,ioutp,incat,invlc,inopc,inpar
  integer:: inpnt,inbnd,invel,inbgr,insam,ibugout
  integer:: isotropic,iadvect,interp,nzone
  integer:: irevz,istream,nres
  integer:: iseed1,iseed2,nran 
  integer:: iloop

  integer:: kstp,kper ! MODFLOW step and stress period from vlc input 

  integer:: inobgr !(if 1 then we read parameter fields)

  integer:: nsam

  real:: dx,dy,dz,dx2,dy2,dz2,ax,ay,az,vol
  real:: dtinit,tinit,dtminimum,dtcntrl,bke
  real:: rmx,rmy,rmz 
  real:: smallxyz(3),dxyz(3),near


  double precision:: xyzl2(3),xyzmax(3),xyzmin(3),xyzc(3)
  double precision:: mass(nspec),netmass(nspec),netmassp(3,nspec),netmassm(3,nspec)
  double precision:: massbc(nbtype+1,nspec)
  double precision:: cres
  double precision:: rm 
  double precision:: curtime,tnextopc,tnextbnd,tnextvel
  double precision:: tnextpnt,tnextsrc,tmax,dt

  character (len=80) fnamevel,fnameopc,fnamepar,fnamesam
  character (len=80) fnamepnt,fnamebnd,fnamebas,fnameout,fnameoutp,string(5)

end module global
!------------------------------------------------------------------------------------------
PROGRAM RWHET
! comment this line to compile on non-windows 
!USE DFPORT
use global
implicit none
!
! dt                                  - current time step
! curtime                             - current time
! tinit                               - initial time 
! tmax                                - maximum time 
! dtvel,dtbnd,dtopc,dtmax             - difference between current time
!                                       and time to update attribute bnd, vel, etc. 
! tnextbnd,tnextopc,tnextvel          - time of next change in BC's, OPC and VEL
! source(maxsource)                   - pointer to bounds for current constant concentration sources
!                                       bounds(i,ibounds), nsource = number of current SOURCE's
! maxsource                           - maximum number of allowed SOURCE's at any one time

! opc(nopc)                           - ouput control opc(iopc) = 1 print output 
! outunit(nopc+1)                     - unit number for output of attribute iopc 
! outfname(nopc)                      - file name for       "                                  "
!
! netxyzm,netxyzp ...                 - breakthrough counters at extents of model
! mass                                - total mass in the system
! netmassm,netmassp                   - mass breakthrough counters at extents of model
! cres                                - minimum concentration for particle splitting
! nres                                - particle resolution at cres
!
type (sample), allocatable:: sam(:)
type (particle), allocatable:: pat(:,:)
type (imparticle), allocatable:: imp(:,:)
type (cell), allocatable::     cat(:,:,:)
type (boundary), allocatable:: bounds(:)
type (recharge), allocatable:: rech(:,:)
type (constant_head),allocatable:: chd(:,:,:)
!.....DARCY VELOCITY
real,allocatable:: vel3(:,:,:,:)
!.....BOUNDARY CONDITIONS AND SOURCE
integer,allocatable:: source(:)
!.....PARAMETERS
real,allocatable:: por(:), ret(:), dlong(:), dtran(:), ddiff(:),decay(:)
!.....OUTPUT CONTROL (OPC)
integer,allocatable:: opc(:), outunit(:)
character (len=80),allocatable:: outfname(:)
!
integer:: iread_error,ierror,idir,inodes,iicat,istat,itype
integer:: iplotm,iplotc,iplotd,iplotb,iplotp,iploti,iplotr,iplotsc,iplotsm,icnf,imt3d
logical:: flexist
!
real:: tarray(2)
!
external movep1,movep2,movep3,movep4
!
!.....output and input unit numbers
!
!istat=dtime(tarray)
!
ioutp=8
ibugout=9
iout=10
inbas=11
inopc=12      
inpar=13           
invlc=14           
inbnd=15          
invel=16
inpnt=17
inbgr=18
insam=19
!
iplotm=21
iplotc=22
iplotd=23
iplotb=24
iplotp=25
iploti=26
iplotr=27
iplotsc=28
iplotsm=29
icnf=30
imt3d=31
!
!.....OPEN MAIN INPUT AND OUTPUT FILES
!
flexist=.false.
do; if(flexist)exit
  print*,' Main Input File Name?'
  read(*,'(a)',err=9993)fnamebas
  inquire (file=fnamebas, exist=flexist)
  if(.not.flexist)print*,' Main input file does not exist'
enddo
open(inbas,file=fnamebas,status='old')
!
do
  print*,' Output File Name?'
  read(*,'(a)',err=9993)fnameout
! do not overwrite main input file
  if(fnameout.eq.fnamebas)then
    print*,' illegal name'
    cycle
  else
    open(iout,file=fnameout,status='unknown')
    exit
  endif
enddo
write(iout,'(a)')'--------------------------------------------------------------------------------'
write(iout,'(a80)')version
write(iout,'(a)')'--------------------------------------------------------------------------------'
!
do
  print*,' Particle Removal Output File Name?'
  read(*,'(a)',err=9993)fnameoutp
! do not overwrite main input file
  if(fnameoutp.eq.fnamebas.or.fnameoutp.eq.fnamebas)then
    print*,' illegal name'
    cycle
  else
!~     open(ioutp,file=fnameoutp,status='unknown',form='binary')
    ! open(ioutp,file=fnameoutp,status='unknown',form='unformatted',access='stream')
    open(ioutp,file=fnameoutp,status='unknown',form='unformatted')
    exit
  endif
enddo
!
!.....BASIC MODEL INPUT
!
call maininput()
string=' ';nline=1
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR VELOCITY ARRAY!';call debug(ibugout,string,nline);endif
allocate (vel3(1:3,0:nx,0:ny,0:nz), STAT = ierror)
if(ierror.ne.0)goto 9994
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR RECHARGE ARRAY!';call debug(ibugout,string,nline);endif
allocate (rech(nx,ny), STAT = ierror)
if(ierror.ne.0)goto 9994
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR CHD ARRAY!';call debug(ibugout,string,nline);endif
allocate (chd(nx,ny,nz), STAT = ierror)
if(ierror.ne.0)goto 9994
!
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR PAT ARRAY!';call debug(ibugout,string,nline);endif
allocate (pat(1,maxnp), STAT = ierror)
if(ierror.ne.0)goto 9994
!
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR CAT ARRAY!';call debug(ibugout,string,nline);endif
allocate (cat(nx,ny,nz), STAT = ierror)
if(ierror.ne.0)goto 9994
!
!.....INITIALIZE SOME VALUES
!
if(ibug.ge.1)then;string(1)=' INITIALIZING VARIABLES!';call debug(ibugout,string,nline);endif
cat(:,:,:)%tc=0.0
do iloop=1,nspec
  cat(:,:,:)%cmass(iloop)=0.0d+0
  cat(:,:,:)%mass_remove(iloop)=0.0d+0
enddo
cat(:,:,:)%np_cell=0
cat(:,:,:)%np_remove=0
cat(:,:,:)%bc_number=0
cat(:,:,:)%zone=0
!
vel3(:,:,:,:)=0.0
np=0
!.model extents
dx2=dx/2.0; dy2=dy/2.0; dz2=dz/2.0
!
dxyz(1)=dx; dxyz(2)=dy; dxyz(3)=dz
xyzmin(:)=dble((ixyzmin(:)-1)*dxyz(:))
xyzmax(:)=dble(ixyzmax(:))*dble(dxyz(:))
xyzl2(:)=(ixyzmax(:)-ixyzmin(:)+1)*dble(dxyz(:))/2.0d+0 ! half total length
xyzc(:)=(xyzmax(:)+xyzmin(:))/2.0d+0 ! center location
!.....specify epsilon as 10**log10(machine precision)/2 
near=10**(0.5*nint(log10(abs(nearest(0.0,-1.0)))))
! small increment in each direction, based on discretization
smallxyz(1)=mx*small*dx
smallxyz(2)=my*small*dy
smallxyz(3)=mz*small*dz
write(iout,1000)near,smallxyz
write(*,1000)near,smallxyz
1000 format(" Epsilon based on platform precision:",e10.5/& 
            "                            small dx:",e10.5/& 
            "                            small dy:",e10.5/& 
            "                            small dz:",e10.5) 
!//    " Epsilon to adjust particle locations (x,y,z): ",3(e10.5,1x)//)
! DAB check these for species expansion:
netmassp=0.0; netmassm=0.0
netxyzm=0; netxyzp=0
mass=0
netsplit=0
!.....cell face areas
ax=dy*dz; ay=dx*dz; az=dx*dy
vol=az*dz ! cell volume
!.....rmx,rmy,rmz = large number if dimension inactive, or else = 1.0
rmx=(dble(large)-dble(mx)*dble(large)+1.0d+0) 
rmy=(dble(large)-dble(my)*dble(large)+1.0d+0) 
rmz=(dble(large)-dble(mz)*dble(large)+1.0d+0) 
!.....time
curtime=tinit
dt=dtinit
!.....OUTPUT CONTROL

if(ibug.ge.1)then;string(1)=' OPENING OUTPUT CONTROL FILE!';call debug(ibugout,string,nline);endif
inquire (file=fnameopc, exist=flexist)
if(.not.flexist)stop ' Output control file does not exist'
open(inopc,file=fnameopc,status="old")
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR OPC ARRAY!';call debug(ibugout,string,nline);endif
allocate (opc(1:nopc), STAT = ierror)
if(ierror.ne.0)goto 9994
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR OUTUNIT ARRAY!';call debug(ibugout,string,nline);endif
allocate (outunit(1:nopc+1), STAT = ierror)
if(ierror.ne.0)goto 9994
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR OUTFNAME ARRAY!';call debug(ibugout,string,nline);endif
allocate (outfname(1:nopc), STAT = ierror)
if(ierror.ne.0)goto 9994
outunit(1)=iplotm
outunit(2)=iplotc
outunit(3)=iplotd
outunit(4)=iplotb
outunit(5)=iplotp
outunit(6)=iploti
outunit(7)=iplotsc
outunit(8)=iplotsm
!
!.....PARAMETERS
!
if(ibug.ge.1)then;string(1)=' OPENING PARAMETER FILE!';call debug(ibugout,string,nline);endif
inquire (file=fnamepar, exist=flexist)
if(.not.flexist)stop ' Parameter file does not exist'
open(inpar,file=fnamepar,status="old")  
!.....read number of zones, if zero use one zone
call skip(inpar)
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR PARAMETER ARRAYS!';call debug(ibugout,string,nline);endif
inobgr=0
read(inpar,*,err=9990)nzone
if(nzone.eq.0)then
  nzone=nx*ny*nz
  inobgr=1
endif
!.....allocate memory
allocate (por(nzone))
if(ierror.ne.0)goto 9994
allocate (ret(nzone))
if(ierror.ne.0)goto 9994
allocate (dlong(nzone))
if(ierror.ne.0)goto 9994
allocate (dtran(nzone))
if(ierror.ne.0)goto 9994
allocate (ddiff(nzone))
if(ierror.ne.0)goto 9994
allocate (decay(nzone))
if(ierror.ne.0)goto 9994
if(ibug.ge.1)then;string(1)=' READING PARAMETERS!';call debug(ibugout,string,nline);endif
call parinput(por,ret,dlong,dtran,ddiff,decay,cat)
!
!.....BOUNDARIES
!
if(ibug.ge.1)then;string(1)=' OPENING BOUNDARY FILE!';call debug(ibugout,string,nline);endif
inquire (file=fnamebnd, exist=flexist)
if(.not.flexist)stop ' Boundary file does not exist'
open(inbnd,file=fnamebnd,status="old")
!.....determine the number of sources and boundaries
if(ibug.ge.1)then;string(1)=' COMPUTING NUMBER OF SOURCES AND BOUNDARIES!';call debug(ibugout,string,nline);endif
maxbnd=0
maxsource=0
do
  call skip(inbnd)
  read(inbnd,*,iostat=iread_error)itype
  if(iread_error.eq.-1)exit
  if(itype.eq.1.or.itype.eq.2.or.itype.eq.6.or.itype.eq.8.or.itype.eq.9.or.itype.eq.10.or.itype.eq.11)maxsource=maxsource+1
  maxbnd=maxbnd+1
enddo
rewind(inbnd)
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR BOUNDS ARRAY!';call debug(ibugout,string,nline);endif
allocate (bounds(1:maxbnd), STAT = ierror)
if(ierror.ne.0)goto 9994
if(ibug.ge.1)then;string(1)=' ALLOCATING MEMORY FOR SOURCE ARRAY!';call debug(ibugout,string,nline);endif
allocate (source(1:maxsource), STAT = ierror)
if(ierror.ne.0)goto 9994
if(ibug.ge.1)then;string(1)=' READING BOUNDARY CONDITIONS!';call debug(ibugout,string,nline);endif
call bndinput(bounds,cat,por,ret)
!
!.....POINT SOURCES
!   
if(ibug.ge.1)then;string(1)=' OPENING POINT SOURCE FILE!';call debug(ibugout,string,nline);endif
inquire (file=fnamepnt, exist=flexist)
if(.not.flexist)stop ' Point source file does not exist'
open(inpnt,file=fnamepnt,status="unknown")
if(ibug.ge.1)then;string(1)=' READING POINT SOURCE INPUT!';call debug(ibugout,string,nline);endif
call skip(inpnt)
tnextpnt=dble(large)
read(inpnt,*,iostat=iread_error)tnextpnt,npntsrc
if(iread_error.eq.-2)goto 9992
!
!.....VELOCITY CONTROL
!
if(ibug.ge.1)then;string(1)=' OPENING VLC FILE!';call debug(ibugout,string,nline);endif
inquire (file=fnamevel, exist=flexist)
if(.not.flexist)then
  write(iout,'(a,a)')' Velocity file does not exist:',fnamevel
  stop ' Velocity file does not exist - see output file'
endif
open(invlc,file=fnamevel,status="old")
!..... MONITORING
if(ibug.ge.1)then;string(1)=' READING MONITORING LOCATIONS!';call debug(ibugout,string,nline);endif
inquire (file=fnamesam, exist=flexist)
if(.not.flexist)stop ' Monitoring file does not exist'
open(insam,file=fnamesam,status="old")
nsam=0
do
  call skip(insam)
  read(insam,*,iostat=iread_error)itype
  if(iread_error.eq.-1)exit
  nsam=nsam+1
enddo
rewind(insam)
allocate (sam(nsam), STAT = ierror)
if(ierror.ne.0)goto 9994
call saminput(sam,cat,por,ret)
!
!.....INITIALIZE SYSTEM
!
!.....initialize output control 
if(ibug.ge.1)then;string(1)=' READY TO CRANK OUT SOLUTION!';call debug(ibugout,string,nline);endif
if(ibug.ge.1)then;string(1)=' INITIALIZING SYSTEM!';call debug(ibugout,string,nline);endif
if(ibug.ge.1)then;string(1)=' OPCUPDT!';call debug(ibugout,string,nline);endif
call opcupdt(opc,outunit,outfname)
!.....initialize Darcy velocities
if(ibug.ge.1)then;string(1)=' VELUPDT!';call debug(ibugout,string,nline);endif
call velupdt(vel3,por,cat,rech,chd)
if(ibug.ge.1)then;string(1)=' MAWUPDT!';call debug(ibugout,string,nline);endif
call mawupdt(bounds,vel3)
!.....initialize time step control field
if(ibug.ge.1)then;string(1)=' TCNTRL!';call debug(ibugout,string,nline);endif
if(dtcntrl.ne.0.0)call tcntrl(vel3,por,ret,dtran,dlong,ddiff,cat)
!.....initialize boundary conditions
if(curtime.eq.tnextbnd)then
  if(ibug.ge.1)then;string(1)=' BNDUPDT!';call debug(ibugout,string,nline);endif
  call bndupdt(bounds,source,cat,pat)
!.....initialize Courant condition for sources 
  call courant(source,bounds,vel3,cat,por)
endif
!.....initial particle distribution (DAB MUST INPUT INITIAL IMMOBILE PARTICLES HERE)
if(curtime.eq.tnextpnt)call pntupdt(pat,imp,cat,vel3)
!.....initialize sources boundary information
if(ibug.ge.1)then;string(1)=' SRCUPDT!';call debug(ibugout,string,nline);endif
call srcupdt(source,bounds,pat,cat)
!.....output at initial time?
if(curtime.eq.tnextopc)then
  if(ibug.ge.1)then;string(1)=' OUTPUT!';call debug(ibugout,string,nline);endif
  call output(opc,outunit,outfname,sam,pat,cat,vel3,bounds,por,ret,decay)
endif
! output mt3d cnf file
call cnf(outfname(2),icnf,nx,ny,nz,dx,dy,dz)
if(np.ne.0)call split(pat,cat,bounds,por,ret)
!
!.....SOLVE PROBLEM (NOW INCLUDING MIXING AND REACTIONS)
!
!.....pass the appropriate function to solve
!.....movep1: advection only
!.....movep2: discrete D, isotropic dispersion
!.....movep3: interpolated D, isotropic dispersion
!.....movep4: interpolated D, anisotropic dispersion
if(iadvect.eq.1)then
  if(ibug.ge.1)then;string(1)=' Using MOVEP1 - advection only';call debug(ibugout,string,nline);endif
  call solve(movep1,sam,cat,pat,rech,chd,vel3,por,ret,&
  dlong,dtran,ddiff,decay,bounds,source,opc,outunit,outfname)
else
  if(isotropic.eq.1.and.interp.eq.0)then                        ! isotropic D
    if(ibug.ge.1)then;string(1)=' Using MOVEP2 - isotropic, no interpolation';call debug(ibugout,string,nline);endif
    call solve(movep2,sam,cat,pat,rech,chd,vel3,por,ret,&
    dlong,dtran,ddiff,decay,bounds,source,opc,outunit,outfname)
  elseif(isotropic.eq.1.and.interp.eq.1)then                    ! isotropic D, interpolated
    if(ibug.ge.1)then;string(1)=' Using MOVEP3 - isotropic, interpolation';call debug(ibugout,string,nline);endif
    call solve(movep3,sam,cat,pat,rech,chd,vel3,por,ret,&
    dlong,dtran,ddiff,decay,bounds,source,opc,outunit,outfname)
  elseif(isotropic.eq.0.and.interp.eq.1)then                    ! anisotropic D, interpolated
    if(ibug.ge.1)then;string(1)=' Using MOVEP4 - anisotropic, interpolation';call debug(ibugout,string,nline);endif
    call solve(movep4,sam,cat,pat,rech,chd,vel3,por,ret,&
    dlong,dtran,ddiff,decay,bounds,source,opc,outunit,outfname)
  else
    stop ' Must use interp=1 for anisotropic dispersion tensor'
  endif
endif
! ELB 7-9-01 call output(opc,outunit,outfname,pat,cat,vel3,bounds,por,ret,decay)
!.....deallocate memory
if(ibug.ge.1)then;string(1)=' DEALLOCATE MEMORY!';call debug(ibugout,string,nline);endif
deallocate (vel3)
deallocate (pat)
deallocate (cat)
deallocate (opc)
deallocate (outunit)
deallocate (outfname)
deallocate (por)
deallocate (ret)
deallocate (dlong)
deallocate (dtran)
deallocate (ddiff)
deallocate (bounds)
deallocate (source)
deallocate (rech)
if(allocated(chd))deallocate (chd)
!
!istat=dtime(tarray)
print*
print*,' Elapsed time:    ',tarray
print*
stop ' Normal Termination'
!.....errors
9990 stop ' Error in parameter input file!'
9991 stop ' Error in boundary input file!'
9992 stop ' Error in point source input file!'
9993 stop ' Error reading file name'
9994 stop ' Error allocating memory'
end
!------------------------------------------------------------
! degbug
!------------------------------------------------------------
subroutine debug(ibugout,string,nline)
character (len=80) string(5)
do iline=1,nline
  write(ibugout,*)string(iline)
!  write(*,*)string(iline)
enddo
return
end

!------------------------------------------------------------
! maininput
!------------------------------------------------------------
subroutine maininput()
! basic information, problem description, allocation of memory, etc.
use global
character (len=80) heading(2),dbgfile
character (len=3) lpack(npack)
! read titles
call skip(inbas)
read(inbas,'(A)',err=9999) heading(1)
call skip(inbas)
read(inbas,'(A)',err=9999) heading(2)
!
write(iout,'(1X,/1X,A)')   heading(1)
write(iout,'(1X,A)')       heading(2)
write(*,'(1X,/1X,A)')      heading(1)
write(*,'(1X,A)')          heading(2)
! read debug option and debug output file
call skip(inbas)
read(inbas,*,err=9999)ibug
call skip(inbas)
read(inbas,'(a40)',err=9999)dbgfile
open(ibugout,file=dbgfile)
! 
if(ibug.ge.1)then;string(1)=' READING MAIN INPUT!';call debug(ibugout,string,1);endif
! read flag to write out mt3d output or not
call skip(inbas)
read(inbas,*,err=9999)nx,ny,nz
call skip(inbas)
read(inbas,*,err=9999)dx,dy,dz
!.minimum a maximum extent of simulation
call skip(inbas)
read(inbas,*,err=9999)ixyzmin(1),ixyzmax(1) ! model domain includes ixmin and ixmax 
call skip(inbas)
read(inbas,*,err=9999)ixyzmin(2),ixyzmax(2)
call skip(inbas)
read(inbas,*,err=9999)ixyzmin(3),ixyzmax(3)
!
if(ixyzmin(1).lt.1.or.ixyzmax(1).gt.nx.or.&
   ixyzmin(2).lt.1.or.ixyzmax(2).gt.ny.or.&
   ixyzmin(3).lt.1.or.ixyzmax(3).gt.nz)then
   write(*,*)' error in ixmin ... izmax '
   write(iout,*)' error in ixmin ... izmax '
   goto 9999
endif
!
call skip(inbas)
read(inbas,*,err=9999)ixyzm(1),ixyzp(1)
call skip(inbas)
read(inbas,*,err=9999)ixyzm(2),ixyzp(2)
call skip(inbas)
read(inbas,*,err=9999)ixyzm(3),ixyzp(3)
if(abs(ixyzm(1)).gt.1.or.abs(ixyzp(1)).gt.1.or.&
   abs(ixyzm(2)).gt.1.or.abs(ixyzp(2)).gt.1.or.&
   abs(ixyzm(3)).gt.1.or.abs(ixyzp(3)).gt.1)then
   write(*,*)' error in ixm ... izp '
   write(iout,*)' error in ixm ... izp '
   goto 9999
endif
call skip(inbas)
read(inbas,*,err=9999)maxnp
call skip(inbas)
read(inbas,*,err=9999)iseed1,iseed2
call skip(inbas)
read(inbas,*,err=9999)interp,bke
if(abs(bke).gt.1)goto 9999
call skip(inbas)
read(inbas,*,err=9999)cres,nres
call skip(inbas)
read(inbas,*,err=9999)tmax,dtinit,tinit,dtcntrl
call skip(inbas)
read(inbas,'(a)',iostat=iostatus)fnameopc
if(iostatus.ne.0)then
  if(iostatus.eq.-1)then
    write(iout,*)' End of file encountered reading output control file name in main input'
    write(*,*)' End of file encountered reading output control file name in main input '
  else
    write(iout,*)' Error reading output control file name in main input'
    write(*,*)' End of file encountered reading output control file name in main input '
  endif
endif  
call skip(inbas)
read(inbas,'(a)',iostat=iostatus)fnamepar
if(iostatus.ne.0)then
  if(iostatus.eq.-1)then
    write(iout,*)' End of file encountered reading parameter file name in main input'
    write(*,*)' End of file encountered reading parameter file name in main input '
  else
    write(iout,*)' Error reading parameter file name in main input'
    write(*,*)' End of file encountered reading parameter file name in main input '
  endif
endif  
call skip(inbas)
read(inbas,'(a)',iostat=iostatus)fnamevel
if(iostatus.ne.0)then
  if(iostatus.eq.-1)then
    write(iout,*)' End of file encountered reading velocity file name in main input'
    write(*,*)' End of file encountered reading velocity file name in main input '
  else
    write(iout,*)' Error reading velocity file name in main input'
    write(*,*)' End of file encountered reading velocity file name in main input '
  endif
endif  
call skip(inbas)
read(inbas,'(a)',iostat=iostatus)fnamebnd
if(iostatus.ne.0)then
  if(iostatus.eq.-1)then
    write(iout,*)' End of file encountered reading boundary file name in main input'
    write(*,*)' End of file encountered reading boundary file name in main input '
  else
    write(iout,*)' Error reading boundary file name in main input'
    write(*,*)' End of file encountered reading boundary file name in main input '
  endif
endif  
call skip(inbas)
read(inbas,'(a)',iostat=iostatus)fnamepnt
if(iostatus.ne.0)then
  if(iostatus.eq.-1)then
    write(iout,*)' End of file encountered reading point file name in main input'
    write(*,*)' End of file encountered reading point file name in main input '
  else
    write(iout,*)' Error reading point file name in main input'
    write(*,*)' End of file encountered reading point file name in main input '
  endif
endif  
call skip(inbas)
read(inbas,'(a)',iostat=iostatus)fnamesam
if(iostatus.ne.0)then
  if(iostatus.eq.-1)then
    write(iout,*)' End of file encountered reading sample file name in main input'
    write(*,*)' End of file encountered reading sample file name in main input '
  else
    write(iout,*)' Error reading sample file name in main input'
    write(*,*)' End of file encountered reading sample file name in main input '
  endif
endif  
! compute total number of nodes
nxyz=nx*ny*nz
! commonly used number
nxy=nx*ny
! indices indicating active dimensions
mx=1; my=1; mz=1
if(nx.eq.1)mx=0; if(ny.eq.1)my=0; if(nz.eq.1)mz=0
! echo main input parameters
write(iout,1001)ibug,dbgfile,nx,ny,nz,nxyz,dx,dy,dz,&
                ixyzmin(1)*dx-dx,ixyzmax(1)*dx,ixyzmin(2)*dy-dy,&
                ixyzmax(2)*dy,ixyzmin(3)*dz-dz,ixyzmax(3)*dz,&
                maxnp,&
                interp,bke,cres,nres,iseed1,iseed2,&
                sngl(tmax),dtinit,tinit,dtcntrl,&
                ixyzm(1),ixyzp(1),ixyzm(2),&
                ixyzp(2),ixyzm(3),ixyzp(3),&
                fnamebas,fnameopc,fnamepar,&
                fnamevel,fnamebnd,fnamepnt,fnamesam
write(*,1001)ibug,dbgfile,nx,ny,nz,nxyz,dx,dy,dz,&
                ixyzmin(1)*dx-dx,ixyzmax(1)*dx,ixyzmin(2)*dy-dy,&
                ixyzmax(2)*dy,ixyzmin(3)*dz-dz,ixyzmax(3)*dz,&
                maxnp,&
                interp,bke,cres,nres,iseed1,iseed2,&
                sngl(tmax),dtinit,tinit,dtcntrl,&
                ixyzm(1),ixyzp(1),ixyzm(2),&
                ixyzp(2),ixyzm(3),ixyzp(3),&
                fnamebas,fnameopc,fnamepar,&
                fnamevel,fnamebnd,fnamepnt,fnamesam
 1001 format('                         G R I D '/&
       ' debug level                     ',20('.'),3x,i15/&
       ' debug file                      ',20('.'),3x,a40/&
       ' number of cells in x direction  ',20('.'),3x,i15/&
       '                    y direction  ',20('.'),3x,i15/&
       '                    z direction  ',20('.'),3x,i15/&
       ' total number of cells           ',20('.'),3x,i15//&

       ' discretization in x direction   ',20('.'),3x,e15.8/&
       '                   y direction   ',20('.'),3x,e15.8/&
       '                   z direction   ',20('.'),3x,e15.8//&

       ' extent of domain: xmin          ',20('.'),3x,e15.8/&
       '                   xmax          ',20('.'),3x,e15.8/&
       '                   ymin          ',20('.'),3x,e15.8/&
       '                   ymax          ',20('.'),3x,e15.8/&
       '                   zmin          ',20('.'),3x,e15.8/&
       '                   zmax          ',20('.'),3x,e15.8//&

       '                  S O L U T I O N  C O N T R O L'/&
       ' maximum number of particles     ',20('.'),3x,i15/&
       ' interpolate D                   ',20('.'),3x,i15/&
       ' backward Kolmogorov''s equation ',20('.'),3x,e15.8/&
       ' minimum conc. for splitting     ',20('.'),3x,e15.8/&
       ' minimum particle resolution     ',20('.'),3x,i15/&
       ' random # seed for displacements ',20('.'),3x,i15/&
       ' random # seed                   ',20('.'),3x,i15//&

       '                            T I M E'/&
       ' simulation time                 ',20('.'),3x,e15.8/&
       ' maximum time step               ',20('.'),3x,e15.8/&
       ' initial time                    ',20('.'),3x,e15.8/&
       ' time step control constant      ',20('.'),3x,e15.8//&

       '               E X T E R N A L  B O U N D A R I E S '/&
       ' reflect at xmin (1) yes (0) no  ',20('.'),3x,i15/&
       '            xmax                 ',20('.'),3x,i15/&
       '            ymin                 ',20('.'),3x,i15/&
       '            ymax                 ',20('.'),3x,i15/&
       '            zmin                 ',20('.'),3x,i15/&
       '            zmax                 ',20('.'),3x,i15//&

       '                      I N P U T   F I L E S '/&
       ' main input       :',3x,a/&
       ' output control   :',3x,a/&
       ' parameters       :',3x,a/&
       ' velocity control :',3x,a/&
       ' boundaries       :',3x,a/&
       ' point sources    :',3x,a/&
       ' monitoring points:',3x,a/)

return
 9999      write(iout,*)' ERROR IN MAIN INPUT FILE'
stop ' ERROR IN MAIN INPUT FILE'
end
!------------------------------------------------------------
! bndinput
!------------------------------------------------------------
subroutine bndinput(bounds,cat,por,ret)
! read boundary conditions
use global
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
real por(nzone),ret(nzone),pmass(nspec),conc(nspec),massrate(nspec),refresh,nptime
integer nitype(nbtype),iostatus,nprint
character(80)::nsprint,fmt1,fmt2,fmt3
write(nsprint,'(i3)')nspec
fmt1=' ('' aq. concentrations    =  '', '//trim(nsprint)//'es12.4)'
fmt2=' ('' particle masses       =  '', '//trim(nsprint)//'es12.4)'
fmt3=' ('' mass loadings (M/T)   =  '', '//trim(nsprint)//'es12.4)'

!print*,trim(nsprint),fmt1
! initialize # of boundaries nbounds
nprint=nspec
nbounds=0
! initialize time to update boundaries
tnextbnd=dble(large)
!
write(iout,'(/a)')&
'                    B O U N D A R Y   C O N D I T I O N S'
write(*,'(/a)')&
'                    B O U N D A R Y   C O N D I T I O N S'
! 
do ibtype=1,nbtype
  nitype(ibtype)=0
enddo

do
  call skip(inbnd)
  read(inbnd,*,iostat=iread_error)itype
  if(iread_error.eq.-1)then
    exit
  elseif(iread_error.eq.-2)then
    goto 9999
  endif
  nbounds=nbounds+1
  if(nbounds.gt.maxbnd)then
    write(iout,*)' Error in boundary file, check formatting'
       stop ' Error in boundary file, check formatting'
  endif
  backspace (inbnd)
!.source
  if(itype.eq.1)then
    read(inbnd,*,iostat=iostatus)itype,i,j,kbot,ktop,npart,conc,&
                                 refresh,tbeg,tend,numbnd
    write(iout,2000,err=9999)nbounds,itype,i,j,kbot,ktop,npart,&
                             refresh,tbeg,tend,numbnd
    write(iout,fmt1,err=9999)conc                                  ! DAB quicky fix
    write(*,fmt1,err=9999)conc                                  ! DAB quicky fix
    nitype(itype)=nitype(itype)+1
!.flux condition - source that does not absorb after each time step
  elseif(itype.eq.2)then
    read(inbnd,*,err=9999)itype,x,y,z,sx,sy,sz,pmass,&
                             massrate,tbeg,tend,numbnd
    write(iout,2001,err=9999)nbounds,itype,x,y,z,sx,sy,sz,&
                             tbeg,tend,numbnd
    write(iout,fmt2,err=9999)pmass
    write(iout,fmt3,err=9999)massrate
    nitype(itype)=nitype(itype)+1
!.absorbing 
  elseif(itype.eq.3)then
    read(inbnd,*,iostat=iostatus)itype,i,j,kbot,ktop,tbeg,tend,numbnd
    write(iout,2002,err=9999)nbounds,itype,i,j,kbot,ktop,tbeg,tend,numbnd
    nitype(itype)=nitype(itype)+1
!.well with specified flux 
  elseif(itype.eq.4)then
    read(inbnd,*,iostat=iostatus)itype,i,j,kbot,ktop,wflux,tbeg,tend,numbnd
    write(iout,2003,err=9999)nbounds,itype,i,j,kbot,ktop,wflux,tbeg,tend,numbnd
    nitype(itype)=nitype(itype)+1
!.well with unknown flux, flux calculated from divergence 
  elseif(itype.eq.5)then
    read(inbnd,*,iostat=iostatus)itype,i,j,kbot,ktop,tbeg,tend,numbnd
    write(iout,2004,err=9999)nbounds,itype,i,j,kbot,ktop,tbeg,tend,numbnd
    nitype(itype)=nitype(itype)+1
!.flux over specified region
  elseif(itype.eq.6)then
    read(inbnd,*,iostat=iostatus)itype,i,j,kbot,ktop,pmass,&
                             massrate,tbeg,tend,numbnd
    write(iout,2006,err=9999)nbounds,itype,i,j,kbot,ktop,&
                             tbeg,tend,numbnd
    write(iout,fmt2,err=9999)pmass
    write(iout,fmt3,err=9999)massrate
    nitype(itype)=nitype(itype)+1
!.well with unknown flux, flux calculated from divergence 
  elseif(itype.eq.7)then
    read(inbnd,*,iostat=iostatus)itype,i,j,kbot,ktop,tbeg,tend,numbnd
    write(iout,2007,err=9999)nbounds,itype,i,j,kbot,ktop,tbeg,tend,numbnd
    nitype(itype)=nitype(itype)+1
!.recharge boundary (must have recharge in the cell-by-cell flux file)
  elseif(itype.eq.8)then
    read(inbnd,*,err=9999)itype,x,y,sx,sy,pmass,&
                             massrate,tbeg,tend,numbnd
    write(iout,2008,err=9999)nbounds,itype,x,y,sx,sy,&
                             tbeg,tend,numbnd
    write(iout,fmt2,err=9999)pmass
    write(iout,fmt3,err=9999)massrate
    nitype(itype)=nitype(itype)+1
!.recharge boundary w/ concentration (must have recharge in the cell-bycell flux file)
  elseif(itype.eq.9)then
    read(inbnd,*,err=9999)itype,im,jm,ip,jp,conc,nptime,tbeg,tend,numbnd
    write(iout,2009,err=9999)nbounds,itype,im,jm,ip,jp,nptime,tbeg,tend,numbnd
    write(iout,fmt1,err=9999)conc                                  ! DAB quicky fix
    nitype(itype)=nitype(itype)+1
  elseif(itype.eq.10)then
    read(inbnd,*,err=9999)itype,im,jm,km,ip,jp,kp,conc,nptime,tbeg,tend,numbnd
    write(iout,2010,err=9999)nbounds,itype,im,jm,km,ip,jp,kp,nptime,tbeg,tend,numbnd
    write(iout,fmt1,err=9999)conc                                  ! DAB quicky fix
    nitype(itype)=nitype(itype)+1
  elseif(itype.eq.11)then
    read(inbnd,*,err=9999)itype,im,jm,km,ip,jp,kp,conc,nptime,tbeg,tend,numbnd
    write(iout,2011,err=9999)nbounds,itype,im,jm,km,ip,jp,kp,nptime,tbeg,tend,numbnd
    write(iout,fmt1,err=9999)conc                                  ! DAB quicky fix
    nitype(itype)=nitype(itype)+1
  endif
!
  if((itype.eq.1.or.itype.eq.3.or.itype.eq.4.or.itype.eq.5.or.itype.eq.6.or.itype.eq.7).and.&
     ((mx.eq.1.and.(i.lt.1.or.i.gt.nx)).or.&
     (my.eq.1.and.(j.lt.1.or.j.gt.ny)).or.&
     (mx.eq.0.and.i.ne.1).or.&
     (my.eq.0.and.j.ne.1)))then
     print*,' Boundary not in domain'
     write(*,1000)nbounds,itype
     write(iout,1000)nbounds,itype
     stop
  elseif((itype.eq.2).and.((x.gt.nx*dx.or.x.lt.0.).or.(y.gt.ny*dy.or.y.lt.0.).or.&
                           (z.gt.nz*dz.or.z.lt.0.)))then
     print*,' Boundary not in domain'
     write(*,1000)nbounds,itype
     write(iout,1000)nbounds,itype
     stop
  elseif((itype.eq.8).and.((x.gt.nx*dx.or.x.lt.0.).or.(y.gt.ny*dy.or.y.lt.0.)))then
     print*,' Boundary not in domain'
     write(*,1000)nbounds,itype
     write(iout,1000)nbounds,itype
     stop
  else
!...location
    if(itype.eq.1)then
      kx=i-1; ky=j-1
      bounds(nbounds)%bc_xyzm(1)=kx*dx
      bounds(nbounds)%bc_xyzm(2)=ky*dy
      bounds(nbounds)%bc_xyzm(3)=(kbot-1)*dz
      bounds(nbounds)%bc_xyzs(1)=dx
      bounds(nbounds)%bc_xyzs(2)=dy
      bounds(nbounds)%bc_xyzs(3)=(ktop-kbot+1)*dz
      bounds(nbounds)%refresh=refresh
      bounds(nbounds)%np=npart
   do iloop=1,nspec
       if(npart.ne.0)then
        bounds(nbounds)%pmass(iloop)=(ktop-kbot+1)*conc(iloop)*&
              vol*por(cat(i,j,kbot)%zone)*ret(cat(i,j,kbot)%zone)/npart
       else
        bounds(nbounds)%pmass(iloop)=0.0
       endif
   enddo
      bounds(nbounds)%ijk(1)=i
      bounds(nbounds)%ijk(2)=j
      bounds(nbounds)%kbot=kbot
      bounds(nbounds)%ktop=ktop
    elseif(itype.eq.2)then
      bounds(nbounds)%bc_xyzm(1)=x-.5*sx
      bounds(nbounds)%bc_xyzm(2)=y-.5*sy
      bounds(nbounds)%bc_xyzm(3)=z-.5*sz
      bounds(nbounds)%bc_xyzs(1)=sx
      bounds(nbounds)%bc_xyzs(2)=sy
      bounds(nbounds)%bc_xyzs(3)=sz
      do iloop=1,nspec
        bounds(nbounds)%massrate(iloop)=massrate(iloop)
        bounds(nbounds)%pmass(iloop)=pmass(iloop)
      enddo
    elseif(itype.eq.3)then
      bounds(nbounds)%ijk(1)=i
      bounds(nbounds)%ijk(2)=j
      bounds(nbounds)%kbot=kbot
      bounds(nbounds)%ktop=ktop
    elseif(itype.eq.4)then
      bounds(nbounds)%ijk(1)=i
      bounds(nbounds)%ijk(2)=j
      bounds(nbounds)%kbot=kbot
      bounds(nbounds)%ktop=ktop
      bounds(nbounds)%flux=wflux  ! DAB check this
    elseif(itype.eq.5)then
      bounds(nbounds)%ijk(1)=i
      bounds(nbounds)%ijk(2)=j
      bounds(nbounds)%kbot=kbot
      bounds(nbounds)%ktop=ktop
    elseif(itype.eq.6)then
      bounds(nbounds)%ijk(1)=i
      bounds(nbounds)%ijk(2)=j
      bounds(nbounds)%kbot=kbot
      bounds(nbounds)%ktop=ktop
      bounds(nbounds)%bc_xyzs(1)=dx
      bounds(nbounds)%bc_xyzs(2)=dy
      bounds(nbounds)%bc_xyzs(3)=dz
      do iloop=1,nspec
        bounds(nbounds)%massrate(iloop)=massrate(iloop)
        bounds(nbounds)%pmass(iloop)=pmass(iloop)
      enddo
    elseif(itype.eq.7)then
      bounds(nbounds)%ijk(1)=i
      bounds(nbounds)%ijk(2)=j
      bounds(nbounds)%kbot=kbot
      bounds(nbounds)%ktop=ktop
    elseif(itype.eq.8)then
      bounds(nbounds)%bc_xyzm(1)=x-.5*sx
      bounds(nbounds)%bc_xyzm(2)=y-.5*sy
      bounds(nbounds)%bc_xyzs(1)=sx
      bounds(nbounds)%bc_xyzs(2)=sy
      do iloop=1,nspec
        bounds(nbounds)%massrate(iloop)=massrate(iloop)
        bounds(nbounds)%pmass(iloop)=pmass(iloop)
      enddo
    elseif(itype.eq.9)then
      bounds(nbounds)%ijkm(1)=im
      bounds(nbounds)%ijkm(2)=jm
      bounds(nbounds)%ijkp(1)=ip
      bounds(nbounds)%ijkp(2)=jp
      do iloop=1,nspec
        bounds(nbounds)%conc(iloop)=conc(iloop)
      enddo
      bounds(nbounds)%nptime=nptime
    elseif(itype.eq.10.or.itype.eq.11)then
      bounds(nbounds)%ijkm(1)=im
      bounds(nbounds)%ijkm(2)=jm
      bounds(nbounds)%ijkm(3)=km
      bounds(nbounds)%ijkp(1)=ip
      bounds(nbounds)%ijkp(2)=jp
      bounds(nbounds)%ijkp(3)=kp
      do iloop=1,nspec
        bounds(nbounds)%conc(iloop)=conc(iloop)
      enddo
      bounds(nbounds)%nptime=nptime ! particles per unit time
    endif
    bounds(nbounds)%bc_type=itype
    bounds(nbounds)%tbeg=tbeg
    bounds(nbounds)%tend=tend
!...zero counters
    bounds(nbounds)%np_remove=0
    bounds(nbounds)%group=numbnd
    do iloop=1,nspec
      bounds(nbounds)%mass_remove(iloop)=0.
      bounds(nbounds)%maw(iloop)=0.
    enddo
!...particle masses
!...particle masses
    tnextbnd=min(dble(tbeg),tnextbnd)
  endif
enddo 
!
nbtotal=0
do ibtype=1,nbtype
  nbtotal=nbtotal+nitype(ibtype)
enddo
write(iout,2005)(nitype(itype),itype=1,nbtype),nbtotal
write(*,2005)(nitype(itype),itype=1,nbtype),nbtotal
!.....initialize mass and particle counters for various boundary types
do ibtype=1,nbtype+1
  npbc(ibtype)=0
  do iloop=1,nspec
    massbc(itype,iloop)=0.
  enddo
enddo
!.....check for multiple boundaries defined at same time and cell location

!do ibounds=1,nbounds
!  itype=bounds(ibounds)%bc_type
!  tbeg=bounds(ibounds)%tbeg
!  tend=bounds(ibounds)%tend
!  i=bounds(ibounds)%ijk(1); j=bounds(ibounds)%ijk(2); k=bounds(ibounds)%ijk(3)
!  l4=i+(j-1)*nx+(k-1)*nxy
!  do ibounds2=1,nbounds
!    itype2=bounds(ibounds2)%bc_type
!    tbeg2=bounds(ibounds2)%tbeg
!    tend2=bounds(ibounds2)%tend
!    i2=bounds(ibounds2)%ijk(1); j2=bounds(ibounds2)%ijk(2); k2=bounds(ibounds2)%ijk(3)
!    l42=i2+(j2-1)*nx+(k2-1)*nxy
!    if(itype.ne.0.and.itype2.ne.0.and.itype.ne.2.and.itype2.ne.2.and.l4.eq.l42)then
!      if((tend2.gt.tbeg.and.tend2.lt.tend).or.(tbeg2.gt.tbeg.and.tbeg2.lt.tend))then
!        write(*,*)bounds(ibounds2)
!        stop ' Multiple boundaries defined at same time and place'
!      endif            
!!    endif
!  enddo
!enddo
return

9999 stop ' Error in boundary file'
1000 format ('Error in boundary condition #',i10,' TYPE = ',i10,', aborting simulation')
2000 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,' SOURCE'/&
       ' location (i,j)        =  ',2i10/&
       ' location (kbot,ktop)  =  ',2i10/&
       ' number of particles   =  ',i10/&
!       ' aq. concentrations    =  ',e10.5/&
       ' refresh rate control  =  ',e10.5/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10)
2001 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,' FLUX'/&
       ' patch location (x,y,z)=  ',3e10.5/&
       ' patch size (sx,sy,sz) =  ',3e10.5/&
!       ' particle masses       =  ',e10.5/&
!       ' mass loadings (M/T)   =  ',100e10.5/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10)
2002 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,' ABSORBING'/&
       ' location (i,j)        =  ',2i10/&
       ' location (kbot,ktop)  =  ',2i10/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10/)
2003 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,&
                                      ' ABSORBING SPECIFIED FLUX'/&
       ' location (i,j)        =  ',2i10/&
       ' location (kbot,ktop)  =  ',2i10/&
       ' flux                  =  ',e10.5/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10/)
2004 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,&
                                       ' ABSORBING COMPUTED FLUX'/&
       ' location (i,j)        =  ',2i10/&
       ' location (kbot,ktop)  =  ',2i10/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10/)
2006 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,&
                                       ' INJECTION WELL'/&
       ' location (i,j)        =  ',2i10/&
       ' location (kbot,ktop)  =  ',2i10/&
!       ' particle masses       =  ',100e10.5/&
!       ' mass loadings (M/T)   =  ',100e10.5/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10)
2007 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,&
                                       ' MAW BOUNDARY'/&
       ' location (i,j)        =  ',2i10/&
       ' location (kbot,ktop)  =  ',2i10/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10/)
2008 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,&
                                       ' RCH1 BOUNDARY'/&
       ' patch location (x,y)  =  ',3e10.5/&
       ' patch size (sx,sy)    =  ',3e10.5/&
!       ' particle masses       =  ',100e10.5/&
!       ' mass loadings (M/T)   =  ',100e10.5/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10)
2009 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,&
                                       ' RCH2 BOUNDARY'/&
       ' rech. i,j start loc.  =  ',2i10/&
       ' rech. i,j end loc.    =  ',2i10/&
!       ' concentration(s)      =  ',100e10.5/&
       ' particles/unit time   =  ',e10.5/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10)
2010 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,&
                                       ' CONSTANT CONC. INFLOW'/&
       ' i,j,k start loc.      =  ',3i10/&
       ' i,j,k end loc.        =  ',3i10/&
!       ' aq. concentration(s)  =  ',100e10.5/&
       ' particles/unit time   =  ',e10.5/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10)
2011 format(49('-'),' boundary # ',i10/&
       ' type                  =  ',i10,&
                                       ' CONSTANT CONC. ON BOUNDARY'/&
       ' i,j,k start loc.      =  ',3i10/&
       ' i,j,k end loc.        =  ',3i10/&
!       ' aq. concentration(s)  =  ',100e10.5/&
       ' particles/unit time   =  ',e10.5/&
       ' beginning time        =  ',e10.5/&
       ' ending time           =  ',e10.5/&
       ' boundary group        =  ',i10)
       
2005 format(' Boundary Type                   Number Read'/&
       ' ------------------------------  -----------'/&
       ' Constant Conc.             (1)',i10/&
       ' Continuous Source          (2)',i10/&
       ' Absorbing                  (3)',i10/&
       ' Absorbing (Specified Flux) (4)',i10/&
       ' Absorbing (Computed Flux)  (5)',i10/&
       ' Injection Well             (6)',i10/&
       ' MAW                        (7)',i10/&
       ' RCH TYPE 1 (mass rate)     (8)',i10/&
       ' RCH TYPE 2 (concentration) (9)',i10/&
       ' Constant Conc. Inflow     (10)',i10/&
       ' Constant Conc. on Bound.  (11)',i10/&
       ' TOTAL ........................',i10/)
end
!------------------------------------------------------------
! saminput
!------------------------------------------------------------
subroutine saminput(sam,cat,por,ret)
! read boundary conditions
use global
implicit none
type (sample)::   sam(nsam)
type (cell)::     cat(nx,ny,nz)
real:: por(nzone),ret(nzone)

real:: rvol,area,radius,xm,ym,zm,xp,yp,zp,x,y,z,zbot,ztop
integer:: itype,isam,i,j,k,kx,ky,kz,kzbot,kzbotm,kztop,izone
integer:: iread_error,iostatus
isam=0
do
  call skip(insam)
  read(insam,*,iostat=iread_error)itype
  if(iread_error.eq.-1)then
    exit
  elseif(iread_error.eq.-2)then
    goto 9999
  endif
  isam=isam+1
  backspace (insam)
!.monitoring point in x,y,zbot,ztop,radius
  if(itype.eq.1)then
    read(insam,*,iostat=iostatus)sam(isam)%itype,x,y,zbot,ztop,radius,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!    write(iout,2000,err=9999)isam,itype,x,y,zbot,ztop,radius
    if((x   .gt.nx*dx.or.x.lt.0.).or.(y   .gt.ny*dy.or.y.lt.0.).or.&
       (zbot.gt.nz*dz.or.zbot.lt.0.).or.(ztop.gt.nz*dz.or.ztop.lt.0.))then
     rvol=0.0
     write(iout,2000,err=9999)isam,itype,x,y,zbot,ztop,radius,vol
     write(iout,*)' monitoring point not in domain'
     print*,' monitoring point not in domain'
     stop
    endif
    kx=ifix(x/dx); ky=ifix(y/dy); kzbotm=ifix(zbot/dz)
    i=kx+1; j=ky+1; kzbot=kzbotm+1
    xm=kx*dx; ym=ky*dy; zm=kzbotm*dz; xp=xm+dx; yp=ym+dy; zp=zm+dz
    area=3.141592*radius*radius
    zm=zbot
    rvol=0.0
    k=kzbot
    do
      if(ztop.lt.zp)exit
      do izone=1,sam(isam)%nzone
        if(cat(i,j,k)%zone.eq.sam(isam)%zone(izone))then
          rvol=rvol+(zp-zm)*area*por(cat(i,j,k)%zone)
          exit
        endif
      enddo
      zm=zm+dz; zp=zp+dz; k=k+1
    enddo
    do izone=1,sam(isam)%nzone
      if(cat(i,j,k)%zone.eq.sam(isam)%zone(izone))then
        rvol=rvol+(ztop-zm)*area*por(cat(i,j,k)%zone)
        exit
      endif
    enddo
    sam(isam)%vol=rvol
    sam(isam)%ijk(1)=i
    sam(isam)%ijk(2)=j
    sam(isam)%ijk(3)=kzbot
    sam(isam)%ijk(4)=k
    sam(isam)%xyz(1)=x
    sam(isam)%xyz(2)=y
    sam(isam)%xyz(3)=zbot
    sam(isam)%xyz(4)=ztop
    sam(isam)%radius=radius
    write(iout,2000,err=9999)isam,sam(isam)%itype,x,y,zbot,ztop,radius,rvol,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!.x plane
  elseif(itype.eq.2)then
    read(insam,*,iostat=iostatus)sam(isam)%itype,i,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
    if(abs(i).gt.nx.or.abs(i).lt.1)then
     write(iout,*)' monitoring point not in domain'
     print*,' monitoring point not in domain'
     stop
    endif
    sam(isam)%ijk(1)=i
    write(iout,2001,err=9999)isam,sam(isam)%itype,i,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!.y plane
  elseif(itype.eq.3)then
    read(insam,*,iostat=iostatus)sam(isam)%itype,j,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
    if(abs(j)   .gt.ny.or.abs(j).lt.1)then
     write(iout,*)' monitoring point not in domain'
     print*,' monitoring point not in domain'
     stop
    endif
    sam(isam)%ijk(2)=j
    write(iout,2002,err=9999)isam,sam(isam)%itype,j,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!.z plane
  elseif(itype.eq.4)then
    read(insam,*,iostat=iostatus)sam(isam)%itype,k,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
    if(abs(k)   .gt.nz.or.abs(k).lt.1)then
     write(iout,*)' monitoring point not in domain'
     print*,' monitoring point not in domain'
     stop
    endif
    sam(isam)%ijk(3)=k
    write(iout,2003,err=9999)isam,sam(isam)%itype,k,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!.i,j,kbot,ktop
  elseif(itype.eq.5)then
    read(insam,*,iostat=iostatus)sam(isam)%itype,i,j,kzbot,kztop,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
    if(kzbot   .gt.nz.or.kzbot.lt.1)then
    if(kztop   .gt.nz.or.kztop.lt.1)then
     write(iout,*)' monitoring point not in domain'
     print*,' monitoring point not in domain',k,1,nz
     stop
    endif
    endif
    sam(isam)%ijk(1)=i
    sam(isam)%ijk(2)=j
    sam(isam)%ijk(3)=kzbot
    sam(isam)%ijk(4)=kztop
    sam(isam)%vol=0.0
    do k=kzbot,kztop
      if(sam(isam)%nzone.ne.0)then
       do izone=1,sam(isam)%nzone
        if(cat(i,j,k)%zone.eq.sam(isam)%zone(izone))then
          sam(isam)%vol=sam(isam)%vol+vol*por(cat(i,j,k)%zone)
          exit
        endif
       enddo
      else
       sam(isam)%vol=sam(isam)%vol+vol*por(cat(i,j,k)%zone)
      endif
    enddo
    write(iout,2004,err=9999)isam,sam(isam)%itype,i,j,kzbot,kztop,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!.x plane
  elseif(itype.eq.6)then
    read(insam,*,iostat=iostatus)sam(isam)%itype,i,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
    if(iabs(i)   .gt.nx)then
     write(iout,*)' monitoring point not in domain'
     print*,' monitoring point not in domain'
     stop
    endif
    sam(isam)%ijk(1)=i
    write(iout,2001,err=9999)isam,sam(isam)%itype,i,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!.y plane
  elseif(itype.eq.7)then
    read(insam,*,iostat=iostatus)sam(isam)%itype,j,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
    if(iabs(j)   .gt.ny)then
     write(iout,*)' monitoring point not in domain'
     print*,' monitoring point not in domain'
     stop
    endif
    sam(isam)%ijk(2)=j
    write(iout,2002,err=9999)isam,sam(isam)%itype,j,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!.z plane
  elseif(itype.eq.8)then
    read(insam,*,iostat=iostatus)sam(isam)%itype,k,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
    if(iabs(k)   .gt.nz)then
     write(iout,*)' monitoring point not in domain'
     print*,' monitoring point not in domain'
     stop
    endif
    sam(isam)%ijk(3)=k
    write(iout,2003,err=9999)isam,sam(isam)%itype,k,&
    sam(isam)%nzone,(sam(isam)%zone(izone),izone=1,sam(isam)%nzone)
!.i,j,kbot,ktop
  else
    goto 9999
  endif
enddo
return

9999 stop ' Error in monitoring file'
2000 format(49('-'),' monitoring interval # ',i10/&
       ' type                         =  ',i10,' x,y,zbot,ztop'/&
       ' location (x,y)               =  ',2(' ',g10.5)/&
       ' location (zbot,ztop)         =  ',2(' ',g10.5)/&
       ' radius                       =  ',g10.5/&
       ' vol                          =  ',g10.5/&
       ' number of zones to sample    =  ',i10/&
   100(' zone                         =  ',i10/))
2001 format(49('-'),' monitoring plane (max. conc. x-plane) # ',i10/&
       ' type                         =  ',i10,' PLANE'/&
       ' location (i)                 =  ',1(' ',i10)/&
   100(' zone                         =  ',i10/))
2002 format(49('-'),' monitoring plane (max. conc., or mass, y-plane) # ',i10/&
       ' type                         =  ',i10,' PLANE'/&
       ' location (j)                 =  ',1(' ',i10)/&
   100(' zone                         =  ',i10/))
2003 format(49('-'),' monitoring plane (max. conc., or mass, z-plane) # ',i10/&
       ' type                         =  ',i10,' PLANE'/&
       ' location (k)                 =  ',1(' ',i10)/&
   100(' zone                         =  ',i10/))
2004 format(49('-'),' monitoring interval for vertical column # ',i10/&
       ' type                         =  ',i10,' CELL'/&
       ' location (i,j,kbot,ktop)     =  ',4(' ',i10)/&
       ' vol                          =  ',g10.5/&
       ' number of zones to sample    =  ',i10/&
   100(' zone                         =  ',i10/))
2006 format(49('-'),' monitoring plane (mass on either side of a x-plane) # ',i10/&
       ' type                         =  ',i10,' PLANE'/&
       ' location (j)                 =  ',1(' ',i10)/&
   100(' zone                         =  ',i10/))
2007 format(49('-'),' monitoring plane (mass on either side of a y-plane) # ',i10/&
       ' type                         =  ',i10,' PLANE'/&
       ' location (k)                 =  ',1(' ',i10)/&
   100(' zone                         =  ',i10/))
2008 format(49('-'),' monitoring plane (mass on either side of a z-plane) # ',i10/&
       ' type                         =  ',i10,' CELL'/&
       ' location (i,j,kbot,ktop)     =  ',4(' ',i10)/&
       ' vol                          =  ',g10.5/&
       ' number of zones to sample    =  ',i10/&
   100(' zone                         =  ',i10/))
end
!----------------------------------------------------------------------
! parinput
!----------------------------------------------------------------------
subroutine parinput(por,ret,dlong,dtran,ddiff,decay,cat)
! parameter input
! parameters include porosity, longitudinal, and lateral dispersivities
! and retardation. Parameters can be indexed to zone values read from a 
! binary grid file (*.bgr).
use global
parameter (npar=6) 
type (cell)::     cat(nx,ny,nz)
real:: por(nzone),ret(nzone),dlong(nzone),dtran(nzone),ddiff(nzone),decay(nzone)
character (len=25) parami(npar),param,temp
character (len=80) bgrfl,zonefl
integer*1, allocatable:: izone1(:,:,:)
integer, allocatable:: izone4(:,:,:)
integer:: idim(3),iparam(npar)
integer:: iread_error
logical:: flexist
!
parami(1)='POROSITY'
parami(2)='RETARDATION'
parami(3)='LONGITUDINAL DISPERSIVITY'
parami(4)='TRANSVERSE DISPERSIVITY'
parami(5)='DIFFUSION COEFFICIENT'
parami(6)='DECAY'
jbeg=1;jend=ny;jstep=1
kbeg=1;kend=nz;kstep=1
if(irevz.eq.1)then
  jbeg=ny;jend=1;jstep=-1
  kbeg=nz;kend=1;kstep=-1
endif
if(inobgr.eq.1)then
!.parameter fields
  do i=1,nx; do j=1,ny; do k=1,nz
    izone=i+(j-1)*nx+(k-1)*nxy
    cat(i,j,k)%zone=izone
  enddo;enddo;enddo
elseif(nzone.eq.1)then
!.constant zone
  cat(:,:,:)%zone=1
else
! lets try several different formats here until isuccess = 1 
  isuccess=0
!.zones specified in bgr file
  read(inpar,'(a)',err=9998)bgrfl
  write(iout,2000)bgrfl,nzone
  write(*,2000)bgrfl,nzone
!.read bgr file            
  inquire (file=bgrfl, exist=flexist)
  if(.not.flexist)then
    write(*,*)' BGR file:',bgrfl
    write(*,*)' does not exist!'
    write(iout,*)' BGR file:',bgrfl
    write(iout,*)' does not exist!'
    goto 9998
  endif
  do iread=0,4
   if(iread.eq.0)then
! ASCII
     open(inbgr,file=bgrfl)
     write(iout,*)' Trying to read TSIM file as default ASCII'
     write(*,*)' Trying to read TSIM file as default ASCII'
   elseif(iread.eq.1.or.iread.eq.2)then
! UNFORMATTED BINARY
     open(inbgr,file=bgrfl,form='unformatted')
     write(iout,*)' Trying to read bgr file as unformatted'
     write(*,*)' Trying to read bgr file as unformatted'
     if(iread.eq.1)then
       write(iout,*)' with izone as integer*1'
       write(*,*)' with izone as integer*1'
     else
       write(iout,*)' with izone as integer*4'
       write(*,*)' with izone as integer*4'
     endif
   else
! BINARY
     ! open(inbgr,file=bgrfl,form='unformatted')
     open(inbgr,file=bgrfl,form='unformatted',access='stream')
     ! open(inbgr,file=bgrfl,form='unformatted')
     write(iout,*)' Trying to read bgr file as binary'
     write(*,*)' Trying to read bgr file as binary'
     if(iread.eq.3)then
       write(iout,*)' with izone as integer*1'
       write(*,*)' with izone as integer*1'
     else
       write(iout,*)' with izone as integer*4'
       write(*,*)' with izone as integer*4'
     endif
   endif    
! READ FIRST LINE
   if(iread.eq.0)then
     read(inbgr,*,iostat=istatus)idim(1),idim(2),idim(3)
     if(istatus.ne.0)then
       close(inbgr)
       write(iout,*)' Error reading bgr file value irank'
       write(*,*)' Error reading bgr file value irank'
       write(iout,*)' Trying new BGR file format '     
       write(*,*)' Trying new BGR file format '     
       cycle
     endif
   elseif(iread.gt.0)then
     read(inbgr,iostat=istatus) irank
     if(istatus.ne.0)then
       close(inbgr)
       write(iout,*)' Error reading bgr file value irank'
       write(*,*)' Error reading bgr file value irank'
       write(iout,*)' Trying new BGR file format '     
       write(*,*)' Trying new BGR file format '     
       cycle
     endif
     write(*,*)' BGR file rank = ',irank
     write(iout,*)' BGR file rank = ',irank
   endif
! READ SECOND LINE
   if(iread.gt.0)then
     read(inbgr,iostat=istatus) idim
     if(istatus.ne.0)then
       close(inbgr)
       write(iout,*)' Error reading bgr file value idim'
       write(*,*)' Error reading bgr file value idim'
       write(iout,*)' Trying new BGR file format '     
       write(*,*)' Trying new BGR file format '     
       cycle
     endif
   endif
   write(*,*)' TSIM file dimensions = ',idim
   write(iout,*)' TSIM file dimensions = ',idim
! CHECK DIMENSIONS
   if(nx.ne.idim(1).or.ny.ne.idim(2).or.nz.ne.idim(3))then
     close(inbgr)
       write(iout,*)' Error reading TSIM file dimensions'
       write(*,*)' Error reading TSIM file dimensions'
       write(iout,*)' Trying new TSIM file format '     
       write(*,*)' Trying new TSIM file format '     
     cycle
   endif
! READ ZONES
   if(iread.eq.0)then
     allocate (izone4(nx,ny,nz))
     read(inbgr,*,iostat=istatus)(((izone4(i,j,k),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
     if(istatus.ne.0)then
       close(inbgr)
       write(iout,*)' Error reading TSIM file indicators'
       write(*,*)' Error reading TSIM file indicators'
       write(iout,*)' Trying new TSIM file format '     
       write(*,*)' Trying new TSIM file format '     
       cycle
     endif
   elseif(iread.eq.1.or.iread.eq.3)then
     allocate (izone1(nx,ny,nz))
     read(inbgr,iostat=istatus)(((izone1(i,j,k),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
     if(istatus.ne.0)then
       close(inbgr)
       write(iout,*)' Error reading bgr file indicators'
       write(*,*)' Error reading bgr file indicators'
       write(iout,*)' Trying new BGR file format '     
       write(*,*)' Trying new BGR file format '     
       cycle
     endif
   else
     allocate (izone4(nx,ny,nz))
     read(inbgr,iostat=istatus)(((izone4(i,j,k),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
     if(istatus.ne.0)then
       close(inbgr)
       write(iout,*)' Error reading bgr file indicators'
       write(*,*)' Error reading bgr file indicators'
       write(iout,*)' Trying new BGR file format '     
       write(*,*)' Trying new BGR file format '     
       cycle
     endif
   endif
   isuccess=1
   do i=1,nx;do j=1,ny;do k=1,nz
     if(iread.eq.1.or.iread.eq.3)then
       ic=iabs(1*izone1(i,j,k))
     else
       ic=iabs(izone4(i,j,k))
     endif
     if(ic.lt.1.or.ic.gt.nzone)then
       write(iout,*)' Error TSIM file value at i,j,k ',i,j,k,' is ',ic
       write(*,*)' Error TSIM file value at i,j,k ',i,j,k,' is ',ic
       write(iout,*)' Number of zones  = ',nzone
       write(*,*)' Number of zones  = ',nzone
       write(iout,*)' Trying new TSIM file format '     
       write(*,*)' Trying new TSIM file format '     
       close(inbgr)
       isuccess=0
       exit
     endif
     cat(i,j,k)%zone=ic
   enddo
   if(isuccess.eq.0)exit
   enddo
   if(isuccess.eq.0)exit
   enddo
   close(inbgr)
   if(isuccess.eq.1)exit
  enddo      
  if(isuccess.eq.0)then
    write(iout,*)' Error in bgr file'
    write(*,*)' Error in bgr file'
    stop ' Exiting' 
  endif
endif        
!.read zone data for each parameter
do ipar=1,npar
  iparam(ipar)=0
enddo
do
  call skip(inpar)
  read(inpar,'(a)',iostat=iread_error)param
  if(iread_error.eq.-1)exit
  if(iread_error.eq.-2)goto 9998
  do ipar=1,npar
    temp=parami(ipar)
    if(param(1:4).eq.temp(1:4))then
      iparam(ipar)=1
      if(inobgr.eq.1)then
        read(inpar,*,iostat=iread_error)iconstant
        if(iconstant.eq.0)then ! constant parameter
          backspace(inpar)
          read(inpar,*,iostat=iread_error)iconstant,val
          do i=1,nx; do j=1,ny; do k=1,nz
            izone=i+(j-1)*nx+(k-1)*nxy
            if(ipar.eq.1)por(izone)=val
            if(ipar.eq.2)ret(izone)=val
            if(ipar.eq.3)dlong(izone)=val
            if(ipar.eq.4)dtran(izone)=val
            if(ipar.eq.5)ddiff(izone)=val
            if(ipar.eq.6)decay(izone)=val
          enddo;enddo;enddo
          write(*,2003)temp,1,nzone,val
          write(iout,2003)temp,1,nzone,val
        elseif(iconstant.eq.1.or.iconstant.eq.-1)then
          backspace(inpar)
          read(inpar,*,iostat=iread_error)iconstant,zonefl
          inquire (file=zonefl, exist=flexist)
          if(.not.flexist)then
            write(*,*)' ZONE file:',zonefl
            write(*,*)' does not exist!'
            write(iout,*)' ZONE file:',zonefl
            write(iout,*)' does not exist!'
            goto 9998
          endif
          if(iconstant.eq.-1)then
            open(inbgr,file=zonefl,form='unformatted')
            if(ipar.eq.1)read(inbgr)(((por  (i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.2)read(inbgr)(((ret  (i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.3)read(inbgr)(((dlong(i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.4)read(inbgr)(((dtran(i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.5)read(inbgr)(((ddiff(i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.6)read(inbgr)(((decay(i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            write(*,2004)zonefl
            write(iout,2004)zonefl
          elseif(iconstant.eq.1)then
            open(inbgr,file=zonefl)
            do izone=1,nzone
            if(ipar.eq.1)read(inbgr,*)(((por  (i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.2)read(inbgr,*)(((ret  (i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.3)read(inbgr,*)(((dlong(i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.4)read(inbgr,*)(((dtran(i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.5)read(inbgr,*)(((ddiff(i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            if(ipar.eq.6)read(inbgr,*)(((decay(i+(j-1)*nx+(k-1)*nxy),i=1,nx),j=jbeg,jend,jstep),k=kbeg,kend,kstep)
            enddo
            write(*,2005)zonefl
            write(iout,2005)zonefl
          endif
        else
          goto 9998
        endif
      else
        do iz=1,nzone
          read(inpar,*,iostat=iread_error)izone,val
          if(iread_error.ne.0)goto 9998      
          if(izone.ne.iz)goto 9998
          if(ipar.eq.1)por(izone)=val
          if(ipar.eq.2)ret(izone)=val
          if(ipar.eq.3)dlong(izone)=val
          if(ipar.eq.4)dtran(izone)=val
          if(ipar.eq.5)ddiff(izone)=val
          if(ipar.eq.6)decay(izone)=val
          if(iz.eq.1)write(*,2001)temp,izone,val
          if(iz.eq.1)write(iout,2001)temp,izone,val
          if(iz.ne.1)write(*,2002)izone,val
          if(iz.ne.1)write(iout,2002)izone,val
        enddo
      endif
      cycle
    endif
  enddo
enddo
close (inpar)
!.did we get all  of the parameters needed?
do ipar=1,npar
  if(iparam(ipar).eq.0)then 
    write(*,*)' Missing input for ',parami(ipar)
    write(iout,*)' Missing input for ',parami(ipar)
    goto 9998
  endif
enddo
!.check to see if simulation is isotropic and if advection only
isotropic=1
iadvect=1
do izone=1,nzone
  if(dtran(izone).ne.dlong(izone))isotropic=0
  if(dlong(izone).ne.0.0)iadvect=0
  if(dtran(izone).ne.0.0)iadvect=0
  if(ddiff(izone).ne.0.0)iadvect=0
enddo
! if anisotropic, make sure that dtran, dlong and ddiff are homogeneous
if(isotropic.eq.0)then
  do izone=1,nzone
  do izone2=1,nzone
    if(dtran(izone).ne.dtran(izone2))then
      write(*,*)' Transverse dispersivity is heterogeneous:'//&
      ' option supported only for isotropic D'
      write(iout,*)' Transverse dispersivity is heterogeneous:'//&
      ' option supported only for isotropic D'
      goto 9998 
    elseif(dlong(izone).ne.dlong(izone2))then
      write(*,*)' Longitudinal dispersivity is heterogeneous:'//&
      ' option supported only for isotropic D'
      write(iout,*)' Longitudinal dispersivity is heterogeneous:'//&
      ' option supported only for isotropic D'
      goto 9998 
      elseif(ddiff(izone).ne.ddiff(izone2))then 
      write(*,*)' Diffusion Coef. is heterogeneous: '//&
      'option supported only for isotropic D'
      write(iout,*)' Diffusion Coef. is heterogeneous: '//&
      'option supported only for isotropic D'
      goto 9998 
    endif
  enddo
  enddo
endif
!write(*,*)' RETARDATION IS NOT ACTIVE IN THIS VERSION'
!write(IOUT,*)' RETARDATION IS NOT ACTIVE IN THIS VERSION'
if(nzone.ne.1.and.inobgr.eq.0)deallocate (izone1)
! write (*,*) 'nzone = ', nzone
! write (*,*) 'inobgr = ', inobgr
return
9998 stop ' Error in parameter input file!'
9999 stop ' Error in BGR file specifying zones!'
2000 format('                           C O E F F I C I E N T S   '/&
            ' bgr file         :',3x,a/&
            ' number of zones                 ',13('.'),8x,i15///&
            ' COEFFICIENT                     ',7(' '),8x,&
            'ZONE #             VALUE'/)
2001 format(a33,7('.'),7x,i3,8x,e15.8)
2002 format(38(' '),6x,i3,8x,e15.8)
2003 format(a33,7('.'),7x,i3,' -',i8,3x,e15.8)
2004 format('Variable Parameter Field, Unformatted File = ',a80)
2005 format('Variable Parameter Field, Formatted File = ',a80)
end
!------------------------------------------------------------
! solve
!------------------------------------------------------------
subroutine solve(movep,sam,cat,pat,imp,rech,chd,vel3,por,ret,dlong,dtran,&
ddiff,decay,bounds,source,opc,outunit,outfname)
use global
implicit none
type (particle):: pat(1,maxnp)
type (imparticle):: imp(1,maxnp)
type (sample):: sam(nsam)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
type (recharge):: rech(nx,ny)
type (constant_head):: chd(nx,ny,nz)
integer:: imawupdt,i,j,k
integer:: source(maxsource),opc(nopc),outunit(nopc+1)
real:: sngldt
real:: vel3(3,0:nx,0:ny,0:nz)
real:: por(nzone),ret(nzone),dlong(nzone),dtran(nzone),ddiff(nzone),decay(nzone)
double precision:: tnextvelo
character (len=80) outfname(nopc)
!.external subroutine to move particles
external movep
!.time loop
dtminimum=dtinit
write(*,*)'  Current Time   # of Particles '
write(*,*)' --------------  -------------- '
do; if(.not.(curtime.lt.tmax))exit
  write(*,'(1x,f15.4,1x,i15,a1,$)')sngl(curtime),np,char(13)
!.calculate dt based on min. time to update velocity, boundaries, output, or tmax
  dt=dtinit
! DAB also need dt based on reaction time ...
  dt=min(dt,tnextvel-curtime,tnextbnd-curtime,tnextopc-curtime,&
  tnextpnt-curtime,tnextsrc-curtime,tmax-curtime)
  curtime=curtime+dt
  sngldt=dt
!.move all existing particles up to curtime
  if(ibug.ge.1)then;string(1)=' CALL MOVEP!';call debug(ibugout,string,nline);endif
  call movep(1,np,sngldt,vel3,por,ret,dlong,dtran,ddiff,decay,pat,cat,bounds)  ! full anisotropic displacement
!.distribute and move particles from continuous source boundary sources
  if(ibug.ge.1)then;string(1)=' CALL SRC!';call debug(ibugout,string,nline);endif
  call src(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,ddiff,decay,bounds,source,sngldt)
!.remove tagged particles that left domain
  if(ibug.ge.1)then;string(1)=' CALL REMOVEP!';call debug(ibugout,string,nline);endif
  call removep(1,pat,cat)
!.check for mass conservation (particle mass should sum to cell-based mass)
  if(ibug.ge.2)call conserve(pat,cat)
!.split particles for resolution of concentration field
  if(ibug.ge.1)then;string(1)=' CALL SPLIT!';call debug(ibugout,string,nline);endif
  if(np.ne.0)call split(pat,cat,bounds,por,ret)
! DAB put in mixing subroutine call here

    nactive = count(pat%active)
    allocate (alive(nactive)) ! maybe preallocate to avoid repeatedly doing this
    alive = pack(indices, mparts%active)
    mparts(alive)%bin = floor(mparts(alive)%loc/dxv) + 1

    call kdtree(

    call mix_particles(imp, pat, nactive, alive, nipart, D, dt, omega)

    call abc(imp,pat)

! DAB put in reaction subroutine here 

!.update velocities
  tnextvelo=tnextvel
  imawupdt=0
  if(curtime.eq.tnextvel)then
    imawupdt=1
    if(ibug.ge.1)then;string(1)=' CALL VELUPDT!';call debug(ibugout,string,nline);endif
    call velupdt(vel3,por,cat,rech,chd)
    if(ibug.ge.1)then;string(1)=' CALL TCNTRL!';call debug(ibugout,string,nline);endif
    if(dtcntrl.ne.0.0)call tcntrl(vel3,por,ret,dtran,dlong,ddiff,cat)
    write(*,*)'  Current Time   # of Particles '
    write(*,*)' --------------  -------------- '
  endif
!.update BCs
  if((curtime.eq.tnextbnd).and.ibug.ge.1)then
    string(1)=' CALL BNDUPT!';call debug(ibugout,string,nline)
  endif
  if(curtime.eq.tnextbnd)call bndupdt(bounds,source,cat,pat)
  if((curtime.eq.tnextbnd.or.curtime.eq.tnextvelo).and.ibug.ge.1)then
    string(1)=' CALL COURANT!';call debug(ibugout,string,nline)
  endif
  if(curtime.eq.tnextbnd.or.curtime.eq.tnextvelo)&   
  call courant(source,bounds,vel3,cat,por)
!.distribute particles from sources
  if((curtime.eq.tnextpnt).and.ibug.ge.1)then
    string(1)=' CALL PNTUPDT!';call debug(ibugout,string,nline)
  endif
  if(curtime.eq.tnextpnt)call pntupdt(pat,imp,cat,vel3)
!.update sources
  if((curtime.eq.tnextsrc).and.ibug.ge.1)then
    string(1)=' CALL SRCUPDT!';call debug(ibugout,string,nline)
  endif
  if(curtime.eq.tnextsrc)call srcupdt(source,bounds,pat,cat)
!.output results
  if((curtime.eq.tnextopc).and.ibug.ge.1)then
    string(1)=' CALL OUTPUT!';call debug(ibugout,string,nline)
  endif
  if(curtime.eq.tnextopc)call output(opc,outunit,outfname,sam,pat,cat,vel3,bounds,por,ret,decay)
  if(imawupdt.eq.1)then
    if(ibug.ge.1)then;string(1)=' CALL MAWUPDT!';call debug(ibugout,string,nline);endif
    call mawupdt(bounds,vel3)
  endif
!.exit if no particles left
  if(np.eq.0)exit
!.end time loop
enddo
return
write(*,'(f15.4,f15.4,i10,a1,$)')sngl(curtime),dt,np,char(13)
return
end
!------------------------------------------------------------
! src
!------------------------------------------------------------
subroutine src(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,&
ddiff,decay,bounds,source,sngldt)
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
type (recharge):: rech(nx,ny)
type (constant_head):: chd(nx,ny,nz)
real sngldt,srctime
integer:: ip,isource,ibounds,itype,npart,iptype,nparttotal,i,j,k,kbot,ktop,nw,iw,iside
integer:: kx,ky,kz,im,jm,km,jp,kp,ijk,ierr
integer:: source(maxsource)
integer:: itrelease,is,nside
real:: dtsrc,nptime,sx,sy,sz,xm,ym,zm,qtotal,q,r,randu01,xmbc,ymbc,zmbc,dxbc,dybc,dzbc
real:: vel3(3,0:nx,0:ny,0:nz),vx4,vy4,vz4,vxyz,v4,u,ux,uy,uz,x,y,z,c(nspec),flow
real:: por(nzone),ret(nzone),dlong(nzone),dtran(nzone),ddiff(nzone),decay(nzone)
real, allocatable:: probc(:),probw(:)
integer, allocatable:: loc(:,:)
double precision pmass(nspec),massrate(nspec),dtsource,trelease,tinitsrc,tnextvelo,pmass2
character (len=80) outfname(nopc)
!.external subroutine to move particles
external movep

  do isource=1,nsource
    ibounds=source(isource)
    itype=bounds(ibounds)%bc_type
    if(itype.eq.2)then
      massrate=abs(bounds(ibounds)%massrate)
      pmass=bounds(ibounds)%pmass
      dtsource=maxval(pmass/massrate)                   ! DAB needed a scalar here! (should all be equal?)
      tinitsrc=curtime-dt ! initial time
      trelease=0.0        ! release time
      itrelease=1
      xm=bounds(ibounds)%bc_xyzm(1)
      ym=bounds(ibounds)%bc_xyzm(2)
      zm=bounds(ibounds)%bc_xyzm(3)
      sx=bounds(ibounds)%bc_xyzs(1)
      sy=bounds(ibounds)%bc_xyzs(2)
      sz=bounds(ibounds)%bc_xyzs(3)
      do; if(trelease.gt.sngldt)exit
        if(dt-trelease.lt.dtsource)pmass=pmass*((dt-trelease)/dtsource) ! distribute remaining fractional mass
        if(maxval(pmass).le.0.0)exit
        mass=mass+pmass     
        bounds(ibounds)%np_remove=bounds(ibounds)%np_remove+1
        bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove+ pmass
        npbc(itype)=npbc(itype)+1
        do iloop=1,nspec
           massbc(itype,iloop)=massbc(itype,iloop)+pmass(iloop)
        enddo
        srctime=sngl(tinitsrc+trelease)
        dtsrc=sngl(dt-trelease)
        iptype=0
! DAB not exactly sure why any of massrate would be (-).  How to break this up?
        if(minval(bounds(ibounds)%massrate).gt.0.)call placeu(srctime,xm,ym,zm,sx,sy,sz,pat,cat,1,pmass,iptype,ierr)
        if(maxval(bounds(ibounds)%massrate).lt.0.)call placeuf(srctime,xm,ym,zm,sx,sy,sz,pat,cat,vel3,1,pmass,iptype)
        call movep(np,np,dtsrc,vel3,por,ret,dlong,dtran,ddiff,decay,pat,cat,bounds)  ! full anisotropic displacement
        trelease=itrelease*dtsource
        itrelease=itrelease+1
      enddo
    elseif(itype.eq.6)then
      massrate=abs(bounds(ibounds)%massrate)
      pmass=bounds(ibounds)%pmass
      i=bounds(ibounds)%ijk(1)
      xm=(i-1)*dx
      j=bounds(ibounds)%ijk(2)
      ym=(j-1)*dy
      ktop=bounds(ibounds)%ktop
      kbot=bounds(ibounds)%kbot
      nw=ktop-kbot+1
      if(.not.allocated(probw))allocate(probw(nw))
      if(.not.allocated(probc))allocate(probc(2*5+(nw-2)*4))
      if(.not.allocated(loc))allocate(loc(2,2*5+(nw-2)*4))
      qtotal=0.
      probc(:)=0.
      iside=0
      do k=kbot,ktop
          iside=iside+1
          if(iside.ne.1)probc(iside)=probc(iside-1)
          if(bke*vel3(1,i-1,j,k).lt.0.)then
            qtotal=qtotal+abs(bke*vel3(1,i-1,j,k)*ax)
            probc(iside)=probc(iside)+abs(bke*vel3(1,i-1,j,k)*ax)
            loc(1,iside)=1
            loc(2,iside)=k
          endif
          iside=iside+1
          probc(iside)=probc(iside-1)
          if(bke*vel3(1,i,j,k).gt.0.)then
            qtotal=qtotal+abs(bke*vel3(1,i,j,k)*ax)
            probc(iside)=probc(iside)+abs(bke*vel3(1,i,j,k)*ax)
            loc(1,iside)=2
            loc(2,iside)=k
          endif
          iside=iside+1
          probc(iside)=probc(iside-1)
          if(bke*vel3(2,i,j-1,k).lt.0.)then
            qtotal=qtotal+abs(bke*vel3(2,i,j-1,k)*ay)
            probc(iside)=probc(iside)+abs(bke*vel3(2,i,j-1,k)*ay)
            loc(1,iside)=3
            loc(2,iside)=k
          endif
          iside=iside+1
          probc(iside)=probc(iside-1)
          if(bke*vel3(2,i,j,k).gt.0.)then
            qtotal=qtotal+abs(bke*vel3(2,i,j,k)*ay)
            probc(iside)=probc(iside)+abs(bke*vel3(2,i,j,k)*ay)
            loc(1,iside)=4
            loc(2,iside)=k
          endif
          if(k.eq.kbot)then
            iside=iside+1
            probc(iside)=probc(iside-1)
            if(bke*vel3(3,i,j,k-1).lt.0.)then
              qtotal=qtotal+abs(bke*vel3(3,i,j,k-1)*az)
              probc(iside)=probc(iside)+abs(bke*vel3(3,i,j,k-1)*az)
              loc(1,iside)=5
              loc(2,iside)=k
            endif
          endif
          if(k.eq.ktop)then
            iside=iside+1
            probc(iside)=probc(iside-1)
            if(bke*vel3(3,i,j,k).gt.0.)then
              qtotal=qtotal+abs(bke*vel3(3,i,j,k)*az)
              probc(iside)=probc(iside)+abs(bke*vel3(3,i,j,k)*az)
              loc(1,iside)=6
              loc(2,iside)=k
            endif
          endif
      enddo
      nside=iside
      if(qtotal.ne.0.)then
        probc=probc/qtotal
        nside=iside
        dtsource=maxval(pmass/massrate)  ! DAB not sure about maxval - some pmass might be zero? 
        tinitsrc=curtime-dt ! initial time
        trelease=0.0        ! release time
        itrelease=1
        iptype=0
        sz=dz
        do; if(trelease.gt.sngldt)exit
          if(dt-trelease.lt.dtsource)pmass=pmass*((dt-trelease)/dtsource) ! distribute remaining fractional mass
          if(maxval(pmass).le.0.0)exit
          mass=mass+pmass
          bounds(ibounds)%np_remove=bounds(ibounds)%np_remove+1
          bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove+ pmass
          npbc(itype)=npbc(itype)+1
         do iloop=1,nspec
            massbc(itype,iloop)=massbc(itype,iloop)+pmass(iloop)
         enddo
! old         massbc(itype)=massbc(itype)+pmass
          srctime=sngl(tinitsrc+trelease)
          dtsrc=sngl(dt-trelease)
          r=randu01()
          do is=1,nside
            if(r.le.probc(is))exit
          enddo
! Find k and iside
          k=loc(2,is) 
          iside=loc(1,is) 
! Place particles and move
          zm=(k-1)*dz
          if(iside.eq.1)then
            xmbc=xm; ymbc=ym+bcdxyz*dy; zmbc=zm+bcdxyz*dz
            dxbc=bcdxyz*dx; dybc=dy-2*bcdxyz*dy; dzbc=dz-2*bcdxyz*dz
          elseif(iside.eq.2)then
            xmbc=xm+dx-bcdxyz*dx; ymbc=ym+bcdxyz*dy; zmbc=zm+bcdxyz*dz
            dxbc=bcdxyz*dx; dybc=dy-2*bcdxyz*dy; dzbc=dz-2*bcdxyz*dz
          elseif(iside.eq.3)then
            xmbc=xm+bcdxyz*dx; ymbc=ym; zmbc=zm+bcdxyz*dz
            dxbc=dx-2*bcdxyz*dx; dybc=bcdxyz*dy; dzbc=dz-2*bcdxyz*dz
          elseif(iside.eq.4)then
            xmbc=xm+bcdxyz*dx; ymbc=ym+dy-bcdxyz*dy; zmbc=zm+bcdxyz*dz
            dxbc=dx-2*bcdxyz*dx; dybc=bcdxyz*dy; dzbc=dz-2*bcdxyz*dz
          elseif(iside.eq.5)then
           xmbc=xm+bcdxyz*dx; ymbc=ym+bcdxyz*dy; zmbc=zm
           dxbc=dx-2*bcdxyz*dx; dybc=dy-2*bcdxyz*dy; dzbc=bcdxyz*dz
          elseif(iside.eq.6)then
           xmbc=xm+bcdxyz*dx; ymbc=ym+bcdxyz*dy; zmbc=zm+dz-bcdxyz*dz
           dxbc=dx-2*bcdxyz*dx; dybc=dy-2*bcdxyz*dy; dzbc=bcdxyz*dz
          endif
          call placeu(srctime,xmbc,ymbc,zmbc,dxbc,dybc,dzbc,pat,cat,1,pmass,iptype,ierr)
          call movep(np,np,dtsrc,vel3,por,ret,dlong,dtran,ddiff,decay,pat,cat,bounds) 
          trelease=itrelease*dtsource
          itrelease=itrelease+1
        enddo
      endif
      if(allocated(probw))deallocate(probw)
      if(allocated(probc))deallocate(probc)
      if(allocated(loc))deallocate(loc)
    elseif(itype.eq.8)then
      massrate=abs(bounds(ibounds)%massrate)
      pmass=bounds(ibounds)%pmass
      xm=bounds(ibounds)%bc_xyzm(1)
      ym=bounds(ibounds)%bc_xyzm(2)
      sx=bounds(ibounds)%bc_xyzs(1)
      sy=bounds(ibounds)%bc_xyzs(2)
      dtsource=maxval(pmass/massrate)   ! DAB not sure, maybe some pmass(nspec) = 0? 
      tinitsrc=curtime-dt ! initial time
      trelease=0.0        ! release time
      itrelease=1
      do; if(trelease.gt.sngldt)exit
        if(dt-trelease.lt.dtsource)pmass=pmass*((dt-trelease)/dtsource) ! distribute remaining fractional mass
        if(minval(pmass).le.0.0)exit      ! DAB check this logic
        mass=mass+pmass
        bounds(ibounds)%np_remove=bounds(ibounds)%np_remove+1
        bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove+ pmass
        npbc(itype)=npbc(itype)+1
        do iloop=1,nspec
           massbc(itype,iloop)=massbc(itype,iloop)+pmass(iloop)
        enddo
!        massbc(itype)=massbc(itype)+pmass
        srctime=sngl(tinitsrc+trelease)
        dtsrc=sngl(dt-trelease)
        iptype=0
        if(minval(bounds(ibounds)%massrate).gt.0.)&
        call placerch1(srctime,xm,ym,sx,sy,pat,cat,1,pmass,iptype,rech)
        call movep(np,np,dtsrc,vel3,por,ret,dlong,dtran,ddiff,decay,pat,cat,bounds)  ! full anisotropic displacement
        trelease=itrelease*dtsource
        itrelease=itrelease+1
      enddo
    elseif(itype.eq.9)then
      c=bounds(ibounds)%conc
      nptime=bounds(ibounds)%nptime ! particles per unit time (real) in each cell
      im=bounds(ibounds)%ijkm(1)
      jm=bounds(ibounds)%ijkm(2)
      ip=bounds(ibounds)%ijkp(1)
      jp=bounds(ibounds)%ijkp(2)
! cycle through all cells
      do i=im,ip;do j=jm,jp
        flow=bke*rech(i,j)%flow
        if(flow.le.0.0)cycle      
        massrate=dble(c)*dble(flow) ! massrate = c*Q(i,j) [M/T], where Q is flow rate
        pmass=massrate/dble(nptime)  ! pmass = massrate [M/T]/nptime [particles/T] in each cell
        dtsource=maxval(pmass/massrate)    ! DAB these should all be equal by definition?
        tinitsrc=curtime-dt ! initial time
        trelease=0.0        ! release time
        itrelease=1
        do; if(trelease.gt.sngldt)exit
          if(dt-trelease.lt.dtsource)pmass=pmass*(dt-trelease)/dtsource ! distribute remaining fractional mass
          if(minval(pmass).le.0.0)exit      ! DAB check this logic
          srctime=sngl(tinitsrc+trelease)
          dtsrc=sngl(dt-trelease)
          mass=mass+pmass
          bounds(ibounds)%np_remove=bounds(ibounds)%np_remove+1
          bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove+ pmass
          npbc(itype)=npbc(itype)+1
        do iloop=1,nspec
           massbc(itype,iloop)=massbc(itype,iloop)+pmass(iloop)
        enddo
!          massbc(itype)=massbc(itype)+pmass
          iptype=0
          call placerch2(srctime,i,j,pat,cat,rech,1,pmass,iptype)
          call movep(np,np,dtsrc,vel3,por,ret,dlong,dtran,ddiff,decay,pat,cat,bounds)  ! full anisotropic displacement
          trelease=itrelease*dtsource
          itrelease=itrelease+1
        enddo
      enddo;enddo
    elseif(itype.eq.10)then
      c=bounds(ibounds)%conc
      nptime=bounds(ibounds)%nptime ! particles per unit time (real) in each cell
      im=bounds(ibounds)%ijkm(1)
      jm=bounds(ibounds)%ijkm(2)
      km=bounds(ibounds)%ijkm(3)
      ip=bounds(ibounds)%ijkp(1)
      jp=bounds(ibounds)%ijkp(2)
      kp=bounds(ibounds)%ijkp(3)  
! cycle through all cells
      do k=km,kp;do j=jm,jp;do i=im,ip
        flow=chd(i,j,k)%flow
        if(bke*flow.le.0.0)cycle      
        massrate=dble(c)*dble(bke*flow) ! massrate = c*Q(i,j) [M/T], where Q is flow rate
        pmass=massrate/dble(nptime)  ! pmass = massrate [M/T]/nptime [particles/T] in each cell
        dtsource=maxval(pmass/massrate)          ! DAB still not sure about maxval 
        tinitsrc=curtime-dt ! initial time
        trelease=0.0        ! release time
        itrelease=1
		xm=(i-1)*dx;ym=(j-1)*dy;zm=(k-1)*dz
        do; if(trelease.gt.sngldt)exit
          if(dt-trelease.lt.dtsource)pmass=pmass*(dt-trelease)/dtsource ! distribute remaining fractional mass
          if(maxval(pmass).le.0.0d+0)exit
          srctime=sngl(tinitsrc+trelease)
          dtsrc=sngl(dt-trelease)
          iptype=0
          do
            call placeu(srctime,xm,ym,zm,dx,dy,dz,pat,cat,1,pmass,iptype,ierr)
            if(ierr.eq.0)exit
            write(iout,*)' WARNING: PLACEU ATTEMPTED TO PLACE A SRC TYPE 10 PARTICLE OUTSIDE OF THE DOMAIN'
          enddo
! particle was on edge or outside of domain; ignore and start over 
          mass=mass+pmass

          bounds(ibounds)%np_remove=bounds(ibounds)%np_remove+1
          bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove+ pmass
          npbc(itype)=npbc(itype)+1
        do iloop=1,nspec
           massbc(itype,iloop)=massbc(itype,iloop)+pmass(iloop)
        enddo
!          massbc(itype)=massbc(itype)+pmass
          call movep(np,np,dtsrc,vel3,por,ret,dlong,dtran,ddiff,decay,pat,cat,bounds)  ! full anisotropic displacement          
          trelease=itrelease*dtsource
          itrelease=itrelease+1
        enddo
      enddo;enddo;enddo
    elseif(itype.eq.11)then
      c=bounds(ibounds)%conc
      nptime=bounds(ibounds)%nptime ! particles per unit time (real) in each cell
      im=bounds(ibounds)%ijkm(1)
      jm=bounds(ibounds)%ijkm(2)
      km=bounds(ibounds)%ijkm(3)
      ip=bounds(ibounds)%ijkp(1)
      jp=bounds(ibounds)%ijkp(2)
      kp=bounds(ibounds)%ijkp(3)  
! cycle through all cells
      do k=km,kp;do j=jm,jp;do i=im,ip
        flow=0.0
        if(k.eq.km)then
          xm=(i-1)*dx;ym=(j-1)*dy;zm=(k-1)*dz
          sx=dx;sy=dy;sz=bcdxyz*dz
          flow=-bke*vel3(3,i,j,k-1)*az
		  call release_type11(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,&
		  ddiff,decay,bounds,source,sngldt,xm,ym,zm,sx,sy,sz,flow,massrate,nptime,c,ibounds,itype)
        endif           
        if(k.eq.kp)then
          xm=(i-1)*dx;ym=(j-1)*dy;zm=(k)*dz
          sx=dx;sy=dy;sz=bcdxyz*dz
          flow=+bke*vel3(3,i,j,k)*az
		  call release_type11(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,&
		  ddiff,decay,bounds,source,sngldt,xm,ym,zm,sx,sy,sz,flow,massrate,nptime,c,ibounds,itype)           
        endif
        if(j.eq.jm)then
          xm=(i-1)*dx;ym=(j-1)*dy;zm=(k-1)*dz
          sx=dx;sy=bcdxyz*dy;sz=dz
          flow=-bke*vel3(2,i,j-1,k)*ay
		  call release_type11(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,&
		  ddiff,decay,bounds,source,sngldt,xm,ym,zm,sx,sy,sz,flow,massrate,nptime,c,ibounds,itype)
        endif           
        if(j.eq.jp)then
          xm=(i-1)*dx;ym=(j)*dy;zm=(k-1)*dz
          sx=dx;sy=bcdxyz*dy;sz=dz
          flow=+bke*vel3(2,i,j,k)*ay  
		  call release_type11(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,&
		  ddiff,decay,bounds,source,sngldt,xm,ym,zm,sx,sy,sz,flow,massrate,nptime,c,ibounds,itype)         
        endif
        if(i.eq.im)then
          xm=(i-1)*dx;ym=(j-1)*dy;zm=(k-1)*dz
          sx=bcdxyz*dx;sy=dy;sz=dz
          flow=-bke*vel3(1,i-1,j,k)*ax
		  call release_type11(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,&
		  ddiff,decay,bounds,source,sngldt,xm,ym,zm,sx,sy,sz,flow,massrate,nptime,c,ibounds,itype)
        endif           
        if(i.eq.ip)then
          xm=(i)*dx;ym=(j-1)*dy;zm=(k-1)*dz
          sx=bcdxyz*dx;sy=dy;sz=dz
          flow=+bke*vel3(1,i,j,k)*ax
		  call release_type11(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,&
		  ddiff,decay,bounds,source,sngldt,xm,ym,zm,sx,sy,sz,flow,massrate,nptime,c,ibounds,itype)          
        endif
      enddo;enddo;enddo
    endif          
  enddo
return
end
!
subroutine release_type11(movep,cat,pat,rech,chd,vel3,por,ret,dlong,dtran,&
ddiff,decay,bounds,source,sngldt,xm,ym,zm,sx,sy,sz,flow,massrate,nptime,c,ibounds,itype)
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
type (recharge):: rech(nx,ny)
type (constant_head):: chd(nx,ny,nz)
real sngldt,srctime
integer:: ibounds,itype,iptype,ierr
integer:: source(maxsource)
integer:: itrelease
real:: dtsrc,nptime,sx,sy,sz,xm,ym,zm
real:: vel3(3,0:nx,0:ny,0:nz),c(nspec),flow
real:: por(nzone),ret(nzone),dlong(nzone),dtran(nzone),ddiff(nzone),decay(nzone)
real, allocatable:: probc(:),probw(:)
double precision pmass(nspec),massrate(nspec),dtsource,trelease,tinitsrc
!
if(flow.le.0.0)return
massrate=dble(c)*dble(flow) ! massrate = c*Q(i,j) [M/T], where Q is flow rate
pmass=massrate/dble(nptime)  ! pmass = massrate [M/T]/nptime [particles/T] in each cell
dtsource=maxval(pmass/massrate)   ! DAB unsure about maxval 
tinitsrc=curtime-dt ! initial time
trelease=0.0        ! release time
itrelease=1
do; if(trelease.gt.sngldt)exit
 if(dt-trelease.lt.dtsource)pmass=pmass*(dt-trelease)/dtsource ! distribute remaining fractional mass
 if(maxval(pmass).le.0.0d+0)exit
 srctime=sngl(tinitsrc+trelease)
 dtsrc=sngl(dt-trelease)
 iptype=0
 do
   call placeu(srctime,xm,ym,zm,sx,sy,sz,pat,cat,1,pmass,iptype,ierr)
   if(ierr.eq.0)exit
   write(iout,*)' WARNING: PLACEU ATTEMPTED TO PLACE A SRC TYPE 11 PARTICLE OUTSIDE OF THE DOMAIN'
 enddo
 ! particle was on edge or outside of domain; ignore and start over 
 mass=mass+pmass
 bounds(ibounds)%np_remove=bounds(ibounds)%np_remove+1
 bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove+ pmass
 npbc(itype)=npbc(itype)+1
   do iloop=1,nspec
      massbc(itype,iloop)=massbc(itype,iloop)+pmass(iloop)
   enddo
! massbc(itype)=massbc(itype)+pmass
 call movep(np,np,dtsrc,vel3,por,ret,dlong,dtran,ddiff,decay,pat,cat,bounds)  ! full anisotropic displacement          
 trelease=itrelease*dtsource
 itrelease=itrelease+1
enddo
return
end
!------------------------------------------------------------------      
! placerch1
!------------------------------------------------------------------
subroutine placerch1(time,xm,ym,sx,sy,pat,cat,npart,pmass,iptype,rech)
! place npart particles uniformly in a block starting at xyzm and ending at sxyz
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (recharge):: rech(nx,ny)
!integer:: maxnp,np,npart,iptype
integer:: k,kx,ky,kz,ip,npart,iptype
real:: x,y,z,xm,ym,zm,sx,sy,sz,time,randu01
double precision:: pmass(nspec)
do ip=1,npart
  if(np+1.gt.maxnp)then
    write(iout,*)' error:too many particles, increase maxnp'
    stop    ' error:too many particles, increase maxnp'
  endif
!
  do ! loop until flow is positive
    x=xm+sx*randu01(); y=ym+sy*randu01() !; z=zm+sz*randu01()
    kx=ifix(x/dx); ky=ifix(y/dy) !; kz=ifix(z/dz)
! determine kz and z from recharge array
    if(rech(kx+1,ky+1)%flow.gt.0.0)exit ! when we have found a positive flow region      
  enddo
  k=rech(kx+1,ky+1)%k
  if(irevz.eq.1)then
    z=float(k)*dz-smallxyz(3)
  else
    z=float(k-1)*dz+smallxyz(3)
  endif
  call addp(time,x,y,z,kx+1,ky+1,k,pmass,pat,cat)
enddo
return
end
!------------------------------------------------------------------      
! placerch2
!------------------------------------------------------------------
subroutine placerch2(time,i,j,pat,cat,rech,npart,pmass,iptype)
! place npart particles uniformly in a cell i,j,k with k determined 
! from the recharge array
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (recharge):: rech(nx,ny)
integer:: i,j,k,ip,npart,iptype
real:: x,y,z,xm,ym,sx,sy,time,randu01
double precision:: pmass(nspec)
do ip=1,npart
    if(np+1.gt.maxnp)then
      write(iout,*)' error:too many particles, increase maxnp'
      stop    ' error:too many particles, increase maxnp'
    endif
    xm=float(i-1)*dx;ym=float(j-1)*dy
    x=xm+dx*randu01(); y=ym+dy*randu01() !; z=zm+sz*randu01()
! determine k and z from recharge array; add particles at upper cell boundary
    k=rech(i,j)%k
    if(irevz.eq.1)then
      z=float(k)*dz-smallxyz(3)
    else
      z=float(k-1)*dz+smallxyz(3)
    endif
    call addp(time,x,y,z,i,j,k,pmass,pat,cat)
enddo
return
end
!------------------------------------------------------------
! movep1
!------------------------------------------------------------
subroutine movep1(ip0,ipn,dtcurrent,vel3,por,ret,&
dlong,dtran,ddiff,decay,pat,cat,bounds)
!.....analytical stream lines on a block-centered finite difference solution
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(*)
integer:: iflag,ip,kx,ky,kz
integer:: ip0,ipn,i,j,k,io,jo,ko
real:: por(*),ret(*)
real:: dlong(*),dtran(*),ddiff(*),decay(nzone) ! included for consistent input to all movep subroutines 
real:: vel3(3,0:nx,0:ny,0:nz)
real:: dtcurrent,dtmin
double precision:: timeleft,totaltime
real:: x,y,z,xm,ym,zm,xp,yp,zp,xold,yold,zold,th4,ret4,retth
real:: xstream,ystream,zstream
data iflag/1/
do ip=ip0,ipn
  if(.not.pat(1,ip)%active)cycle
!
  timeleft=dble(dtcurrent)
  totaltime=0.0
  do; if(.not.sngl(timeleft).gt.0.0)exit
!...particle location
    x=pat(1,ip)%xyz(1); y=pat(1,ip)%xyz(2); z=pat(1,ip)%xyz(3) 
!
    i=pat(1,ip)%ijk(1); j=pat(1,ip)%ijk(2); k=pat(1,ip)%ijk(3) 
    io=pat(1,ip)%ijk(1); jo=pat(1,ip)%ijk(2); ko=pat(1,ip)%ijk(3) 
    kx=i-1; ky=j-1; kz=k-1
!...parameters
    th4=por(cat(i,j,k)%zone)
    ret4=ret(cat(i,j,k)%zone)
    retth=ret4*th4
!...subtime step
    dtmin=min(sngl(timeleft),dtcurrent)
!...streamlines
    call stream(vel3,x,y,z,retth,kx,ky,kz,i,j,k,dtmin,xstream,ystream,zstream)
    totaltime=totaltime+dble(dtmin)
    timeleft=dble(dtcurrent)-totaltime
    x=xstream
    y=ystream
    z=zstream
    xold=x
    yold=y
    zold=z
!...decay m=m0*exp(k*dt)
!  DAB either handle here with nspec decay constants or do separately in reaction subroutine
!  DAB I opt for the latter at this point.
    if(decay(cat(i,j,k)%zone).ne.0.0)&
    pat(1,ip)%pmass=pat(1,ip)%pmass*dexp(dble(decay(cat(i,j,k)%zone)*dtmin))
!...reflect particles at boundaries as needed or tag to remove
    call reflect(sngl(timeleft),pat,cat,x,y,z,ip)
    if(.not.pat(1,ip)%active)exit
!...stream returns new cell location, 
!...compute new cell location only if reflected
    if(x.ne.xold.or.y.ne.yold.or.z.ne.zold)then
      kx=ifix(x/dx); ky=ifix(y/dy); kz=ifix(z/dz)
      i=kx+1; j=ky+1; k=kz+1
      xm=kx*dx; ym=ky*dy; zm=kz*dz; xp=xm+dx; yp=ym+dy; zp=zm+dz
      call icell_correct(kx,ky,kz,i,j,k,x,y,z,xm,ym,zm,xp,yp,zp,vel3)
    endif
!...absorbing boundaries
    if(cat(i,j,k)%bc_number.ne.0)call absorb(pat,cat,bounds,vel3,por,ret,&
    kx,ky,kz,i,j,k,ip,dtmin,sngl(timeleft),x,y,z)
!...update particle and cell attributes
    call updtp(io,jo,ko,i,j,k,ip,x,y,z,pat,cat)
    if(.not.pat(1,ip)%active)exit
  enddo
enddo
return
end   
!------------------------------------------------------------
! movep2
!------------------------------------------------------------
subroutine movep2(ip0,ipn,dtcurrent,vel3,por,ret,&
dlong,dtran,ddiff,decay,pat,cat,bounds)
! variable porosity, isotropic dispersion displacement
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(*)
integer:: ip,kx,ky,kz,ii,jj,kk,kxt,kyt,kzt
integer:: ip0,ipn,l1,l2,l3,i,j,k,it,jt,kt,io,jo,ko
real:: vel3(3,0:nx,0:ny,0:nz)
real:: por(*),ret(*),dlong(*),dtran(*),ddiff(*),decay(nzone)
real:: dtcurrent,dtmin
double precision:: timeleft,totaltime
real:: x,y,z,xo,yo,zo,xd,yd,zd,xt,yt,zt,xm,ym,zm,xp,yp,zp,th4,ret4,retth,sqrtretth
real:: w1,w2,w3,vx4,vy4,vz4,v4,vx4t,vy4t,vz4t,v4t,th4t
real:: vxm,vym,vzm,vxp,vyp,vzp
real:: at,diff,ratv4d,ratv4dr2dt,dxd,dyd,dzd,r2dt
real:: xstream,ystream,zstream
!
do ip=ip0,ipn
  if(.not.pat(1,ip)%active)cycle
        
!
  timeleft=dble(dtcurrent)
  totaltime=0.0d+0
  do; if(.not.sngl(timeleft).gt.0.0)exit
    call random(w1,w2,w3)
!...particle location
    x=pat(1,ip)%xyz(1); y=pat(1,ip)%xyz(2); z=pat(1,ip)%xyz(3) 
    xo=pat(1,ip)%xyz(1); yo=pat(1,ip)%xyz(2); zo=pat(1,ip)%xyz(3) 
!
    i=pat(1,ip)%ijk(1); j=pat(1,ip)%ijk(2); k=pat(1,ip)%ijk(3) 
    io=pat(1,ip)%ijk(1); jo=pat(1,ip)%ijk(2); ko=pat(1,ip)%ijk(3) 
    kx=i-1; ky=j-1; kz=k-1
!...interpolate velocity          
    vx4=vel3(1,kx,j,k)+(vel3(1,i,j,k)-vel3(1,kx,j,k))*(x/dx-float(kx))
    vy4=vel3(2,i,ky,k)+(vel3(2,i,j,k)-vel3(2,i,ky,k))*(y/dy-float(ky))
    vz4=vel3(3,i,j,kz)+(vel3(3,i,j,k)-vel3(3,i,j,kz))*(z/dz-float(kz))
    v4=sqrt(vx4*vx4+vy4*vy4+vz4*vz4)
!...parameters
    th4=por(cat(i,j,k)%zone)
    at=dtran(cat(i,j,k)%zone)
    diff=ddiff(cat(i,j,k)%zone)
    ret4=ret(cat(i,j,k)%zone)
    retth=ret4*th4
    sqrtretth=sqrt(ret4*th4)
!...time-step control
    if(dtcntrl.ne.0.)then
      ii=mx*(mod(ifix(x/dx2),2)-min(kx,1))
      jj=my*(mod(ifix(y/dy2),2)-min(ky,1))
      kk=mz*(mod(ifix(z/dz2),2)-min(kz,1))
      dtmin=min(large,cat(i+ii,j+jj,k+kk)%tc,sngl(timeleft)) 
    else
      dtmin=dtcurrent
    endif
!...streamlines
    call stream(vel3,x,y,z,retth,kx,ky,kz,i,j,k,dtmin,xstream,ystream,zstream)
    totaltime=totaltime+dble(dtmin)
    timeleft=dble(dtcurrent)-totaltime
!...diagonal-tensor displacement
    ratv4d=sqrt(at*v4/retth+diff)
    r2dt=sqrt(2.0*dtmin)
    ratv4dr2dt=ratv4d*r2dt
    xt=pat(1,ip)%xyz(1)+mx*ratv4dr2dt*w1
    yt=pat(1,ip)%xyz(2)+my*ratv4dr2dt*w2
    zt=pat(1,ip)%xyz(3)+mz*ratv4dr2dt*w3
!...reflect particles at boundaries as needed or tag to remove
    call reflect(sngl(timeleft),pat,cat,xt,yt,zt,ip)     
    if(.not.pat(1,ip)%active)exit
!...compute new cell location from xt,yt,zt 
    kxt=ifix(xt/dx); kyt=ifix(yt/dy); kzt=ifix(zt/dz)
    it=kxt+1; jt=kyt+1; kt=kzt+1
!...interpolate velocity
    vx4t=vel3(1,kxt,jt,kt)+(vel3(1,it,jt,kt)-vel3(1,kxt,jt,kt))*(xt/dx-float(kxt))
    vy4t=vel3(2,it,kyt,kt)+(vel3(2,it,jt,kt)-vel3(2,it,kyt,kt))*(yt/dy-float(kyt))
    vz4t=vel3(3,it,jt,kzt)+(vel3(3,it,jt,kt)-vel3(3,it,jt,kzt))*(zt/dz-float(kzt))
    v4t=sqrt(vx4t*vx4t+vy4t*vy4t+vz4t*vz4t)
!...parameters at new location
    th4t=por(cat(it,jt,kt)%zone)
    at=dtran(cat(it,jt,kt)%zone)
    diff=ddiff(cat(it,jt,kt)%zone)
!...
    ratv4dr2dt= r2dt*sqrt(at*v4t+diff*th4t)/sqrt(th4)
!print*,ip,dt,diff,ratv4dr2dt
    x=xstream+mx*ratv4dr2dt*w1
    y=ystream+my*ratv4dr2dt*w2
    z=zstream+mz*ratv4dr2dt*w3
!...decay m=m0*exp(k*dt)
    if(decay(cat(i,j,k)%zone).ne.0.0)&
    pat(1,ip)%pmass=pat(1,ip)%pmass*dexp(dble(decay(cat(i,j,k)%zone)*dtmin))
!...reflect particles at boundaries as needed or tag to remove
    call reflect(sngl(timeleft),pat,cat,x,y,z,ip)
    if(.not.pat(1,ip)%active)exit
!...compute cell location
    kx=ifix(x/dx); ky=ifix(y/dy); kz=ifix(z/dz)
    i=kx+1; j=ky+1; k=kz+1
    xd=x-xo;yd=y-yo;zd=z-zo
!    if(dtmin.lt.near)then   
    if(sqrt(xd*xd+yd*yd+zd*zd).lt.near)then
      xm=kx*dx; ym=ky*dy; zm=kz*dz; xp=xm+dx; yp=ym+dy; zp=zm+dz
      call icell_correct(kx,ky,kz,i,j,k,x,y,z,xm,ym,zm,xp,yp,zp,vel3)
    endif
!...absorbing boundaries
    if(cat(i,j,k)%bc_number.ne.0)&
    call absorb(pat,cat,bounds,vel3,por,ret,kx,ky,kz,i,j,k,ip,dtmin,sngl(timeleft),x,y,z)
!...update particle and cell attributes
    call updtp(io,jo,ko,i,j,k,ip,x,y,z,pat,cat)
    if(.not.pat(1,ip)%active)exit
  enddo
enddo
return
end   
!------------------------------------------------------------
! movep3
!------------------------------------------------------------
subroutine movep3(ip0,ipn,dtcurrent,vel3,por,ret,&
dlong,dtran,ddiff,decay,pat,cat,bounds)
! variable porosity, isotropic dispersion displacement
! with bilinear interpolation of the dispersion tensor
!
use global
implicit none
!
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
integer:: ip,ip0,ipn
integer:: kx,ky,kz,i,j,k,ii,jj,kk,iii,jjj,kkk,imm,jmm,kmm
integer:: immii,jmmjj,kmmkk,ith,io,jo,ko,ix,iy,iz
real:: por(nzone),ret(nzone),dlong(nzone),dtran(nzone),ddiff(nzone),decay(nzone)
real:: vel3(3,0:nx,0:ny,0:nz)
real:: th(8),rt(8),vx(8),vy(8),vz(8)
!
real:: dtcurrent,dtmin
double precision:: timeleft,totaltime
real:: x,y,z,xo,yo,zo,xd,yd,zd,xm,ym,zm,xp,yp,zp,x1,y1,z1,x1n,y1n,z1n
real:: xstream,ystream,zstream
real:: fx,fy,fz,fx1,fy1,fz1,fxn,fyn,fzn,fx1n,fy1n,fz1n
real:: at,diff,r2dt,thi,th4,alt,al,retth,sqrtretth,reti,ret4,dtminreti
real:: t1,t2,t3,t4,t5,t6,t7,r1,r2,r3,r4,r5,r6,r7
real:: vx1,vx2,vx3,vx4,vx5,vx6,vx7,vy1,vy2,vy3,vy4,vy5,vy6,vy7
real:: vz1,vz2,vz3,vz4,vz5,vz6,vz7,v4,vv4
real:: e1,e2,e3,w1,w2,w3
real:: dispx,dispy,dispz,dzvz,dxvx,dyvy,dxx,dyy,dzz,diffatv4
real:: diffxx,diffyy,diffzz,adx,ady,adz,dxth,dyth,dzth
real:: fzfx,fz1fx,fzfx1,fz1fx1,fyfx,fy1fx,fyfx1,fy1fx1
real:: fzfy,fz1fy,fzfy1,fz1fy1,fzfxn,fz1fxn,fzfx1n,fz1fx1n
real:: fyfxn,fy1fxn,fyfx1n,fy1fx1n,fzfyn,fz1fyn,fzfy1n,fz1fy1n
real:: fynfx,fy1nfx,fynfx1,fy1nfx1,fznfx,fz1nfx,fznfx1,fz1nfx1
real:: fznfy,fz1nfy,fznfy1,fz1nfy1
!
do ip=ip0,ipn
  if(.not.pat(1,ip)%active)cycle
!
  at=dtran(cat(1,1,1)%zone)
  diff=ddiff(cat(1,1,1)%zone)
  timeleft=dble(dtcurrent)
  totaltime=0.0
  do; if(.not.sngl(timeleft).gt.0.0)exit
    call random(w1,w2,w3)
!...particle location
    x=pat(1,ip)%xyz(1); y=pat(1,ip)%xyz(2); z=pat(1,ip)%xyz(3) 
    xo=pat(1,ip)%xyz(1); yo=pat(1,ip)%xyz(2); zo=pat(1,ip)%xyz(3) 
!
    i=pat(1,ip)%ijk(1); j=pat(1,ip)%ijk(2); k=pat(1,ip)%ijk(3) 
    io=pat(1,ip)%ijk(1); jo=pat(1,ip)%ijk(2); ko=pat(1,ip)%ijk(3) 
    kx=i-1; ky=j-1; kz=k-1
!...bilinearly interpolate on cells centered on the corners of the grid blocks
!
!...we use these indices to find in which corner the particle is located, e.g.,
!...ii=1 if the partcle is located in the second half of the cell in the x direction
!...ii=0 if "                            " first half "                          "
    ii=mod(ifix(x/dx2),2); jj=mod(ifix(y/dy2),2); kk=mod(ifix(z/dz2),2)
!...get parameters from eight cells surrounding corner where partice is located
    iii=i+ii; jjj=j+jj; kkk=k+kk
    ith=1
    do iz=kkk-1,kkk; do iy=jjj-1,jjj; do ix=iii-1,iii
!.....limit index from 1 to nx
      imm=max(1,min(ix,nx)); jmm=max(1,min(iy,ny)); kmm=max(1,min(iz,nz))
!.....limit index in direction of velocity component to a minimum of 1
      kmmkk=kmm-kk; jmmjj=jmm-jj; immii=imm-ii
      th(ith)=por(cat(imm,jmm,kmm)%zone)
      rt(ith)=ret(cat(imm,jmm,kmm)%zone)
      vx(ith)=vel3(1,immii,jmm,kmm)/por(cat(i,jmm,kmm)%zone) 
      vy(ith)=vel3(2,imm,jmmjj,kmm)/por(cat(imm,j,kmm)%zone) 
      vz(ith)=vel3(3,imm,jmm,kmmkk)/por(cat(imm,jmm,k)%zone) 
      ith=ith+1
    enddo; enddo; enddo
!...location of lower faces of the cell over which interolation occurs
    x1=float(kx+ii)*dx-dx/2.0; y1=float(ky+jj)*dy-dy/2.0; z1=float(kz+kk)*dz-dz/2.0
    x1n=float(kx)*dx; y1n=float(ky)*dy; z1n=float(kz)*dz
!...coefficients for interpolation
    fx=(x-x1)/dx; fy=(y-y1)/dy; fz=(z-z1)/dz
    fx1=(1.0-fx); fy1=(1.0-fy); fz1=(1.0-fz)
!
    fxn=(x-x1n)/dx; fx1n=(1.0-fxn); fyn=(y-y1n)/dy
    fy1n=(1.0-fyn); fzn=(z-z1n)/dz; fz1n=(1.0-fzn)
!
    fzfx=fz*fx; fz1fx=fz1*fx; fzfx1=fz*fx1; fz1fx1=fz1*fx1
    fzfy=fz*fy; fz1fy=fz1*fy; fzfy1=fz*fy1; fz1fy1=fz1*fy1
    fyfx=fy*fx; fy1fx=fy1*fx; fyfx1=fy*fx1; fy1fx1=fy1*fx1
    fznfx=fzn*fx; fz1nfx=fz1n*fx; fznfx1=fzn*fx1; fz1nfx1=fz1n*fx1
    fznfy=fzn*fy; fz1nfy=fz1n*fy; fznfy1=fzn*fy1; fz1nfy1=fz1n*fy1
    fzfxn=fz*fxn; fz1fxn=fz1*fxn; fzfx1n=fz*fx1n; fz1fx1n=fz1*fx1n
    fzfyn=fz*fyn; fz1fyn=fz1*fyn; fzfy1n=fz*fy1n; fz1fy1n=fz1*fy1n
    fyfxn=fy*fxn; fy1fxn=fy1*fxn; fyfx1n=fy*fx1n; fy1fx1n=fy1*fx1n
    fynfx=fyn*fx; fy1nfx=fy1n*fx; fynfx1=fyn*fx1; fy1nfx1=fy1n*fx1
!...interpolate porosty to particle location
    t2=fzfx1*th(5)+fzfx*th(6)+fz1fx1*th(1)+fz1fx*th(2) !-y face
    t6=fzfx1*th(7)+fzfx*th(8)+fz1fx1*th(3)+fz1fx*th(4) !+y face
    t3=fzfy1*th(5)+fzfy*th(7)+fz1fy1*th(1)+fz1fy*th(3) !-x face
    t5=fzfy1*th(6)+fzfy*th(8)+fz1fy1*th(2)+fz1fy*th(4) !+x face
    t1=fyfx1*th(3)+fyfx*th(4)+fy1fx1*th(1)+fy1fx*th(2) !-z face
    t7=fyfx1*th(7)+fyfx*th(8)+fy1fx1*th(5)+fy1fx*th(6) !+z face
    thi=fz1*t1+fz*t7 
    dxth=(t5-t3)/dx; dyth=(t6-t2)/dy; dzth=(t7-t1)/dz
    th4=por(cat(i,j,k)%zone) ! porosity in current cell
!...get eight retardation coefficients
    r2=fzfx1*rt(5)+fzfx*rt(6)+fz1fx1*rt(1)+fz1fx*rt(2) !-y face
    r6=fzfx1*rt(7)+fzfx*rt(8)+fz1fx1*rt(3)+fz1fx*rt(4) !+y face
    r3=fzfy1*rt(5)+fzfy*rt(7)+fz1fy1*rt(1)+fz1fy*rt(3) !-x face
    r5=fzfy1*rt(6)+fzfy*rt(8)+fz1fy1*rt(2)+fz1fy*rt(4) !+x face
    r1=fyfx1*rt(3)+fyfx*rt(4)+fy1fx1*rt(1)+fy1fx*rt(2) !-z face
    r7=fyfx1*rt(7)+fyfx*rt(8)+fy1fx1*rt(5)+fy1fx*rt(6) !+z face
    reti=fz1*r1+fz*r7    ! retardation coef. at particle location 
    ret4=ret(cat(i,j,k)%zone) ! ret in current cell
!
    retth=ret4*th4
    sqrtretth=sqrt(ret4*th4)
!...x-velocities; on faces of block
    vx2=fzfx1n*vx(5)+fzfxn*vx(6)+fz1fx1n*vx(1)+fz1fxn*vx(2) !-y face
    vx6=fzfx1n*vx(7)+fzfxn*vx(8)+fz1fx1n*vx(3)+fz1fxn*vx(4) !+y face
    vx3=fzfy1 *vx(5)+ fzfy*vx(7)+fz1fy1 *vx(1)+ fz1fy*vx(3) !-x face
    vx5=fzfy1 *vx(6)+ fzfy*vx(8)+fz1fy1 *vx(2)+ fz1fy*vx(4) !+x face
    vx1=fyfx1n*vx(3)+fyfxn*vx(4)+fy1fx1n*vx(1)+fy1fxn*vx(2) !-z face
    vx7=fyfx1n*vx(7)+fyfxn*vx(8)+fy1fx1n*vx(5)+fy1fxn*vx(6) !+z face
!...at particle location
    vx4=fz1*vx1+fz*vx7
!...y-velocities
    vy2=fzfx1 *vy(5)+ fzfx*vy(6)+fz1fx1 *vy(1)+ fz1fx*vy(2) !-y face
    vy6=fzfx1 *vy(7)+ fzfx*vy(8)+fz1fx1 *vy(3)+ fz1fx*vy(4) !+y face
    vy3=fzfy1n*vy(5)+fzfyn*vy(7)+fz1fy1n*vy(1)+fz1fyn*vy(3) !-x face
    vy5=fzfy1n*vy(6)+fzfyn*vy(8)+fz1fy1n*vy(2)+fz1fyn*vy(4) !+x face
    vy1=fynfx1*vy(3)+fynfx*vy(4)+fy1nfx1*vy(1)+fy1nfx*vy(2) !-z face
    vy7=fynfx1*vy(7)+fynfx*vy(8)+fy1nfx1*vy(5)+fy1nfx*vy(6) !+z face
!...at particle location
    vy4=fz1*vy1+fz*vy7
!...z-velocites
    vz2=fznfx1*vz(5)+fznfx*vz(6)+fz1nfx1*vz(1)+fz1nfx*vz(2) !-y face
    vz6=fznfx1*vz(7)+fznfx*vz(8)+fz1nfx1*vz(3)+fz1nfx*vz(4) !+y face
    vz3=fznfy1*vz(5)+fznfy*vz(7)+fz1nfy1*vz(1)+fz1nfy*vz(3) !-x face
    vz5=fznfy1*vz(6)+fznfy*vz(8)+fz1nfy1*vz(2)+fz1nfy*vz(4) !+x face
    vz1=fyfx1 *vz(3) +fyfx*vz(4)+fy1fx1 *vz(1) +fy1fx*vz(2) !-z face
    vz7=fyfx1 *vz(7) +fyfx*vz(8)+fy1fx1 *vz(5) +fy1fx*vz(6) !+z face
!...interpolated vel
    vz4=fz1n*vz1+fzn*vz7
!...dispersive-diffusive step w/ symmetric square root
    vv4=vx4*vx4+vy4*vy4+vz4*vz4
    v4=sqrt(vv4)
    diffatv4=diff+at*v4
    e1=sqrt(diffatv4)
    dispx=e1*w1; dispy=e1*w2; dispz=e1*w3
    if(vv4.lt.100.0*near)then
      adx=diffatv4*dxth/thi; ady=diffatv4*dyth/thi; adz=diffatv4*dzth/thi
    else
!.....derivative of D =  dj(Dij) = dj(at*v4I)=0.5*at*(vv4**-1/2)di(vv4)=0.5*at*divv4/v4
!.....vv4 = v4*v4, dj(vv4)/v4 = dj(sumi(vi*vi))/v4 = 2vidjvi/v4
!.....Therefore dj(Dij) = 0.5*at*divv4/v4 = at*vidjvi/v4
      dxx=at*(vx4*(vx5-vx3)/dx+vy4*(vy5-vy3)/dx+vz4*(vz5-vz3)/dx)/v4
      dyy=at*(vx4*(vx6-vx2)/dy+vy4*(vy6-vy2)/dy+vz4*(vz6-vz2)/dy)/v4
      dzz=at*(vx4*(vx7-vx1)/dz+vy4*(vy7-vy1)/dz+vz4*(vz7-vz1)/dz)/v4
!.....gradient drift term
!.....We need (1/(R theta)) grad(D theta)= 
!.....D grad (theta)/(R theta) + grad D/R
!.....grad theta term, (Dij/(R*theta)) gradj (theta)
      adx=dxx+(diffatv4*dxth)/thi   
      ady=dyy+(diffatv4*dyth)/thi
      adz=dzz+(diffatv4*dzth)/thi
    endif
!...summary of the dispersion tensor and velocity
    if(ibug.ge.2)then
      write(ibugout,*)'---------------  location  ---------------'
      write(ibugout,*)'x,y,z ',x,y,z
      write(ibugout,*)'i,j,k ',kx+1,ky+1,kz+1
      write(ibugout,*)'---------------  velocity  ---------------'
      write(ibugout,*)'vx,vy,vz,magnitude(v) ',vx4,vy4,vz4,v4
!
      write(ibugout,*)'---------- B/sqrt(2) =sqrt(D) ------------'
      write(ibugout,*)'b11,b12,b13       ',e1,0.,0.
      write(ibugout,*)'b12,b22,b23       ',0.,e1,0.
      write(ibugout,*)'b13,b23,b33       ',0.,0.,e1
!.....B dot B Transpose/2
      write(ibugout,*)'------------- D=B dot B/2 ----------------'
      write(ibugout,*)'bb11,bb12,bb13 ',(e1*e1),0.,0.
      write(ibugout,*)'bb12,bb22,bb23 ',0.,(e1*e1),0.
      write(ibugout,*)'bb13,bb23,bb33 ',0.,0.,(e1*e1)
!.....Dispersion Tensor
!      write(ibugout,*)'--------------     D       ---------------'
!      write(ibugout,*)'D11,D12,D13       ',diffxx,0.,0.
!      write(ibugout,*)'D12,D22,D23       ',0.,diffyy,0.
!      write(ibugout,*)'D13,D22,D13       ',0.,0.,diffzz
!.....Grad D
      write(ibugout,*)'-------------  Grad D    -----------------'
      write(ibugout,*)'dD11/dx,dD12/dy,dD13/dz ',dxx,0.,0.
      write(ibugout,*)'dD21/dx,dD22/dy,dD23/dz ',0.,dyy,0.
      write(ibugout,*)'dD31/dx,dD32/dy,dD33/dz ',0.,0.,dzz
    endif
!...time-step control
    if(dtcntrl.ne.0.0)then
      ii=mx*(ii-min(kx,1))
      jj=my*(jj-min(ky,1))
      kk=mz*(kk-min(kz,1))
      dtmin=min(large,cat(i+ii,j+jj,k+kk)%tc,sngl(timeleft))
    else
      dtmin=dtcurrent
    endif
!...move particle
    call stream(vel3,x,y,z,retth,kx,ky,kz,i,j,k,dtmin,xstream,ystream,zstream)
    r2dt=sqrt(2.0*dtmin/reti)
    dtminreti=dtmin/reti
    x=xstream+mx*(adx*dtminreti+dispx*r2dt)
    y=ystream+my*(ady*dtminreti+dispy*r2dt)
    z=zstream+mz*(adz*dtminreti+dispz*r2dt)
    totaltime=totaltime+dble(dtmin)
    timeleft=dble(dtcurrent)-totaltime
!...decay m=m0*exp(k*dt)
    if(decay(cat(i,j,k)%zone).ne.0.0)&
    pat(1,ip)%pmass=pat(1,ip)%pmass*dexp(dble(decay(cat(i,j,k)%zone)*dtmin))
!...reflect particles at boundaries as needed or tag to remove
    call reflect(sngl(timeleft),pat,cat,x,y,z,ip)
    if(.not.pat(1,ip)%active)exit
!...compute cell location
    kx=ifix(x/dx); ky=ifix(y/dy); kz=ifix(z/dz)
    i=kx+1; j=ky+1; k=kz+1
    xd=x-xo;yd=y-yo;zd=z-zo
!    if(dtmin.lt.near)then   
    if(sqrt(xd*xd+yd*yd+zd*zd).lt.near)then
      xm=kx*dx; ym=ky*dy; zm=kz*dz; xp=xm+dx; yp=ym+dy; zp=zm+dz
!...correct cell location for odd case where particle lands on cell boundary          
      call icell_correct(kx,ky,kz,i,j,k,x,y,z,xm,ym,zm,xp,yp,zp,vel3)
    endif
!...absorbing boundaries
    if(cat(i,j,k)%bc_number.ne.0)call absorb(pat,cat,bounds,vel3,por,ret,&
    kx,ky,kz,i,j,k,ip,dtmin,sngl(timeleft),x,y,z)
!...update particle and cell attributes
    call updtp(io,jo,ko,i,j,k,ip,x,y,z,pat,cat)
    if(.not.pat(1,ip)%active)exit
  enddo
enddo
return
end   
!------------------------------------------------------------
! movep4
!------------------------------------------------------------
subroutine movep4(ip0,ipn,dtcurrent,vel3,por,ret,&
dlong,dtran,ddiff,decay,pat,cat,bounds)
! variable porosity, anisotropic dispersion displacement
! with bilinear interpolation of the dispersion tensor
!
! Block-Centered FD velocity field for a 2-D 3 x 3 system.
! o - particle location.
! X - Region over which bilinear interpolation occurs for this particle.
! Velocities associated with the interpolation of this
! particle have arrow heads shown.
!
!               ------------- ------------ -------------
!              |             |            |             |
!              |             |            |             |
!              |             |            |             |
!              |            --->         --->           |
!              |             |            |             |
!              |             |      ^     |      ^      |
!              |             |      |     |      |      |
!               -------|------------|------------|------
!              |             |      |XXXXX|      |      |
!              |             |       XXXXo|             |
!              |             |       XXXXX|             |
!              |            --->         --->           |
!              |             |            |             |
!              |             |      ^     |      ^      |
!              |             |      |     |      |      |
!               -------|------------|------------|------
!              |             |      |     |      |      |
!              |             |            |             |
!              |             |            |             |   
!              |             -            -             |
!              |             |            |             |
!              |             |            |             |
!              |             |            |             |
!              -------------- ------------ --------------
!
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
integer:: ip,ip0,ipn
integer:: kx,ky,kz,i,j,k,ii,jj,kk,iii,jjj,kkk,imm,jmm,kmm
integer:: immii,jmmjj,kmmkk,ith,ix,iy,iz,io,jo,ko
real:: por(nzone),ret(nzone),dlong(nzone),dtran(nzone),ddiff(nzone),decay(nzone)
real:: th(8),rt(8),vx(8),vy(8),vz(8)
!
real:: vel3(3,0:nx,0:ny,0:nz)
real:: dtcurrent,dtmin
double precision:: timeleft,totaltime
real:: x,y,z,xo,yo,zo,xd,yd,zd,xm,ym,zm,xp,yp,zp,xold,yold,zold,x1,y1,z1,x1n,y1n,z1n
real:: fx,fy,fz,fx1,fy1,fz1,fxn,fyn,fzn,fx1n,fy1n,fz1n
real:: xstream,ystream,zstream
real:: thi,th4,alt,al,retth,sqrtretth,reti,ret4,dtminreti
real:: t1,t2,t3,t4,t5,t6,t7,r1,r2,r3,r4,r5,r6,r7
real:: vx1,vx2,vx3,vx4,vx5,vx6,vx7,vy1,vy2,vy3,vy4,vy5,vy6,vy7
real:: vz1,vz2,vz3,vz4,vz5,vz6,vz7,v4,vv4,vxx,vyy,vzz,vxy,vxz,vyz
real:: e1,e2,e3,w1,w2,w3
real:: vxpzxpz,vxpz,vyxpz,vxxpz,vxypz,vzxpz,vxzxz,vyzyz,vxymyz
real:: vyyvzxpz,vyyvxxpz,beta1,beta2,beta3,vyybeta2,vyxpzbeta2
real:: b11,b12,b13,b22,b23,b33,dispx,dispy,dispz
real:: dispxt,dispyt,dispzt
real:: dyvz,dzvz,dxvx,dyvx,dzvx,dxvy,dyvy,dzvy,dxvz,v42
real:: vx4dxvx,vx4dyvx,vx4dzvx,vy4dxvy,vy4dyvy,vy4dzvy
real:: vz4dxvz,vz4dyvz,vz4dzvz,dxvv4,dyvv4,dzvv4,altv4
real:: dxvv4v4,dyvv4v4,dzvv4v4,at2
real:: dxx,dxy,dxz,dyx,dyy,dyz,dzx,dzy,dzz,diffatv4
real:: diffxx,diffxy,diffxz,diffyy,diffyz,diffzz
real:: dgradthx,dgradthy,dgradthz,adx,ady,adz
real:: at,diff,ratv4d,dxd,dyd,dzd,dxth,dyth,dzth,r2dt
!
real:: fzfx,fz1fx,fzfx1,fz1fx1,fyfx,fy1fx,fyfx1,fy1fx1
real:: fzfy,fz1fy,fzfy1,fz1fy1,fzfxn,fz1fxn,fzfx1n,fz1fx1n
real:: fyfxn,fy1fxn,fyfx1n,fy1fx1n,fzfyn,fz1fyn,fzfy1n,fz1fy1n
real:: fynfx,fy1nfx,fynfx1,fy1nfx1,fznfx,fz1nfx,fznfx1,fz1nfx1
real:: fznfy,fz1nfy,fznfy1,fz1nfy1
!
do ip=ip0,ipn
  if(.not.pat(1,ip)%active)cycle
!
  al=dlong(cat(1,1,1)%zone)
  at=dtran(cat(1,1,1)%zone)
  diff=ddiff(cat(1,1,1)%zone)
  alt=al-at
  timeleft=dble(dtcurrent)
  totaltime=0.0
  do; if(.not.sngl(timeleft).gt.0.0)exit
    call random(w1,w2,w3)
!...particle location
    x=pat(1,ip)%xyz(1); y=pat(1,ip)%xyz(2); z=pat(1,ip)%xyz(3) 
    xo=pat(1,ip)%xyz(1); yo=pat(1,ip)%xyz(2); zo=pat(1,ip)%xyz(3) 
!
    i=pat(1,ip)%ijk(1); j=pat(1,ip)%ijk(2); k=pat(1,ip)%ijk(3) 
    io=pat(1,ip)%ijk(1); jo=pat(1,ip)%ijk(2); ko=pat(1,ip)%ijk(3) 
    kx=i-1; ky=j-1; kz=k-1
!...bilinearly interpolate on cells centered on the corners of the grid blocks
!
!...we use these indices to find in which corner the particle is located
!...ii=1 if the partcle is located in the second half of the cell in the x direction
!...=0 if "                            " first half "                          "
    ii=mod(ifix(x/dx2),2); jj=mod(ifix(y/dy2),2); kk=mod(ifix(z/dz2),2)
!...get parameters from eight cells surrounding corner where partice is located
    iii=i+ii; jjj=j+jj; kkk=k+kk
    ith=1
    do iz=kkk-1,kkk; do iy=jjj-1,jjj; do ix=iii-1,iii
      imm=max(1,min(ix,nx)); jmm=max(1,min(iy,ny)); kmm=max(1,min(iz,nz))
      kmmkk=kmm-kk; jmmjj=jmm-jj; immii=imm-ii
      th(ith)=por(cat(imm,jmm,kmm)%zone)
      rt(ith)=ret(cat(imm,jmm,kmm)%zone)
      vx(ith)=vel3(1,immii,jmm,kmm)/por(cat(i,jmm,kmm)%zone) 
      vy(ith)=vel3(2,imm,jmmjj,kmm)/por(cat(imm,j,kmm)%zone) 
      vz(ith)=vel3(3,imm,jmm,kmmkk)/por(cat(imm,jmm,k)%zone) 
      ith=ith+1
    enddo; enddo; enddo
!.........location of lower faces of the cell over which interolation occurs
    x1=float(kx+ii)*dx-dx/2.0; y1=float(ky+jj)*dy-dy/2.0; z1=float(kz+kk)*dz-dz/2.0
    x1n=float(kx)*dx; y1n=float(ky)*dy; z1n=float(kz)*dz
!.........coefficients for interpolation
    fx=(x-x1)/dx; fy=(y-y1)/dy; fz=(z-z1)/dz
    fx1=(1.0-fx); fy1=(1.0-fy); fz1=(1.0-fz)
!
    fxn=(x-x1n)/dx; fx1n=(1.0-fxn); fyn=(y-y1n)/dy
    fy1n=(1.0-fyn); fzn=(z-z1n)/dz; fz1n=(1.0-fzn)
!
    fzfx=fz*fx; fz1fx=fz1*fx; fzfx1=fz*fx1; fz1fx1=fz1*fx1
    fzfy=fz*fy; fz1fy=fz1*fy; fzfy1=fz*fy1; fz1fy1=fz1*fy1
    fyfx=fy*fx; fy1fx=fy1*fx; fyfx1=fy*fx1; fy1fx1=fy1*fx1
    fznfx=fzn*fx; fz1nfx=fz1n*fx; fznfx1=fzn*fx1; fz1nfx1=fz1n*fx1
    fznfy=fzn*fy; fz1nfy=fz1n*fy; fznfy1=fzn*fy1; fz1nfy1=fz1n*fy1
    fzfxn=fz*fxn; fz1fxn=fz1*fxn; fzfx1n=fz*fx1n; fz1fx1n=fz1*fx1n
    fzfyn=fz*fyn; fz1fyn=fz1*fyn; fzfy1n=fz*fy1n; fz1fy1n=fz1*fy1n
    fyfxn=fy*fxn; fy1fxn=fy1*fxn; fyfx1n=fy*fx1n; fy1fx1n=fy1*fx1n
    fynfx=fyn*fx; fy1nfx=fy1n*fx; fynfx1 =fyn*fx1; fy1nfx1=fy1n*fx1
!....interpolate porosty to particle location
    t2=fzfx1*th(5)+fzfx*th(6)+fz1fx1*th(1)+fz1fx*th(2) !-y face
    t6=fzfx1*th(7)+fzfx*th(8)+fz1fx1*th(3)+fz1fx*th(4) !+y face
    t3=fzfy1*th(5)+fzfy*th(7)+fz1fy1*th(1)+fz1fy*th(3) !-x face
    t5=fzfy1*th(6)+fzfy*th(8)+fz1fy1*th(2)+fz1fy*th(4) !+x face
    t1=fyfx1*th(3)+fyfx*th(4)+fy1fx1*th(1)+fy1fx*th(2) !-z face
    t7=fyfx1*th(7)+fyfx*th(8)+fy1fx1*th(5)+fy1fx*th(6) !+z face
    thi=fz1*t1+fz*t7 
    dxth=(t5-t3)/dx; dyth=(t6-t2)/dy; dzth=(t7-t1)/dz
    th4=por(cat(i,j,k)%zone) ! porosity in current cell
!...get eight retardation coefficients
    r2=fzfx1*rt(5)+fzfx*rt(6)+fz1fx1*rt(1)+fz1fx*rt(2) !-y face
    r6=fzfx1*rt(7)+fzfx*rt(8)+fz1fx1*rt(3)+fz1fx*rt(4) !+y face
    r3=fzfy1*rt(5)+fzfy*rt(7)+fz1fy1*rt(1)+fz1fy*rt(3) !-x face
    r5=fzfy1*rt(6)+fzfy*rt(8)+fz1fy1*rt(2)+fz1fy*rt(4) !+x face
    r1=fyfx1*rt(3)+fyfx*rt(4)+fy1fx1*rt(1)+fy1fx*rt(2) !-z face
    r7=fyfx1*rt(7)+fyfx*rt(8)+fy1fx1*rt(5)+fy1fx*rt(6) !+z face
    reti=fz1*r1+fz*r7    ! retardation coef. at particle location 
    ret4=ret(cat(i,j,k)%zone) ! ret in current cell
!
    retth=ret4*th4
    retth=ret4*th4
    sqrtretth=sqrt(ret4*th4)
!...x-velocities; on faces of block
    vx2=fzfx1n*vx(5)+fzfxn*vx(6)+fz1fx1n*vx(1)+fz1fxn*vx(2) !-y face
    vx6=fzfx1n*vx(7)+fzfxn*vx(8)+fz1fx1n*vx(3)+fz1fxn*vx(4) !+y face
    vx3=fzfy1 *vx(5)+ fzfy*vx(7)+fz1fy1 *vx(1)+ fz1fy*vx(3) !-x face
    vx5=fzfy1 *vx(6)+ fzfy*vx(8)+fz1fy1 *vx(2)+ fz1fy*vx(4) !+x face
    vx1=fyfx1n*vx(3)+fyfxn*vx(4)+fy1fx1n*vx(1)+fy1fxn*vx(2) !-z face
    vx7=fyfx1n*vx(7)+fyfxn*vx(8)+fy1fx1n*vx(5)+fy1fxn*vx(6) !+z face
!...at particle location
    vx4=fz1*vx1+fz*vx7
!...y-velocities
    vy2=fzfx1 *vy(5)+ fzfx*vy(6)+fz1fx1 *vy(1)+ fz1fx*vy(2) !-y face
    vy6=fzfx1 *vy(7)+ fzfx*vy(8)+fz1fx1 *vy(3)+ fz1fx*vy(4) !+y face
    vy3=fzfy1n*vy(5)+fzfyn*vy(7)+fz1fy1n*vy(1)+fz1fyn*vy(3) !-x face
    vy5=fzfy1n*vy(6)+fzfyn*vy(8)+fz1fy1n*vy(2)+fz1fyn*vy(4) !+x face
    vy1=fynfx1*vy(3)+fynfx*vy(4)+fy1nfx1*vy(1)+fy1nfx*vy(2) !-z face
    vy7=fynfx1*vy(7)+fynfx*vy(8)+fy1nfx1*vy(5)+fy1nfx*vy(6) !+z face
!...at particle location
    vy4=fz1*vy1+fz*vy7
!...z-velocites
    vz2=fznfx1*vz(5)+fznfx*vz(6)+fz1nfx1*vz(1)+fz1nfx*vz(2) !-y face
    vz6=fznfx1*vz(7)+fznfx*vz(8)+fz1nfx1*vz(3)+fz1nfx*vz(4) !+y face
    vz3=fznfy1*vz(5)+fznfy*vz(7)+fz1nfy1*vz(1)+fz1nfy*vz(3) !-x face
    vz5=fznfy1*vz(6)+fznfy*vz(8)+fz1nfy1*vz(2)+fz1nfy*vz(4) !+x face
    vz1=fyfx1 *vz(3) +fyfx*vz(4)+fy1fx1 *vz(1) +fy1fx*vz(2) !-z face
    vz7=fyfx1 *vz(7) +fyfx*vz(8)+fy1fx1 *vz(5) +fy1fx*vz(6) !+z face
!...interpolated vel
    vz4=fz1n*vz1+fzn*vz7
!...dispersive-diffusive step w/ symmetric square root
    vxx=vx4*vx4; vyy=vy4*vy4; vzz=vz4*vz4 
    vv4=vxx+vyy+vzz
    v4=sqrt(vv4)  
    vxy=vx4*vy4; vxz=vx4*vz4; vyz=vy4*vz4
    e1=sqrt(al*v4+diff); e2=sqrt(at*v4+diff)
    if(vv4.lt.100.0*near)then
      adx=0.; ady=0.; adz=0.
      dispx=e1*w1; dispy=e1*w2; dispz=e1*w3
    else
      vxpz=vx4+vz4
      vxpzxpz=vxpz*vxpz
      vxxpz=vx4*vxpz; vyxpz=vy4*vxpz; vzxpz=vz4*vxpz; vxzxz=vxz*vxz; vyzyz=vyz*vyz
      vxymyz=vxy-vyz
      vyyvzxpz=vyy+vzxpz
      vyyvxxpz=vyy+vxxpz
      beta1=vv4/e1; beta2=(2.0*vyy+vxpzxpz)/e2; beta3=beta1*beta2*e1
      vyybeta2=vyy/beta2; vyxpzbeta2=vyxpz/beta2
!.....3 eigenvectors: vx,     vy,         vz
!                    -vy,  vx+vz,        -vy
!           -vyy-vz*vxpz,vxy-vyz,vx*vxpz+vyy
      b11=vxx/beta1+     vyybeta2+vyyvzxpz*vyyvzxpz/beta3
      b12=vxy/beta1-   vyxpzbeta2-  vyyvzxpz*vxymyz/beta3
      b13=vxz/beta1+     vyybeta2-vyyvzxpz*vyyvxxpz/beta3
      b22=vyy/beta1+vxpzxpz/beta2+    vxymyz*vxymyz/beta3
      b23=vyz/beta1-   vyxpzbeta2+  vyyvxxpz*vxymyz/beta3
      b33=vzz/beta1+     vyybeta2+vyyvxxpz*vyyvxxpz/beta3
      dispx=b11*w1+b12*w2+b13*w3; dispy=b12*w1+b22*w2+b23*w3; dispz=b13*w1+b23*w2+b33*w3
      dispxt=b11+b12+b13; dispyt=b12+b22+b23; dispzt=b13+b23+b33
!.....gradient advective step (bilinear interpolation)
      dxvx=(vx5-vx3)/dx; dyvx=(vx6-vx2)/dy; dzvx=(vx7-vx1)/dz 
      dxvy=(vy5-vy3)/dx; dyvy=(vy6-vy2)/dy; dzvy=(vy7-vy1)/dz
      dxvz=(vz5-vz3)/dx; dyvz=(vz6-vz2)/dy; dzvz=(vz7-vz1)/dz
!.....vv4 = v4*v4, dj(vv4)/v4 = dj(sumi(vi*vi))/v4 = 2vidjvi/v4
      vx4dxvx=vx4*dxvx; vx4dyvx=vx4*dyvx; vx4dzvx=vx4*dzvx
      vy4dxvy=vy4*dxvy; vy4dyvy=vy4*dyvy; vy4dzvy=vy4*dzvy
      vz4dxvz=vz4*dxvz; vz4dyvz=vz4*dyvz; vz4dzvz=vz4*dzvz
      v42=v4/2.0
      dxvv4=(vx4dxvx+vy4dxvy+vz4dxvz)/v42 
      dyvv4=(vx4dyvx+vy4dyvy+vz4dyvz)/v42 
      dzvv4=(vx4dzvx+vy4dzvy+vz4dzvz)/v42
!.....derivative of D =  dj(Dij) = dj(alt(vivj/v4)+at*v4I)
!.....=(alt/v4)(vjdj(vi)+vidj(vj))-.5*(alt)vivj(vv4**-3/2)dj(vv4)+
!.....0.5*at*(vv4**-1/2)di(vv4)
!.....=(alt/v4)(vjdj(vi)+vidj(vj))-.5*(alt)vivj(1/vv4)djvv4/v4+
!.....0.5*at*divv4/v4
!
!.....We need (1/(R theta)) grad(D theta)= 
!.....D grad (theta)/(R theta) + grad D/R
      altv4=alt/v4
      dxvv4v4=.5*dxvv4/v4; dyvv4v4=.5*dyvv4/v4; dzvv4v4=.5*dzvv4/v4
      at2=at/2.0
!
      dxx=altv4*(vx4dxvx+vx4dxvx-vxx*dxvv4v4)+at2*dxvv4
      dxy=altv4*(vy4*dyvx+vx4*dyvy-vxy*dyvv4v4)
      dxz=altv4*(vz4*dzvx+vx4*dzvz-vxz*dzvv4v4)
!
      dyx=altv4*(vx4*dxvy+vy4*dxvx-vxy*dxvv4v4)
      dyy=altv4*(vy4dyvy+vy4dyvy-vyy*dyvv4v4)+at2*dyvv4
      dyz=altv4*(vz4*dzvy+vy4*dzvz-vyz*dzvv4v4)
!                                                                       
      dzx=altv4*(vx4*dxvz+vz4*dxvx-vxz*dxvv4v4)
      dzy=altv4*(vy4*dyvz+vz4*dyvy-vyz*dyvv4v4)
      dzz=altv4*(vz4dzvz+vz4dzvz-vzz*dzvv4v4)+at2*dzvv4
!.....grad theta term, (Dij/(R*theta)) gradj (theta)
      diffatv4=diff+at*v4
      diffxx=diffatv4+vxx*altv4; diffxy=vxy*altv4;          diffxz=vxz*altv4
                                 diffyy=diffatv4+vyy*altv4; diffyz=vyz*altv4
                                                            diffzz=diffatv4+vzz*altv4
      dgradthx=(diffxx*dxth+diffxy*dyth+diffxz*dzth)/thi
      dgradthy=(diffxy*dxth+diffyy*dyth+diffyz*dzth)/thi
      dgradthz=(diffxz*dxth+diffyz*dyth+diffzz*dzth)/thi
!.....drift coefficients
      adx=((dxx+dxy+dxz)+dgradthx); ady=((dyx+dyy+dyz)+dgradthy); adz=((dzx+dzy+dzz)+dgradthz)
    endif
    if(ibug.ge.2)then
!.....summary of the dispersion tensor and velocity
      write(ibugout,*)'---------------  location  ---------------'
      write(ibugout,*)'x,y,z ',x,y,z
      write(ibugout,*)'i,j,k ',kx+1,ky+1,kz+1
      write(ibugout,*)'---------------  velocity  ---------------'
      write(ibugout,*)'vx,vy,vz,magnitude(v) ',vx4,vy4,vz4,v4
!
      write(ibugout,*)'---------- B/sqrt(2) =sqrt(D) ------------'
      write(ibugout,*)'b11,b12,b13       ',b11,b12,b13
      write(ibugout,*)'b12,b22,b23       ',b12,b22,b23
            write(ibugout,*)'b13,b23,b33       ',b13,b23,b33
!.....B dot B Transpose/2
      write(ibugout,*)'------------- D=B dot B/2 ----------------'
      write(ibugout,*)'bb11,bb12,bb13 ',(b11*b11+b12*b12+b13*b13),&
                                 (b11*b12+b12*b22+b13*b23),&
                                 (b11*b13+b12*b23+b13*b33)
      write(ibugout,*)'bb12,bb22,bb23 ',(b12*b11+b22*b12+b23*b13),&
                                 (b12*b12+b22*b22+b23*b23),&
                                 (b12*b13+b22*b23+b23*b33)
      write(ibugout,*)'bb13,bb23,bb33 ',(b13*b11+b23*b12+b33*b13),&
                                 (b13*b12+b23*b22+b33*b23),&
                                 (b13*b13+b23*b23+b33*b33)
!.....Dispersion Tensor
      write(ibugout,*)'--------------     D       ---------------'
      write(ibugout,*)'D11,D12,D13       ',diffxx,diffxy,diffxz
      write(ibugout,*)'D12,D22,D23       ',diffxy,diffyy,diffyz
      write(ibugout,*)'D13,D22,D13       ',diffxz,diffyz,diffzz
!.....Grad D
      write(ibugout,*)'-------------  Grad D    -----------------'
      write(ibugout,*)'dD11/dx,dD12/dy,dD13/dz ',dxx,dxy,dxz
      write(ibugout,*)'dD21/dx,dD22/dy,dD23/dz ',dyx,dyy,dyz
      write(ibugout,*)'dD31/dx,dD32/dy,dD33/dz ',dzx,dzy,dzz
    endif
!...time-step control
    if(dtcntrl.ne.0.0)then 
      ii=mx*(ii-min(kx,1))
      jj=my*(jj-min(ky,1))
      kk=mz*(kk-min(kz,1))
      dtmin=min(large,cat(i+ii,j+jj,k+kk)%tc,sngl(timeleft))
    else
      dtmin=dtcurrent
    endif
!...move particle
    call stream(vel3,x,y,z,retth,kx,ky,kz,i,j,k,dtmin,xstream,ystream,zstream)
    r2dt=sqrt(2.0*dtmin/reti)
    dtminreti=dtmin/reti
    x=xstream+mx*(adx*dtminreti+dispx*r2dt)
    y=ystream+my*(ady*dtminreti+dispy*r2dt)
    z=zstream+mz*(adz*dtminreti+dispz*r2dt)
    totaltime=totaltime+dble(dtmin)
    timeleft=dble(dtcurrent)-totaltime
!...decay m=m0*exp(k*dt)
    if(decay(cat(i,j,k)%zone).ne.0.0)&
    pat(1,ip)%pmass=pat(1,ip)%pmass*dexp(dble(decay(cat(i,j,k)%zone)*dtmin))
!...reflect particles at boundaries as needed or tag to remove
    call reflect(sngl(timeleft),pat,cat,x,y,z,ip)
    if(.not.pat(1,ip)%active)exit
!...compute cell location
    kx=ifix(x/dx); ky=ifix(y/dy); kz=ifix(z/dz)
    i=kx+1; j=ky+1; k=kz+1
    xd=x-xo;yd=y-yo;zd=z-zo
!    if(dtmin.lt.near)then   
    if(sqrt(xd*xd+yd*yd+zd*zd).lt.near)then
      xm=kx*dx; ym=ky*dy; zm=kz*dz; xp=xm+dx; yp=ym+dy; zp=zm+dz
!.....correct cell location for odd case where particle lands on cell boundary          
      call icell_correct(kx,ky,kz,i,j,k,x,y,z,xm,ym,zm,xp,yp,zp,vel3)
    endif
!...absorbing boundaries
    if(cat(i,j,k)%bc_number.ne.0)call absorb(pat,cat,bounds,vel3,por,ret,&
                kx,ky,kz,i,j,k,ip,dtmin,sngl(timeleft),x,y,z)
!...update particle and cell attributes
    call updtp(io,jo,ko,i,j,k,ip,x,y,z,pat,cat)
    if(.not.pat(1,ip)%active)exit
  enddo
enddo
return
end   
!------------------------------------------------------------
! tcntrl
!------------------------------------------------------------
subroutine tcntrl(vel3,por,ret,dtran,dlong,ddiff,cat)
!
! Develop a static time-step control as a real number cell attribute
! based on the maximum dispersive displacement. 
use global
type (cell)::     cat(nx,ny,nz)
real:: por(nzone),ret(nzone),dtran(nzone),dlong(nzone),ddiff(nzone)
real:: vel3(3,0:nx,0:ny,0:nz),diff
!
dtminimum=large
dtmaximum=0
imin=0
jmin=0
kmin=0
cat(:,:,:)%tc=dtinit
do k=1,nz; do j=1,ny; do i=1,nx
!.......get all 12 component velocities near corner of cells
! e.g., x-diection i-1,j,k; i,j,k; i+1,j,k  
!                  i-1,j+1,k; i,j+1,k; i+1,j+1,k  
!                  i-1,j+1,k+1; i,j+1,k+1; i+1,j+1,k+1
!                  i-1,j+1,k+1; i,j+1,k+1; i+1,j+1,k+1
  im=1; ip=1
  jm=1; jp=1
  km=1; kp=1
  if(i.eq.1)im=0; if(i.eq.nx)ip=0
  if(j.eq.1)jm=0; if(j.eq.ny)jp=0
  if(k.eq.1)km=0; if(k.eq.nz)kp=0
!
  diffxyz=0.0
!.find max representative dispersion coefficient
  do kk=0,kp*mz; do jj=0,jp*my; do ii=0,ip*mx
    th=por(cat(i+ii,j+jj,k+kk)%zone)         
    diff=ddiff(cat(i+ii,j+jj,k+kk)%zone)         
    at=dtran(cat(i+ii,j+jj,k+kk)%zone)
    al=dlong(cat(i+ii,j+jj,k+kk)%zone)
    r=ret(cat(i+ii,j+jj,k+kk)%zone)
!...average to cell centers 
    vx=(vel3(1,i-mx*im+ii,j      +jj,k      +kk)+&
        vel3(1,i      +ii,j      +jj,k      +kk))/2.0
    vy=(vel3(2,i      +ii,j-my*jm+jj,k      +kk)+&
        vel3(2,i      +ii,j      +jj,k      +kk))/2.0
    vz=(vel3(3,i      +ii,j      +jj,k-mz*km+kk)+&
        vel3(3,i      +ii,j      +jj,k      +kk))/2.0
!...magnitude 
    v4=sqrt(vx**2+vy**2+vz**2)
!...compute representative dispersion coefficient 
    diffxyz=max(diffxyz,(max(at,al)/r)*v4+diff*th/r)
  enddo; enddo; enddo
! minimum cell size
  rmindxyz2=min(rmx*dx*dx,rmy*dy*dy,rmz*dz*dz)
  dtmin=sngl(dtinit)
  if(diffxyz.ne.0.0)then
    cat(i,j,k)%tc=min(dtmin,dtcntrl*rmindxyz2/(48.0*diffxyz))
    if(dtcntrl*rmindxyz2/(48.0*diffxyz).lt.dtminimum)then
      dtminimum=dtcntrl*rmindxyz2/(48.0*diffxyz)
      imin=i
      jmin=j
      kmin=k
    endif
    dtmaximum=max(dtcntrl*rmindxyz2/(48.0*diffxyz),dtmaximum)
  endif
enddo; enddo; enddo
write(*,'(4(a,e10.5))')'CURTIME:',sngl(curtime),' DTINIT =',dtinit,' DTMINIMUM =',dtminimum ,' DTMAXIMUM =',dtmaximum 
write(*,'((a,3i10))')'I,J,K LOCATION OF MINIMUM = ',imin,jmin,kmin
write(iout,'(4(a,e10.5))')'CURTIME:',sngl(curtime),' DTINIT =',dtinit,' DTMINIMUM =',dtminimum ,' DTMAXIMUM =',dtmaximum 
write(iout,'((a,3i10))')'I,J,K LOCATION OF MINIMUM = ',imin,jmin,kmin
return
end               
!------------------------------------------------------------
! courant
!------------------------------------------------------------
subroutine courant(source,bounds,vel3,cat,por)
! determine time step over which to apply type 1 BC
use global
implicit none
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
integer:: source(maxsource),ibounds,itype,kx,ky,kz,isource
integer:: i,j,k,im,jm,km,ii,jj,kk,l4,kbot,ktop
real:: vel3(3,0:nx,0:ny,0:nz),por(nzone),th4,tc
real:: tc1,tc2,tc3,tc4,tc5,tc6,tc7,tc8,vxm,vym,vzm,vxp,vyp,vzp      
!
do isource=1,nsource
  ibounds=source(isource)
  itype=bounds(ibounds)%bc_type
  if(itype.eq.1)then
    i=bounds(ibounds)%ijk(1); j=bounds(ibounds)%ijk(2); kbot=bounds(ibounds)%kbot; ktop=bounds(ibounds)%ktop
    bounds(ibounds)%dt_type1=sngl(dt)
    tc=0.0
    do k=kbot,ktop
      kx=i-1; ky=j-1; kz=k-1
!...porosity
      th4=por(cat(i,j,k)%zone)
!...time increment from current time to update boundary
      if(max(abs(vel3(1,kx,j,k)),abs(vel3(1,i,j,k))).ne.0.0)&
      bounds(ibounds)%dt_type1=th4*rmx*dx2/max(abs(vel3(1,kx,j,k)),abs(vel3(1,i,j,k)))
      if(max(abs(vel3(2,i,ky,k)),abs(vel3(2,i,j,k))).ne.0.0)&
      bounds(ibounds)%dt_type1=min(bounds(ibounds)%dt_type1,th4*rmy*dy2/&
      max(abs(vel3(2,i,ky,k)),abs(vel3(2,i,j,k))))
      if(max(abs(vel3(3,i,j,kz)),abs(vel3(3,i,j,k))).ne.0.0)&
      bounds(ibounds)%dt_type1=min(bounds(ibounds)%dt_type1,th4*rmz*dz2/&
      max(abs(vel3(3,i,j,kz)),abs(vel3(3,i,j,k))))
!...neighboring time step control for constant conc. cell 
      im=0; jm=0; km=0
      if(i.ne.1)im=1; if(j.ne.1)jm=1; if(k.ne.1)km=1
      do ii=0,mx; do jj=0,my; do kk=0,mz
        tc=max(tc,cat(i-im+ii,j-jm+jj,k-km+kk)%tc)
      enddo; enddo; enddo
!...minimum based on advection as well as dispersion (through tcntrl)
      if(iadvect.ne.1)bounds(ibounds)%dt_type1=min(bounds(ibounds)%dt_type1,tc)
    enddo
!...finally, allow user to control the time step through refresh
    bounds(ibounds)%dt_type1=bounds(ibounds)%dt_type1*bounds(ibounds)%refresh     
  endif
enddo
return
end
!------------------------------------------------------------
! split
!------------------------------------------------------------
subroutine split(pat,cat,bounds,por,ret)
! split particles to achieve a specifed mass resolution
! Two controls:   cres:  minimum concentration of interest
!                 nres:  minimum particle resolution within cell
!  
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
integer:: ip,itype,ibounds,ip_res
integer:: isplit,nsplit,iptype,i,j,k
real:: por(nzone),ret(nzone),x,y,z
real:: time
double precision:: pmass_old(nspec),pmass(nspec),pmass_min(nspec),cmass(nspec),conc(nspec)
!
if(nres.le.1)return
do ip=1,np
  pmass_old=pat(1,ip)%pmass
  i=pat(1,ip)%ijk(1); j=pat(1,ip)%ijk(2); k=pat(1,ip)%ijk(3)
! don't split particles with mass less than pmass_min:
  pmass_min=cres*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone)*vol/float(nres)   
! old  if(pmass_old.ge.pmass_min)then                            ! particles with pmass > pmass_min
   if(maxval(pmass_old-pmass_min).gt.0.0)then                    ! DAB particles with ALL pmass > pmass_min
    cmass=cat(i,j,k)%cmass                                          ! mass in cell
    conc=cmass/dble(por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone)*vol)  ! current concentration
!  DAB don't penalize a low-mass species on a particle:
    if(maxval(conc).ge.cres)then                                    ! cell mass >= minimum c of interest
      itype=0
      ibounds=cat(i,j,k)%bc_number
      if(ibounds.ne.0)itype=bounds(ibounds)%bc_type
      if(itype.ne.3)then                               ! particles outside of absorbing boundary type 3
!........force all particles to a mass resolution (pmass/cmass) of at least 1/nres
! DAB using maxval in case all mass turns to one species (don't force split):
        ip_res=nint(maxval(cat(i,j,k)%cmass/pmass_old)) 
        if(ip_res.lt.nres)then
!.........split particle into enough equal parts to force the above condition false
          nsplit=nres/ip_res+1
          x=pat(1,ip)%xyz(1); y=pat(1,ip)%xyz(2); z=pat(1,ip)%xyz(3)
          pmass=pmass_old/dble(nsplit)
          cat(i,j,k)%cmass=cat(i,j,k)%cmass-pmass_old+pmass
          pat(1,ip)%pmass=pmass
          time=pat(1,ip)%birth_day
          do isplit=1,nsplit-1
            call addp(time,x,y,z,i,j,k,pmass,pat,cat)
!...........assign birth date for split particle to that of parent                    
            pat(1,np)%birth_day=pat(1,ip)%birth_day
          enddo
          netsplit=netsplit+nsplit-1
        endif
      endif
    endif
  endif
enddo
return
end



!-----------------------------------------------------------------
! stream: particle displacement along semianalytical stream lines
!-----------------------------------------------------------------
subroutine stream(vel3,x,y,z,retth,kx,ky,kz,i,j,k,dtmin,xstream,ystream,zstream)
!
! user enters computes particle location (xstream, ystream, and zstream) and time to edge
! dtmin (if particle reaches edge in specified dtmin as input), new locations along the stream
! line a time step tstep away from x,y,z where tstep is constrained by 
! the time to the edge of the cell
!
!
!   IN
!       vel3    - velocity
!       x,y,z   - current particle location
!       retth   - porosity time retardation
!       kx,ky,kz- i-1,j-1,k-1
!       i,j,k   - current cell
!       dtmin   - time step
!   OUT
!       xstream,ystream,zstream - new particle location
!   INOUT
!       dtmin   - min(time to edge of cell, dtmin)
!       i,j,k   - new cell indices if particle reaches cell boundary
!       kx,ky,kz- new cell indices if particle reaches cell boundary
!       
!
! MODFIED: 10/27/01 ELB
!   
use global
implicit none
integer:: kx,ky,kz,ixend,iyend,izend,ix,iy,iz,i,j,k
real:: x,y,z,retth,dtmin,xstream,ystream,zstream
real:: vel3(3,0:nx,0:ny,0:nz)
double precision:: vxm,vym,vzm,vxp,vyp,vzp,xm,ym,zm,xp,yp,zp,xxm,yym,zzm
double precision:: xx,yy,zz,ddx,ddy,ddz,sx,sy,sz,dxstream,dystream,dzstream
double precision:: tstep,xyzstream2
intent(in):: vel3,x,y,z,retth
intent(out):: xstream,ystream,zstream
intent(inout):: i,j,k,kx,ky,kz,dtmin
!
tstep=dtmin
dxstream=x; dystream=y; dzstream=z
ixend=0; iyend=0; izend=0
! compute new location dxstream,dystream,dzstream
! if particle hits edge of cell, time step tstep is returned as time to edge
if(mx.eq.1)then
  ix=0
  vxm=bke*vel3(1,kx,j,k)/retth; vxp=bke*vel3(1,i,j,k)/retth
  xm=kx*dx; xp=xm+dx
  xxm=x-xm
  ddx=dx
  xx=x
  call xyzstream1(dxstream,xx,ddx,xm,xp,xxm,kx,nx,&
  vxm,vxp,tstep,sx,ixend,ix,1)
endif
if(my.eq.1)then
  iy=0
  vym=bke*vel3(2,i,ky,k)/retth; vyp=bke*vel3(2,i,j,k)/retth
  ym=ky*dy; yp=ym+dy
  yym=y-ym
  ddy=dy
  yy=y
  call xyzstream1(dystream,yy,ddy,ym,yp,yym,ky,ny,&
  vym,vyp,tstep,sy,iyend,iy,2)
endif
if(mz.eq.1)then
  iz=0
  vzm=bke*vel3(3,i,j,kz)/retth; vzp=bke*vel3(3,i,j,k)/retth
  zm=kz*dz; zp=zm+dz
  zzm=z-zm
  ddz=dz
  zz=z
  call xyzstream1(dzstream,zz,ddz,zm,zp,zzm,kz,nz,&
  vzm,vzp,tstep,sz,izend,iz,3)
endif
!
! if particle hit edge of cell recalculate position for modified tstep
!
! if it passed the z face, recompute the x and y locations, and new cell
if(izend.ne.0)then
  k=k+iz
  kz=k-1
  if(mx.eq.1)dxstream=xyzstream2(xx,xm,xxm,vxm,tstep,sx)
  if(my.eq.1)dystream=xyzstream2(yy,ym,yym,vym,tstep,sy)
! if it hit the y face, recompute the x location and new cell
elseif(iyend.ne.0)then
  j=j+iy
  ky=j-1
  if(mx.eq.1)dxstream=xyzstream2(xx,xm,xxm,vxm,tstep,sx)
! if it hit the x face, recompute new cell
elseif(ixend.ne.0)then
  i=i+ix
  kx=i-1
endif
dtmin=tstep
xstream=dxstream
ystream=dystream
zstream=dzstream
return
end
!------------------------------------------------------------
subroutine xyzstream1(dxstream,x,dx,xm,xp,xxm,kx,nx,vxm,vxp,&
tstep,sx,ixend,ix,idir)
implicit none
integer:: kx,ixend,ix,nx,idir
double precision:: x,dx,xm,xp,vxm,vxp,z1,xxm
double precision:: sx,tstep,dxstream,xyzstream2
intent(in):: xxm,dx,xm,xp,kx,nx,vxm,vxp,idir
intent(out):: dxstream,sx
intent(inout):: tstep,ixend,ix,x
sx=(vxp-vxm)/dx      
! compute location
dxstream=xyzstream2(x,xm,xxm,vxm,tstep,sx)
! if particle hit edge, compute time to edge
if(dxstream.gt.xp.and.vxp.gt.0.0)then
  ixend=1
  dxstream=xp
! In computing ix, don't let it leave the active grid.
! If it wants to leave, on the next round through it will start from the edge 
! and dtmin will be zero. In this case, icell_correct will bump the particle 
! outside of the domain, where it is either absorbed or reflected.
  if(kx+1.ne.nx)ix=1
  if(sx.ne.0.0)then
    tstep=dlog((sx*dx+vxm)/(sx*xxm+vxm))/sx
  else
    tstep=(real(xp)-real(x))/vxm
  endif
elseif(dxstream.lt.xm.and.vxm.lt.0.0)then
  ixend=1
  dxstream=(xm)
  if(kx.ne.0)ix=-1
  if(sx.ne.0.0)then
    tstep=dlog(vxm/(sx*xxm+vxm))/sx
  else
    tstep=-xxm/vxm
  endif
endif
return
end
!------------------------------------------------------------
double precision function xyzstream2(x,xm,xxm,vxm,tstep,sx)
implicit none
double precision, intent(in):: x,xm,vxm,xxm,sx
double precision, intent(inout):: tstep
double precision:: dxstream,vxmsx
!
if(sx.ne.0.0)then
! the following is hard wired to limit the maximum value of tstep*sx.
! If tstep*sx is too large, dexp(sx*step) cannot be evaluated. 
! These routines will still function is we change tstep here to a reasonable value.
  tstep=min(tstep,600.0d+00/dabs(sx))
  vxmsx=vxm/sx
  xyzstream2=(xxm+vxmsx)*dexp(sx*tstep)-vxmsx+xm
else
  xyzstream2=vxm*tstep+x
endif
return
end



!----------------------------------------------------------------------
!     update velocities t=tnextvel
!     unit = invlc
!
! MODFLOW x-velocity field for a 2-D 3 x 2 system
!
!              ------------------------------------------
!              |             |            |              |
!              |             |            |              |
!              |             | (1,2)      | (2,2)        |
!              |            --->         --->            |
!              |             |            |              |
!              |             |            |              |
!              |             |            |              |
!              ------------------------------------------
!              |             |            |              |
!              |             |            |              |
!              |             | (1,1)      | (2,1)        |  (3,1) = 0
!              |            --->         --->           --->
!              |             |            |              |
!              |             |            |              |
!              |             |            |              |
!              ------------------------------------------
!
! Note: Components (3,j) are all zero because midside velocities
!       cannot be computed at these locations
!----------------------------------------------------------------------
subroutine velupdt(vel3,por,cat,rech,chd)
use global
implicit none
type (cell)::     cat(nx,ny,nz)
type (recharge):: rech(nx,ny)
type (constant_head):: chd(nx,ny,nz)
real vel3(3,0:nx,0:ny,0:nz),por(nzone)
real,allocatable:: temp(:,:,:),cr(:,:,:),cc(:,:,:),cv(:,:,:),head(:,:,:)
character (len=16):: text,textx,texty,textz,textr,textc
character (len=80):: vfile,vfileo,vfile_type
integer checkx,checky,checkz,checkr,checkc,irecharge,ichd
integer, allocatable:: ibound(:,:,:)
integer:: i,j,k,n,kpero,kstpo,kstpin,kperin,kstpino,kperino,iread_error,ierror,ivtype
integer:: ncol,nrow,nlay,kbeg,kend,kstep,jbeg,jend,jstep,iz,idir,k2,j2,&
          imax1,jmax1,kmax1, nflux,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10, &
          imax2,jmax2,kmax2,imax3,jmax3,kmax3, cc1,cc2,cc3,cc4,cc5, &
          cc6,cc7,cc8,cc9,cc10
real:: pert,divflux,divmax,fluxmax,vx,vy,vz,totalt,hdiff,f,f1,f2,f3, &
        maxdivflux,tfm, maxfluxmax, sumfluxmax, divmaxcell, divmaxmaxf, &
        divmaxavgf, divcell, olddiv, divinc, divinc2, ctot,d1,d2,d3
logical flexist
data vfile,vfileo,kpero,kstpo,irecharge,ichd/' ',' ',0,0,0,0/
save 
!.....strings identifying x,y ad z component fluxes in MODFLOW cell-by-cell flow file
textx='FLOW RIGHT FACE'; texty='FLOW FRONT FACE'; textz='FLOW LOWER FACE';textr='        RECHARGE'
textc='   CONSTANT HEAD'
checkx=0; checky=0; checkz=0; checkr=0; checkc=0
vfileo=vfile
call skip(invlc)
!.....tnextvel:      next time to update velocities
!.....kstp:          time step of MODFLOW output, kstp=0 for conductance file
!.....kper:          stress-period of MODFLOW output, kper=0 for conductance file
!.....irevz:         =1 reverse coordinate system to z upwrd (y- and z-directions) of MODFLOW velocities
!.....ivtype:        1=binary cbc file,2=binary head file,-1=unformatted cbc, -2=unformatted head
!.....vfile:         MODFLOW cell-by-cell flow file name
1 continue
read(invlc,*,iostat=iread_error)tnextvel,kstp,kper,irevz,ivtype,vfile
if(iread_error.eq.-1)then
! end of velocity input, make tnextvel large number
  tnextvel=dble(large)
  write(*,*)' '
  return
elseif(iread_error.eq.-2)then
  goto 9999
endif
! set file type
vfile_type='binary'
if(ivtype.lt.0)then
  vfile_type='unformatted'
endif
!
if(tnextvel.eq.tmax)tnextvel=dble(large)             
!
if(iabs(ivtype).eq.2)then
!.compute velocities from heads
  if(kstp.eq.0)then
!...read the conductance file to compute velocities from heads
    write(iout,1000)vfile
    write(   *,1000)vfile
    if(.not.allocated(cr))allocate (ibound(1:nx,1:ny,1:nz),STAT = ierror)
    if(.not.allocated(cr))allocate (cr(1:nx,1:ny,1:nz),STAT = ierror)
    if(.not.allocated(cc))allocate (cc(1:nx,1:ny,1:nz),STAT = ierror)
    if(.not.allocated(cv))allocate (cv(1:nx,1:ny,1:nz),STAT = ierror)
!...open conductance file
    inquire(file=vfile,exist=flexist)
    if(flexist)then
      open(invel,file=vfile,status='old',form=vfile_type)
    else
      write(iout,*)' error: conductance file does not exist'
      stop ' error: conductance file does not exist'
    endif
    call readcc(iout,invel,ibound,cr,cc,cv,nx,ny,nz,iread_error)
    if(iread_error.lt.0)goto 9998
!...now go back and read the velocity file info
    goto 1
  else
!...echo head file info
    write(iout,1500)vfile,kstp,kper,abs(sngl(tnextvel)),irevz
    write(   *,1500)vfile,kstp,kper,abs(sngl(tnextvel)),irevz
  endif
else
  write(iout,2000)vfile,kstp,kper,abs(sngl(tnextvel)),irevz
  write(   *,2000)vfile,kstp,kper,abs(sngl(tnextvel)),irevz
endif
!
kbeg=1
kend=nz
kstep=1
jbeg=1
jend=ny
jstep=1
if(irevz.eq.1)then
  kbeg=nz
  kend=1
  kstep=-1
  jbeg=ny
  jend=1
  jstep=-1
endif
!.......open new velocity or head files?
if(kper.lt.kpero.or.vfileo.ne.vfile.or.((kper.eq.kpero).and.(kstp.le.kstpo)))then
  irecharge=0
  close(invel)
  inquire(file=vfile,exist=flexist)
  if(flexist)then
    open(invel,file=vfile,status='old',form=vfile_type)
  else
    write(iout,*)' error: velocity file does not exist'
    stop ' error: velocity file does not exist'
  endif
!.check to see if file contains recharge or chd fluxes
! set checkr=1 if recharge fluxes are not present
  if(iabs(ivtype).ne.2)then  ! if flux file
    if(.not.allocated(temp))allocate (temp(1:nx,1:ny,1:nz),STAT = ierror)
    i=0
    do 
      call readcbc(iout,invel,nx,ny,nz,kstpin,text,kperin,ncol,nrow,nlay,temp,iread_error)
      if(iread_error.eq.-1)then
        exit
      elseif(iread_error.eq.-2)then
        goto 9998
      endif             
      if(i.ne.0)then
        if(kstpino.ne.kstpin)exit
        if(kperino.ne.kperin)exit
      endif
      kstpino=kstpin
      kperino=kperin
      if(text.eq.textr)then
        irecharge=1
      endif
      if(text.eq.textc)then
        ichd=1
      endif
      i=i+1
    enddo
    rewind(invel) ! rewind file for reading
    if(irecharge.eq.0)then
      checkr=1 ! don't check to make sure recharge is read
      write(*,*)' No recharge fluxes in cbc file'
      write(iout,*)' No recharge fluxes in cbc file'
    else
      write(*,*)' Found recharge fluxes in cbc file'
      write(iout,*)' Found recharge fluxes in cbc file'
    endif
    if(ichd.eq.0)then
      checkc=1 ! don't check to make sure chd is read
      write(*,*)' No chd fluxes in cbc file'
      write(iout,*)' No chd fluxes in cbc file'
    else
      write(*,*)' Found chd fluxes in cbc file'
      write(iout,*)' Found chd fluxes in cbc file'
      checkc=1 ! don't check to make sure chd is read
    endif
  endif        
endif
vfileo=vfile
kpero=kper
kstpo=kstp
!.....read flux
if(.not.allocated(temp))allocate (temp(1:nx,1:ny,1:nz),STAT = ierror)
if(iabs(ivtype).eq.2)then
!.read the head file
  do k=1,nz
   read(invel,iostat=iread_error)kstpin,kperin,pert,totalt,text,ncol,nrow,iz
   if(ncol.ne.nx.or.nrow.ne.ny)then
    write(*,1001)ncol,nrow,nz
    write(iout,1001)ncol,nrow,nz
    stop
   endif
   if(iread_error.eq.-1)then
     exit
   elseif(iread_error.eq.-2)then
     goto 9998
   endif             
   if(.not.allocated(head))allocate (head(nx,ny,nz))
   write(*,1002)iz,nx,ny,kstp,kper,char(13)
   1002   format('Reading layer ',i3,' NCOL =',I3,' NROW =',I3,&
                 ' kstp ',i3,' kper ',i3,a1,$)
   read(invel,iostat=iread_error)((head(i,j,iz),i=1,nx),j=1,ny)
   if(iread_error.eq.-1)then
     exit
   elseif(iread_error.eq.-2)then
     goto 9998
   endif             
  enddo
endif
idir=1
do
  if(iabs(ivtype).ne.2)then
!.read fluxes
    call readcbc(iout,invel,nx,ny,nz,kstpin,text,kperin,ncol,nrow,nlay,temp,iread_error)
    if(iread_error.eq.-1)then
      exit
    elseif(iread_error.eq.-2)then
      goto 9998
    endif             
  else
!...compute fluxes from heads
    if(idir.eq.1)then
      temp=0.0
      do k=1,nz
      do j=1,ny
      do i=1,nx-1
      if((ibound(i,j,k).ne.0) .AND. (ibound(i+1,j,k).ne.0))then
        hdiff=head(i,j,k)-head(i+1,j,k)
        temp(i,j,k)=hdiff*cr(i,j,k)
      endif
      enddo
      enddo
      enddo
      text=textx
    elseif(idir.eq.2)then
      do k=1,nz
      do j=1,ny-1
      do i=1,nx
      if((ibound(i,j,k).ne.0) .AND. (ibound(i,j+1,k).ne.0))then
        hdiff=head(i,j,k)-head(i,j+1,k)
        temp(i,j,k)=hdiff*cc(i,j,k)
      endif
      enddo
      enddo
      enddo
      text=texty
    elseif(idir.eq.3)then
      temp=0.0
      do k=1,nz-1
      do j=1,ny
      do i=1,nx
      if((ibound(i,j,k).ne.0) .AND. (ibound(i,j,k+1).ne.0))then
        hdiff=head(i,j,k)-head(i,j,k+1)
        temp(i,j,k)=hdiff*cv(i,j,k)
      endif
      enddo
      enddo
      enddo
      text=textz
    endif
!.....
  endif
  if(text.eq.textx.and.kstp.eq.kstpin.and.kper.eq.kperin.and.mx.ne.0)then ! check for flux component, time and stress period
    write(*,'(/a)')'................Assigning x-velocities' !,char(13)
    k2=0
    do k=kbeg,kend,kstep
    k2=k2+1
    j2=0
    do j=jbeg,jend,jstep
    j2=j2+1
    do i=1,nx
      if(irevz.eq.1)vel3(1,i,j2,k2)=temp(i,j,k)
      if(irevz.eq.0)vel3(1,i,j,k)=temp(i,j,k)
    enddo
    enddo
    enddo
    checkx=1
  elseif(text.eq.texty.and.kstp.eq.kstpin.and.kper.eq.kperin.and.my.ne.0)then
    write(*,'(/a)')'................Assigning y-velocities' !,char(13)
    k2=0
    do k=kbeg,kend,kstep
    k2=k2+1
    j2=0
    do j=jbeg,jend,jstep
    j2=j2+1
    do i=1,nx
      if(j.ne.1.and.irevz.eq.1)vel3(2,i,j2,k2)=-temp(i,j-1,k)
      if(j.eq.1.and.irevz.eq.1)vel3(2,i,j2,k2)=0.0
      if(irevz.eq.0)vel3(2,i,j,k)=temp(i,j,k)
    enddo
    enddo
    enddo
    checky=1
  elseif(text.eq.textz.and.kstp.eq.kstpin.and.kper.eq.kperin.and.mz.ne.0)then
    write(*,'(/a)')'................Assigning z-velocities' !,char(13)
    k2=0
    do k=kbeg,kend,kstep
    k2=k2+1
    j2=0
    do j=jbeg,jend,jstep
    j2=j2+1
    do i=1,nx
      if(k.ne.1.and.irevz.eq.1)vel3(3,i,j2,k2)=-temp(i,j,k-1)
      if(k.eq.1.and.irevz.eq.1)vel3(3,i,j2,k2)=0.0
      if(irevz.eq.0)vel3(3,i,j,k)=temp(i,j,k)
    enddo
    enddo
    enddo
    checkz=1
  elseif(irecharge.eq.1.and.text.eq.textr.and.kstp.eq.kstpin.and.kper.eq.kperin.and.mz.ne.0)then
    write(*,'(/a)')'..............Assigning rch-velocities' !,char(13)
    do i=1,nx
    j2=0
    do j=jbeg,jend,jstep
    j2=j2+1
    k2=0
    do k=kbeg,kend,kstep
      k2=k2+1
      if(temp(i,j,k).ne.0)then
        rech(i,j2)%flow=temp(i,j,k)
        rech(i,j2)%k=k2
        exit
      else
        rech(i,j2)%flow=0.0
        rech(i,j2)%k=1
      endif
    enddo
    enddo
    enddo
    checkr=1
  elseif(ichd.eq.1.and.text.eq.textc.and.kstp.eq.kstpin.and.kper.eq.kperin.and.mz.ne.0)then
    write(*,'(/a)')'..............Assigning chd-fluxes' !,char(13)
    k2=0
    do k=kbeg,kend,kstep
    k2=k2+1
    j2=0
    do j=jbeg,jend,jstep
    j2=j2+1
    do i=1,nx
      chd(i,j2,k2)%flow=temp(i,j,k)
    enddo
    enddo
    enddo
    checkc=1
  endif
  idir=idir+1
  if((checkx.eq.0.and.mx.eq.1).or.(checky.eq.0.and.my.eq.1).or.(checkz.eq.0.and.mz.eq.1).or.&
     (irecharge.eq.1.and.mz.eq.1.and.checkr.eq.0).or.(ichd.eq.1.and.checkc.eq.0))then     
     cycle
  else
     exit
  endif
enddo
!.....some additional error checking
if((checkx.eq.0.and.mx.eq.1).or.(checky.eq.0.and.my.eq.1).or.(checkz.eq.0.and.mz.eq.1))then
  write(*,2003)
  write(iout,2003)
  stop
endif
deallocate(temp)
!......compute divergence
if(ibug.ge.1)then
  write(ibugout,*)'             Divergence of the Flux'     
  write(ibugout,*)'        i          j          k         div'
  maxfluxmax=0.0     
  maxdivflux=0.0     
  sumfluxmax=0.0     
  divmaxcell=0.0
  nflux=0
  do k=1+mz,nz-mz; do j=1+my,ny-my; do i=1+mx,nx-mx
!    if(i.gt.mx.and.i.lt.nx.and.j.gt.1.and.j.lt.ny.and.k.gt.1.and.k.lt.nz)then
      nflux=nflux+1
	  divflux=vel3(1,i-1,j,k)-vel3(1,i,j,k)+vel3(2,i,j-1,k)-&
              vel3(2,i,j,k)+vel3(3,i,j,k-1)-vel3(3,i,j,k)
      fluxmax=max(abs(vel3(1,i-1,j,k)),abs(vel3(1,i,j,k)),&
                  abs(vel3(2,i,j-1,k)),abs(vel3(2,i,j,k)),&
                  abs(vel3(3,i,j,k-1)),abs(vel3(3,i,j,k)))
      if(ibug.ge.3.and.fluxmax.ne.0.0)write(ibugout,'(3(i5,1x),e11.5,1x,e11.5)')i,j,k,divflux,divflux/fluxmax
      maxfluxmax=max(maxfluxmax,fluxmax)
      maxdivflux=max(maxdivflux,abs(divflux))
      sumfluxmax=sumfluxmax+fluxmax
      if(fluxmax.ne.0.0)divmaxcell=max(divmaxcell,abs(divflux)/fluxmax)
!    endif
  enddo; enddo; enddo
  if(maxfluxmax.ne.0.0)divmaxmaxf=maxdivflux/maxfluxmax
  if(sumfluxmax.ne.0.0)divmaxavgf=maxdivflux/(sumfluxmax/float(nflux))
  write(ibugout,'(a)')'----------------------------------------------------------------------------'
  write(ibugout,'(a)')'                                  Normalization by'
  write(ibugout,'(a)')'                                  Max. Flux   Max. Flux  By Avg. Flux'
  write(ibugout,'(a)')'                                   in Cell    in Domain   in Domain  '
  write(ibugout,'(a)')'                                  ----------- ----------- -----------'
  write(ibugout,'(a,3(1x,e11.5))')' Maximum normalized divergence = ',divmaxcell,divmaxmaxf,divmaxavgf
  write(ibugout,'(a)')'----------------------------------------------------------------------------'
  write(*      ,'(a)')'----------------------------------------------------------------------------'
  write(*      ,'(a)')'                                  Normalization by'
  write(*      ,'(a)')'                                  Max. Flux   Max. Flux  By Avg. Flux'
  write(*      ,'(a)')'                                   in Cell    in Domain   in Domain  '
  write(*      ,'(a)')'                                  ----------- ----------- -----------'
  write(*,'(a,3(1x,e11.5))')' Maximum normalized divergence = ',divmaxcell,divmaxmaxf,divmaxavgf
  write(*      ,'(a)')'----------------------------------------------------------------------------'
  endif
  !.....
  do k=1,nz; do j=1,ny; do i=1,nx
!.......convert flux to Darcy velocity
    vx=vel3(1,i,j,k)/ax; vy=vel3(2,i,j,k)/ay; vz=vel3(3,i,j,k)/az
!.......set minimum velocities for anisotropic dispersion tensor
    vel3(1,i,j,k)=vx; vel3(2,i,j,k)=vy; vel3(3,i,j,k)=vz 
  enddo; enddo; enddo
!..adjust velocities for recharge
  if(irecharge.eq.1)then
    do j=1,ny; do i=1,nx
! adjust all of the way from the water table to the top of the model
      if(irevz.eq.1)then
        do k=rech(i,j)%k,nz
          vel3(3,i,j,k)=vel3(3,i,j,k)-rech(i,j)%flow/az
        enddo
      else    
        do k=1,rech(i,j)%k-1
          vel3(3,i,j,k)=vel3(3,i,j,k)+rech(i,j)%flow/az
        enddo
      endif
    enddo; enddo
  endif
!.....some additional error checking
  if((checkx.eq.0.and.mx.eq.1).or.(checky.eq.0.and.my.eq.1).or.(checkz.eq.0.and.mz.eq.1))then
    write(*,2003)
    write(iout,2003)
    stop
  endif
return
 1000 format(/' MODFLOW conductance file        ',20('.'),3x,a//)
 1001 format(/1x,70('-')/1x,'ERROR IN CONDUCTANCE FILE DIMENSIONS '/&
                   1x,'NCOL,NROW,NLAY:',i5,i5,i5/70('-'))
 1500 format(/' MODFLOW head file               ',20('.'),3x,a/&
       '                  time step      ',20('.'),3x,i15/&
       '                  stress period  ',20('.'),3x,i15/&
       '                  tnextvel       ',20('.'),3x,e15.8/&
       '                  irevz          ',20('.'),3x,i15//)
 2000 format(/' MODFLOW velocity file           ',20('.'),3x,a/&
       '                  time step      ',20('.'),3x,i15/&
       '                  stress period  ',20('.'),3x,i15/&
       '                  tnextvel       ',20('.'),3x,e15.8/&
       '                  irevz          ',20('.'),3x,i15//)
 2001 format(/1x,70('-')/1x,'ERROR IN VELOCITY FILE DIMENSIONS '/&
                   1x,'NCOL,NROW,NLAY:',i5,i5,i5/70('-'))
 2002 format(/1x,70('-')/1x,'ERROR: CELL-BY-CELL FLOWS IN ',a1,'-DIRECTION '//&
          ' DIMENSION:',i5/1x,70('-'))
 2003 format(/1x,70('-')/1x,'ERROR READING VELOCITIES'/70('-'))
 2004 format(/1x,70('-')/1x,'INDICATOR SIMULATION FILE',a50/&
   ' DOES NOT EXIST, ASSIGNING CONSTANT POROSITY'/70('-'))
 2005 format(/1x,70('-')/1x,'MORE CATEGORIES THAN POROSITY'//&
                      ' VALUES'/70('-'))
 2006 format(/1x,70('-')/1x,'DIMENSIONS OF INDICATOR FILE   ',3I10/&
                      'DO NOT MATCH MODEL DIMENSIONS ',3I10/70('-'))
 9999 stop 'Error in Velocity Input File'
 9998 stop 'Error in file containing velocities'
end

!------------------------------------------------------------
! read conductances
!------------------------------------------------------------
subroutine readcc(iout,invel,ibound,cr,cc,cv,nx,ny,nz,iread_error)
implicit none
integer:: ncol,nrow,nlay,iread_error,nx,ny,nz,i,j,k,iout,invel,ibound(nx,ny,nz)
character (len=16) text
real:: cr(nx,ny,nz),cc(nx,ny,nz),cv(nx,ny,nz)

!...read ibound
read(invel,iostat=iread_error)text,ncol,nrow,nlay
if(iread_error.lt.0)return
write(*,'(a,a,a1,$)')' READING:    ',text,char(13)
if(ncol.ne.nx.or.nrow.ne.ny.or.nlay.ne.nz)then
  write(*,1001)ncol,nrow,nlay
  write(iout,1001)ncol,nrow,nlay
  stop
endif
read(invel,iostat=iread_error)(((ibound(i,j,k),i=1,ncol),j=1,nrow),k=1,nlay)
if(iread_error.lt.0)return
!...read cr
read(invel,iostat=iread_error)text,ncol,nrow,nlay
if(iread_error.lt.0)return
write(*,'(a,a,a1,$)')' CONDUCTANCE FOR:    ',text,char(13)
if(ncol.ne.nx.or.nrow.ne.ny.or.nlay.ne.nz)then
  write(*,1001)ncol,nrow,nlay
  write(iout,1001)ncol,nrow,nlay
  stop
endif
read(invel,iostat=iread_error)(((cr(i,j,k),i=1,ncol),j=1,nrow),k=1,nlay)
if(iread_error.lt.0)return
!...read cc
read(invel,iostat=iread_error)text,ncol,nrow,nlay
if(iread_error.lt.0)return
write(*,'(a,a,a1,$)')' CONDUCTANCE FOR:    ',text,char(13)
if(ncol.ne.nx.or.nrow.ne.ny.or.nlay.ne.nz)then
  write(*,1001)ncol,nrow,nlay
  write(iout,1001)ncol,nrow,nlay
  stop
endif
read(invel,iostat=iread_error)(((cc(i,j,k),i=1,ncol),j=1,nrow),k=1,nlay)
if(iread_error.lt.0)return
!...read cv
read(invel,iostat=iread_error)text,ncol,nrow,nlay
if(iread_error.lt.0)return
write(*,'(a,a,a1,$)')' CONDUCTANCE FOR:    ',text,char(13)
if(ncol.ne.nx.or.nrow.ne.ny.or.nlay.ne.nz)then
  write(*,1001)ncol,nrow,nlay
  write(iout,1001)ncol,nrow,nlay
  stop
endif
read(invel,iostat=iread_error)(((cv(i,j,k),i=1,ncol),j=1,nrow),k=1,nlay)
if(iread_error.lt.0)return
1001 format(/1x,70('-')/1x,'ERROR IN CONDUCTANCE FILE DIMENSIONS '/&
                   1x,'NCOL,NROW,NLAY:',i5,i5,i5/70('-'))
return
end

!------------------------------------------------------------
! read cell by cell flux terms
!------------------------------------------------------------
subroutine readcbc(iout,invel,nx,ny,nz,kstpin,text,kperin,ncol,nrow,nlay,temp,iread_error)
implicit none
integer:: kstpin,kperin,ncol,nrow,nlay,iread_error,nx,ny,nz,i,j,k,n,iout,invel
character (len=16) text
real:: temp(nx,ny,nz),f,f1,f2,f3
read(invel,iostat=iread_error)kstpin,kperin,text,ncol,nrow,nlay
if(iread_error.lt.0)return
! MODFLOW 2000 fix; may need further modification in some cases 
if(nlay.lt.0)then ! read additional data from MODFLOW 2000
  nlay=-nlay
  read(invel,iostat=iread_error) i,f1,f2,f3
endif
!
write(*,'(a,a,i6,i6,a1,$)')' CELL-BY-CELL FLOW:    ',text,kstpin,kperin,char(13)
write(iout,'(a,a,i6,i6)')' CELL-BY-CELL FLOW:    ',text,kstpin,kperin
if(ncol.ne.nx.or.nrow.ne.ny.or.nlay.ne.nz)then
  write(*,2001)ncol,nrow,nlay
  write(iout,2001)ncol,nrow,nlay
  stop
endif
! MODFLOW 2000 fix; may need further modification in some cases 
if(i.eq.2)then
  read(invel,iostat=iread_error)n
  if(iread_error.lt.0)return
  do i=1,n
    read(invel,iostat=iread_error)f
    if(iread_error.lt.0)return
  enddo
else
  read(invel,iostat=iread_error)(((temp(i,j,k),i=1,ncol),j=1,nrow),k=1,nlay)
  if(iread_error.lt.0)return
endif
2001 format(/1x,70('-')/1x,'ERROR IN VELOCITY FILE DIMENSIONS '/&
                   1x,'NCOL,NROW,NLAY:',i5,i5,i5/70('-'))
return
end
!------------------------------------------------------------
! update boundary conditions
!------------------------------------------------------------
subroutine bndupdt(bounds,source,cat,pat)
! update BCs
use global
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
integer source(maxsource)
double precision pmass(nspec),tbeg,tend
!
nsource=0
tnextbnd=dble(large)
do ibounds=1,nbounds
  itype=bounds(ibounds)%bc_type
  if(itype.gt.0)then
    tbeg=bounds(ibounds)%tbeg
    tend=bounds(ibounds)%tend
    if(itype.ne.2.and.itype.ne.8.and.itype.ne.9.and.itype.ne.10.and.itype.ne.11)then
      i=bounds(ibounds)%ijk(1)
      j=bounds(ibounds)%ijk(2)
      kbot=bounds(ibounds)%kbot
      ktop=bounds(ibounds)%ktop
    endif
!...if boundary is active set it by referencing in icat
    if(curtime.eq.tbeg)then
     if(itype.ne.2.and.itype.ne.8.and.itype.ne.9.and.itype.ne.10.and.itype.ne.11)then
        do k=kbot,ktop
          cat(i,j,k)%bc_number=ibounds
        enddo
     endif
    endif
!...if curtime=tend, set itype=0 and cat(i,j,k)%bc_number = 0
    if(curtime.eq.tend)then
      if(itype.ne.2.and.itype.ne.8.and.itype.ne.9.and.itype.ne.10.and.itype.ne.11)then
        do k=kbot,ktop
          cat(i,j,k)%bc_number=0
        enddo
      endif
!     if(itype.ne.2)cat(i,j,k)%bc_number=0
      bounds(ibounds)%bc_type=-bounds(ibounds)%bc_type
      if(itype.eq.1)then
!.........distribute constant source one last time before turning off
          npart=bounds(ibounds)%np
          bounds(ibounds)%np_remove=bounds(ibounds)%np_remove+npart
          pmass=bounds(ibounds)%pmass
          mass=mass+npart*pmass
          npbc(itype)=npbc(itype)+npart
        do iloop=1,nspec
           massbc(itype,iloop)=massbc(itype,iloop)+pmass(iloop)*npart
        enddo
!  old        massbc(itype)=massbc(itype)+pmass*npart
          xm=bounds(ibounds)%bc_xyzm(1)
          ym=bounds(ibounds)%bc_xyzm(2)
          zm=bounds(ibounds)%bc_xyzm(3)
          sx=bounds(ibounds)%bc_xyzs(1)
          sy=bounds(ibounds)%bc_xyzs(2)
          sz=bounds(ibounds)%bc_xyzs(3)
          iptype=0
          snglcurtime=curtime
          call placeu(snglcurtime,xm,ym,zm,sx,sy,sz,pat,cat,npart,pmass,iptype,ierr)
      endif
    endif
!...set source/concentration BCs
    if((itype.eq.1.or.itype.eq.2.or.itype.eq.6.or.itype.eq.8.or.itype.eq.9.or.itype.eq.10.or.itype.eq.11).and.&
    (curtime.ge.tbeg.and.curtime.lt.tend))then
      nsource=nsource+1
      source(nsource)=ibounds
    endif
!...find next time to update bounds 
    if(itype.gt.0.and.curtime.ge.tbeg.and.curtime.ne.tend)tnextbnd=min(tend,tnextbnd)
    if(itype.gt.0.and.curtime.lt.tbeg)tnextbnd=min(tbeg,tnextbnd)
  endif
enddo
return
end
!------------------------------------------------------------
! update fluxin and fluxout for maw conditions
!------------------------------------------------------------
subroutine mawupdt(bounds,vel3)
! compute fluxes for MAW boundaries after each call to update velocity field
use global
implicit none
type (boundary):: bounds(1:maxbnd)
double precision pmass(nspec),tbeg,tend
real:: vel3(3,0:nx,0:ny,0:nz)
integer:: numbndo,numbnd,i,j,k,ibounds,itype,itypeo,ibou,numb,kbot,ktop
real:: fluxin,fluxout,q   ! DAB make sure fluxes are only water fluxes
!
numbndo=-9999
fluxin=0.0
fluxout=0.0
itype=0
itypeo=-9999
do ibounds=1,nbounds
  itype=bounds(ibounds)%bc_type
  numbnd=bounds(ibounds)%group
  if(ibounds.ne.1)then  
    if(numbnd.ne.numbndo.and.itypeo.eq.7)then
!     assign fluxin and fluxout to all previous boundaries in this group
      ibou=ibounds-1 
      numb=bounds(ibou)%group 
      do
        bounds(ibou)%fluxin=fluxin
        bounds(ibou)%fluxout=fluxout
        if(ibou.eq.1)exit
        ibou=ibou-1
        if(bounds(ibou)%group.ne.numb)exit
      enddo       
!     reinitialize fluxin and fluxout
      fluxin=0.0
      fluxout=0.0
    endif
  endif
  if(itype.eq.7)then
    i=bounds(ibounds)%ijk(1)
    j=bounds(ibounds)%ijk(2)
    kbot=bounds(ibounds)%kbot
    ktop=bounds(ibounds)%ktop
!   compute cell flux
    do k=kbot,ktop
      q=bke*((vel3(1,i-1,j,k)-vel3(1,i,j,k))*ax+&
             (vel3(2,i,j-1,k)-vel3(2,i,j,k))*ay+&
             (vel3(3,i,j,k-1)-vel3(3,i,j,k))*az)
      if(q.gt.0)fluxin=fluxin+q
      if(q.lt.0)fluxout=fluxout+q
    enddo
    if(ibounds.eq.nbounds)then
!     assign fluxin and fluxout to all previous boundaries in this group
      ibou=ibounds 
      numb=bounds(ibou)%group 
      do
        bounds(ibou)%fluxin=fluxin
        bounds(ibou)%fluxout=fluxout
        if(ibug.ge.1)write(ibugout,*)' Boundry #',ibounds,' Boundary Type 7 fluxin = ',fluxin
        if(ibug.ge.1)write(ibugout,*)' Boundry #',ibounds,' Boundary Type 7 fluxout = ',fluxout
        ibou=ibou-1
        if(ibou.le.1)exit
        if(bounds(ibou)%group.ne.numb)exit
      enddo       
    endif
  endif
  numbndo=numbnd  
  itypeo=itype
enddo
return
end
!------------------------------------------------------------
! source
!------------------------------------------------------------
subroutine srcupdt(source,bounds,pat,cat)
use global
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
! distribute particles in source's
integer source(maxsource)
double precision pmass(nspec)

!.determine next time to update constant conc. bounds
tnextsrc=dble(large)
do isource=1,nsource
  ibounds=source(isource)
  itype=bounds(ibounds)%bc_type
!.determine extent of patch source
  if(itype.eq.1)then
!...update all constant conc. sources at minimum courant condition
    tnextsrc=min(tnextsrc,curtime+dble(bounds(ibounds)%dt_type1),dble(bounds(ibounds)%tend))
    i=bounds(ibounds)%ijk(1); j=bounds(ibounds)%ijk(2); kbot=bounds(ibounds)%kbot; ktop=bounds(ibounds)%ktop
    npart=bounds(ibounds)%np
    pmass=bounds(ibounds)%pmass
    mass=mass+npart*pmass
    bounds(ibounds)%np_remove=bounds(ibounds)%np_remove+npart
    bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove+pmass*npart
    npbc(itype)=npbc(itype)+npart
     do iloop=1,nspec
        massbc(itype,iloop)=massbc(itype,iloop)+pmass(iloop)*npart
     enddo
! DAB old line  massbc(itype)=massbc(itype)+pmass*npart
          xm=bounds(ibounds)%bc_xyzm(1)
          ym=bounds(ibounds)%bc_xyzm(2)
          zm=bounds(ibounds)%bc_xyzm(3)
          sx=bounds(ibounds)%bc_xyzs(1)
          sy=bounds(ibounds)%bc_xyzs(2)
          sz=bounds(ibounds)%bc_xyzs(3)
    iptype=0
    call placeu(sngl(curtime),xm,ym,zm,sx,sy,sz, pat,cat,npart,pmass,iptype,ierr)          
  endif
enddo
return
end
!------------------------------------------------------------
! pntupdt
!------------------------------------------------------------
subroutine pntupdt(pat,imp,cat,vel3)
use global
implicit none
type (particle):: pat(1,maxnp)
type (imparticle):: imp(1,maxnp)
type (cell)::     cat(nx,ny,nz)
!.....read point source information
integer:: ipntsrc,indomain,iread_error,kx,ky,kz,i,j,k
integer:: inppnt,nppnt,iptype
integer:: ptype
real:: x,y,z,xm,ym,zm,xp,yp,zp,xmax,ymax,zmax,xmin,ymin,zmin
real:: vel3(3,0:nx,0:ny,0:nz)
double precision:: pmass(nspec)
!
do ipntsrc=1,npntsrc

  call skip(inpnt)
  read(inpnt,*,err=9999)xmin,xmax,ymin,ymax,zmin,zmax,nppnt,pmass,ptype  

!.check if particle is in domain
  if(x.lt.0.0.or.xmax.gt.nx*dx.or.& 
     y.lt.0.0.or.ymax.gt.ny*dy.or.&
     z.lt.0.0.or.zmax.gt.nz*dz)goto 9998

     if(ptype.eq.0) mass=mass+nppnt*pmass  ! DAB count mobile mass

     npbc(nbtype+1)=npbc(nbtype+1)+nppnt
     do iloop=1,nspec
        massbc(nbtype+1,iloop)=massbc(nbtype+1,iloop)+nppnt*pmass(iloop)
     enddo

! DAB Do this point by point to add random positions between x and xmax, etc.
   do inppnt=1,nppnt
   x=xmin+rand()*(xmax-xmin)
   y=ymin+rand()*(ymax-ymin)
   z=zmin+rand()*(zmax-zmin)

!.compute cell from x,y,z
  kx=ifix(x/dx); ky=ifix(y/dy); kz=ifix(z/dz)
  i=kx+1; j=ky+1; k=kz+1
  xm=kx*dx; ym=ky*dy; zm=kz*dz; xp=xm+dx; yp=ym+dy; zp=zm+dz
  call icell_correct(kx,ky,kz,i,j,k,x,y,z,xm,ym,zm,xp,yp,zp,vel3)
  if(i.lt.1.or.j.lt.1.or.k.lt.1.or.i.gt.nx.or.j.gt.ny.or.k.gt.nz)then
    print*,' PARTICLE OUTSIDE OF DOMAIN'
    write(iout,*)' PARTICLE OUTSIDE OF DOMAIN'
    print*,' xmin,ymin,zmin',0.0,0.0,0.0
    print*,' xmax,ymax,zmax',dx*nx,dy*ny,dz*nz
    print*,'    x,   y,   z',x,y,z
    write(iout,*)' PARTICLE OUTSIDE OF DOMAIN'
    cycle
  endif
!
!  massbc(nbtype+1)=massbc(nbtype+1)+nppnt*pmass
! DAB moved up this loop start ... do inppnt=1,nppnt
      if(np+1.gt.maxnp)then
        write(iout,*)' error: maximum # of particles exceeded in pntupdt'
        stop ' error: maximum # of particles exceeded in pntupdt'
      endif
      if(ptype.eq.q)call addimp(sngl(curtime),x,y,z,i,j,k,pmass,imp,cat)
      if(ptype.eq.0)call   addp(sngl(curtime),x,y,z,i,j,k,pmass,pat,cat)
  enddo
enddo
call skip(inpnt)
read(inpnt,*,iostat=iread_error)tnextpnt,npntsrc
if(iread_error.eq.0)then
  return
elseif(iread_error.eq.-1)then
  tnextpnt=dble(large)
  return
else
  goto 9999
endif
9999 write(*,*)' Error reading point source '
     write(iout,*)' Error reading point source '
     stop ' Error in point source input file'
9998 write(*,*)' Point source out of domain '
     write(iout,*)' Specified point source out of domain '
     stop ' Error in point source input file'
     end
!----------------------------------------------------------------------
!     output control input
!     unit = inopc
!----------------------------------------------------------------------
subroutine opcupdt(opc,outunit,outfname)
use global
integer opc(nopc),outunit(nopc+1),iread_error
character (len=80) outfname(nopc),flname
character (len=3) llopc(nopc)
data iflag/1/
if(iflag.eq.1)then 
  do iopc=1,nopc 
    call skip(inopc)
    read(inopc,'(a)',err=9999)outfname(iopc)
    if(iopc.eq.7)then
!     open two monitoring files
      open(outunit(iopc),file=outfname(iopc),status='unknown')
      flname=outfname(iopc)
      flname=flname(1:index(flname,' ')-1)//'m'
      open(outunit(iopc+1),file=flname,status='unknown')
    elseif(iopc.eq.5)then
!     prt file is unformatted
      open(outunit(iopc),file=outfname(iopc),form='unformatted',status='unknown')      
    elseif(iopc.eq.2)then
!     concentration file is unformatted
!     DAB 10-2-17 Not any more!
      open(outunit(iopc),file=outfname(iopc),status='unknown')
!      open(outunit(iopc),file=outfname(iopc),form='unformatted',status='unknown')
    else
!     open other files
      open(outunit(iopc),file=outfname(iopc),status='unknown')
    endif
! ELB5-9-08    if(iopc.ne.2)open(outunit(iopc),file=outfname(iopc),status='unknown')
! ELB5-9-08    if(iopc.eq.2)open(outunit(iopc),file=outfname(iopc),form='unformatted',status='unknown')
  enddo
  write(iout,1000)(OUTFNAME(i),i=1,nopc)
  write(*,1000)(OUTFNAME(i),i=1,nopc)
  iflag=0
endif
!
call skip(inopc)
do
  read(inopc,*,iostat=iread_error)tnextopc,(opc(iopc),iopc=1,nopc)
  if(iread_error.eq.-1)then
    tnextopc=dble(large)
    return
  elseif(iread_error.eq.-2)then
    goto 9999
  endif
  ioutput=0
  do iopc=1,nopc
    if(opc(iopc).gt.0)then
      ioutput=1
      llopc(iopc)='YES'
    else
      llopc(iopc)=' NO'
    endif
  enddo
  if(ioutput.ne.0)exit
enddo
return
1000  format('                           O U T P U T   F I L E S '/&
             ' moments                    :',3x,a/&
             ' concentrations             :',3x,a/&
             ' macro-dispersion tensor    :',3x,a/&
             ' breakthrough counters      :',3x,a/&
             ' particle locations         :',3x,a/&
             ' internal boundary counters :',3x,a/&
             ' monitoring                 :',3x,a//)

return
9999 write(iout,*)' Error reading output control file '
write(*,*)' Error in output control file '
stop
end








!------------------------------------------------------------------      
! placeu
!------------------------------------------------------------------
subroutine placeu(time,xm,ym,zm,sx,sy,sz,pat,cat,npart,pmass,iptype,ierr)
! place npart particles uniformly in a block starting at xyzm and ending at sxyz
use global
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
double precision pmass(nspec)
ierr=0
do ip=1,npart
    if(np+1.gt.maxnp)then
      write(iout,*)' error:too many particles, increase maxnp'
      stop    ' error:too many particles, increase maxnp'
    endif
    x=xm+sx*randu01(); y=ym+sy*randu01(); z=zm+sz*randu01()
    kx=ifix(x/dx); ky=ifix(y/dy); kz=ifix(z/dz)
    if(kx.lt.0.or.ky.lt.0.or.kz.lt.0.or.kx.gt.nx-1.or.ky.gt.ny-1.or.kz.gt.nz-1)then
      ierr=-1
      return
    endif
    call addp(time,x,y,z,kx+1,ky+1,kz+1,pmass,pat,cat)
enddo
return
end
!------------------------------------------------------------------      
! placeuf
!------------------------------------------------------------------
subroutine placeuf(time,xm,ym,zm,sx,sy,sz,pat,cat,vel3,npart,pmass,iptype)
! place npart particles in a block starting at xyzm and ending at sxyz
! allocating then by flux over the vertical (specialized for perchlorate problem)
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
real:: vel3(3,0:nx,0:ny,0:nz),sx,sy,sz,xm,ym,zm,qtotal,q,r,randu01,sznew,zmnew,x,y,z,time
real, allocatable:: probw(:)
double precision pmass(nspec)
integer:: iw,nw,kxm,im,kym,jm,kzm,km,i,j,k,kxp,kyp,kzp,ip,kp,kx,ky,kz,npart,iptype
kxm=ifix(xm/dx); kym=ifix(ym/dy); kzm=ifix(zm/dz)
im=kxm+1; jm=kym+1; km=kzm+1
kxp=ifix((xm+sx)/dx); kyp=ifix((ym+sy)/dy); kzp=ifix((zm+sz)/dz)
kp=kzp+1
nw=kp-km+1
allocate(probw(nw))
probw(:)=0
! compute total flux leaving cell
iw=0
qtotal=0.
do k=km,km+nw-1
  kz=k-1
  iw=iw+1
  q=-bke*((vel3(1,kxm,jm,k)-vel3(1,im,jm,k))*ax+&
          (vel3(2,im,kym,k)-vel3(2,im,jm,k))*ay+&
          (vel3(3,im,jm,kz)-vel3(3,im,jm,k))*az)
  probw(iw)=0.
  if(q.gt.0.)probw(iw)=q
  if(q.gt.0.)qtotal=qtotal+q
enddo
if(qtotal.ne.0)probw(:)=probw(:)/qtotal ! normalize
! make into cumulative histogram
do iw=2,nw
  probw(iw)=probw(iw-1)+probw(iw)
enddo
do ip=1,npart
! determine which interval to release particle
    r=randu01()
    do iw=1,nw
      if(r.le.probw(iw))exit
    enddo
    k=km+iw-1
    sznew=dz
    zmnew=(k-1)*dz
    if(np+1.gt.maxnp)then
      write(iout,*)' error:too many particles, increase maxnp'
      stop    ' error:too many particles, increase maxnp'
    endif
    x=xm+sx*randu01(); y=ym+sy*randu01(); z=zmnew+sznew*randu01()
    kx=ifix(x/dx); ky=ifix(y/dy); kz=ifix(z/dz)
    call addp(time,x,y,z,kx+1,ky+1,kz+1,pmass,pat,cat)
enddo
return
end



!---------------------------------------------------------------------
! output 
!---------------------------------------------------------------------
subroutine output(opc,outunit,outfname,sam,pat,cat,vel3,bounds,por,ret,decay)
use global
type (particle):: pat(1,maxnp)
type (sample):: sam(nsam)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
real por(nzone),ret(nzone),decay(nzone)
real vel3(3,0:nx,0:ny,0:nz)
integer opc(nopc),opcnow(nopc),outunit(nopc+1)
character (len=80) outfname(nopc)
double precision masstotal(nspec)
! opc(1)      -MOMENTS                  
! opc(2)      -CONCENTRATIONS           
! opc(3)      -BREAKTHROUGH LOCATIONS  
! opc(4)      -BREAKTHROUGH COUNTERS     
! opc(5)      -PARTICLE LOCATIONS
! opc(6)      -INTERNAL BOUNDARY BREAKTHROUGH COUNTERS
if(opc(1).ge.1)call plotm(pat,outunit(1),outfname(1),cat)
if(opc(2).ge.1)call plotc(cat,outunit(2),opc,por,ret)
if(opc(2).ge.3)call plotmt3d(outfname(2),cat,imt3d,opc,por,ret) ! mt3d output
if(opc(3).eq.1)call plotd(cat,outunit(3))
if(opc(4).ge.1)call plotb(outunit(4),opc)
if(opc(5).eq.1)call plotp(pat,outunit(5))
if(opc(6).ge.1)call ploti(bounds,outunit(6),opc)
if(opc(7).ge.1)call plots(sam,pat,cat,vel3,por,ret,outunit(7),outunit(8))
!.store current output control for output to summary file
do iopc=1,nopc
  opcnow(iopc)=opc(iopc)
enddo
!.update output control
call opcupdt(opc,outunit,outfname)
!.print summary statistics
nptotal=0
masstotal=0.0
do ibtype=1,nbtype+1
  nptotal=nptotal+npbc(ibtype)
  do iloop=1,nspec
    masstotal(iloop)=masstotal(iloop)+massbc(ibtype,iloop)
  enddo
enddo
do idir=1,3
  nptotal=nptotal+netxyzm(idir)
  nptotal=nptotal+netxyzp(idir)
  do iloop=1,nspec
    masstotal(iloop)=masstotal(iloop)+netmassm(idir,iloop)
    masstotal(iloop)=masstotal(iloop)+netmassp(idir,iloop)
  enddo
enddo
nptotal=nptotal+netsplit
! recompute total mass if decay is active
idecay=0
do izone=1,nzone
  if(decay(izone).ne.0.0)idecay=1
enddo
if(idecay.eq.1)then
  mass=0.0
  do ip=1,np
    mass=mass+pat(1,ip)%pmass
  enddo
endif
write(iout,1000)sngl(curtime),(npbc(ibtype),ibtype=1,nbtype+1),&
(netxyzm(idir),netxyzp(idir),idir=1,3),netsplit,nptotal,np  

do iloop=1,nspec
  write(iout,1002)iloop,(sngl(massbc(ibtype,iloop)),ibtype=1,nbtype+1),&
  (sngl(netmassm(idir,iloop)),sngl(netmassp(idir,iloop)),idir=1,3),&
  sngl(masstotal(iloop)),sngl(mass(iloop))
enddo

write(iout,1004)(opcnow(iopc),iopc=1,nopc),abs(sngl(tnextvel)),sngl(tnextbnd),&
sngl(tnextpnt),sngl(tnextopc),sngl(tnextsrc)

if(nran.ne.0)write(iout,1005)sngl(rm)/float(nran)
! do not report correlation coef.
!     &sngl(cor(1)/dfloat(nran/3)),sngl(cor(2)/dfloat(nran/3)),
!     &sngl(cor(3)/dfloat(nran/3)),sngl(cor(2)/dfloat(nran/3)),
!     &sngl(cor(4)/dfloat(nran/3)),sngl(cor(5)/dfloat(nran/3)),
!     &sngl(cor(3)/dfloat(nran/3)),sngl(cor(5)/dfloat(nran/3)),
!     &sngl(cor(6)/dfloat(nran/3))
 1000 format(68('-')/&
       '               C U R R E N T   S T A T U S'/&
       68('-')/&
       ' time                            ',20('.'),3x,e15.8//&

       ' particles entering bc type 1    ',20('.'),3x,i15/&
       ' particles entering bc type 2    ',20('.'),3x,i15/&
       ' particles entering bc type 3    ',20('.'),3x,i15/&
       ' particles entering bc type 4    ',20('.'),3x,i15/&
       ' particles entering bc type 5    ',20('.'),3x,i15/&
       ' particles entering bc type 6    ',20('.'),3x,i15/&
       ' particles entering bc type 7    ',20('.'),3x,i15/&
       ' particles entering bc type 8    ',20('.'),3x,i15/&
       ' particles entering bc type 9    ',20('.'),3x,i15/&
       ' particles entering bc type 10   ',20('.'),3x,i15/&
       ' particles entering bc type 11   ',20('.'),3x,i15/&
       ' particles entering point sources',20('.'),3x,i15/&
       ' particles entering -x face      ',20('.'),3x,i15/&
       ' particles entering +x face      ',20('.'),3x,i15/&
       ' particles entering -y face      ',20('.'),3x,i15/&
       ' particles entering +y face      ',20('.'),3x,i15/&
       ' particles entering -z face      ',20('.'),3x,i15/&
       ' particles entering +z face      ',20('.'),3x,i15/&
       ' particles added due to splitting',20('.'),3x,i15/&
       ' TOTAL                           ',20('.'),3x,i15/&
       ' particles on record             ',20('.'),3x,i15//)

 1002 format(' ** Species number ',i2,' **'/&
       ' mass entering bc type 1         ',20('.'),3x,e15.8/&
       ' mass entering bc type 2         ',20('.'),3x,e15.8/&
       ' mass entering bc type 3         ',20('.'),3x,e15.8/&
       ' mass entering bc type 4         ',20('.'),3x,e15.8/&
       ' mass entering bc type 5         ',20('.'),3x,e15.8/&
       ' mass entering bc type 6         ',20('.'),3x,e15.8/&
       ' mass entering bc type 7         ',20('.'),3x,e15.8/&
       ' mass entering bc type 8         ',20('.'),3x,e15.8/&
       ' mass entering bc type 9         ',20('.'),3x,e15.8/&
       ' mass entering bc type 10        ',20('.'),3x,e15.8/&
       ' mass entering bc type 11        ',20('.'),3x,e15.8/&
       ' mass entering point sources     ',20('.'),3x,e15.8/&
       ' mass entering -x face           ',20('.'),3x,e15.8/&
       ' mass entering +x face           ',20('.'),3x,e15.8/&
       ' mass entering -y face           ',20('.'),3x,e15.8/&
       ' mass entering +y face           ',20('.'),3x,e15.8/&
       ' mass entering -z face           ',20('.'),3x,e15.8/&
       ' mass entering +z face           ',20('.'),3x,e15.8/&
       ' TOTAL                           ',20('.'),3x,e15.8/&
       ' mass on record                  ',20('.'),3x,e15.8//)

 1004 format(' output: moments                 ',20('.'),3x,i15/&
       '         concentrations          ',20('.'),3x,i15/&
       '         macrodispersion         ',20('.'),3x,i15/&
       '         breakthrough            ',20('.'),3x,i15/&
       '         particle attributes     ',20('.'),3x,i15/&
       '         internal counters       ',20('.'),3x,i15//&
       '         monitoring              ',20('.'),3x,i15//&

       ' time of next velocity update    ',20('.'),3x,e15.8/&
       ' time of next boundary update    ',20('.'),3x,e15.8/&
       ' time of next point source input ',20('.'),3x,e15.8/&
       ' time of next output             ',20('.'),3x,e15.8/&
       ' time of next type 1 BC update   ',20('.'),3x,e15.8//)
 1005 format('         Random Number Statistics'/&
       '         AVG          ',5('.'),1x,f15.12/)
!     &       '         COV11,12,13  ',5('.'),3(1x,f15.12)/
!     &       '         COV21,22,23  ',5('.'),3(1x,f15.12)/
!     &       '         COV31,32,33  ',5('.'),3(1x,f15.12)/)
return
end
!---------------------------------------------------------------------
! ploti 
!---------------------------------------------------------------------
subroutine ploti(bounds,iunit,opc)
!.....Report breakthrough for all internal boundary groups of the model. 
use global
type (boundary):: bounds(1:maxbnd)
integer:: opc(nopc)
double precision:: massin(nspec),maw(nspec)
real:: fluxin,fluxout
if(nbounds.eq.0)return
!.....boundary type
itypeo=bounds(1)%bc_type  
!.....initial/previous boundary group number
numbndo=bounds(1)%group
!.....initialize counter for breathrough
nin=0
massin=0.0
mawin=0.0
! initialize counter for number of boundaries in a group
inumbnd=0
!.loop through all boundaries
do ibounds=1,nbounds
  inumbnd=inumbnd+1
! boundary type
  itype=bounds(ibounds)%bc_type
! cell location
  i=0;j=0;kbot=0;ktop=0
  if(itype.eq.1.or.(itype.ge.3.and.itype.le.7))then
      i=bounds(ibounds)%ijk(1); j=bounds(ibounds)%ijk(2); kbot=bounds(ibounds)%kbot; ktop=bounds(ibounds)%ktop
  endif
!.boundary group number
  numbnd=bounds(ibounds)%group
!.
  fluxin=bounds(ibounds)%fluxin;fluxout=bounds(ibounds)%fluxout

!.if we are considering a new boundary group number then summarize statistics for the previous group
  if(numbnd.ne.numbndo)then
!...if there is only one boundary cell in this group write its ijk location
    if(iabs(itypeo).ne.7)then
      fluxino=0.0
      fluxouto=0.0
    endif
    if(inumbnd-1.eq.1)then

      if(opc(6).eq.1)write(iunit,'(f15.4,9(",",i8))')sngl(curtime),itypeo,numbndo,nin,inumbnd-1,io,jo,kboto,ktopo
!      if(opc(6).eq.2)write(iunit,'(f15.4,2(",",i8),4(",",e16.10),5(",",i8))')&
!      sngl(curtime),itypeo,numbndo,sngl(massin),sngl(maw),fluxino,fluxouto,inumbnd-1,io,jo,kboto,ktopo
      if(opc(6).eq.2)write(iunit,*)&
      sngl(curtime),itypeo,numbndo,sngl(massin),sngl(maw),fluxino,fluxouto,inumbnd-1,io,jo,kboto,ktopo

    else
      if(opc(6).eq.1)write(iunit,'(f15.4,9(",",i8))')sngl(curtime),itypeo,numbndo,nin,inumbnd-1
! maw - mass transferred internally within a MAW boundary
! fluxin - total flux into the well
! fluxout - flux out of the well back into the aquifer
! DAB old      if(opc(6).eq.2)write(iunit,'(f15.4,2(",",i8),4(",",e16.10),4(",",i8))')&
!      sngl(curtime),itypeo,numbndo,sngl(massin),sngl(maw),fluxino,fluxouto,inumbnd-1
      if(opc(6).eq.2)write(iunit,*)&
      sngl(curtime),itypeo,numbndo,sngl(massin),sngl(maw),fluxino,fluxouto,inumbnd-1
    endif
!...initialize counter for the number of bounds
    inumbnd=1
!...initialize nin for new boundary group
    nin=bounds(ibounds)%np_remove
    massin=bounds(ibounds)%mass_remove
    maw=bounds(ibounds)%maw
  else
!...increment counter for breathrough
    nin=nin+bounds(ibounds)%np_remove
    massin=massin+bounds(ibounds)%mass_remove
    maw=maw+bounds(ibounds)%maw
  endif
!.zero counter
  bounds(ibounds)%np_remove=0
  bounds(ibounds)%mass_remove=0.0
  bounds(ibounds)%maw=0.0
  numbndo=numbnd
  fluxino=fluxin; fluxouto=fluxout
  io=i;jo=j;ktopo=ktop;kboto=kbot
  itypeo=itype         
enddo
!.write output for the last boundary
if(iabs(itypeo).ne.7)then
  fluxino=0.0
  fluxouto=0.0
endif
if(inumbnd.eq.1)then
  if(opc(6).eq.1)write(iunit,'(f15.4,9(",",i8))')sngl(curtime),itypeo,numbndo,nin,inumbnd,&
  io,jo,kboto,ktopo
  if(opc(6).eq.2)write(iunit,*)&
  sngl(curtime),itypeo,numbndo,sngl(massin),sngl(maw),fluxino,fluxouto,inumbnd,io,jo,kboto,ktopo
! DAB old  if(opc(6).eq.2)write(iunit,'(f15.4,2(",",i8),4(",",e16.10),5(",",i8))')&
!  sngl(curtime),itypeo,numbndo,sngl(massin),sngl(maw),fluxino,fluxouto,inumbnd,io,jo,kboto,ktopo
else
  if(opc(6).eq.1)write(iunit,'(f15.4,9(",",i8))')sngl(curtime),itypeo,numbndo,nin,inumbnd
  if(opc(6).eq.2)write(iunit,*)&
  sngl(curtime),itypeo,numbndo,sngl(massin),sngl(maw),fluxino,fluxouto,inumbnd
! DAB old  if(opc(6).eq.2)write(iunit,'(f15.4,2(",",i8),4(",",e16.10),4(",",i8))')&
!  sngl(curtime),itypeo,numbndo,sngl(massin),sngl(maw),fluxino,fluxouto,inumbnd
endif
return
end


!---------------------------------------------------------------------
! MT3D FORMATTED OUTPUT
!---------------------------------------------------------------------
subroutine plotmt3d(flname,cat,iunit,opc,por,ret)
use global
implicit none
character (len=16) text
character (len=80) flname
character(80)::blah
integer:: i,j,k,iflag,ntrans,iunit
type (cell)::     cat(nx,ny,nz)
integer:: opc(nopc)
real:: por(nzone),ret(nzone)
data iflag,ntrans/1,0/
save ntrans,iflag
! if flag is set in main .dat file, print ucn file
do iloop = 1,nspec    ! now need a species loop
  write(blah,'(i3.3)')iloop
if(iflag.eq.1)then
  text='CONCENTRATION'
  FLNAME=FLNAME(1:index(FLNAME,'.')-1)//trim(blah)//'._001.ucn'
!~   open(iunit,file=flname,form='binary')
  open(iunit,file=flname,form='unformatted',access='stream')
  iflag=0
endif
ntrans=ntrans+1
if(irevz.eq.1)then
  do k=nz,1,-1
    write(iunit) ntrans,kstp,kper,sngl(curtime),text,nx,ny,nz-k+1
    write(iunit) ((sngl(cat(i,j,k)%cmass(iloop)/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone))),i=1,nx),j=ny,1,-1)
  enddo
else
  do k=1,nz
    write(iunit) ntrans,kstp,kper,sngl(curtime),text,nx,ny,k
    write(iunit) ((sngl(cat(i,j,k)%cmass(iloop)/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone))),i=1,nx),j=1,ny)
  enddo
endif
enddo   ! species loop
return
end
!---------------------------------------------------------------------
! plotc
!---------------------------------------------------------------------
subroutine plotc(cat,iunit,opc,por,ret)
! plot particle density and concentration 
use global
type (cell)::     cat(nx,ny,nz)
integer opc(nopc)
real:: por(nzone),ret(nzone)
integer,allocatable:: ncell2(:)
integer,allocatable:: ncell1(:),ncellt(:)                    
real,allocatable:: cell2(:,:)  
integer:: i,j,k,ijk                  
data iflag/1/
!.cell number and np file
if(iflag.eq.1)then
! DAB formatted output
  write(iunit,*)nx,ny,nz,dx,dy,dz,nspec
  iflag=0
endif
ncl=0 
allocate (ncellt(1:nxyz), STAT = ierror)
do i=1,nx; do j=1,ny; do k=1,nz
   if(cat(i,j,k)%np_cell.ne.0)then
     ncl=ncl+1
     ijk=i+(j-1)*nx+(k-1)*nxy
     ncellt(ncl)=ijk
   endif
enddo; enddo; enddo
allocate (ncell1(1:ncl), STAT = ierror)
if(ierror.ne.0)stop ' error allocating memory in plotc'
if(opc(2).eq.1)allocate (ncell2(1:ncl), STAT = ierror)
if(opc(2).eq.2)allocate (cell2(1:ncl,1:nspec), STAT = ierror)
if(ierror.ne.0)stop ' error allocating memory in plotc'
n=1
do icl=1,ncl
  ijk=ncellt(icl)
  k=(ijk-1)/nxy+1
  j=(ijk-1-(k-1)*nxy)/nx+1
  i= ijk-(k-1)*nxy-(j-1)*nx  
  ncell1(n)=ijk  
  if(opc(2).eq.1)ncell2(n)=cat(i,j,k)%np_cell
  do iloop=1,nspec
  if(opc(2).eq.2)cell2(n,iloop)=cat(i,j,k)%cmass(iloop)/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone))
  enddo
  n=n+1
enddo
! DAB formatted output
write(iunit,*)sngl(curtime),ncl
write(iunit,*)ncell1
if(opc(2).eq.1)write(iunit,*)ncell2
do iloop=1,nspec
write(iunit,*)iloop
if(opc(2).eq.2)write(iunit,*)cell2(:,iloop)
enddo

! clean house
deallocate (ncellt)
deallocate (ncell1)
if(opc(2).eq.1)deallocate (ncell2)
if(opc(2).eq.2)deallocate (cell2)
return
end
!---------------------------------------------------------------------
! plotb (after SLIM written by AFB Tompson, LLNL)
!---------------------------------------------------------------------
subroutine plotb(iunit,opc)
!.breakthrough counters
use global
integer::opc(nopc),netxyzm_last(3),netxyzp_last(3)
double precision :: netmassm_last(3,nspec),netmassp_last(3,nspec)
double precision :: mass_last(nspec)
data mass_last /nspec*0.0/
data netxyzm_last,netxyzp_last,np_last/0,0,0,0,0,0,0/
data netmassm_last,netmassp_last/nspec*0.,nspec*0.,nspec*0.,nspec*0.,nspec*0.,nspec*0./
save
!
if(opc(4).eq.1)write(iunit,10)sngl(curtime),np,np-np_last,&
(netxyzp(idir)-netxyzp_last(idir),&
 netxyzm(idir)-netxyzm_last(idir),idir=1,3)   !netxp,netxm,netyp,netym,netzp,netzm

if(opc(4).eq.2)write(iunit,11)sngl(curtime),mass,mass-mass_last, &   !netmassxp,netmassxm,netmassyp,
(netmassp(idir,:)-netmassp_last(idir,:),&
netmassm(idir,:)-netmassm_last(idir,:),idir=1,3) !netxp,netxm,netyp,netym,netzp,netzm np_last=np mass_last=mass
np_last=np
mass_last=mass
do idir=1,3
 netxyzm_last(idir)=netxyzm(idir)
 netxyzp_last(idir)=netxyzp(idir)
 netmassm_last(idir,:)=netmassm(idir,:)
 netmassp_last(idir,:)=netmassp(idir,:)
enddo
10 format(1x,f15.4,8(1x,i8))
11 format(1x,f15.4,50(1x,e11.5))
return
end
!---------------------------------------------------------------------
! plotd
!---------------------------------------------------------------------
subroutine plotd(cat,iunit)
! breakthrough for particles removed from internal cells 
use global
type (cell)::     cat(nx,ny,nz)
nout=0
do l4=1,nxyz
    kz=(l4-1)/nxy
    ky=(l4-1-kz*nxy)/nx
    kx= l4-1-kz*nxy-ky*nx
    i=kx+1; j=ky+1; k=kz+1
    if(cat(i,j,k)%np_remove.ne.0)nout=nout+1
enddo
write(iunit,'(1x,f15.4,8(1x,i8))')sngl(curtime),nout
if(nout.gt.0)then
  do l4=1,nxyz
    kz=(l4-1)/nxy
    ky=(l4-1-kz*nxy)/nx
    kx= l4-1-kz*nxy-ky*nx
    i=kx+1; j=ky+1; k=kz+1
    if(cat(i,j,k)%np_remove.ne.0)write(iunit,'(3(i10,1x),i10)')i,j,k,cat(i,j,k)%np_remove
    cat(i,j,k)%np_remove=0
  enddo
endif
return
end
!---------------------------------------------------------------------
! plotp (after SLIM written by AFB Tompson, LLNL)
!---------------------------------------------------------------------
subroutine plotp(pat,iunit)
!     particle attributes
use global
type (particle):: pat(1,maxnp)
!
write(iunit) np,sngl(curtime)
do ip=1,np	
  write(iunit)&
  pat(1,ip)%pnumber,&  
!  pat(1,ip)%birth_place(1),pat(1,ip)%birth_place(2),pat(1,ip)%birth_place(3),&
  pat(1,ip)%xyz(1),pat(1,ip)%xyz(2),pat(1,ip)%xyz(3),&
  sngl(curtime)-pat(1,ip)%birth_day,sngl(pat(1,ip)%pmass)
enddo
return
! skip formatted output

write(iunit,*) np,sngl(curtime)
do ip=1,np	
  write(iunit,'(1x,i15,1x,4(E15.8,1x),50D15.8)')&
  pat(1,ip)%pnumber,&  
!  pat(1,ip)%birth_place(1),pat(1,ip)%birth_place(2),pat(1,ip)%birth_place(3),&
  pat(1,ip)%xyz(1),pat(1,ip)%xyz(2),pat(1,ip)%xyz(3),&
  curtime-pat(1,ip)%birth_day,pat(1,ip)%pmass
enddo
return
end
!---------------------------------------------------------------------
! plots - output concentrations at monitoring locations
!---------------------------------------------------------------------
subroutine plots(sam,pat,cat,vel3,por,ret,iunit1,iunit2)
use global
implicit none
type (particle):: pat(1,maxnp)
type (sample):: sam(nsam)
type (cell)::   cat(nx,ny,nz)
real:: por(nzone),ret(nzone)
real vel3(3,0:nx,0:ny,0:nz),vmag,vmagtotal

integer:: isam,itime,ip,iunit1,iunit2,i,j,k,kx,ky,kz,izone
integer:: ibeg,iend,jbeg,jend,kbeg,kend
real:: zdiff,thick,dxy,cmax(nspec),msamplet1(nspec),msamplet2(nspec),&
       msamplet3(nspec),msamplet4(nspec),msamplet5(nspec),mtotal
real, allocatable:: msample(:,:)
data itime/0/
save itime
!
if(itime.eq.0)then
  write(iunit1,1000)"TIME,MTOTAL,MSAMPLE",(",LOCATION ",isam,isam=1,nsam)  
  write(iunit2,1000)"TIME,MTOTAL,MSAMPLE",(",LOCATION ",isam,isam=1,nsam)  
  do isam=1,nsam
    sam(isam)%conc=0.0
    sam(isam)%mass=0.0
  enddo
endif
itime=itime+1
if(.not.allocated(msample))allocate (msample(nsam,nspec))
msamplet1=0.0
msamplet3=0.0
do isam=1,nsam
  sam(isam)%conc=0.0
  sam(isam)%mass=0.0
  msample(isam,:)=0.0
  if(sam(isam)%itype.eq.1)then 
    do ip=1,np
      dxy=sqrt((pat(1,ip)%xyz(1)-sam(isam)%xyz(1))**2+(pat(1,ip)%xyz(2)-sam(isam)%xyz(2))**2)
      zdiff=abs(pat(1,ip)%xyz(3)-(sam(isam)%xyz(3)+sam(isam)%xyz(4))/2.0)
      thick=(sam(isam)%xyz(4)-sam(isam)%xyz(3))/2.0
      if(dxy.le.sam(isam)%radius.and.zdiff.le.thick)then
        kx=ifix(pat(1,ip)%xyz(1)/dx); ky=ifix(pat(1,ip)%xyz(2)/dy); kz=ifix(pat(1,ip)%xyz(3)/dz)
        i=kx+1; j=ky+1; k=kz+1
        do izone=1,sam(isam)%nzone
          if(cat(i,j,k)%zone.eq.sam(isam)%zone(izone))then
  	        sam(isam)%conc=sam(isam)%conc+sngl(pat(1,ip)%pmass)/sam(isam)%vol
	        sam(isam)%mass=sam(isam)%mass+sngl(pat(1,ip)%pmass)
     do iloop=1,nspec
	        msample(isam,iloop)=msample(isam,iloop)+sngl(pat(1,ip)%pmass(iloop))
     enddo
	        msamplet1=msamplet1+sngl(pat(1,ip)%pmass)
            exit
          endif
        enddo
	  endif
    enddo
  elseif(sam(isam)%itype.eq.2)then 
    i=abs(sam(isam)%ijk(1))
    cmax=0.0
    mtotal=0.0
    vmagtotal=0.0
    do j=1,ny;do k=1,nz
      if(sam(isam)%nzone.ne.0)then
        do izone=1,sam(isam)%nzone
            if(cat(i,j,k)%zone.eq.abs(sam(isam)%zone(izone)))then
                if(sam(isam)%zone(izone).lt.0)then
                  cmax=cmax+cat(i,j,k)%cmass
                else
                  cmax=max(cmax,cat(i,j,k)%cmass/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone)))
                endif
     do iloop=1,nspec
	        msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
!	            msample(isam)=msample(isam)+cat(i,j,k)%cmass
	            msamplet2=msamplet2+cat(i,j,k)%cmass
                exit
            endif
        enddo
      else
        if(sam(isam)%ijk(1).gt.0)then ! no flux weighting
          cmax=cmax+cat(i,j,k)%cmass
     do iloop=1,nspec
	      msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
	      msamplet2=msamplet2+cat(i,j,k)%cmass
        else 
          vmag=sqrt(((vel3(1,i-1,j,k)+vel3(1,i,j,k))/2.0)**2+&
                    ((vel3(2,i,j-1,k)+vel3(2,i,j,k))/2.0)**2+&
                    ((vel3(3,i,j,k-1)+vel3(3,i,j,k))/2.0)**2)
          cmax=cmax+cat(i,j,k)%cmass*vmag
          vmagtotal=vmag+vmagtotal
     do iloop=1,nspec
	        msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
!	      msample(isam)=msample(isam)+cat(i,j,k)%cmass
	      msamplet2=msamplet2+cat(i,j,k)%cmass          
        endif
      endif
    enddo;enddo
	if(sam(isam)%ijk(1).gt.0)then
      sam(isam)%conc=cmax    
    else
      if(vmagtotal.ne.0)then
        sam(isam)%conc=cmax/vmagtotal    
      else
        sam(isam)%conc=0.0    
      endif
    endif
  elseif(sam(isam)%itype.eq.3)then 
    j=abs(sam(isam)%ijk(2))
    cmax=0.0
    mtotal=0.0
    vmagtotal=0
    do i=1,nx;do k=1,nz
      if(sam(isam)%nzone.ne.0)then
        do izone=1,sam(isam)%nzone
            if(cat(i,j,k)%zone.eq.abs(sam(isam)%zone(izone)))then
              if(sam(isam)%zone(izone).lt.0)then
                cmax=cmax+cat(i,j,k)%cmass
              else
                cmax=max(cmax,cat(i,j,k)%cmass/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone)))
              endif
     do iloop=1,nspec
	        msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
!    	      msample(isam)=msample(isam)+cat(i,j,k)%cmass
  	          msamplet3=msamplet3+cat(i,j,k)%cmass
              exit
            endif
        enddo
      else
        if(sam(isam)%ijk(2).gt.0)then ! no flux weighting
          cmax=cmax+cat(i,j,k)%cmass
     do iloop=1,nspec
	        msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
!	      msample(isam)=msample(isam)+cat(i,j,k)%cmass
	      msamplet2=msamplet2+cat(i,j,k)%cmass
        else 
          vmag=sqrt(((vel3(1,i-1,j,k)+vel3(1,i,j,k))/2.0)**2+&
                    ((vel3(2,i,j-1,k)+vel3(2,i,j,k))/2.0)**2+&
                    ((vel3(3,i,j,k-1)+vel3(3,i,j,k))/2.0)**2)
          cmax=cmax+cat(i,j,k)%cmass*vmag
          vmagtotal=vmag+vmagtotal
     do iloop=1,nspec
	        msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
!	      msample(isam)=msample(isam)+cat(i,j,k)%cmass
	      msamplet2=msamplet2+cat(i,j,k)%cmass          
        endif
      endif
    enddo;enddo
	if(sam(isam)%ijk(2).gt.0)then
      sam(isam)%conc=cmax    
    else
      if(vmagtotal.ne.0)then
        sam(isam)%conc=cmax/vmagtotal    
      else
        sam(isam)%conc=0.0    
      endif
    endif
  elseif(sam(isam)%itype.eq.4)then 
    k=abs(sam(isam)%ijk(3))
    cmax=0.0
    mtotal=0.0
    vmagtotal=0
    do i=1,nx;do j=1,ny
      if(sam(isam)%nzone.ne.0)then
        do izone=1,sam(isam)%nzone
            if(cat(i,j,k)%zone.eq.abs(sam(isam)%zone(izone)))then
              if(sam(isam)%zone(izone).lt.0)then
                cmax=cmax+cat(i,j,k)%cmass
              else
                cmax=max(cmax,cat(i,j,k)%cmass/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone)))
              endif
     do iloop=1,nspec
      	      msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
	          msamplet4=msamplet4+cat(i,j,k)%cmass
              exit
            endif
        enddo
      else
        if(sam(isam)%ijk(3).gt.0)then ! no flux weighting
          cmax=cmax+cat(i,j,k)%cmass
     do iloop=1,nspec
      	      msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
!	      msample(isam)=msample(isam)+cat(i,j,k)%cmass
	      msamplet2=msamplet2+cat(i,j,k)%cmass
        else 
          vmag=sqrt(((vel3(1,i-1,j,k)+vel3(1,i,j,k))/2.0)**2+&
                    ((vel3(2,i,j-1,k)+vel3(2,i,j,k))/2.0)**2+&
                    ((vel3(3,i,j,k-1)+vel3(3,i,j,k))/2.0)**2)
          cmax=cmax+cat(i,j,k)%cmass*vmag
          vmagtotal=vmag+vmagtotal
     do iloop=1,nspec
      	      msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
!	      msample(isam)=msample(isam)+cat(i,j,k)%cmass
	      msamplet2=msamplet2+cat(i,j,k)%cmass          
        endif
      endif
    enddo;enddo
	if(sam(isam)%ijk(3).gt.0)then
      sam(isam)%conc=cmax    
    else
      if(vmagtotal.ne.0)then
        sam(isam)%conc=cmax/vmagtotal    
      else
        sam(isam)%conc=0.0    
      endif
    endif
  elseif(sam(isam)%itype.eq.5)then 
    do ip=1,np
      kx=ifix(pat(1,ip)%xyz(1)/dx); ky=ifix(pat(1,ip)%xyz(2)/dy); kz=ifix(pat(1,ip)%xyz(3)/dz)
      i=kx+1; j=ky+1; k=kz+1
      if(i.eq.sam(isam)%ijk(1).and.j.eq.sam(isam)%ijk(2).and.&
        (k.ge.sam(isam)%ijk(3).and.k.le.sam(isam)%ijk(4)))then
        if(sam(isam)%nzone.ne.0)then
           do izone=1,sam(isam)%nzone
             if(cat(i,j,k)%zone.eq.sam(isam)%zone(izone))then
  	           sam(isam)%conc=sam(isam)%conc+sngl(pat(1,ip)%pmass)/sam(isam)%vol
	           sam(isam)%mass=sam(isam)%mass+sngl(pat(1,ip)%pmass)
     do iloop=1,nspec
	           msample(isam,iloop)=msample(isam,iloop)+sngl(pat(1,ip)%pmass(iloop))
     enddo
	           msamplet5=msamplet5+sngl(pat(1,ip)%pmass)
               exit
             endif
           enddo
        else
  	       sam(isam)%conc=sam(isam)%conc+sngl(pat(1,ip)%pmass)/sam(isam)%vol
           sam(isam)%mass=sam(isam)%mass+sngl(pat(1,ip)%pmass)
     do iloop=1,nspec
	       msample(isam,iloop)=msample(isam,iloop)+sngl(pat(1,ip)%pmass(iloop))
     enddo
	       msamplet5=msamplet5+sngl(pat(1,ip)%pmass)
        endif
	  endif
    enddo
  elseif(sam(isam)%itype.eq.6)then 
    i=sam(isam)%ijk(1)
    cmax=0.0
    mtotal=0.0
    if(i.lt.0)then
      ibeg=1;iend=-i
    else
      ibeg=i;iend=nx
    endif
    do i=ibeg,iend;do j=1,ny;do k=1,nz
      if(sam(isam)%nzone.ne.0)then
        do izone=1,sam(isam)%nzone
            if(cat(i,j,k)%zone.eq.abs(sam(isam)%zone(izone)))then
              if(sam(isam)%zone(izone).lt.0)then
                cmax=cmax+cat(i,j,k)%cmass
              else
                cmax=max(cmax,cat(i,j,k)%cmass/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone)))
              endif
     do iloop=1,nspec
	          msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
	          msamplet2=msamplet2+cat(i,j,k)%cmass
              exit
            endif
        enddo
      else
        cmax=cmax+cat(i,j,k)%cmass
     do iloop=1,nspec
	    msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
	    msamplet2=msamplet2+cat(i,j,k)%cmass
      endif
    enddo;enddo;enddo
	sam(isam)%conc=cmax    
  elseif(sam(isam)%itype.eq.7)then 
    j=sam(isam)%ijk(2)
    cmax=0.0
    mtotal=0.0
    if(j.lt.0)then
      jbeg=1;jend=-j
    else
      jbeg=j;jend=ny
    endif
    do i=1,nx;do j=jbeg,jend;do k=1,nz
      if(sam(isam)%nzone.ne.0)then
        do izone=1,sam(isam)%nzone
            if(cat(i,j,k)%zone.eq.abs(sam(isam)%zone(izone)))then
              if(sam(isam)%zone(izone).lt.0)then
                cmax=cmax+cat(i,j,k)%cmass
              else
                cmax=max(cmax,cat(i,j,k)%cmass/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone)))
              endif
     do iloop=1,nspec
    	      msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
  	          msamplet3=msamplet3+cat(i,j,k)%cmass
              exit
            endif
        enddo
      else
        cmax=cmax+cat(i,j,k)%cmass
     do iloop=1,nspec
	    msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
	    msamplet2=msamplet2+cat(i,j,k)%cmass
      endif
    enddo;enddo;enddo
	sam(isam)%conc=cmax    
  elseif(sam(isam)%itype.eq.8)then 
    k=sam(isam)%ijk(3)
    cmax=0.0
    mtotal=0.0
    if(k.lt.0)then
      kbeg=1;kend=-k
    else
      kbeg=k;kend=nz
    endif
    do i=1,nx;do j=1,ny;do k=kbeg,kend
      if(sam(isam)%nzone.ne.0)then
        do izone=1,sam(isam)%nzone
            if(cat(i,j,k)%zone.eq.abs(sam(isam)%zone(izone)))then
              if(sam(isam)%zone(izone).lt.0)then
                cmax=cmax+cat(i,j,k)%cmass
              else
                cmax=max(cmax,cat(i,j,k)%cmass/dble(vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone)))
              endif
     do iloop=1,nspec
      	      msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
	          msamplet4=msamplet4+cat(i,j,k)%cmass
              exit
            endif
        enddo
      else
        cmax=cmax+cat(i,j,k)%cmass
     do iloop=1,nspec
	    msample(isam,iloop)=msample(isam,iloop)+cat(i,j,k)%cmass(iloop)
     enddo
	    msamplet2=msamplet2+cat(i,j,k)%cmass
      endif
    enddo;enddo;enddo
	sam(isam)%conc=cmax    
  endif
enddo
write(iunit1,'(g15.10,1000('','',g15.10))')curtime,mass,msamplet1,(sam(isam)%conc,isam=1,nsam)
write(iunit2,'(g15.10,1000('','',g15.10))')curtime,mass,msamplet1,(msample(isam,:),isam=1,nsam)
1000 format(a,1000(a,i5))
return
end
!---------------------------------------------------------------------
! plotm (after SLIM written by AFB Tompson, LLNL)
!---------------------------------------------------------------------
subroutine plotm(pat,iunit,fname,cat)
! 
use global
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
character (len=80) fname,fname2
double precision xmean,ymean,zmean,xd,yd,zd,pmass,totmas
double precision sigxx,sigxy,sigxz,sigyy,sigyz,sigzz,th4d
data iflag,iflag2/1,1/
! 
totmas=0.
xmean=0.; ymean=0.; zmean=0.
sigxx=0.; sigxy=0.; sigxz=0.; sigyy=0.; sigyz=0.; sigzz=0.
if(iflag.eq.1)then
  iflag=0
  write(iunit,*)' moments'
  write(iunit,*)10
  write(iunit,*)' time'
  write(iunit,*)' species #'
  write(iunit,*)' first spatial moment x-direction'
  write(iunit,*)' first spatial moment y-direction'
  write(iunit,*)' first spatial moment z-direction'
  write(iunit,*)' covariance xx'
  write(iunit,*)' covariance xy'
  write(iunit,*)' covariance xz'
  write(iunit,*)' covariance yy'
  write(iunit,*)' covariance yz'
  write(iunit,*)' covariance zz'
endif
! 
if(np.gt.0)then
  rmaxx=0.0; rmaxy=0.0; rmaxz=0.0     
  rminx=large; rminy=large; rminz=large  
  do iloop=1,nspec   

    do ip=1,np
      xd=pat(1,ip)%xyz(1); yd=pat(1,ip)%xyz(2); zd=pat(1,ip)%xyz(3)
      pmass=dble(pat(1,ip)%pmass(iloop))
!      pmass=dble(pat(1,ip)%pmass)
      rmaxx=max(rmaxx,sngl(xd))
      rmaxy=max(rmaxy,sngl(yd))
      rmaxz=max(rmaxz,sngl(zd))
      rminx=min(rminx,sngl(xd))
      rminy=min(rminy,sngl(yd))
      rminz=min(rminz,sngl(zd))
      totmas=totmas+pmass
      xmean=xmean+xd*pmass; ymean=ymean+yd*pmass; zmean=zmean+zd*pmass
      sigxx=sigxx+xd*xd*pmass; sigxy=sigxy+xd*yd*pmass; sigxz=sigxz+xd*zd*pmass
      sigyy=sigyy+yd*yd*pmass; sigyz=sigyz+yd*zd*pmass
      sigzz=sigzz+zd*zd*pmass
    enddo
    xmean=xmean/totmas; ymean=ymean/totmas; zmean=zmean/totmas
    sigxx=mx*(sigxx/totmas-xmean*xmean) 
    sigxy=mx*my*(sigxy/totmas-xmean*ymean)
    sigxz=mx*mz*(sigxz/totmas-xmean*zmean)
    sigyy=my*(sigyy/totmas-ymean*ymean)
    sigyz=my*mz*(sigyz/totmas-ymean*zmean)
    sigzz=mz*(sigzz/totmas-zmean*zmean)
    write(iunit,20) sngl(curtime),iloop,xmean,ymean,zmean,sigxx,sigxy,sigxz,sigyy,sigyz,sigzz
  enddo
endif

20 format(1x,f15.6,1x,i3,1x,16(f15.6,1x)) 
return
end

!-------------------------------------------------------------------
! skip
!-------------------------------------------------------------------
subroutine skip(iunit)
!.....skip a commented line of input
implicit none
integer:: iunit,iread_error
character (len=1) c
!
do
  read(iunit,'(a)',iostat=iread_error) c
  if(iread_error.eq.-1)return
  if(c.eq.'!' .or.c.eq.'#') then
    cycle
  else
    backspace(iunit)
    exit
  endif
enddo
return
end

!-------------------------------------------------------------------
! icell_correct
!-------------------------------------------------------------------
subroutine icell_correct(kx,ky,kz,i,j,k,x,y,z,xm,ym,zm,xp,yp,zp,vel3)
!.....correct cell location for odd case where particle lands on cell boundary
use global
implicit none
integer kx,ky,kz,i,j,k
real:: x,y,z,xm,ym,zm,xp,yp,zp,vxm,vym,vzm,vxp,vyp,vzp
real:: vel3(3,0:nx,0:ny,0:nz)
if(abs(x-xm).lt.smallxyz(1))then
if(bke*vel3(1,kx,j,k).lt.0.0.or.bke*vel3(1,kx,j,k).eq.0.0)then
  kx=kx-mx
  i=i-mx
  x=xm-mx*smallxyz(1)
  return
endif
endif
if(abs(y-ym).lt.smallxyz(2))then
if(bke*vel3(2,i,ky,k).lt.0.0.or.bke*vel3(2,i,ky,k).eq.0.0)then
  ky=ky-my
  j=j-my
  y=ym-my*smallxyz(2)
  return
endif
endif
if(abs(z-zm).lt.smallxyz(3))then
if(bke*vel3(3,i,j,kz).lt.0.0.or.bke*vel3(3,i,j,kz).eq.0.0)then
  kz=kz-mz
  k=k-mz
  z=zm-mz*smallxyz(3)
  return
endif
endif
if(abs(x-xp).lt.smallxyz(1))then
if(bke*vel3(1,i,j,k).gt.0.0.or.bke*vel3(1,i,j,k).eq.0.0)then
  kx=kx+mx
  i=i+mx
  x=xp+mx*smallxyz(1)
  return
endif
endif
if(abs(y-yp).lt.smallxyz(2))then
if(bke*vel3(2,i,j,k).gt.0.0.or.bke*vel3(2,i,j,k).eq.0.0)then
  ky=ky+my
  j=j+my
  y=yp+my*smallxyz(2)
  return
endif
endif
if(abs(z-zp).lt.smallxyz(3))then
if(bke*vel3(3,i,j,k).gt.0.0.or.bke*vel3(3,i,j,k).eq.0.0)then
  kz=kz+mz
  k=k+mz
  z=zp+mz*smallxyz(3)
  return
endif
endif
return
end

!-------------------------------------------------------------------
! random1
!-------------------------------------------------------------------
subroutine random(z1,z2,z3)
! reliable and fast uniform random number generator with
! mean zero and variance 1 (after SLIM written by AFB Tompson, LLNL)
use global
implicit none
integer:: ll,mm,nn
real:: z1,z2,z3
double precision:: zz1,zz2,zz3,tsq3 
data mm/1048576/,ll/1027/
data tsq3/3.46410161513775d+0/
nn=ll*iseed1
iseed1=mod(nn,mm)
z1=float(iseed1)/float(mm)
!      zz1=dble(iseed1)/dble(mm)
nn=ll*iseed1
iseed1=mod(nn,mm)
z2=float(iseed1)/float(mm)
!      zz2=dble(iseed1)/dble(mm)
nn=ll*iseed1
iseed1=mod(nn,mm)
z3=float(iseed1)/float(mm)
!      zz3=dble(iseed1)/dble(mm)
nran=nran+3
rm=rm+z1+z2+z3-1.5
z1=(z1-.5)*tsq3
z2=(z2-.5)*tsq3
z3=(z3-.5)*tsq3
!      rm=rm+zz1+zz2+zz3
!      zz1=(zz1-.5d+0)*tsq3
!      zz2=(zz2-.5d+0)*tsq3
!      zz3=(zz3-.5d+0)*tsq3
! don't waste time computing covariance
!      cor(1)=cor(1)+zz1*zz1
!      cor(2)=cor(2)+zz1*zz2
!      cor(3)=cor(3)+zz1*zz3
!      cor(4)=cor(4)+zz2*zz2
!      cor(5)=cor(5)+zz2*zz3
!      cor(6)=cor(6)+zz3*zz3
!      z1=sngl(zz1)
!      z2=sngl(zz2)
!      z3=sngl(zz3)
if(mod(nran,200000).eq.0)iseed1=iseed1+2
return
end


!-------------------------------------------------------------------
! randu01
!-------------------------------------------------------------------
real function randu01()
! reliable and fast uniform random number generator (between 0 and 1)
! (after SLIM written by AFB Tompson, LLNL)      
use global
implicit none
integer ll,mm,nn
data mm/1048576/,ll/1027/
!
nn=ll*iseed2
iseed2=mod(nn,mm)
randu01=float(iseed2)/float(mm)
return
end

!-------------------------------------------------------------------
! random2
!-------------------------------------------------------------------
subroutine random2(z1,z2,z3)
! often overly robust, but slow, uniform random number generator 
! with mean zero and varance 1
use global
implicit none
real:: ran3,z1,z2,z3
double precision zz1,zz2,zz3,tsq3 
data tsq3/3.46410161513775d+0/
!      zz1=dble(ran3(iseed1))
!      zz2=dble(ran3(iseed1))
!      zz3=dble(ran3(iseed1))
!      rm=rm+zz1+zz2+zz3
!      zz1=(zz1-.5d+0)*tsq3
!      zz2=(zz2-.5d+0)*tsq3
!      zz3=(zz3-.5d+0)*tsq3
! don't waste time computing covariance
!      cor(1)=cor(1)+zz1*zz1
!      cor(2)=cor(2)+zz1*zz2
!      cor(3)=cor(3)+zz1*zz3
!      cor(4)=cor(4)+zz2*zz2
!      cor(5)=cor(5)+zz2*zz3
!      cor(6)=cor(6)+zz3*zz3
!      z1=sngl(zz1)
!      z2=sngl(zz2)
!      z3=sngl(zz3)
! single precision
z1=ran3(iseed1)
z2=ran3(iseed1)
z3=ran3(iseed1)
nran=nran+3
rm=rm+z1+z2+z3-1.5
z1=(z1-.5)*tsq3
z2=(z2-.5)*tsq3
z3=(z3-.5)*tsq3
return
end

!-------------------------------------------------------------------
REAL FUNCTION ran3(idum)
implicit none
INTEGER idum
INTEGER MBIG,MSEED,MZ
REAL FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
if(idum.lt.0.or.iff.eq.0)then
  iff=1
  mj=MSEED-iabs(idum)
  mj=mod(mj,MBIG)
  ma(55)=mj
  mk=1
  do i=1,54
    ii=mod(21*i,55)
    ma(ii)=mk
    mk=mj-mk
    if(mk.lt.MZ)mk=mk+MBIG
    mj=ma(ii)
  end do
  do k=1,4
    do i=1,55
      ma(i)=ma(i)-ma(1+mod(i+30,55))
      if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
    enddo
  enddo
  inext=0
  inextp=31
  idum=1
endif
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
return
END


!------------------------------------------------------------
! addp
!------------------------------------------------------------
subroutine addp(time,x,y,z,i,j,k,pmass,pat,cat)
!.....add a particle at location x,y,x,in cell, with masses pmass
use global
implicit none
type (particle):: pat(1,maxnp)
!type (imparticle):: imp(1,maxnp)
type (cell)::     cat(nx,ny,nz)
integer:: i,j,k,ptype
real:: time,x,y,z
double precision:: pmass(nspec)
!

np=np+1
pnumber=pnumber+1
if(np.le.maxnp)then
!.initialize new particle
  pat(1,np)%xyz(1)=x; pat(1,np)%xyz(2)=y; pat(1,np)%xyz(3)=z
!  pat(1,np)%birth_place(1)=x; pat(1,np)%birth_place(2)=y; pat(1,np)%birth_place(3)=z
  pat(1,np)%birth_day=time
  pat(1,np)%pmass=pmass
  pat(1,np)%pnumber=pnumber
  pat(1,np)%ijk(1)=i; pat(1,np)%ijk(2)=j; pat(1,np)%ijk(3)=k
  pat(1,np)%active=.true.                  
  pat(1,np)%death_day=0.0  
!
  cat(i,j,k)%np_cell=cat(i,j,k)%np_cell+1       ! number of particles in cell
!
  cat(i,j,k)%cmass=cat(i,j,k)%cmass+pmass   ! mass in cell
! log the bithplace in outp and tag with a NEGATIVE ip
  write(ioutp)-pat(1,np)%pnumber  
  write(ioutp)x,y,z  
else
  write(iout,*)' maxnp exceded'
  write(*,*)' maxnp exceded'
  stop
endif
return
end
!------------------------------------------------------------
! addimp
!------------------------------------------------------------
subroutine addimp(time,x,y,z,i,j,k,pmass,imp,cat)
!.....add a particle at location x,y,x,in cell, with masses pmass
use global
implicit none
!type (particle):: pat(1,maxnp)
type (imparticle):: imp(1,maxnp)
type (cell)::     cat(nx,ny,nz)
integer:: i,j,k,ptype
real:: time,x,y,z
double precision:: pmass(inspec)
!

nimp=nimp+1
impnumber=impnumber+1
if(nimp.le.maxnp)then
!.initialize new IMMOBILE particle
  imp(1,nimp)%xyz(1)=x; imp(1,np)%xyz(2)=y; imp(1,np)%xyz(3)=z
!  pat(1,nimp)%birth_place(1)=x; pat(1,np)%birth_place(2)=y; pat(1,np)%birth_place(3)=z
  imp(1,nimp)%birth_day=time
  imp(1,nimp)%pmass=pmass
  imp(1,nimp)%pnumber=impnumber
  imp(1,nimp)%ijk(1)=i; imp(1,np)%ijk(2)=j; imp(1,np)%ijk(3)=k
  imp(1,nimp)%active=.true.                  
  imp(1,nimp)%death_day=0.0  
!
  cat(i,j,k)%nimp_cell=cat(i,j,k)%nimp_cell+1       ! number of particles in cell
!
  cat(i,j,k)%cimmass=cat(i,j,k)%cmass+pmass   ! mass in cell
! log the bithplace in outp and tag with a NEGATIVE ip
  write(ioutp)-imp(1,np)%pnumber  
  write(ioutp)x,y,z  
else
  write(iout,*)' maxnp exceded'
  write(*,*)' maxnp exceded'
  stop
endif
return
end
        
!------------------------------------------------------------
! updtp
!------------------------------------------------------------
subroutine updtp(io,jo,ko,i,j,k,ip,x,y,z,pat,cat)
!.....update particle and cell attributes
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
integer:: i,j,k,io,jo,ko,ip
real:: x,y,z
!
cat(io,jo,ko)%np_cell=cat(io,jo,ko)%np_cell-1
cat(i,j,k)%np_cell=cat(i,j,k)%np_cell+1
cat(io,jo,ko)%cmass=cat(io,jo,ko)%cmass-pat(1,ip)%pmass
cat(i,j,k)%cmass=cat(i,j,k)%cmass+pat(1,ip)%pmass
!
pat(1,ip)%xyz(1)=x; pat(1,ip)%xyz(2)=y; pat(1,ip)%xyz(3)=z
pat(1,ip)%ijk(1)=i; pat(1,ip)%ijk(2)=j; pat(1,ip)%ijk(3)=k
return
end
!------------------------------------------------------------
! conserve
!------------------------------------------------------------
subroutine conserve(pat,cat)
!.....check for cell-based particle and mass conservation 
use global
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
integer i,j,k
allocatable mcell(:,:,:,:),pcell(:,:,:)
integer pcell
double precision mcell
!
allocate (mcell(nx,ny,nz,nspec),pcell(nx,ny,nz))
mcell=0.
pcell=0
do ip=1,np
 i=pat(1,ip)%ijk(1); j=pat(1,ip)%ijk(2); k=pat(1,ip)%ijk(3)
 do iloop=1,nspec
   mcell(i,j,k,iloop)=mcell(i,j,k,iloop)+pat(1,ip)%pmass(iloop)
 enddo
 pcell(i,j,k)=pcell(i,j,k)+1
enddo
do i=1,nx; do j=1,ny; do k=1,nz
  do iloop=1,nspec
    if(abs(cat(i,j,k)%cmass(iloop)-mcell(i,j,k,iloop)).gt.small)then
      write(ibugout,*)' Non-zero Cell Mass-Balance Error '
      write(ibugout,*)' time,i,j,k,species,tracked mass,current mass',&
      sngl(curtime),i,j,k,cat(i,j,k)%cmass(iloop),mcell(i,j,k,iloop)
    endif
  enddo
  if(cat(i,j,k)%np_cell.ne.pcell(i,j,k))then
    write(ibugout,*)' Non-zero Cell Particle-Balance Error '
    write(ibugout,*)' time,i,j,k,tracked particle count,current count',&
    sngl(curtime),i,j,k,cat(i,j,k)%np_cell,pcell(i,j,k)
  endif
enddo; enddo; enddo
deallocate (mcell,pcell)
return
end
!------------------------------------------------------------
! removep (after SLIM written by AFB Tompson, LLNL)
!------------------------------------------------------------
subroutine removep(istart,pat,cat)
! remove tagged particles
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
integer:: i,j,k,ip,istart,i1,jj
double precision:: pmass   ! DAB unneeded
!.....remove those particles at end of array first
jj=0
do ip=np,istart,-1
  if(pat(1,ip)%active) then
    exit
  else
    jj=jj+1
!.........update particle and cell attributes
    i=pat(1,ip)%ijk(1)
    j=pat(1,ip)%ijk(2)
    k=pat(1,ip)%ijk(3)
    cat(i,j,k)%np_cell=cat(i,j,k)%np_cell-1
    cat(i,j,k)%cmass=cat(i,j,k)%cmass-pat(1,ip)%pmass
    mass=mass-pat(1,ip)%pmass
!.........log particle removal for output in plotd
    cat(i,j,k)%np_remove=cat(i,j,k)%np_remove+1
    cat(i,j,k)%mass_remove=cat(i,j,k)%mass_remove+pat(1,ip)%pmass
!.........log particle removal in output file
    write(ioutp)pat(1,ip)%pnumber !,'(1x,i15,1x,8(G15.5,1x),10G15.5,3(1x,i6))')&  
    write(ioutp)& !,'(1x,i15,1x,8(G15.5,1x),G15.5,3(1x,i6))')&
    pat(1,ip)%birth_day,pat(1,ip)%death_day,&  
!    pat(1,ip)%birth_place(1),pat(1,ip)%birth_place(2),pat(1,ip)%birth_place(3),&
    pat(1,ip)%xyz(1),pat(1,ip)%xyz(2),pat(1,ip)%xyz(3),&
    pat(1,ip)%pmass,i,j,k
!    write(*,*)& !,'(1x,i15,1x,8(G15.5,1x),G15.5,3(1x,i6))')&
!    pat(1,ip)%pnumber,pat(1,ip)%death_day,i,j,k
  endif
enddo
np=np-jj
jj=0
i1=np
!.....now remove the rest
do ip=np,istart,-1
  if(.not.pat(1,ip)%active)then
!....update particle and cell attributes
    i=pat(1,ip)%ijk(1)
    j=pat(1,ip)%ijk(2)
    k=pat(1,ip)%ijk(3)
    cat(i,j,k)%np_cell=cat(i,j,k)%np_cell-1
    cat(i,j,k)%cmass=cat(i,j,k)%cmass-pat(1,ip)%pmass
!.........log particle removal for output in plotd
    cat(i,j,k)%np_remove=cat(i,j,k)%np_remove+1
    cat(i,j,k)%mass_remove=cat(i,j,k)%mass_remove+pat(1,ip)%pmass
    mass=mass-pat(1,ip)%pmass
!.........log particle removal in output file
    write(ioutp)pat(1,ip)%pnumber !,'(1x,i15,1x,8(G15.5,1x),10G15.5,3(1x,i6))')&  
    write(ioutp) & !,'(1x,i15,1x,8(G15.5,1x),G15.5,3(1x,i6))')&
    pat(1,ip)%birth_day,pat(1,ip)%death_day,&  
!    pat(1,ip)%birth_place(1),pat(1,ip)%birth_place(2),pat(1,ip)%birth_place(3),&
    pat(1,ip)%xyz(1),pat(1,ip)%xyz(2),pat(1,ip)%xyz(3),&
    pat(1,ip)%pmass,i,j,k
!    write(*,*)& !,'(1x,i15,1x,8(G15.5,1x),G15.5,3(1x,i6))')&
!    pat(1,ip)%pnumber,pat(1,ip)%death_day,i,j,k
!
    pat(1,ip)=pat(1,i1)
!
    i1=i1-1
    jj=jj+1
  endif
enddo
np=np-jj
return
end

!------------------------------------------------------------
! reflect
!------------------------------------------------------------
subroutine reflects(timeleft,pat,cat,x,y,z,ip)
! Single precision
! reflect partcles at edges of domain
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
integer:: idir,ip
real:: x,y,z,timeleft
real:: xd,yd,zd 
real:: xyz(3),xyzd(3)
!double precision:: xd,yd,zd 
double precision:: pmass(nspec) !,xyz(3),xyzd(3)

!
!xyz(1)=dble(x)
!xyz(2)=dble(y)
!xyz(3)=dble(z)
xyz(1)=x
xyz(2)=y
xyz(3)=z
do idir=1,3
!  xyzd(idir)=dabs(dble(xyz(idir))-xyzc(idir))
  xyzd(idir)=abs(xyz(idir)-sngl(xyzc(idir)))
!.is.particle within shell of cells sourounding boundary, or outside of domain?
  if(xyzd(idir).ge.sngl(xyzl2(idir)))then
!...if it's on the domain, nudge it out 
    if(xyz(idir).eq.sngl(xyzmax(idir)))xyz(idir)=xyz(idir)+smallxyz(idir)
    if(xyz(idir).eq.sngl(xyzmin(idir)))xyz(idir)=xyz(idir)-smallxyz(idir)
!...reflect particles if required
    if(xyz(idir).ge.sngl(xyzmax(idir)))&
    xyz(idir)=xyz(idir)-2.*ixyzp(idir)*(xyz(idir)-sngl(xyzmax(idir)))
    if(xyz(idir).le.sngl(xyzmin(idir)))&
    xyz(idir)=xyz(idir)+2.*ixyzm(idir)*(sngl(xyzmin(idir))-xyz(idir))
    xyzd(idir)=abs(xyz(idir)-sngl(xyzc(idir)))
!...Now that we've reflected where required, if particle is still outside of domain remove it!
!...Is particle outside of active domain?
    if(xyzd(idir).ge.sngl(xyzl2(idir)))then
!.....check for out of bounds, update breakthrough counters, and tag partcle for removal
      pat(1,ip)%active=.false.    ! tag particle for removal
      pat(1,ip)%death_day=curtime-timeleft  
      pmass=pat(1,ip)%pmass 
      if(xyz(idir).ge.sngl(xyzmax(idir))) then
        netxyzp(idir)=netxyzp(idir)-1
  do iloop=1,nspec
        netmassp(idir,iloop)=netmassp(idir,iloop)-pmass(iloop)
  enddo
      else if(xyz(idir).le.sngl(xyzmin(idir))) then
        netxyzm(idir)=netxyzm(idir)-1
 do iloop=1,nspec
        netmassm(idir,iloop)=netmassm(idir,iloop)-pmass(iloop)
 enddo
      end if 
      return
    endif
  endif
enddo
x=xyz(1)
y=xyz(2)
z=xyz(3)
return
end
!------------------------------------------------------------
! reflect
!------------------------------------------------------------
subroutine reflect(timeleft,pat,cat,x,y,z,ip)
!Double precision
! reflect partcles at edges of domain
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
integer:: idir,ip
real:: x,y,z,timeleft
double precision:: xd,yd,zd 
double precision:: pmass(nspec),xyz(3),xyzd(3)
!
xyz(1)=dble(x)
xyz(2)=dble(y)
xyz(3)=dble(z)
do idir=1,3
  xyzd(idir)=dabs(dble(xyz(idir))-xyzc(idir))
!.is.particle within shell of cells sourounding boundary, or outside of domain?
  if(xyzd(idir).ge.xyzl2(idir))then
!...if it's on the domain, nudge it out 
    if(xyz(idir).eq.xyzmax(idir))xyz(idir)=xyz(idir)+smallxyz(idir)
    if(xyz(idir).eq.xyzmin(idir))xyz(idir)=xyz(idir)-smallxyz(idir)
!...reflect particles if required
    if(xyz(idir).ge.xyzmax(idir))xyz(idir)=xyz(idir)-2.*ixyzp(idir)*(xyz(idir)-xyzmax(idir))
    if(xyz(idir).le.xyzmin(idir))xyz(idir)=xyz(idir)+2.*ixyzm(idir)*(xyzmin(idir)-xyz(idir))
    xyzd(idir)=dabs(dble(xyz(idir))-xyzc(idir))
!...Now that we've reflected where required, if particle is still outside of domain remove it!
!...Is particle outside of active domain?
    if(xyzd(idir).ge.xyzl2(idir))then
!.....check for out of bounds, update breakthrough counters, and tag partcle for removal
      pat(1,ip)%active=.false.    ! tag particle for removal
      pat(1,ip)%death_day=curtime-timeleft  
      pmass=pat(1,ip)%pmass 
      if(xyz(idir).ge.xyzmax(idir)) then
        netxyzp(idir)=netxyzp(idir)-1
  do iloop=1,nspec
        netmassp(idir,iloop)=netmassp(idir,iloop)-pmass(iloop)
  enddo
!        netmassp(idir)=netmassp(idir)-pmass
      else if(xyz(idir).le.xyzmin(idir)) then
        netxyzm(idir)=netxyzm(idir)-1
 do iloop=1,nspec
        netmassm(idir,iloop)=netmassm(idir,iloop)-pmass(iloop)
 enddo
!        netmassm(idir)=netmassm(idir)-pmass
      end if 
      return
    endif
  endif
enddo
x=sngl(xyz(1))
y=sngl(xyz(2))
z=sngl(xyz(3))
return
end
!------------------------------------------------------------
! absorb
!------------------------------------------------------------
subroutine absorb(pat,cat,bounds,vel3,por,ret,kx,ky,kz,i,j,k,ip,dtmin,timeleft,&
x,y,z)
use global
implicit none
type (particle):: pat(1,maxnp)
type (cell)::     cat(nx,ny,nz)
type (boundary):: bounds(1:maxbnd)
! tag particles in absorbing boundaries for removal
integer:: ibounds,ip,itype,kx,ky,kz,i,j,k,kk,numbnd,ibou,istart_ibou,ktop,kbot,kkz,iside
real:: z1,q,vxm,vym,vzm,vxp,vyp,vzp,dtmin,timeleft,prob,randu01
real:: vel3(3,0:nx,0:ny,0:nz),por(nzone),ret(nzone)
real:: fluxin,fluxout,r,probc(6),fmassin,fmassout,qtotal,x,y,z
real:: xmbc,xm,ymbc,ym,zmbc,zm,dxbc,dybc,dzbc
!.....boundary condition #
ibounds=cat(i,j,k)%bc_number
!.....if no boundary, return
!      if(ibounds.eq.0)return
itype=bounds(ibounds)%bc_type
!.....constant conc. boundaries (updatad only at tnextsrc)
if((itype.eq.1).and.(curtime.eq.tnextsrc).and.(timeleft.eq.0.))then
  pat(1,ip)%active=.false.
  pat(1,ip)%death_day=curtime-timeleft  
!.......removal counters
  bounds(ibounds)%np_remove=bounds(ibounds)%np_remove-1
  bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove-pat(1,ip)%pmass
  npbc(itype)=npbc(itype)-1
  do iloop=1,nspec
     massbc(itype,iloop)=massbc(itype,iloop)-pat(1,ip)%pmass(iloop)
  enddo
!  massbc(itype)=massbc(itype)-pat(1,ip)%pmass
!.absorbing boundary
elseif(itype.eq.3)then
  pat(1,ip)%active=.false.
  pat(1,ip)%death_day=curtime-timeleft  
!.removal counters
  bounds(ibounds)%np_remove=bounds(ibounds)%np_remove-1
  bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove-pat(1,ip)%pmass
  npbc(itype)=npbc(itype)-1
  do iloop=1,nspec
     massbc(itype,iloop)=massbc(itype,iloop)-pat(1,ip)%pmass(iloop)
  enddo
!  massbc(itype)=massbc(itype)-pat(1,ip)%pmass
!.absorbing boundary with computed or specified flux
elseif(itype.eq.4.or.itype.eq.5)then
  if(itype.eq.4)then
    q=bounds(ibounds)%flux
  elseif(itype.eq.5)then
    q=bke*((vel3(1,kx,j,k)-vel3(1,i,j,k))*ax+&
           (vel3(2,i,ky,k)-vel3(2,i,j,k))*ay+&
           (vel3(3,i,j,kz)-vel3(3,i,j,k))*az)
!...no absorbing below epsilon
!    if(q.lt.small*vol*por(cat(i,j,k)%zone))q=0.0
  endif
!.probability of absorbing based on dtmin*flux/volume of voids in cell
  prob=dtmin*q/(dtmin*q+vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone))
  if(randu01().lt.prob)then
    pat(1,ip)%active=.false.
    pat(1,ip)%death_day=curtime-timeleft  
!...removal counters
    bounds(ibounds)%np_remove=bounds(ibounds)%np_remove-1
    bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove-pat(1,ip)%pmass
    npbc(itype)=npbc(itype)-1
  do iloop=1,nspec
     massbc(itype,iloop)=massbc(itype,iloop)-pat(1,ip)%pmass(iloop)
  enddo
!    massbc(itype)=massbc(itype)-pat(1,ip)%pmass
  endif
elseif(itype.eq.7)then
  q=bke*((vel3(1,kx,j,k)-vel3(1,i,j,k))*ax+&
         (vel3(2,i,ky,k)-vel3(2,i,j,k))*ay+&
         (vel3(3,i,j,kz)-vel3(3,i,j,k))*az)
!.probability of absorbing based on dtmin*flux/volume of voids in cell
  prob=dtmin*q/(dtmin*q+vol*por(cat(i,j,k)%zone)*ret(cat(i,j,k)%zone))
  if(randu01().lt.prob)then
    fluxout=bounds(ibounds)%fluxout
    if(fluxout.ge.0.0)then
!.....MAW boundary has no flow into aquifer - just remove particle 
      pat(1,ip)%active=.false.
      pat(1,ip)%death_day=curtime-timeleft  
!.....removal counters
      bounds(ibounds)%np_remove=bounds(ibounds)%np_remove-1
      bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove-pat(1,ip)%pmass
      npbc(itype)=npbc(itype)-1
  do iloop=1,nspec
     massbc(itype,iloop)=massbc(itype,iloop)-pat(1,ip)%pmass(iloop)
  enddo
!      massbc(itype)=massbc(itype)-pat(1,ip)%pmass
    else
!.....MAW boundary flowing into system - remove mass from particle and move particle in well
!
!.....compute fraction to remove
      fluxin=bounds(ibounds)%fluxin
!.....fluxout could be greater than fluxin due to error in solution
!.....assume mass completely mixed in inflow
      fmassout=min(abs(fluxout),fluxin)/fluxin ! fraction transferred
      fmassin=1.0-fmassout ! fraction absorbed
!.....account for mass leaving well - this is the mass that leaves the domain 
      bounds(ibounds)%mass_remove=bounds(ibounds)%mass_remove-dble(fmassin)*pat(1,ip)%pmass
  do iloop=1,nspec
     massbc(itype,iloop)=massbc(itype,iloop)-dble(fmassin)*pat(1,ip)%pmass(iloop)
  enddo
!      massbc(itype)=massbc(itype)-dble(fmassin)*pat(1,ip)%pmass
!.....change particle mass
      pat(1,ip)%pmass=pat(1,ip)%pmass*dble(fmassout) ! updated particle mass transferred to new location
!.....Find beginning of MAW boundary
      ibou=ibounds
      numbnd=bounds(ibounds)%group
      do 
        ibou=ibou-1
        if(ibou.lt.1)exit
        if(bounds(ibou)%group.ne.numbnd)exit
      enddo
      istart_ibou=ibou+1
!.....choose cell to transfer particle
      r=randu01()
      prob=0.0
      do ibou=istart_ibou,nbounds
        if(bounds(ibou)%group.ne.numbnd)exit              
        ktop=bounds(ibou)%ktop
        kbot=bounds(ibou)%kbot
        do kk=kbot,ktop
          kkz=kk-1
          q=((vel3(1,kx,j,kk)-vel3(1,i,j,kk))*ax+&
             (vel3(2,i,ky,kk)-vel3(2,i,j,kk))*ay+&
             (vel3(3,i,j,kkz)-vel3(3,i,j,kk))*az)
          if(q.lt.0.0)then
            prob=prob+q/fluxout
            if(r.le.prob)exit
          endif
        enddo
        if(r.le.prob)exit
      enddo
      if(ibou.gt.nbounds)ibou=nbounds      
      if(kk.gt.ktop)kk=ktop
!.....ibou = boundary number
!.....kk = layer in new cell
!.....choose location
!.....decide which cell edge to release particle from
      probc(:)=0.
      qtotal=0.
      iside=1
      if(bke*vel3(1,i-1,j,kk).lt.0.)then
        qtotal=qtotal+abs(bke*vel3(1,i-1,j,kk)*ax)
        probc(iside)=abs(bke*vel3(1,i-1,j,kk)*ax)
      endif
      iside=iside+1
      probc(iside)=probc(iside-1)
      if(bke*vel3(1,i,j,kk).gt.0.)then
        qtotal=qtotal+abs(bke*vel3(1,i,j,kk)*ax)
        probc(iside)=probc(iside)+abs(bke*vel3(1,i,j,kk)*ax)
      endif
      iside=iside+1
      probc(iside)=probc(iside-1)
      if(bke*vel3(2,i,j-1,kk).lt.0.)then
        qtotal=qtotal+abs(bke*vel3(2,i,j-1,kk)*ay)
        probc(iside)=probc(iside)+abs(bke*vel3(2,i,j-1,kk)*ay)
      endif
      iside=iside+1
      probc(iside)=probc(iside-1)
      if(bke*vel3(2,i,j,kk).gt.0.)then
        qtotal=qtotal+abs(bke*vel3(2,i,j,kk)*ay)
        probc(iside)=probc(iside)+abs(bke*vel3(2,i,j,kk)*ay)
      endif
      iside=iside+1
      probc(iside)=probc(iside-1)
      if(bke*vel3(3,i,j,kkz).lt.0.)then
        qtotal=qtotal+abs(bke*vel3(3,i,j,kkz)*az)
        probc(iside)=probc(iside)+abs(bke*vel3(3,i,j,kkz)*az)
      endif
      iside=iside+1
      probc(iside)=probc(iside-1)
      if(bke*vel3(3,i,j,kk).gt.0.)then
        qtotal=qtotal+abs(bke*vel3(3,i,j,kk)*az)
        probc(iside)=probc(iside)+abs(bke*vel3(3,i,j,kk)*az)
      endif
      probc=probc/qtotal
      r=randu01()
      do iside=1,6
        if(r.le.probc(iside))exit
      enddo
      xm=(i-1)*dx; ym=(j-1)*dy;zm=(kk-1)*dz
      if(iside.eq.1)then
        xmbc=xm; ymbc=ym+bcdxyz*dy; zmbc=zm+bcdxyz*dz
        dxbc=bcdxyz*dx; dybc=dy-2*bcdxyz*dy; dzbc=dz-2*bcdxyz*dz
      elseif(iside.eq.2)then
        xmbc=xm+dx-bcdxyz*dx; ymbc=ym+bcdxyz*dy; zmbc=zm+bcdxyz*dz
        dxbc=bcdxyz*dx; dybc=dy-2*bcdxyz*dy; dzbc=dz-2*bcdxyz*dz
      elseif(iside.eq.3)then
        xmbc=xm+bcdxyz*dx; ymbc=ym; zmbc=zm+bcdxyz*dz
        dxbc=dx-2*bcdxyz*dx; dybc=bcdxyz*dy; dzbc=dz-2*bcdxyz*dz
      elseif(iside.eq.4)then
        xmbc=xm+bcdxyz*dx; ymbc=ym+dy-bcdxyz*dy; zmbc=zm+bcdxyz*dz
        dxbc=dx-2*bcdxyz*dx; dybc=bcdxyz*dy; dzbc=dz-2*bcdxyz*dz
      elseif(iside.eq.5)then
        xmbc=xm+bcdxyz*dx; ymbc=ym+bcdxyz*dy; zmbc=zm
        dxbc=dx-2*bcdxyz*dx; dybc=dy-2*bcdxyz*dy; dzbc=bcdxyz*dz
      elseif(iside.eq.6)then
        xmbc=xm+bcdxyz*dx; ymbc=ym+bcdxyz*dy; zmbc=zm+dz-bcdxyz*dz
        dxbc=dx-2*bcdxyz*dx; dybc=dy-2*bcdxyz*dy; dzbc=bcdxyz*dz
      endif
!.....reposition particle
      kz=kkz
      k=kk
      x=xmbc+dxbc*randu01(); y=ymbc+dybc*randu01(); z=zmbc+dzbc*randu01()
!.....account for mass entering boundary
      bounds(ibou)%maw=bounds(ibou)%maw+pat(1,ip)%pmass
    endif
  endif
endif
return
end



!---------------------------------------------------------------------
! MT3D FORMATTED OUTPUT
!---------------------------------------------------------------------
subroutine cnf(FLNAME,icnf,nx,ny,nz,dx,dy,dz)
! MT3D CNF FILE
CHARACTER (len=80) FLNAME
CHARACTER (len=80) FLNAMECNF
CINACT=-999.0
CDRY=CINACT
FLNAMECNF=FLNAME(1:index(FLNAME,'.'))//'cnf'
! if flag is set in main .dat file, print cnf file
OPEN(ICNF,FILE=FLNAMECNF)
WRITE(ICNF,*) nz,ny,nx !NLAY,NROW,NCOL
WRITE(ICNF,*) (dx,i=1,nx) !(DELR(J),J=1,NCOL)
WRITE(ICNF,*) (dy,j=1,ny)
WRITE(ICNF,*) ((0.0,i=1,nx),j=1,ny)
WRITE(ICNF,*) (((dz,i=1,nx),j=1,ny),k=1,nz)
WRITE(ICNF,*) CINACT,CDRY
CLOSE(ICNF)
RETURN
END
