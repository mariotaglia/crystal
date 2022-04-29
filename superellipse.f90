subroutine update_matrix_superellipse(flag)
use system
use channel
use ematrix
use MPI
use const
use chainsdat
use molecules
use rotchain
implicit none

real*8 rchannel2, rchannelL2, rchannelS2
real*8, external :: rands
integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
logical flag
integer j,ix,iy,iz
real pnumber
real*8 area
real*8 sumpolseg 
real*8 sstemp,vvtemp, maxss
real*8 cutarea
real*8 temp
real*8 temp2
real*8 sumvoleps1, sumvolprot1, sumvolq1, sumvolx1
integer ncha1
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
integer i
real*8 volxx1(dimx,dimy,dimz)
real*8 volxx(dimx,dimy,dimz)
integer nbands


cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
sumpolseg = 0.0

rchannel2 = rchannel**2
rchannelL2 = (rchannel + 3*delta)**2
rchannelS2 = (rchannel - delta)**2

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

! channel center in x, y plane

 originc(1) = float(dimx)*delta/2.0 
 originc(2) = float(dimy)*delta/2.0 

 npoints = 50

 flag = .false.


 call integrate_superellipse(rchannelL2,RdimZ, originc ,npoints, voleps1 , sumvoleps1, flag)

 flag = .false. ! not a problem if eps lays outside boundaries

 call integrate_superellipse(rchannel2, RdimZ, originc,npoints, volprot1, sumvolprot1, flag)

 call integrate_superellipse(rchannelS2,RdimZ, originc,npoints, volq1, sumvolq1, flag)

 call newintegrateg_superellipse(rchannel2,RdimZ,originc,npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1)

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eepsc

! epstype

select case (epstype)

case (1)
nbands = dimz/8
do iz = 1, dimz
if (mod(int((iz-1)/nbands),2).eq.1) then
voleps1(:,:,iz) = 0.0
endif
enddo

endselect


!! charge
 volq1 = volprot1-volq1
 temp = sumvolprot1-sumvolq1
 volq1 = volq1/temp*echargec/(delta**3) ! sum(volq) is echarge

!! grafting

 area = 2.0*pi*rchannel*float(dimz)*delta

 temp2 = maxval(volx1)

 where(volx1<temp2*cutarea) ! remove cells with very little area
 volx1 = 0.0
 end where 

 do i = 1, ncha1
 volx1(i) = volx1(i)/sumvolx1*area*(sigmac+sigmar*(rands(seed)-0.5))
 volxx1(p1(i,1),p1(i,2),p1(i,3)) = & 
     volxx1(p1(i,1),p1(i,2),p1(i,3))/sumvolx1*area*(sigmac+sigmar*(rands(seed)-0.5))
 enddo

 maxss = 1.0d100
 sumpolseg = sumpolseg + area*sigmac*long

!! volume  
 volprot1 = volprot1 * 0.9999
 volprot = volprot+volprot1

! CHECK COLLISION HERE...
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true. 
 endif
 
 voleps = voleps + voleps1
 volq = volq + volq1 

! add com1 and volx to list

 volxx = volxx + volxx1

ncha = ncha1
do i = 1, ncha
volx(i)=volx1(i)
com(i,:)=com1(i,:)
p0(i,:)=p1(i,:)
rotangle(i) = atan2(com1(i,1)-originc(1), com1(i,2)-originc(2))
enddo

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

if (verbose.ge.2) then
temp = pi*rchannel2*float(dimz)*delta

if (rank.eq.0) then
write(stdout,*) 'superellipse:', 'update_matrix: Total nanosuperellipse volumen real space= ', temp
write(stdout,*) 'superellipse:', 'update_matrix: Total discretized volumen =', (sum(volprot))*delta**3
write(stdout,*) 'superellipse:', 'number of monomers in system =', sumpolseg 
endif
endif

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)

title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)

end subroutine

subroutine newintegrateg_superellipse(rchannel2,RdimZ,originc, npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)
use system
use transform
use chainsdat
use ematrix
use const
implicit none
real*8 sumvolx1
integer npoints
integer indexvolx(dimx,dimy,dimz)
integer listvolx(ncha,3)
real*8 radio
real*8 rchannel, rchannel2, originc(2)
real*8 phi, dphi, tetha,dtetha, as, ds
integer mphi, mtetha
integer ix,iy,iz,jx,jy,jz
real*8 x(3), v(3)
integer i,j
integer ncount
real*8 comshift ! how far from the surface of the sphere the grafting point is
integer ncha1 ! count for current sphere
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
real*8 volxx1(dimx,dimy,dimz)
integer flagin
integer dims(3), is(3), js(3)
integer jjjz, jjjt, npointz, npointt
integer RdimZ

pi=acos(-1.0)

dims(1) = dimx
dims(2) = dimy
dims(3) = dimz

rchannel = sqrt(rchannel2)

indexvolx = 0
ncha1 = 0
volx1 = 0.0
sumvolx1 = 0.0 ! total volumen, including that outside system
com1 = 0.0
p1 = 0
volxx1 = 0.0

! This routine determines the surface coverage and grafting positions only for cylinder
!

npointz = npoints*dimz
npointt = int(2.0*pi*rchannel/delta)*npoints

do jjjz = 1, npointz-1
do jjjt = 1, npointt

x(1) = cos(float(jjjt)/float(npointt)*2.0*pi)*rchannel
x(2) = sin(float(jjjt)/float(npointt)*2.0*pi)*rchannel
x(3) = float(jjjz)/float(npointz)*float(dimz)*delta

x(1) = x(1) + originc(1)
x(2) = x(2) + originc(2)

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
if(indexvolx(jx,jy,jz).eq.0) then

 if(ncha1.eq.maxvolx) then
   write(stdout,*) 'channel: increase maxvolx'
   stop
 endif

 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz
endif

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
com1(indexvolx(jx,jy,jz),:) = com1(indexvolx(jx,jy,jz),:) + x(:)
sumvolx1 = sumvolx1 + 1.0

enddo ! jjjt
enddo ! jjjz

do i = 1, ncha1
com1(i,:) = com1(i,:)/volx1(i)

! Moves the position of the first segment lseg/2 away from the surface to prevent collision due to round errors.

com1(i,1) = com1(i,1) + 0.5*lseg*((com1(i,1)-originc(1)))/rchannel 
com1(i,2) = com1(i,2) + 0.5*lseg*((com1(i,2)-originc(2)))/rchannel 
enddo
end

subroutine integrate_superellipse(rchannel2,RdimZ, originc, npoints,volprot,sumvolprot, flag)
use system
use transform

implicit none
real*8 sumvolprot
integer npoints
real*8 rchannel2, originc(2)
real*8 volprot(dimx,dimy,dimz)
real*8 dr(3), dxr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect
logical flagin, flagout
real*8 intcell_superellipse
real*8 mmmult
integer jx,jy, jz
logical flag
integer RdimZ

real*8 box(4)
real*8 x(3), v(3)
integer xmin,xmax,ymin,ymax,zmin,zmax
integer i,j

logical flagsym
real*8 voltemp

volprot = 0.0
sumvolprot = 0.0 ! total volumen, including that outside the system

! scan over all cells

do ix = 1, dimx
do iy = 1, dimy
do iz = RdimZ+1, dimz-RdimZ

flagin = .false.
flagout = .false.

do ax = 0,1
do ay = 0,1
do az = 0,1

! v in transformed space
v(1) = float(ax+ix-1)*delta
v(2) = float(ay+iy-1)*delta
v(3) = 0.0

! x in real space, v in transformed space
    x = MATMUL(IMAT,v)

x(1) = x(1) - originc(1)
x(2) = x(2) - originc(2)

if((x(1)**2+x(2)**2).lt.rchannel2)flagin=.true. ! inside the channel
if((x(1)**2+x(2)**2).gt.rchannel2)flagout=.true. ! outside the channel

enddo
enddo
enddo

if((flagin.eqv..true.).and.(flagout.eqv..false.)) then ! cell all inside channel
    voltemp = 1.0
endif
if((flagin.eqv..false.).and.(flagout.eqv..true.)) then ! cell all outside channel
    voltemp = 0.0
endif
if((flagin.eqv..true.).and.(flagout.eqv..true.)) then ! cell part inside annd outside channel
    voltemp = intcell_superellipse(rchannel2, originc,ix,iy,iz, npoints)
endif

sumvolprot = sumvolprot + voltemp
volprot(ix,iy,iz) = voltemp

enddo ! ix
enddo ! iy
enddo ! iz


end subroutine

double precision function intcell_superellipse(rchannel2,originc,ix,iy,iz,n)
use system
use transform

implicit none
real*8 rchannel2
real*8 originc(2)
integer ix,iy,iz,ax,ay,az
integer cc
real*8 vect
integer n
real*8 mmmult
real*8 dr(3), dxr(3)

cc = 0
do ax = 1, n
do ay = 1, n
do az = 1, n

dr(1) = ix*delta-(ax)*delta/float(n) 
dr(2) = iy*delta-(ay)*delta/float(n) 
dr(3) = iz*delta-(az)*delta/float(n) 

! dr in transformed space
dxr = MATMUL(IMAT, dr)

dxr(1) = dxr(1)-originc(1)
dxr(2) = dxr(2)-originc(2)

vect = dxr(1)**2+dxr(2)**2
if(vect.lt.rchannel2)cc=cc+1 ! outside channel, integrate

enddo
enddo
enddo

intcell_superellipse = float(cc)/(float(n)**3)
end function


