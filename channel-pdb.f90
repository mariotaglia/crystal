subroutine update_matrix_70(flag)
use system
use ellipsoid
use channel
use ematrix
use MPI
use const
use chainsdat
use molecules
use channel
use transform, only : MAT, IMAT
use rotchain
use pdb
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
real*8 x(3), v(3), hcyl
integer nbands
real*8 AAApdb(3,3), Aellpdb(3)

call make_ellipsoid ! update matrixes for all particles

cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
sumpolseg = 0.0
rchannel2 = rchannel**2

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADD CHANNEL AND POLYMERS ON CHANNEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! channel center in x, y plane

 originc(1) = float(dimx)*delta/2.0 
 originc(2) = float(dimy)*delta/2.0 

 npoints = 50

 flag = .false.

 call integrate_c(rchannel2,RdimZ, originc,npoints, volprot1, sumvolprot1, flag)
 call newintegrateg_c_4(rchannel2,RdimZ,originc,npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1, NBRUSH)

!! grafting

v(1) = 0.0
v(2) = 0.0
v(3) = float(dimz-2*RdimZ)*delta

! v in transformed space, x in real space
! only work for gam = 90, cdiva any value

x = MATMUL(IMAT,v)

hcyl = x(3) ! height of the cylinder

area = 2.0*pi*rchannel*hcyl

! add com1 and volx to list

 volxx = volxx1

ncha = ncha1
do i = 1, ncha
volx(i)=volx1(i)
com(i,:)=com1(i,:)
p0(i,:)=p1(i,:)
rotangle(i) = atan2(com1(i,1)-originc(1), com1(i,2)-originc(2))
enddo

 volprot = volprot+volprot1*0.99 ! sum channel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADD PROTEIN FROM PDB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call pdbfromfile ! read protein from pdb



volprot1 = 0.0

do j = 1, naa ! loop over all beads
! integrate volume of aminoacid j

! radius

 AAApdb = 0.0
 do i = 1,3
 AAApdb(i,i) = 1.0/(radiuspdb(j)**2)
 Aellpdb(i) = radiuspdb(j)
 enddo

 npoints = 5
 call integrate(AAApdb(:,:),Aellpdb(:), aapos(:,j),npoints, volpdb, sumvolpdb, flag)

 volprot1 = volprot1 + volpdb

enddo ! j

where (volprot1 > 1.0) volprot1 = 1.0 ! protein has a maximum volume fraction of 1.

sumvolprot1 = sum(volprot1)

!! volume

volprot1 = volprot1 * 0.99
volprot = volprot+volprot1 ! sum particle to channel


! CHECK COLLISION HERE...
if(maxval(volprot).gt.1.0) then ! collision
   flag=.true.
endif

if (rank.eq.0) then
title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)
title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)
title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)
title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)
endif

sumpolseg = ncha

if (rank.eq.0) then
write(stdout,*) 'channel-part:', 'update_matrix: Total discretized volumen =', (dimx*dimy*dimz-sum(volprot))*delta**3
endif

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

end subroutine


