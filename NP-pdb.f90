subroutine update_matrix_80(flag)
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
! ADD NP AND POLYMERS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j = 1, 1

! rotate ellipsoid matrixes according to current rotation matrix

 call rotv(AAA(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAS(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAL(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAX(:,:,j), rotmatrix(:,:,j))

 call rotvo(orient(:,j), rotmatrix(:,:,j))

 npoints = 50

 flag = .false.


 call integrate(AAAL(:,:,j),AellL(:,j), Rell(:,j),npoints, voleps1 , sumvoleps1,flag)
 flag = .false. ! not a problem if eps lays outside boundaries
 call integrate(AAA(:,:,j),Aell(:,j), Rell(:,j),npoints, volprot1, sumvolprot1,flag)
 call integrate(AAAS(:,:,j),AellS(:,j), Rell(:,j),npoints, volq1, sumvolq1,flag)

 npoints = 100000000
 call newintegrateg(Aell(:,j),Rell(:,j),npoints,volx1,sumvolx1, com1, p1, ncha1,volxx1)

!! volume
 temp = 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)/(sumvolprot1*delta**3) !
!rescales volume
 volprot1 = volprot1*temp       
                                          ! OJO:
 !transformation should mantain cell volumen
 sumvolprot1 = sumvolprot1*temp

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eeps(j)

!! charge
 volq1 = volprot1-volq1
 temp = sumvolprot1-sumvolq1
 volq1 = volq1/temp*echarge(j)/(delta**3) ! sum(volq) is echarge

!! grafting
 pnumber = 1.6075

 area = (Aell(1,j)*Aell(2,j))**pnumber
 area = area+(Aell(1,j)*Aell(3,j))**pnumber
 area = area+(Aell(2,j)*Aell(3,j))**pnumber
 area = 4.0*pi*(area/3.0)**(1.0/pnumber) ! approximate (< 1% error) area of elipsoid, see wikipedia

 temp2 = maxval(volx1)

 where(volx1<temp2*cutarea) ! remove cells with very little area
 volx1 = 0.0
 end where

 volx1 = volx1/sumvolx1*area*sigma(j)
 volxx1 = volxx1/sumvolx1*area*sigma(j)

 maxss = 1.0d100
sumpolseg = sumpolseg + area*sigma(j)*long

!! volume  
 volprot1 = volprot1 * 0.99
 volprot = volprot+volprot1

! CHECK COLLISION HERE...
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true.
   exit
 endif

 voleps = voleps + voleps1
 volq = volq + volq1

! add com1 and volx to list

 volxx = volxx + volxx1


do i = 1, ncha1
ncha = ncha+1
volx(ncha)=volx1(i)
com(ncha,:)=com1(i,:)
p0(ncha,:)=p1(i,:)
enddo

enddo





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

 npoints = 20
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
write(stdout,*) 'NP-PDB:', 'update_matrix: Total discretized volumen =', (dimx*dimy*dimz-sum(volprot))*delta**3
endif

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

end subroutine


