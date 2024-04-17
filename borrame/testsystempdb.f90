integer function testsystempdb (x)
use system
use ellipsoid
use transform
use pdb
use ematrix
implicit none
real*8 x(3), xx(3), v(3), maxx(3)
integer j, i, ix,iy,iz
real*8 vect
real*8 mmmult
real, external :: PBCSYMR, PBCREFR
real, external :: PBCSYMI, PBCREFI
real*8 dims(3)
real*8 cutoffphi


cutoffphi = 0.5 ! max cutoff value for PDB - chain collision 

dims(1) = delta*dimx
dims(2) = delta*dimy
dims(3) = delta*dimz


maxx(1) = float(dimx)*delta
maxx(2) = float(dimy)*delta
maxx(3) = float(dimz)*delta

v = MATMUL(MAT,x) ! to transformed space

! collision with walls and out of system

testsystempdb = 0

if (v(1).le.0.0) then ! system is checked in transformed space
 if (PBC(1).eq.2)testsystempdb = -1
 if (PBC(1).eq.0)testsystempdb = -2
endif

if (v(1).gt.(float(dimx)*delta)) then
 if (PBC(2).eq.2)testsystempdb = -1
 if (PBC(2).eq.0)testsystempdb = -2
endif

if (v(2).le.0.0) then
 if (PBC(3).eq.2)testsystempdb = -1
 if (PBC(3).eq.0)testsystempdb = -2
endif

if (v(2).gt.(float(dimy)*delta)) then
 if (PBC(4).eq.2)testsystempdb = -1
 if (PBC(4).eq.0)testsystempdb = -2
endif

if (v(3).le.0.0) then
 if (PBC(5).eq.2)testsystempdb = -1
 if (PBC(5).eq.0)testsystempdb = -2
endif

if (v(3).gt.(float(dimz)*delta)) then
 if (PBC(6).eq.2)testsystempdb = -1
 if (PBC(6).eq.0)testsystempdb = -2
endif

if (testsystempdb.eq.0) then ! saves some time

! put chain in lattice

ix = floor(v(1)/delta) + 1
if(ix.lt.1) then
    if(PBC(1).eq.1)ix = PBCSYMI(ix,dimx)
    if(PBC(1).eq.3)ix = PBCREFI(ix,dimx)
endif
if(ix.gt.dimx) then
    if(PBC(2).eq.1)ix = PBCSYMI(ix,dimx)
    if(PBC(2).eq.3)ix = PBCREFI(ix,dimx)
endif

iy = floor(v(2)/delta) + 1
if(iy.lt.1) then
    if(PBC(3).eq.1)iy = PBCSYMI(iy,dimy)
    if(PBC(3).eq.3)iy = PBCREFI(iy,dimy)
endif
if(iy.gt.dimy) then
    if(PBC(4).eq.1)iy = PBCSYMI(iy,dimy)
    if(PBC(4).eq.3)iy = PBCREFI(iy,dimy)
endif

iz = floor(v(3)/delta) + 1
if(iz.lt.1) then
    if(PBC(5).eq.1)iz = PBCSYMI(iz,dimz)
    if(PBC(5).eq.3)iz = PBCREFI(iz,dimz)
endif
if(iz.gt.dimz) then
    if(PBC(6).eq.1)iz = PBCSYMI(iz,dimz)
    if(PBC(6).eq.3)iz = PBCREFI(iz,dimz)
endif
 
! check volprot 
! collision with pdb

 if(volprot(ix,iy,iz).gt.cutoffphi) then 
  testsystempdb  = -1
 endif

endif ! testsystempdb = 0

return

end function
