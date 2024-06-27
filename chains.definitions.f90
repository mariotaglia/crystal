subroutine chains_definitions

use chainsdat
use MPI
use branches
implicit none
integer i, ii

! Read ligand sequence 
ALLOCATE (segtype(long))

open(file="sequence.in",unit=333)
do i=1,long
 read(333,*) segtype(i)
enddo

close(333)

! Read solvent sequence
ALLOCATE (segtypesv(longsv))

open(file="sequencesv.in",unit=333)
do i=1,longsv
 read(333,*) segtypesv(i)
enddo

close(333)

!if(branched.eq.1) then
!   segtype(1:longbb+longb(1)) = 2 ! backbone and first branch is hydrophobic
!endif

!if(branched.ne.1) segtype = 2

end
