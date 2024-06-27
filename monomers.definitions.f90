subroutine monomer_definitions

use MPI
use chainsdat
use mparameters_monomer

implicit none
integer i, j

N_monomer = N_poorsol + 1 ! 0: hydrophilic, 1 < x < N_poorsol 

ALLOCATE (st_matrix(0:N_poorsol, 0:N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent

! ELECTRO
!ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
!ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Cross interactions as Berthelot Rules
open(file="epsilon.in",unit=333)
st_matrix(0,0)=0.0
do i = 1, N_poorsol
  read(333, *) st_matrix(i,i)
enddo

do i = 1, N_poorsol
  do j = 1, N_poorsol
    st_matrix(i,j) = sqrt(st_matrix(i,i)*st_matrix(j,j))
  enddo
enddo

! Segment type 1 for NPC, positive base, hydrophilic

i = 1
hydroph(i) = 1

! ELECTRO
!zpol(i) = 0
!pKa(i) = 7.0


!i = 2
!zpol(i) = 0
!hydroph(i) = 1
!pKa(i) = 11.0

!! Segment type 2 for NPC, negative , hydrophilic
!
!zpol(2) = -1
!hydroph(2) = 0
!pKa(2) = 5.0
!
!! Segment type 3 for NPC, neutral , hydrophilic
!
!zpol(3) = 0
!hydroph(3) = 0
!pKa(3) = 1 ! set any number if zpol = 0...
!
!! Segment type 4 for NPC , neutral, hydrophobic, 1
!
!zpol(4) = 0
!hydroph(4) = 1
!pKa(4) = 1 ! set any number if zpol = 0....
!
end

