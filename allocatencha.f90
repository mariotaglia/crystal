subroutine allocatencha

use fields_fkfun
use chainsdat
use conformations
use rotchain
use solventchains
implicit none

! fields_fkfun
ALLOCATE(q(ncha))
ALLOCATE(sumgauche(ncha))
ALLOCATE(ngauche(cuantas,ncha))

! chainsdat
allocate(posicion(ncha,3))
allocate(ngpol(ncha))
allocate(newcuantas(ncha))

! solventchains
ALLOCATE(ngauchesv(cuantassv))

end subroutine
