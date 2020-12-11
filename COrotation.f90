subroutine COrotation(rotmatrixCO,lcubeCO,loctaCO)

use ellipsoid
use COrotMod

implicit none

real*8 rotmatrixCO(3,3)
real*8 lcubeCO, loctaCO

! Plane x + y + z = k
plane1 = (/1.0,1.0,1.0/)
plane1 = MATMUL(rotmatrixCO,plane1)
! Plane - x + y + z = k
plane2 = (/-1.0,1.0,1.0/)
plane2 = MATMUL(rotmatrixCO,plane2)
! Plane x - y + z = k
plane3 = (/1.0,-1.0,1.0/)
plane3 = MATMUL(rotmatrixCO,plane3)
! Plane - x - y + z = k
plane4 = (/-1.0,-1.0,1.0/)
plane4 = MATMUL(rotmatrixCO,plane4)
! Plane x = k
plane5 = (/1.0,0.0,0.0/)
plane5 = MATMUL(rotmatrixCO,plane5)
! Plane y = k
plane6 = (/0.0,1.0,0.0/)
plane6 = MATMUL(rotmatrixCO,plane6)
! Plane z = k
plane7 = (/0.0,0.0,1.0/)
plane7 = MATMUL(rotmatrixCO,plane7)

! Plane constants
klocta1v(1) = 0.0
klocta1v(2) = 0.0
klocta1v(3) = loctaCO/2

klocta2v(1) = 0.0
klocta2v(2) = 0.0
klocta2v(3) = -loctaCO/2

klocta1 = DOT_PRODUCT(MATMUL(rotmatrixCO,klocta1v),plane1)
klocta1b = DOT_PRODUCT(MATMUL(rotmatrixCO,klocta2v),plane1)
klocta2 = DOT_PRODUCT(MATMUL(rotmatrixCO,klocta1v),plane2)
klocta2b = DOT_PRODUCT(MATMUL(rotmatrixCO,klocta2v),plane2)
klocta3 = DOT_PRODUCT(MATMUL(rotmatrixCO,klocta1v),plane3)
klocta3b = DOT_PRODUCT(MATMUL(rotmatrixCO,klocta2v),plane3)
klocta4 = DOT_PRODUCT(MATMUL(rotmatrixCO,klocta1v),plane4)
klocta4b = DOT_PRODUCT(MATMUL(rotmatrixCO,klocta2v),plane4)

klcubex1v(1) = lcubeCO/2
klcubex1v(2) = 0.0
klcubex1v(3) = 0.0
klcubex1 = DOT_PRODUCT(MATMUL(rotmatrixCO,klcubex1v),plane5)

klcubex2v(1) = -lcubeCO/2
klcubex2v(2) = 0.0
klcubex2v(3) = 0.0
klcubex2 = DOT_PRODUCT(MATMUL(rotmatrixCO,klcubex2v),plane5)

klcubey1v(1) = 0.0
klcubey1v(2) = lcubeCO/2
klcubey1v(3) = 0.0
klcubey1 = DOT_PRODUCT(MATMUL(rotmatrixCO,klcubey1v),plane6)

klcubey2v(1) = 0.0
klcubey2v(2) = -lcubeCO/2
klcubey2v(3) = 0.0
klcubey2 = DOT_PRODUCT(MATMUL(rotmatrixCO,klcubey2v),plane6)

klcubez1v(1) = 0.0
klcubez1v(2) = 0.0
klcubez1v(3) = lcubeCO/2
klcubez1 = DOT_PRODUCT(MATMUL(rotmatrixCO,klcubez1v),plane7)

klcubez2v(1) = 0.0
klcubez2v(2) = 0.0
klcubez2v(3) = -lcubeCO/2
klcubez2 = DOT_PRODUCT(MATMUL(rotmatrixCO,klcubez2v),plane7)

end subroutine COrotation
