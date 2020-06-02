subroutine update_matrix_cuboctahedron(flag)
use system
use cuboctahedron !cuidado agregar channel quizas
use ematrix
use MPI
use const
use chainsdat
use molecules
use channel
use transform, only : MAT, IMAT
use rotchain

implicit none

integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
logical flag
integer j,ix,iy,iz
real*8 lcubeL, lcubeS, loctaL, loctaS
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


cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
sumpolseg = 0.0
lcubeL = lcube + 2*delta
lcubeS = lcube - 2*delta
loctaL = locta + 2*delta
loctaS = locta - 2*delta

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

 call integrate_cuboctahedron(lcubeL,loctaL,center,npoints,voleps1,sumvoleps1,flag)

 flag = .false. ! not a problem if eps lays outside boundaries
 
 call integrate_cuboctahedron(lcube,locta,center,npoints,volprot1,sumvolprot1,flag)

 call integrate_cuboctahedron(lcubeS,loctaS,center,npoints,volq1,sumvolq1,flag)

 npoints = 100
 call newintegrateg_cuboctahedron(lcube,locta,center,npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eepsc

!! charge
 volq1 = volprot1-volq1
 temp = sum(volq1)
 volq1 = volq1/temp*echargec/(delta**3) ! sum(volq) is echarge

! area = 6.0*lcube**2   !no seeee

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

 volxx = volxx1

ncha = ncha1
do i = 1, ncha
volx(i)=volx1(i)
com(i,:)=com1(i,:)
p0(i,:)=p1(i,:)
rotangle(i) = 0.0
enddo

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

stop

sumpolseg = ncha

if (verbose.ge.2) then
temp = lcube**3 ! no seeee
!do j = 1, NNN
!temp = temp + 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)
!enddo
if (rank.eq.0) then
write(stdout,*) 'cuboctahedron:', 'update_matrix: Total nanocuboct volumen real space= ', temp
write(stdout,*) 'cuboctahedron:', 'update_matrix: Total discretized volumen =', (sum(volprot))*delta**3
write(stdout,*) 'cuboctahedron:', 'number of polymers in system =', sumpolseg
write(stdout,*) 'cuboctahedron:', 'surface area =', area
write(stdout,*) 'cuboctahedron:', 'surface density =', sumpolseg/area
endif
endif

title = 'aveps'
counter = 1
!call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
!call savetodisk(volq, title, counter)

title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)

end subroutine

subroutine integrate_cuboctahedron(lcube,locta,center,npoints,volprot,sumvolprot,flag)
use system
use transform

implicit none
real*8 sumvolprot
integer npoints
real*8 lcube, locta
real*8 center(3)
real*8 volprot(dimx,dimy,dimz)
real*8 dr(3), dxr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect
logical flagin, flagout
real*8 intcell_cuboctahedron
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
do iz = 1, dimz

flagin = .false.
flagout = .false.

do ax = 0,1
do ay = 0,1
do az = 0,1

! v in transformed space
v(1) = float(ax+ix-1)*delta
v(2) = float(ay+iy-1)*delta
v(3) = float(az+iz-1)*delta

! x in real space, v in transformed space
    x = MATMUL(IMAT,v)

x(1) = x(1) - center(1)
x(2) = x(2) - center(2)
x(3) = x(3) - center(3)

if(((x(1)+x(2)+x(3)).gt.(-locta/2)).and.((x(1)+x(2)+x(3)).lt.(locta/2)))then
   if(((-x(1)+x(2)+x(3)).gt.(-locta/2)).and.((-x(1)+x(2)+x(3)).lt.(locta/2)))then
      if(((x(1)-x(2)+x(3)).gt.(-locta/2)).and.((x(1)-x(2)+x(3)).lt.(locta/2)))then
         if(((-x(1)-x(2)+x(3)).gt.(-locta/2)).and.((-x(1)-x(2)+x(3)).lt.(locta/2)))then
            if(((abs(x(1)).lt.(lcube/2)).and.(abs(x(2)).lt.(lcube/2))).and.(abs(x(3)).lt.(lcube/2)))then
               flagin=.true.
            else
               flagout=.true.
            endif
         else
            flagout=.true.
         endif
      else
         flagout=.true.
      endif
   else
      flagout=.true.
   endif
else
   flagout=.true.
endif

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
    voltemp = intcell_cuboctahedron(lcube,locta,center,ix,iy,iz,npoints)
endif

sumvolprot = sumvolprot + voltemp
volprot(ix,iy,iz) = voltemp

enddo ! ix
enddo ! iy
enddo ! iz

end subroutine

double precision function intcell_cuboctahedron(lcube,locta,center,ix,iy,iz,n)
use system
use transform

implicit none
real*8 lcube,locta
real*8 center(3)
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

dxr(1) = dxr(1) - center(1)
dxr(2) = dxr(2) - center(2)
dxr(3) = dxr(3) - center(3)

!if (((abs(dxr(1)).lt.(l_cube/2)).and.(abs(dxr(2)).lt.(l_cube/2))).and.(abs(dxr(3)).lt.(l_cube/2)))cc=cc+1 ! integra dentro del cubo

if(((dxr(1)+dxr(2)+dxr(3)).gt.(-locta/2)).and.((dxr(1)+dxr(2)+dxr(3)).lt.(locta/2)))then
   if(((-dxr(1)+dxr(2)+dxr(3)).gt.(-locta/2)).and.((-dxr(1)+dxr(2)+dxr(3)).lt.(locta/2)))then
      if(((dxr(1)-dxr(2)+dxr(3)).gt.(-locta/2)).and.((dxr(1)-dxr(2)+dxr(3)).lt.(locta/2)))then
         if(((-dxr(1)-dxr(2)+dxr(3)).gt.(-locta/2)).and.((-dxr(1)-dxr(2)+dxr(3)).lt.(locta/2)))then
            if (((abs(dxr(1)).lt.(lcube/2)).and.(abs(dxr(2)).lt.(lcube/2))).and.(abs(dxr(3)).lt.(lcube/2)))then
            cc=cc+1
            endif
         endif
      endif
   endif
endif

enddo
enddo
enddo

intcell_cuboctahedron  = float(cc)/(float(n)**3)
end function


subroutine newintegrateg_cuboctahedron(lcube,locta,center,npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)
use system
use transform
use chainsdat
use ematrix
use const

implicit none
real*8 sumvolx1
integer npoints
!real*8 AAA(3,3), AAAX(3,3)
integer indexvolx(dimx,dimy,dimz)
integer listvolx(ncha,3)
real*8 Rell(3), Aell(3)
real*8 radio
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
real*8 lcuber, pasoc, pasoo
real*8 xx, yy, zz
real*8 vector(3)


indexvolx = 0
ncha1 = 0
volx1 = 0.0
sumvolx1 = 0.0 ! total volumen, including that outside system
com1 = 0.0
p1 = 0
volxx1 = 0.0

! This routine determines the surface coverage and grafting positions only for spheres
!
lcuber = lcube/locta

pasoc = delta/float(npoints)/locta

pasoo = pasoc*(1/3)^(1/4) ! % different integration steps are needed to have the same area element
                          ! % Lets R(u,v) = x(u,v)i + y(u,v)j + z(u,v)k
                          ! % then dA = dR/du x dR/dv * du*dv
                          ! % For octahedro x = u, y = v, z = u+v-1
                          ! % For cube x = v-1, y = 0, z = 0  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%% OCTAHEDRO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx = 0.0

do while (xx < 1.0)
    yy = 0.0
    do while (yy < (1.0-xx))
    zz = -xx-yy+1
    
    if((abs(xx).lt.lcuber).and.(abs(yy).lt.lcuber).and.(abs(zz).lt.lcuber)) then 

       x(1) = xx
       x(2) = yy
       x(3) = -xx-yy+1
       call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

       x(1) = -xx
       x(2) = -yy
       x(3) = -xx-yy+1
       call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = -xx
       x(2) = yy
       x(3) = -xx-yy+1
       call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = xx
       x(2) = -yy
       x(3) = -xx-yy+1
       call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

       x(1) = xx
       x(2) = yy
       x(3) = xx+yy-1
       call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = -xx
       x(2) = -yy
       x(3) = +xx+yy-1
       call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = -xx
       x(2) = yy
       x(3) = +xx+yy-1
       call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = xx
       x(2) = -yy
       x(3) = +xx+yy-1 
       call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
     
       endif    
        
        yy = yy + pasoo
    enddo    
    xx = xx + pasoo
enddo


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%% CUBO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xx = -lcuber
do while (xx < lcuber)
  yy = -lcuber
  do while (yy < lcuber)

  zz = lcuber    
        
  if (((xx+yy+zz).gt.-1.0).and.((xx+yy+zz).lt.1.0)) then
    if(((-xx+yy+xx).gt.-1.0).and.((-xx+yy+zz).lt.1.0)) then
       if(((xx-yy+zz).gt.-1.0).and.((xx-yy+zz).lt.1.0)) then 
          if(((-xx-yy+zz).gt.-1.0).and.((-xx-yy+zz).lt.1.0)) then
                
        x(1) = xx
        x(2) = yy
        x(3) = lcuber
        call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = xx
        x(2) = yy
        x(3) = -lcuber
        call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
        x(1) = lcuber
        x(2) = xx
        x(3) = yy
        call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = -lcuber
        x(2) = xx
        x(3) = yy
        call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = xx
        x(2) = lcuber
        x(3) = yy
        call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = xx
        x(2) = -lcuber
        x(3) = yy
        call integrar_matrices(x,centro,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

          endif 
       endif
    endif
  endif
        
  yy = yy + pasoc
  enddo  
  xx = xx + pasoc
enddo 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aca termina de tirar puntos en la superficie del cubooctahedro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i = 1, ncha1
com1(i,:) = com1(i,:)/volx1(i)
! Moves the position of the first segment lseg/2 away from the surface to prevent collision due to round errors.
vector(:) = com1(i,:)-centro(:)
vector(:) = vector(:)/norm2(vector)
com1(i,:) = com1(i,:) + 1.5*lseg*vector(:)
enddo

end


subroutine integrar_matrices(x, centro, locta, indexvolx, ncha1, p1, volxx1, volx1, com1, sumvolx1)
implicit none
use system
real*8, x(3), centro(3), locta
integer flagin
integer j
integer is(3), js(3), dims(3)
integer, external :: PBCSYMI
integer jx,jy,jz
integer ncha1
integer indexvolx(dimx,dimy,dimz)
integer p1(maxvolx,3)
real*8 volxx1(dimx,dimy,dimz)
real*8 com1(maxvolx,3)
real*8 volx1(maxvolx)

x(:) = x(:)*locta/2.0 + centro(:)

dims(1) = dimx
dims(2) = dimy
dims(3) = dimz

! x in real space, v in transformed space
    v = MATMUL(MAT,x)

! PBC

flagin = 1

do j = 1,3

    is(j) = floor(v(j)/delta)+1
    js(j) = is(j)

select case (PBC((j-1)*2+1))
  case (0 , 2)
    if(is(j).lt.1) then
    write(stdout,*) 'Error in newintegrateg: out of boundary'
    endif
  case (1)
    js(j)=PBCSYMI(is(j), dims(j)) 
  case (3)
    if(v(j).lt.0.0)flagin=0
endselect

select case (PBC((j-1)*2+2))
  case (0 , 2)
    if(is(j).gt.dims(j)) then
    write(stdout,*) 'Error in newintegrateg: out of boundary'
    endif
  case (1)
    js(j)=PBCSYMI(is(j), dims(j)) 
  case (3)
    if(v(j).gt.float(dims(j))*delta)flagin=0
endselect
enddo

jx = js(1)
jy = js(2)
jz = js(3)

if(flagin.eq.1) then

! increase counter
if(indexvolx(jx,jy,jz).eq.0) then

 if(ncha1.eq.maxvolx) then
   write(stdout,*) 'ellipsoid: increase maxvolx'
   stop
 endif

 ncha1 = ncha1 + 1
 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz
endif

! agrega el punto
volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
com1(indexvolx(jx,jy,jz),:) = com1(indexvolx(jx,jy,jz),:) + x(:)
endif

sumvolx1 = sumvolx1 + 1.0

end

