integer function testsystem_cuboctahedron(x)
use system
use transform
use ellipsoid

implicit none
real*8 x(3), xx(3), v(3), maxx(3)
integer j, i
real*8 vect
real*8 mmmult
real, external :: PBCSYMR, PBCREFR
real*8 dims(3)

dims(1) = delta*dimx
dims(2) = delta*dimy
dims(3) = delta*dimz

maxx(1) = float(dimx)*delta
maxx(2) = float(dimy)*delta
maxx(3) = float(dimz)*delta

! collision with walls and out of system

testsystem_cuboctahedron = 0

v = MATMUL(MAT,x) ! to transformed space

! out-of-boundaries check is performed in transformed space

if (v(1).le.0.0) then
 if (PBC(1).eq.2)testsystem_cuboctahedron = -1
 if (PBC(1).eq.0)testsystem_cuboctahedron = -2
endif

if (v(1).gt.(float(dimx)*delta)) then
 if (PBC(2).eq.2)testsystem_cuboctahedron = -1
 if (PBC(2).eq.0)testsystem_cuboctahedron = -2
endif

if (v(2).le.0.0) then
 if (PBC(3).eq.2)testsystem_cuboctahedron = -1
 if (PBC(3).eq.0)testsystem_cuboctahedron = -2
endif

if (v(2).gt.(float(dimy)*delta)) then
 if (PBC(4).eq.2)testsystem_cuboctahedron = -1
 if (PBC(4).eq.0)testsystem_cuboctahedron = -2
endif

if (v(3).le.0.0) then
 if (PBC(5).eq.2)testsystem_cuboctahedron = -1
 if (PBC(5).eq.0)testsystem_cuboctahedron = -2
endif

if (v(3).gt.(float(dimz)*delta)) then
 if (PBC(6).eq.2)testsystem_cuboctahedron = -1
 if (PBC(6).eq.0)testsystem_cuboctahedron = -2
endif



if (testsystem_cuboctahedron.eq.0) then ! saves some time

do i = 1, 3 ! put into cell
if(v(i).lt.0.0) then
   if(PBC(2*i-1).eq.1)v(i) = PBCSYMR(v(i),maxx(i))
   if(PBC(2*i-1).eq.3)v(i) = PBCREFR(v(i),maxx(i))
endif
if(v(i).gt.dims(i)) then
   if(PBC(2*i).eq.1)v(i) = PBCSYMR(v(i),maxx(i))
   if(PBC(2*i).eq.3)v(i) = PBCREFR(v(i),maxx(i))
endif
enddo


! collision with particles


do j = 1, NNN

xx(:) = x(:) - Rellf(:,j) ! distance to center of CO

v = xx
xx = MATMUL(MAT,v) ! to transformed space

! Looks for near neighbor, only important for PBC
if(PBC(1).eq.1)xx(1) = xx(1) - nint(xx(1) / (float(dimx)*delta)) * float(dimx)*delta
if(PBC(3).eq.1)xx(2) = xx(2) - nint(xx(2) / (float(dimy)*delta)) * float(dimy)*delta
if(PBC(5).eq.1)xx(3) = xx(3) - nint(xx(3) / (float(dimz)*delta)) * float(dimz)*delta
xx(:) = MATMUL(IMAT,xx) ! to real space


! Check in real space for collision with CO

if(((xx(1)+xx(2)+xx(3)).gt.(-Loctall(j)/2)).and.((xx(1)+xx(2)+xx(3)).lt.(Loctall(j)/2)))then
   if(((-xx(1)+xx(2)+xx(3)).gt.(-Loctall(j)/2)).and.((-xx(1)+xx(2)+xx(3)).lt.(Loctall(j)/2)))then
      if(((xx(1)-xx(2)+xx(3)).gt.(-Loctall(j)/2)).and.((xx(1)-xx(2)+xx(3)).lt.(Loctall(j)/2)))then
         if(((-xx(1)-xx(2)+xx(3)).gt.(-Loctall(j)/2)).and.((-xx(1)-xx(2)+xx(3)).lt.(Loctall(j)/2)))then
            if(((abs(xx(1)).lt.(Lcubell(j)/2))).and.(abs(xx(2)).lt.(Lcubell(j)/2)).and.(abs(xx(3)).lt.(Lcubell(j)/2))) then
                 testsystem_cuboctahedron = -1
                 return
            endif
         endif
      endif
   endif
endif

enddo

endif ! testsystem = 0

return

end function

