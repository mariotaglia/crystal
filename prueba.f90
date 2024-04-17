nx=4
ny=4
nz=4

do ix = 0,nx-1
 do iy = 0,ny-1
  do iz = 1,nz-2

  dimz=nz-2
  dimy=ny
  dimx=nx

  print*, ix,iy,iz, iz+dimz*iy+dimz*dimy*ix - 1

enddo
enddo
enddo




end
