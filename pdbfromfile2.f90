subroutine pdbfromfile
use pdb
use ellipsoid
use system
implicit none
integer i
integer ix,iy,iz
integer NN
!aaID = 0.0 ! numero de aminoacido en ix,iy,iz

! read aa from disk
open(file='pdb.in', unit=3333)
read(3333,*)naa ! number of beads

! allocate amino acid properties
allocate(aapos(3,naa)) ! position of aminoacid naa
allocate(aal(naa))     ! string with aminoacid ID (one-letter code)
allocate(aan(naa))     ! number of aminoacid
allocate(xxpdb(naa)) ! keeps info of position of original aminoacids
allocate(yypdb(naa)) ! keeps info of position of original aminoacids
allocate(zzpdb(naa)) ! keeps info of position of original aminoacids
allocate(fdispdb(naa)) ! keeps info of position of original aminoacids

! allocate discretization list
!allocate(maxelement_list(naa))
!allocate(coords_list(naa,3,maxel))
!allocate(vol_list(naa,maxel))

!maxelement_list = 0
!coords_list = 0
!vol_list = 0.0

! read aa pos from file
do i = 1, naa
read(3333,*)aapos(1,i),aapos(2,i),aapos(3,i),aal(i),aan(i)
enddo

call assign_aa
if (systemtype.eq.90) NN = 1
if (systemtype.eq.80) NN = 2

! translate  
do i = 1, naa
aapos(1,i) = aapos(1,i) + Rell(1,NN) ! translate
aapos(2,i) = aapos(2,i) + Rell(2,NN)
aapos(3,i) = aapos(3,i) + Rell(3,NN)
enddo

write(*,*)'Prot pos', Rell(1,NN), Rell(2,NN), Rell(3,NN)


! print real coordinates to file
open(file='real_coords.txt',unit=1111)
do i = 1, naa
write(1111,*)i, aan(i), aal(i), aapos(1,i), aapos(2,i), aapos(3,i)
enddo
close(1111) 

! generate amino-acid discretization and generate lists
!volprot = 0.0

!print*, 'Generating aa discretization'

do i = 1, naa
!center(:) = aapos(i,:)
!print*, i
!call sphere(radiuspdb(i), center, protn)
!call matrixtolist(protn,i)
!volprot(:,:,:) = volprot(:,:,:) + protn(:,:,:)/(delta**3)

! PROJECTS TO THE LATTICE

ix=int(aapos(1,i)/delta)+1
iy=int(aapos(2,i)/delta)+1
iz=int(aapos(3,i)/delta)+1

xxpdb(i) = ix
yypdb(i) = iy
zzpdb(i) = iz

!aaID(ix,iy,iz) = aan(i)
enddo ! loop over number of aa, i

! evite que la fraccion de volumen sea mayor que 1 (ej. dos aminoacidos en una misma celda
!where (volprot > 1.0) volprot = 1.0

!title = 'avpro'
!counter = 1
!call savetodisk(volprot, title, counter)

!title = 'aaID_'
!counter = 1
!call savetodisk(aaID, title, counter)

close(3333)

end subroutine



subroutine assign_aa
use pdb
use const
use bulk
use molecules
implicit none

      integer i

      ALLOCATE (zpdb(naa))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (pKapdb(naa), Kapdb(naa), K0pdb(naa), radiuspdb(naa))

      do i = 1, naa
      select case (aal(i))

! Chain definition following Biophys Journal 107 1393-1402
! kai integration routine also changed

      case('A')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 28.6 !!!! Molar volumes (cm^3/mol), transformed to nm at the end of subroutine

      case('I')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 75.8

      case('L')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 75.8

      case('F')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 89.8

      case('W')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 112.2

      case('Y')
      zpdb(i) = -1
      pKapdb(i) = 9.57 ! 10.5
      radiuspdb(i) = 91.9

      case('K')
      zpdb(i) = 1
      pKapdb(i) = 9.5 !11.43 ! 10.54
      radiuspdb(i) = 77.3

      case('R')
      zpdb(i) = 1
      pKapdb(i) = 13.37 ! 12.48
      radiuspdb(i) = 94.6

      case('N')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 45.9

      case('Q')
      zpdb(i) =  0
      pKapdb(i) = 0.0
      radiuspdb(i) = 62.2

      case('M')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) =  73.4  ! 45 corre

      case('P')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 51.13

      case('S')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 28.5

      case('T')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 45.1

      case('V')
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 59.6

      case('D')
      zpdb(i) = -1
      pKapdb(i) = 2.97 !3.9
      radiuspdb(i) = 42.0

      case('E')
      zpdb(i) = -1
      pKapdb(i) = 3.13 ! 4.07
      radiuspdb(i) = 56.42

      case('C')    ! Reduced Cysteine !
      zpdb(i) = -1
      pKapdb(i) = 8.37
      radiuspdb(i) = 41.7

      case('X')    ! Oxidized Cysteine !
      zpdb(i) = 0
      pKapdb(i) = 0.0
      radiuspdb(i) = 41.7

      case('H')
      zpdb(i) = 1
      pKapdb(i) =  6.98 !6.04
      radiuspdb(i) = 67.1
      
      case('Hm') !CASE HISTIDINE COORDINATED TO HEME C
      zpdb(i) = 0
      radiuspdb(i) = 67.1
      pKapdb(i) = 14.0
    
      case('Z')    ! CASE FLUOROPHORE OF GFP !
      zpdb(i) = 0
      pKapdb(i) = 0
      radiuspdb(i) = 67.1  ! SAME VOLUME AS HIS !

      case('B', 'G')
      radiuspdb(i) = 31.7
      zpdb(i) = 0
      pKapdb(i) = 0.0
        if(aan(i).eq.1) then ! N terminal
          zpdb(i) = 0 ! 1 
          pKapdb(i) = 10.4333 !9.5
        endif
        if(aan(i).eq.aan(naa)) then ! N terminal
          zpdb(i) = -1
          pKapdb(i) = 3.57 !4.5
        endif

      case('Fe2')
      zpdb(i) = 2
      pKapdb(i) = 14 !cheuqear
      radiuspdb(i) = 1.1  !chequear 1.1
      case('Fe3')
      zpdb(i) = 3.0
      pKapdb(i) = 1 !chequear
      radiuspdb(i) = 0.83 !chequear  0.83  ok
      case('Pi')
      zpdb(i) = 0 !chequear
      pKapdb(i) = 14 !chequear
      radiuspdb(i) = 67.1 !chequear  67.1 con 20 corre con 4 pirroles
      case('Pm')
      zpdb(i) = 0 !chequear
      pKapdb(i) = 14 !chequear
      radiuspdb(i) = 20.0 !chequear  67.1 con 20 corre con 4 pirroles

      case('Pr')
      zpdb(i) = -1 !chequear
      pKapdb(i) = 4.5 !chequear
      radiuspdb(i) = 56.4 !chequear
      case('Cs')
      zpdb(i) = 0 !chequear
      pKapdb(i) = 14 !chequear
      radiuspdb(i) = 30.1 !45.1  chequear  45.1





      case default
        print*, 'aminoacid not recognized. stop'
        write(*,*) aan(i)

        stop

      endselect

      radiuspdb(i) = (radiuspdb(i)*(1.0d21/6.02d23)/(4.0/3.0*pi))**(1.0/3.0)

      Kapdb(i)=10**(-pKapdb(i)) ! calculate thermodynamic equilibrium constants





      enddo ! loop sobre naa

      end

