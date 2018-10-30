module ImpactFEL
  use ConstDataModule
  use mathModule
  use pGrid2DModule
  use CompDomModule
  use RadiationModule
  use eBeamModule
  use ElementModule
  use UndulatorModule
  use QuadModule
  use DriftModule
  use PropagatorModule
  use DipoleModule
  use DiagnoElemModule
  use fileIO_Module
  
  ! type, private :: utilities
    ! contains
      ! procedure
  ! end type
  ! type(utilities), public :: until
contains
  
  
  subroutine write_undulator_vtk(lattice,beam,fID)
    implicit none
    type(arrElem),intent(in) :: lattice(:)
    type(eBeam),  intent(in) :: beam
    integer,intent(in) :: fID
    integer :: nElem,i,nSec,iUnit
    real*8 :: z,zin,zout
    character(7) :: cfID
    logical :: flagIn
    nElem = size(lattice)
    
    zin=0
    iUnit = fileIO%get_free_unit()
    write(cfID,'(I0)') fID
    open(iUnit,file='undulator.'//trim(adjustl(cfID))//'.vtk',action='write')    
    write(iUnit,*) '# vtk DataFile Version 2.0'
    write(iUnit,*) 'Lattice'
    write(iUnit,*) 'ASCII'
    write(iUnit,*) 'DATASET POLYDATA'
    
    flagIn = .false.
    nSec = 0
    do i=1,nElem
      select type(elemp => lattice(i)%op)
        type is (Undulator)
          if(.not. flagIn) nSec = nSec+1
          flagIn = .true.
        class default
          if(lattice(i)%op%L > 0) then
            flagIn = .false.
          endif
      end select
    enddo
    
    write(iUnit,*) 'POINTS', nSec*4, 'float'
    z =0 
    zin = 0
    flagIn = .false.
    do i=1,nElem
      select type(elemp => lattice(i)%op)
        type is (Undulator)
          if(.not. flagIn) zin = z
          flagIn = .true.
        class default
          if(lattice(i)%op%L > 0 .and. flagIn) then
            zout = z + lattice(i)%op%L
            write(iUnit,*) 1.2*beam%cDom%rng(1,1),beam%cDom%rng(2,2),zin*beam%ku/beam%ks
            write(iUnit,*) 1.2*beam%cDom%rng(1,2),beam%cDom%rng(2,2),zin*beam%ku/beam%ks
            write(iUnit,*) 1.2*beam%cDom%rng(1,1),beam%cDom%rng(2,2),zout*beam%ku/beam%ks
            write(iUnit,*) 1.2*beam%cDom%rng(1,2),beam%cDom%rng(2,2),zout*beam%ku/beam%ks
            flagIn = .false.
          endif
      end select
      if(i==nElem .and. flagIn) then
        zout = z + lattice(i)%op%L
        write(iUnit,*) 1.2*beam%cDom%rng(1,1),beam%cDom%rng(2,2),zin*beam%ku/beam%ks
        write(iUnit,*) 1.2*beam%cDom%rng(1,2),beam%cDom%rng(2,2),zin*beam%ku/beam%ks
        write(iUnit,*) 1.2*beam%cDom%rng(1,1),beam%cDom%rng(2,2),zout*beam%ku/beam%ks
        write(iUnit,*) 1.2*beam%cDom%rng(1,2),beam%cDom%rng(2,2),zout*beam%ku/beam%ks
        flagIn = .false.
      endif
      z = z + lattice(i)%op%L
    enddo
    write(iUnit,*) 'POLYGONS',nSec,'8'
    do i=1,nSec
      write(iUnit,*) '4',i*4-4,i*4-3,i*4-2,i*4-1
    enddo
    close(iUnit)
  end subroutine write_undulator_vtk
  
  
  
  ! subroutine find_periodic_matching_condition(betx,alfx,psix,bety,alfy,psiy,lattice,gamma)
  ! !=============================================================
  ! ! input
  ! !   type(arrElem) :: lattice(:) 
  ! !   real(8) :: gamma = reference energy
  ! ! output
  ! !   real(8) :: betx, bety, psix, alfx, alfy, psiy
  ! !              = twiss beta,alpha, and phase advance(in unit of 2pi)
  ! !-------------------------------------------------------------
    ! implicit none
    ! type(arrElem), intent(in) :: lattice(:)
    ! real(8), intent(in)  :: gamma
    ! real(8), intent(out) :: betx, bety, alfx, alfy, psix, psiy
    ! integer :: i
    ! real(8), dimension(4,4) :: M, M1
    
    ! M = 0d0
    ! do i=1,4
      ! M(i,i) = 1d0
    ! enddo
    ! do i=1,size(lattice)
      ! call lattice(i)%op%get_M(M1,gamma)
      ! M = matmul(M1,M)
    ! enddo
    ! psix = acos((M(1,1)+M(2,2))/2d0)
    ! psiy = acos((M(3,3)+M(4,4))/2d0)
  ! end subroutine find_periodic_matching_condition
  
  
end module ImpactFEL
