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
