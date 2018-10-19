module PropagatorModule
  use ElementModule
  use ConstDataModule
  use eBeamModule
  use RadiationModule
  implicit none
  private
  type, public, extends(Element) :: Propagator
    contains
      procedure :: act_Rad => propagate
      !========== dummy overloadings ==========
      !procedure :: get_M => dummy_get_M
      procedure :: act_Beam => dummy_Beam
      procedure :: act_BeamRad => dummy_BeamRad
      !======= end of dummy overloading =======
  end type Propagator
  
  interface Propagator
    procedure Propagator_constructor
  end interface

contains

function Propagator_constructor(L,nStep)
  class(Propagator),pointer :: Propagator_constructor
  integer,intent(in) :: nStep
  real*8, intent(in) :: L
  allocate(Propagator_constructor)
  Propagator_constructor%L = L
  Propagator_constructor%nStep = nStep
  Propagator_constructor%dz = L/dble(nStep)
  Propagator_constructor%flagBeam = .false.
  Propagator_constructor%flagRad = .true.
end function

subroutine propagate(self,Rad,dz)
  implicit none
  class(Propagator),intent(in) :: self
  type(Radiation),intent(inout):: Rad
  real*8, optional,intent(in) :: dz
  integer :: i,nStep
  real*8 :: dz1
  
  if(present(dz)) then
    nStep = 1
    dz1 = dz
  else
    nStep = self%nStep
    dz1 = self%dz
  endif
  Rad%Src = (0d0, 0d0)
  do i=1,nStep
    call Rad%FldSolve(dz1)
  enddo
  call Rad%Slip(dz1*nStep)
  
end subroutine propagate

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMM dummy deffered subroutines MMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! subroutine dummy_get_M(self,M,gamma)
! ! get transverse transfer matrix
  ! implicit none
  ! class(Propagator), intent(in) :: self
  ! real(8), intent(in) :: gamma
  ! real(8), intent(out):: M(4,4)
  ! M = 0d0
  ! do i=1,4
    ! M(i,i) = 1.0
  ! enddo
  ! return
! end subroutine dummy_get_M

subroutine dummy_Beam(self,beam,dz)
  implicit none
  class(Propagator),intent(in) :: self
  type(eBeam),intent(inout)    :: beam
  real*8, optional,intent(in)  :: dz
end subroutine dummy_Beam

subroutine dummy_BeamRad(self,beam,rad,dz)
  implicit none
  class(Propagator),intent(in) :: self
  type(eBeam),intent(inout)    :: beam
  type(Radiation),intent(inout):: rad
  real*8, optional,intent(in)  :: dz
end subroutine dummy_BeamRad

end module PropagatorModule

