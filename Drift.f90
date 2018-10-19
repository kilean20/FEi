module DriftModule
  use ElementModule
  use ConstDataModule
  use eBeamModule
  use RadiationModule
  implicit none
  private
  type, public, extends(Element) :: Drift
    contains
      procedure :: act_Beam => drift_track
      !procedure :: get_M => drift_get_M
      !======== dummy overloadings =========
      procedure :: act_Rad => dummy_Rad
      procedure :: act_BeamRad => dummy_BeamRad
      !====== end of dummy overloading =====
  end type Drift
  
  interface Drift
    procedure Drift_constructor
  end interface

contains

function Drift_constructor(L,nStep)
  class(Drift),pointer :: Drift_constructor
  integer,intent(in) :: nStep
  real*8, intent(in) :: L
  allocate(Drift_constructor)
  Drift_constructor%L     = L
  Drift_constructor%nStep = nStep
  Drift_constructor%dz = L/dble(nStep)
  Drift_constructor%flagBeam = .true.
  Drift_constructor%flagRad = .false.
end function

subroutine drift_track(self,Beam,dz)
  implicit none
  class(Drift),intent(in)   :: self
  type(eBeam),intent(inout):: Beam
  real*8, optional,intent(in)   :: dz
  real*8 :: L
  
  if(present(dz)) then
    L = dz
  else
    L = self%L
  endif
  call driftmap(self,Beam,L)
end subroutine drift_track

subroutine driftmap(self,beam,dz)
  class(Drift),intent(in)   :: self
  real*8,      intent(in)   :: dz
  type(eBeam), intent(inout):: beam
  real*8 :: dummy(beam%npt)
  
  dummy = sqrt(beam%pData(:,g_) *beam%pData(:,g_) -1d0 &
              -beam%pData(:,px_)*beam%pData(:,px_)-beam%pData(:,py_)*beam%pData(:,py_))
  beam%pData(:,x_)=beam%pData(:,x_)+beam%pData(:,px_)/dummy*dz
  beam%pData(:,y_)=beam%pData(:,y_)+beam%pData(:,py_)/dummy*dz
  beam%pData(:,t_)=beam%pData(:,t_)+(beam%ku+beam%ks*(1d0-beam%pData(:,g_)/dummy))*dz
end subroutine driftmap

! subroutine drift_get_M(self,M,gamma)
! ! get transverse transfer matrix
  ! implicit none
  ! class(Drift), intent(in) :: self
  ! real(8), intent(in) :: gamma
  ! real(8), intent(out):: M(4,4)
  ! real(8) :: bg, Mx(2,2)
  ! bg = sqrt(gamma**2-1d0)
  ! Mx = reshape([1d0,self%L/bg,0d0,1d0],[2,2],order=[2,1])
  ! M = 0d0
  ! M(1:2,1:2) = Mx
  ! M(3:4,3:4) = Mx
  ! return
! end subroutine drift_get_M

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMM dummy deffered subroutines MMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine dummy_Rad(self,rad,dz)
  implicit none
  class(Drift),intent(in)      :: self
  type(Radiation),intent(inout):: rad
  real*8,optional,intent(in)   :: dz
end subroutine dummy_Rad

subroutine dummy_BeamRad(self,beam,rad,dz)
  implicit none
  class(Drift),intent(in)      :: self
  type(eBeam),intent(inout)    :: beam
  type(Radiation),intent(inout):: rad
  real*8,optional,intent(in)   :: dz
end subroutine dummy_BeamRad

end module DriftModule

