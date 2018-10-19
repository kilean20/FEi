module DipoleModule
  use ElementModule
  use ConstDataModule
  use eBeamModule
  use RadiationModule
  implicit none
  private
  type, public, extends(Element) :: Dipole
    private
    real*8 :: gamma0 ! normalized energy of reference particle whose orbit is on-axis in reference frame defeind by h
    real*8 :: h      ! =1/rho inverse of geometric curvature defining reference orbit
    contains
      procedure :: act_Beam => dipole_track
      !procedure :: get_M => dipole_get_M
      !======== dummy overloadings =========
      procedure :: act_Rad => dummy_Rad
      procedure :: act_BeamRad => dummy_BeamRad
      !====== end of dummy overloading =====
  end type Dipole
  
  interface Dipole
    procedure Dipole_constructor
  end interface
  
contains

function Dipole_constructor(L,nStep,h,gamma0)
  class(Dipole),pointer :: Dipole_constructor
  integer,intent(in) :: nStep
  real*8, intent(in) :: L,h,gamma0
  allocate(Dipole_constructor)
  Dipole_constructor%L     = L
  Dipole_constructor%nStep = nStep
  Dipole_constructor%dz = L/dble(nStep)
  Dipole_constructor%h = h
  Dipole_constructor%gamma0 = gamma0
  Dipole_constructor%flagBeam = .true.
  Dipole_constructor%flagRad = .false.
end function Dipole_constructor

subroutine dipole_track(self,Beam,dz)
  implicit none
  class(Dipole),    intent(in)   :: self
  type(eBeam), intent(inout):: Beam
  real*8, optional,intent(in)   :: dz
  integer :: i,nStep
  real*8 :: dz1
  
  if(present(dz)) then
    nStep = 1
    dz1 = dz
  else
    nStep = self%nStep
    dz1 = self%dz
  endif
  do i=1,nStep
    call dipoleMap(self,Beam,dz1)
  enddo

end subroutine dipole_track


subroutine dipoleMap(self,beam,dz)
  class(Dipole),intent(in)   :: self
  real*8,     intent(in)   :: dz
  type(eBeam),intent(inout):: beam
  real*8 :: dummy(beam%npt), bg0

  dummy = sqrt(beam%pData(:,g_) *beam%pData(:,g_) -1d0 &  
              -beam%pData(:,px_)*beam%pData(:,px_)-beam%pData(:,py_)*beam%pData(:,py_))
  beam%pData(:,x_)=beam%pData(:,x_)+beam%pData(:,px_)/dummy*0.5d0*dz
  beam%pData(:,y_)=beam%pData(:,y_)+beam%pData(:,py_)/dummy*0.5d0*dz
  beam%pData(:,t_)=beam%pData(:,t_)+(beam%ku + beam%ks*(1d0-beam%pData(:,g_)/dummy))*0.5d0*dz
  
  dummy = sqrt(beam%pData(:,g_)*beam%pData(:,g_)-1d0)
  bg0 = sqrt(self%gamma0-1d0)
  beam%pData(:,px_) = beam%pData(:,px_) + (self%h*dz)*(dummy-bg0) &
                                        - (bg0*self%h*self%h*dz)*beam%pData(:,x_)
  beam%pData(:,t_ ) = beam%pData(:,t_ ) + beam%ku*dz &
                                        + beam%ks*dz*(1d0 -self%h*Beam%pData(:,x_)*Beam%pData(:,g_)/dummy)
  
  dummy = sqrt(beam%pData(:,g_) *beam%pData(:,g_) -1d0 &  
              -beam%pData(:,px_)*beam%pData(:,px_)-beam%pData(:,py_)*beam%pData(:,py_))
  beam%pData(:,x_)=beam%pData(:,x_)+beam%pData(:,px_)/dummy*0.5d0*dz
  beam%pData(:,y_)=beam%pData(:,y_)+beam%pData(:,py_)/dummy*0.5d0*dz
  beam%pData(:,t_)=beam%pData(:,t_)+(beam%ku + beam%ks*(1d0-beam%pData(:,g_)/dummy))*0.5d0*dz
end subroutine dipoleMap


! subroutine dipole_get_M(self,M,gamma)
! ! get transverse transfer matrix
  ! implicit none
  ! class(Dipole), intent(in) :: self
  ! real(8), intent(in) :: gamma
  ! real(8), intent(out):: M(4,4)
  ! real(8) :: bg, Gx, c,s,ch,sh, Mx(2,2),My(2,2)
  ! bg = sqrt(gamma**2-1d0)
  ! if(abs(self%G)<=esp) then
    ! Mx = reshape([1d0,self%L/bg,0d0,1d0],[2,2],order=[2,1])
    ! My = reshape([1d0,self%L/bg,0d0,1d0],[2,2],order=[2,1])
  ! endif
  ! M = 0d0
  ! M(1:2,1:2) = Mx
  ! M(3,3) = 1d0
  ! M(4,4) = 1d0
  ! return
! end subroutine dipole_get_M

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMM dummy deffered subroutines MMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine dummy_Rad(self,rad,dz)
  implicit none
  class(Dipole),intent(in)   :: self
  type(Radiation),intent(inout):: rad
  real*8, optional,intent(in)   :: dz
end subroutine dummy_Rad

subroutine dummy_BeamRad(self,beam,rad,dz)
  implicit none
  class(Dipole),intent(in)      :: self
  type(eBeam),intent(inout)    :: beam
  type(Radiation),intent(inout):: rad
  real*8, optional,intent(in)  :: dz
end subroutine dummy_BeamRad

end module DipoleModule

