module QuadModule
  use ElementModule
  use ConstDataModule
  use eBeamModule
  use RadiationModule
  implicit none
  private
  type, public, extends(Element) :: Quad
    private
    real*8 :: B1 ! field gradient [Tesla / meter]
    real*8 :: G  ! B1/mc
    contains
      procedure :: act_Beam => quad_track  ! kick-drift-kick
      !procedure :: get_M => quad_get_M
      procedure :: get_B1, set_B1
      !======== dummy overloadings =========
      procedure :: act_Rad => dummy_Rad
      procedure :: act_BeamRad => dummy_BeamRad
      !====== end of dummy overloading =====
  end type Quad
  
  interface Quad
    procedure Quad_constructor
  end interface
  
contains

function Quad_constructor(L,nStep,B1)
  class(Quad),pointer :: Quad_constructor
  integer,intent(in) :: nStep
  real*8, intent(in) :: L,B1
  allocate(Quad_constructor)
  Quad_constructor%L     = L
  Quad_constructor%nStep = nStep
  Quad_constructor%dz = L/dble(nStep)
  Quad_constructor%B1 = B1
  Quad_constructor%G  = B1*physConst%c/physConst%Me
  Quad_constructor%flagBeam = .true.
  Quad_constructor%flagRad = .false.
end function Quad_constructor

subroutine quad_track(self,Beam,dz)
  implicit none
  class(Quad),    intent(in)   :: self
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
    call quadMap(self,Beam,dz1)
  enddo

end subroutine quad_track


subroutine quadMap(self,beam,dz) ! kick-drift-kick
  class(Quad),intent(in)   :: self
  real*8,     intent(in)   :: dz
  type(eBeam),intent(inout):: beam
  real*8 :: dummy(beam%npt)
  ! x'' = G*x/(beta*gamma)
  beam%pData(:,px_) = beam%pData(:,px_) - self%G*beam%pData(:,x_)*0.5d0*dz
  beam%pData(:,py_) = beam%pData(:,py_) + self%G*beam%pData(:,y_)*0.5d0*dz
  
  dummy = sqrt(beam%pData(:,g_) *beam%pData(:,g_) -1d0 &  
              -beam%pData(:,px_)*beam%pData(:,px_)-beam%pData(:,py_)*beam%pData(:,py_))
  beam%pData(:,x_)=beam%pData(:,x_)+beam%pData(:,px_)/dummy*dz
  beam%pData(:,y_)=beam%pData(:,y_)+beam%pData(:,py_)/dummy*dz
  beam%pData(:,t_)=beam%pData(:,t_)+(beam%ku + beam%ks*(1d0-beam%pData(:,g_)/dummy))*dz
  
  beam%pData(:,px_) = beam%pData(:,px_) - self%G*beam%pData(:,x_)*0.5d0*dz
  beam%pData(:,py_) = beam%pData(:,py_) + self%G*beam%pData(:,y_)*0.5d0*dz
end subroutine quadMap

real*8 function get_B1(self)
  class(Quad) :: self
  get_B1 = self%B1
end function

subroutine set_B1(self,B1)
  class(Quad) :: self
  real*8 :: B1
  self%B1 = B1
  self%G  = B1*physConst%c/physConst%Me
end subroutine

! subroutine quad_get_M(self,M,gamma)
! ! get transverse transfer matrix
  ! implicit none
  ! class(Quad), intent(in) :: self
  ! real(8), intent(in) :: gamma
  ! real(8), intent(out):: M(4,4)
  ! real(8) :: bg, esp, Gx, c,s,ch,sh, Mx(2,2),My(2,2)
  ! bg = sqrt(gamma**2-1d0)
  ! esp = tiny(self%G/bg)
  ! if(abs(self%G/bg)<=esp) then
    ! Mx = reshape([1d0,self%L/bg,0d0,1d0],[2,2],order=[2,1])
    ! My = reshape([1d0,self%L/bg,0d0,1d0],[2,2],order=[2,1])
  ! elseif(self%G/bg > esp) then
    ! Gx = sqrt(self%G/bg)
    ! c = cos(Gx*self%L)
    ! s = sin(Gx*self%L)
    ! ch= cosh(Gx*self%L)
    ! sh= sinh(Gx*self%L)
    ! Mx = reshape([c,s/Gx,-Gx*s,c],[2,2],order=[2,1])
    ! My = reshape([ch,sh/Gx,Gx*sh,ch],[2,2],order=[2,1])
  ! else
    ! Gx = sqrt(-self%G/bg)
    ! c = cos(Gx*self%L)
    ! s = sin(Gx*self%L)
    ! ch= cosh(Gx*self%L)
    ! sh= sinh(Gx*self%L)
    ! Mx = reshape([ch,sh/Gx,Gx*sh,ch],[2,2],order=[2,1])
    ! My = reshape([c,s/Gx,-Gx*s,c],[2,2],order=[2,1])
  ! endif
  ! M = 0d0
  ! M(1:2,1:2) = Mx
  ! M(3:4,3:4) = My
  ! return
! end subroutine quad_get_M

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMM dummy deffered subroutines MMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine dummy_Rad(self,rad,dz)
  implicit none
  class(Quad),intent(in)   :: self
  type(Radiation),intent(inout):: rad
  real*8, optional,intent(in)   :: dz
end subroutine dummy_Rad

subroutine dummy_BeamRad(self,beam,rad,dz)
  implicit none
  class(Quad),intent(in)      :: self
  type(eBeam),intent(inout)    :: beam
  type(Radiation),intent(inout):: rad
  real*8, optional,intent(in)  :: dz
end subroutine dummy_BeamRad

end module QuadModule

