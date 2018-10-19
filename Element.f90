module ElementModule
  use RadiationModule
  use eBeamModule
  implicit none
  public
  type, abstract :: Element
    integer:: nStep
    real*8 :: L,dz
    logical:: flagBeam, flagRad
    contains
      procedure :: act
      procedure(virtual_act_Beam), deferred :: act_Beam
      procedure(virtual_act_Rad),  deferred :: act_Rad
      procedure(virtual_act_BeamRad),  deferred :: act_BeamRad
      !procedure(virtual_get_M), deferred :: get_M
  end type Element
  
  type arrElem
    class(Element), pointer :: op
  end type arrElem
  
  private :: virtual_act_Beam,virtual_act_Rad,virtual_act_BeamRad,act!virtual_get_M
  
  interface
    subroutine virtual_act_BeamRad(self,beam,rad,dz)
      import Element, eBeam, Radiation
      class(Element),  intent(in)   :: self
      real*8, optional,intent(in)   :: dz
      type(eBeam),    intent(inout) :: beam
      type(Radiation),intent(inout) :: rad
    end subroutine virtual_act_BeamRad
    
    subroutine virtual_act_Beam(self,beam,dz)
      import Element, eBeam
      class(Element),  intent(in)   :: self
      real*8, optional,intent(in)   :: dz
      type(eBeam),    intent(inout) :: beam
    end subroutine virtual_act_Beam
    
    subroutine virtual_act_Rad(self,rad,dz)
      import Element, Radiation
      class(Element),  intent(in)   :: self
      real*8, optional,intent(in)   :: dz
      type(Radiation),intent(inout) :: rad
    end subroutine virtual_act_Rad
    
    ! subroutine virtual_get_M(self,M,gamma)
    ! ! get transverse transfer matrix
      ! import Element
      ! implicit none
      ! class(Element), intent(in) :: self
      ! real(8), intent(in) :: gamma
      ! real(8), intent(out):: M(4,4)
    ! end subroutine virtual_get_M
  end interface

  ! === example use ==================================
  ! type(arrElem) :: lattice(4)
  ! lattice(1) -> quad(L=0.2d0,nStep=2,DB= 7.86435d0)
  ! lattice(2) -> drift(L=0.1d0,nStep=1)
  ! lattice(3) -> quad(L=0.2d0,nStep=2,DB= -7.86435d0)
  ! lattice(4) -> drift(L=0.1d0,nStep=1)
  ! do i=1,4
  !   call lattice(i)%op%act(Beam,Rad)
  ! enddo
  ! ==================================================
contains

  subroutine act(self,beam,rad,dz)
    implicit none
    class(Element),  intent(in)   :: self
    real*8, optional,intent(in)   :: dz
    type(eBeam),intent(inout) :: beam
    type(Radiation),intent(inout) :: rad
    
    if(self%flagBeam .and. self%flagRad) then
      call self%act_BeamRad(beam,rad,dz)
    elseif(self%flagBeam) then
      call self%act_Beam(beam,dz)
    elseif(self%flagRad)  then
      call self%act_Rad(rad,dz)
    endif
  end subroutine
  
  
end module ElementModule

