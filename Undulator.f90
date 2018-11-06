module undulatorModule
  use ElementModule
  use ConstDataModule
  use eBeamModule
  use RadiationModule
  implicit none
  private 
  type, public, extends(Element) :: Undulator
    integer :: nPeriod, nHarm
    integer,allocatable :: harm_tab(:)
    real*8 :: K,kx,ky,ku
    real*8 :: lu
    real*8 :: K2,kx2,ky2,ku2
    ! integration parameters
    real*8, allocatable :: I0(:,:),IC(:,:)
    integer :: HeffID, RKorder
    ! HeffID = 0  -fastest, effective Hamiltonian includes the leading order terms
    !        = 1  -slower,  effective Hamiltonian includes the next leading order terms
    ! RKorder = 4(default), Runge-Kutta order 
    contains
      procedure :: destruct => undulator_destructor
      procedure :: xyMap => undulator_xyMap
      procedure :: tMap0 => undulator_tMap0
      procedure :: tMap1 => undulator_tMap1
      procedure :: act_BeamRad => track_undulator
      !procedure :: track_steady => undulator_track_steady
      procedure :: track_leapfrog => undulator_track_leapfrog
      !procedure :: track_steady_leapfrog => undulator_track_steady_leapfrog
      procedure :: track_schmitt => undulator_track_schmitt
      procedure :: track_slipsplit => undulator_track_slipsplit
      !procedure :: get_M => undulator_get_M
      !======== dummy overloadings =========
      procedure :: act_Rad => dummy_rad
      procedure :: act_Beam => dummy_beam
      !====== end of dummy overloading =====
  end type Undulator

  interface Undulator
    procedure undulator_constructor
  end interface


contains

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMM   Setting   MMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
function undulator_constructor(nPeriod,nStep,K,kx,ky,ku,harm_tab,nHarm,HeffID,RKorder) result(self)
  implicit none
  class(Undulator), pointer :: self
  integer, intent(in) :: nStep, nPeriod, nHarm, harm_tab(nHarm)
  real*8 , intent(in) :: K,kx,ky,ku
  integer, optional, intent(in) :: HeffID,RKorder
  
  allocate(self)
  self%nStep = nStep
  self%nPeriod = nPeriod
  self%K = K
  self%kx = kx
  self%ky = ky
  self%ku = ku
  self%kx2 = kx*kx
  self%ky2 = ky*ky
  self%ku2 = ku*ku
  self%lu = twopi/ku
  self%L = self%lu*dble(nPeriod)
  self%dz = self%L/dble(nStep)
  self%flagBeam = .true.
  self%flagRad  = .true.
  
  if(present(HeffID)) then
    self%HeffID = HeffID
  else
    self%HeffID = 0
  endif
  
  if(present(RKorder)) then
    self%RKorder = RKorder
  else
    self%RKorder = 4
  endif
  
  self%nHarm = nHarm
  allocate(self%harm_tab(nHarm))
  self%harm_tab = harm_tab
  
  if(self%HeffID==0) then
    call calculate_IntParam0(self)
  elseif(self%HeffID==1) then
    call calculate_IntParam1(self)
  endif
end function undulator_constructor

subroutine undulator_destructor(self)
  class(Undulator), intent(inout) :: self
  if(allocated(self%I0))  deallocate(self%I0)
  if(allocated(self%IC))  deallocate(self%IC)
  if(allocated(self%harm_tab))  deallocate(self%harm_tab)
end subroutine

function calculate_J_tab(xi,harm_tab,nHarm,maxJ) result(J_tab)
  implicit none
  integer, intent(in) :: nHarm, maxJ
  integer, dimension(nHarm), intent(in) :: harm_tab
  real*8, intent(in) :: xi
  
  real*8, dimension(-maxJ:maxJ,nHarm) :: J_tab
  integer :: n,h,l
  
  do n=1,nHarm
    h=harm_tab(n)
    J_tab(0:maxJ,n) = BESSEL_JN(0, maxJ, h*xi)
  enddo      
  do l=1,maxJ,2
    J_tab(-l,:) = -J_tab(l,:)
  enddo
  do l=2,maxJ,2
    J_tab(-l,:) = J_tab(l,:)
  enddo
end function calculate_J_tab

subroutine calculate_IntParam0(self)
  implicit none
  type(undulator), intent(inout) :: self
  double precision :: xi,lu
  double precision, dimension(-self%nHarm/2-1:self%nHarm/2+1,self%nHarm) :: J_tab
  double precision :: ku
  integer :: n,h

  lu = self%lu
  xi = self%K*self%K/(4d0+2d0*self%K*self%K)
  ku = self%ku
  J_tab = calculate_J_tab(xi,self%harm_tab,self%nHarm,self%nHarm/2+1) 

  allocate(self%IC(1,self%nHarm))
  allocate(self%I0(1,self%nHarm))
  do n=1,self%nHarm
    h=self%harm_tab(n)
    if(mod(h,2)==0) then
      self%IC(1,n) = 0.5d0*(J_tab(-(h+2)/2,n)-J_tab(-(h-2)/2,n))
      self%I0(1,n) = J_tab(-h/2,n)
    else
      self%IC(1,n) = 0.5d0*(J_tab(-(h-1)/2,n)+J_tab(-(h+1)/2,n))
      self%I0(1,n) = 0d0
    endif
  enddo
end subroutine calculate_IntParam0

subroutine calculate_IntParam1(self)
  !==================================================
  !  pre-calculate integration parameter, or 
  !  so called 'coupling factor' in focus of
  !  off resonant condition. (energy jitter, seed freq jitter, etc)
  !--------------------------------------------------
  implicit none
  type(undulator), intent(inout) :: self
  integer, parameter :: maxJ=10
  real*8 :: xi
  real*8,  dimension(-maxJ:maxJ,self%nHarm) :: J_tab
  integer :: i,n,h
  
  xi = self%K*self%K/(4d0+2d0*self%K*self%K)
  J_tab = calculate_J_tab(xi,self%harm_tab,self%nHarm,maxJ)
      
  allocate(self%IC(3,self%nHarm))
  allocate(self%I0(1,self%nHarm))
  self%IC = 0d0
  self%I0 = 0d0
  do n=1,self%nHarm
    h=self%harm_tab(n)
    if(mod(h,2)==0) then
      self%IC(1,n) = 0.5d0*(J_tab(-(h+2)/2,n)-J_tab(-(h-2)/2,n))
      self%I0(1,n) = J_tab(-h/2,n)
    else
      self%IC(1,n) = 0.5d0*(J_tab(-(h-1)/2,n)+J_tab(-(h+1)/2,n))
      if(h==1) then
        ! IC(2,n)
        do i=-maxJ,-2
          self%IC(2,n) = self%IC(2,n) + J_tab(i,n)*(2*i+1)/dble(i*i+i)
        enddo
        do i=1,maxJ
          self%IC(2,n) = self%IC(2,n) + J_tab(i,n)*(2*i+1)/dble(i*i+i)
        enddo
        ! case i=-1 and 0
        self%IC(2,n) = self%IC(2,n) + J_tab(0,n) - J_tab(-1,n)
        self%IC(2,n) = 0.5d0*self%IC(2,n)
        ! IC(3,n)
        self%IC(3,n) = J_tab(1,n)*(1d0-1d0/xi) + J_tab(0,n)
      endif
    endif
  enddo
end subroutine calculate_IntParam1


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMM   Track   MMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine undulator_xyMap(self,Beam,dz)
!==============================================================
!Transverse map under undulator field
!==============================================================
  implicit none
  class(Undulator), intent(in) :: self
  type(eBeam) :: Beam
  real*8, intent(in) :: dz  ! half of the step size
  real*8, dimension(Beam%npt) :: iG

  iG = 1d0/Beam%pData(:,g_)
  !==== half drift ===============================================
  Beam%pData(:,x_) = Beam%pData(:,x_) + dz*Beam%pData(:,px_)*iG 
  Beam%pData(:,y_) = Beam%pData(:,y_) + dz*Beam%pData(:,py_)*iG
  !===============================================================
  !==== kick ====================================================
  Beam%pData(:,px_) = Beam%pData(:,px_) -dz*self%kx2*self%K*self%K*Beam%pData(:,x_)*iG
  Beam%pData(:,py_) = Beam%pData(:,py_) -dz*self%ky2*self%K*self%K*Beam%pData(:,y_)*iG
  !===============================================================
  !==== half drift ===============================================
  Beam%pData(:,x_) = Beam%pData(:,x_) + dz*Beam%pData(:,px_)*iG 
  Beam%pData(:,y_) = Beam%pData(:,y_) + dz*Beam%pData(:,py_)*iG
  !===============================================================
end subroutine undulator_xyMap

subroutine undulator_tMap0(self,Beam,Rad,dz)
  implicit none
  class(Undulator),intent(in) :: self
  type(Radiation), intent(in) :: Rad
  real*8,          intent(in) :: dz
  type(eBeam),     intent(inout) :: Beam
  complex(8) :: fld(Beam%npt/Beam%nBin,self%nHarm)
  
  if(Beam%npt==0) return
  call Beam%scatter(Rad,fld)
  
  select case(self%RKorder)
    case (1)
      call undulator_tMap0_RK1(self,Beam,fld,dz)
    case (2)
      call undulator_tMap0_RK2(self,Beam,fld,dz)
    case default
      call undulator_tMap0_RK4(self,Beam,fld,dz)
  end select
  
end subroutine undulator_tMap0

subroutine undulator_tMap0_RK4(self,Beam,fld,dz)
!=================================================================
! Longitudinal map under undulator and radiation field
!-----------------------------------------------------------------
  implicit none
  type(undulator),intent(in) :: self
  type(eBeam),    intent(inout) :: Beam
  real*8,         intent(in) :: dz
  complex(8),     intent(in) :: fld(Beam%npt/Beam%nBin,self%nHarm)
  
  integer :: h,n,i
  real*8  :: gR,zetaR
  real*8,    dimension(Beam%npt)  :: Keff,th0,th1,delta
  complex(8),dimension(Beam%npt)  :: B
  real*8,    dimension(Beam%npt)  :: Veff
  real*8,    dimension(Beam%npt,2):: f1,f2,tmp
  
  associate(x =>Beam%pData(:,x_), px =>Beam%pData(:,px_),&
            y =>Beam%pData(:,y_), py =>Beam%pData(:,py_),&
            gamma=>Beam%pData(:,g_),&
            K=>self%K,  ks=>Beam%ks)
  
    gR = sqrt(0.5d0*ks/self%ku*(1d0+0.5d0*K*K))
    zetaR = 2d0*K/(1d0+0.5d0*K*K)
    tmp(:,1) = self%kx2*x*x +self%ky2*y*y
    tmp(:,2) = px*px+py*py

    th0   = -(0.5d0*ks/(gR*gR))*( tmp(:,2) +0.5d0*K*K*tmp(:,1)*(1d0+0.25d0*tmp(:,1)))
    Keff  =  K*(1d0+0.5d0*tmp(:,1))
    th1   =  (0.5d0*ks/(gR*gR))*(1d0+tmp(:,2)+0.5d0*Keff*Keff)
    !===============================================================
    !==== Longitudinal RK4 =========================================
    !f1-------------------------------------------------------------
    delta = gamma/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))*gamma
    f1(:,2) = 0d0
    delta = 1d0+delta
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -((self%IC(1,n)*0.5d0*h*ks)*(Keff*zetaR/(delta*delta)) &
                                 +self%I0(1,n))*px
        f1(:,1) = f1(:,1) + real(B)*((self%IC(1,n)*ks*h*zetaR/gR)*px/delta**3)
      else
        Veff = -(self%IC(1,n)*ks)*Keff
      endif
      Veff = Veff/gamma
      f1(:,1) = f1(:,1) -  Veff*real(B)
      f1(:,2) = f1(:,2) +h*Veff*aimag(B)
    enddo
    f1(:,1) = f1(:,1)/gamma
    !f2-------------------------------------------------------------
    tmp = Beam%pData(:,t_:g_) + f1*0.5d0*dz
    delta = tmp(:,2)/gR - 1d0
    f2(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))*tmp(:,2)
    f2(:,2) = 0d0
    delta = 1d0+delta
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -((self%IC(1,n)*0.5d0*h*ks)*(Keff*zetaR/(delta*delta)) &
                                 +self%I0(1,n))*px
        f2(:,1) = f2(:,1) + real(B)*((self%IC(1,n)*ks*h*zetaR/gR)*px/delta**3)
      else
        Veff = -(self%IC(1,n)*ks)*Keff
      endif
      Veff = Veff/tmp(:,2)
      f2(:,1) = f2(:,1) -  Veff*real(B)
      f2(:,2) = f2(:,2) +h*Veff*aimag(B)
    enddo
    f2(:,1) = f2(:,1)/tmp(:,2)
    !k3-------------------------------------------------------------
    tmp = Beam%pData(:,t_:g_) + f2*0.5d0*dz
    f2 = f1+2*f2
    delta = tmp(:,2)/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))*tmp(:,2)
    f1(:,2) = 0d0
    delta = 1d0+delta
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -((self%IC(1,n)*0.5d0*h*ks)*(Keff*zetaR/(delta*delta)) &
                                 +self%I0(1,n))*px
        f1(:,1) = f1(:,1) + real(B)*((self%IC(1,n)*ks*h*zetaR/gR)*px/delta**3)
      else
        Veff = -(self%IC(1,n)*ks)*Keff
      endif
      Veff = Veff/tmp(:,2)
      f1(:,1) = f1(:,1) -  Veff*real(B)
      f1(:,2) = f1(:,2) +h*Veff*aimag(B)
    enddo
    f1(:,1) = f1(:,1)/tmp(:,2)
    !k4-------------------------------------------------------------
    tmp = Beam%pData(:,t_:g_) + f1*dz
    f2 = f2+2*f1
    delta = tmp(:,2)/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))*tmp(:,2)
    f1(:,2) = 0d0
    delta = 1d0+delta
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -((self%IC(1,n)*0.5d0*h*ks)*(Keff*zetaR/(delta*delta)) &
                                 +self%I0(1,n))*px
        f1(:,1) = f1(:,1) + real(B)*((self%IC(1,n)*ks*h*zetaR/gR)*px/delta**3)
      else
        Veff = -(self%IC(1,n)*ks)*Keff
      endif
      Veff = Veff/tmp(:,2)
      f1(:,1) = f1(:,1) -  Veff*real(B)
      f1(:,2) = f1(:,2) +h*Veff*aimag(B)
    enddo
    f1(:,1) = f1(:,1)/tmp(:,2)
    !sum(f1,f2,k3,k4)-----------------------------------------------
    Beam%pData(:,t_:g_) = Beam%pData(:,t_:g_) + dz*(f1 + f2)/6d0
    !===============================================================
  end associate
end subroutine undulator_tMap0_RK4

subroutine undulator_tMap0_RK2(self,Beam,fld,dz)
!=================================================================
! Longitudinal map under undulator and radiation field
!-----------------------------------------------------------------
  implicit none
  type(undulator),intent(in) :: self
  type(eBeam),    intent(inout) :: Beam
  real*8,         intent(in) :: dz
  complex(8),     intent(in) :: fld(Beam%npt/Beam%nBin,self%nHarm)
  
  integer :: h,n,i
  real*8  :: gR,zetaR
  real*8,    dimension(Beam%npt)  :: Keff,th0,th1,delta
  complex(8),dimension(Beam%npt)  :: B
  real*8,    dimension(Beam%npt)  :: Veff
  real*8,    dimension(Beam%npt,2):: f1,tmp  
  
  associate(x =>Beam%pData(:,x_), px =>Beam%pData(:,px_),&
            y =>Beam%pData(:,y_), py =>Beam%pData(:,py_),&
            gamma=>Beam%pData(:,g_),&
            K=>self%K,  ks=>Beam%ks)
  
    gR = sqrt(0.5d0*ks/self%ku*(1d0+0.5d0*K*K))
    zetaR = 2d0*K/(1d0+0.5d0*K*K)
    tmp(:,1) = self%kx2*x*x +self%ky2*y*y
    tmp(:,2) = px*px+py*py

    th0   = -(0.5d0*ks/(gR*gR))*( tmp(:,2) +0.5d0*K*K*tmp(:,1)*(1d0+0.25d0*tmp(:,1)))
    Keff  =  K*(1d0+0.5d0*tmp(:,1))
    th1   =  (0.5d0*ks/(gR*gR))*(1d0+tmp(:,2)+0.5d0*Keff*Keff)
    !===============================================================
    !==== Longitudinal RK2 =========================================
    !f1-------------------------------------------------------------
    delta = gamma/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))*gamma
    f1(:,2) = 0d0
    delta = 1d0+delta
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -((self%IC(1,n)*0.5d0*h*ks)*(Keff*zetaR/(delta*delta)) &
                                 +self%I0(1,n))*px
        f1(:,1) = f1(:,1) + real(B)*((self%IC(1,n)*ks*h*zetaR/gR)*px/delta**3)
      else
        Veff = -(self%IC(1,n)*ks)*Keff
      endif
      Veff = Veff/gamma
      f1(:,1) = f1(:,1) -  Veff*real(B)
      f1(:,2) = f1(:,2) +h*Veff*aimag(B)
    enddo
    f1(:,1) = f1(:,1)/gamma
    !f2-------------------------------------------------------------
    tmp = Beam%pData(:,t_:g_) + f1*dz
    Beam%pData(:,t_:g_) = Beam%pData(:,t_:g_) + f1*0.5d0*dz
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))*tmp(:,2)
    f1(:,2) = 0d0
    delta = 1d0+delta
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -((self%IC(1,n)*0.5d0*h*ks)*(Keff*zetaR/(delta*delta)) &
                                 +self%I0(1,n))*px
        f1(:,1) = f1(:,1) + real(B)*((self%IC(1,n)*ks*h*zetaR/gR)*px/delta**3)
      else
        Veff = -(self%IC(1,n)*ks)*Keff
      endif
      Veff = Veff/tmp(:,2)
      f1(:,1) = f1(:,1) -  Veff*real(B)
      f1(:,2) = f1(:,2) +h*Veff*aimag(B)
    enddo
    f1(:,1) = f1(:,1)/tmp(:,2)
    !sum(f1,f2,k3,k4)-----------------------------------------------
    Beam%pData(:,t_:g_) = Beam%pData(:,t_:g_) + 0.5d0*dz*f1
    !===============================================================
  end associate
end subroutine undulator_tMap0_RK2

subroutine undulator_tMap0_RK1(self,Beam,fld,dz)
!=================================================================
! Longitudinal map under undulator and radiation field
!-----------------------------------------------------------------
  implicit none
  type(undulator),intent(in) :: self
  type(eBeam),    intent(inout) :: Beam
  real*8,         intent(in) :: dz
  complex(8),     intent(in) :: fld(Beam%npt/Beam%nBin,self%nHarm)
  
  integer :: h,n,i
  real*8  :: gR,zetaR
  real*8,    dimension(Beam%npt)  :: Keff,th0,th1,delta
  complex(8),dimension(Beam%npt)  :: B
  real*8,    dimension(Beam%npt)  :: Veff
  real*8,    dimension(Beam%npt,2):: tmp
 
  associate(x =>Beam%pData(:,x_), px =>Beam%pData(:,px_),&
            y =>Beam%pData(:,y_), py =>Beam%pData(:,py_),&
            gamma=>Beam%pData(:,g_),&
            K=>self%K,  ks=>Beam%ks)
  
    gR = sqrt(0.5d0*ks/self%ku*(1d0+0.5d0*K*K))
    zetaR = 2d0*K/(1d0+0.5d0*K*K)
    tmp(:,1) = self%kx2*x*x +self%ky2*y*y
    tmp(:,2) = px*px+py*py

    th0   = -(0.5d0*ks/(gR*gR))*( tmp(:,2) +0.5d0*K*K*tmp(:,1)*(1d0+0.25d0*tmp(:,1)))
    Keff  =  K*(1d0+0.5d0*tmp(:,1))
    th1   =  (0.5d0*ks/(gR*gR))*(1d0+tmp(:,2)+0.5d0*Keff*Keff)
    !===============================================================
    !==== Longitudinal RK1 =========================================
    !f1-------------------------------------------------------------
    delta = gamma/gR - 1d0
    tmp(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))*gamma
    tmp(:,2) = 0d0
    delta = 1d0+delta
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -((self%IC(1,n)*0.5d0*h*ks)*(Keff*zetaR/(delta*delta)) &
                                 +self%I0(1,n))*px
        tmp(:,1) = tmp(:,1) + real(B)*((self%IC(1,n)*ks*h*zetaR/gR)*px/delta**3)
      else
        Veff = -(self%IC(1,n)*ks)*Keff
      endif
      Veff = Veff/gamma
      tmp(:,1) = tmp(:,1) -  Veff*real(B)
      tmp(:,2) = tmp(:,2) +h*Veff*aimag(B)
    enddo
    tmp(:,1) = tmp(:,1)/gamma
    !sum(f1,f2,k3,k4)-----------------------------------------------
    Beam%pData(:,t_:g_) = Beam%pData(:,t_:g_) + dz*tmp
    !===============================================================
  end associate
end subroutine undulator_tMap0_RK1


subroutine undulator_tMap1(self,Beam,Rad,fldTmp,dz)
  implicit none
  class(Undulator),intent(in) :: self
  type(Radiation), intent(in) :: Rad
  real*8,          intent(in) :: dz
  type(eBeam),     intent(inout) :: Beam
  complex(8),allocatable, intent(inout) :: fldTmp(:,:)
  complex(8),dimension(Beam%npt/Beam%nBin,self%nHarm) :: fld,fldz

  if(Beam%npt==0) return
  call Beam%scatter(Rad,fld)
  fldz = (fld-fldTmp)/dz
  deallocate(fldTmp)
  select case(self%RKorder)
    case (1)
      call undulator_tMap1_RK1(self,Beam,fld,fldz,dz)
    case (2)
      call undulator_tMap1_RK2(self,Beam,fld,fldz,dz)
    case default
      call undulator_tMap1_RK4(self,Beam,fld,fldz,dz)
  end select
  
end subroutine undulator_tMap1

subroutine undulator_tMap1_RK4(self,Beam,fld,fldz,dz)
!=================================================================
! Longitudinal map under undulator and radiation field
!-----------------------------------------------------------------
  implicit none
  type(undulator),intent(in) :: self
  type(eBeam),    intent(inout) :: Beam
  real*8,         intent(in) :: dz
  complex(8),     intent(in),&
    dimension(Beam%npt/Beam%nBin,self%nHarm) :: fld, fldz

  integer :: h,n,i
  real*8  :: gR,xiR,zetaR
  real*8,    dimension(Beam%npt)  :: Keff,th0,th1,xi1,delta
  complex(8)                      :: IIc
  complex(8),dimension(Beam%npt)  :: Veff,B,II_th0
  real*8,    dimension(Beam%npt,2):: f1,f2,tmp
  
  associate(x =>Beam%pData(:,x_), px =>Beam%pData(:,px_),&
            y =>Beam%pData(:,y_), py =>Beam%pData(:,py_),&
            gamma=>Beam%pData(:,g_),  ks =>Beam%ks ,&
            K=>self%K, kx2=>self%kx2, ky2=>self%ky2, ku=>self%ku)
  
    gR = sqrt(0.5d0*ks/ku*(1d0+0.5d0*K*K))
    xiR = 0.125d0*ks*K*K/(ku*gR*gR)
    zetaR = ks*K/(ku*gR*gR)
    tmp(:,1) = Kx2*x*x + Ky2*y*y
    tmp(:,2) = px*px+py*py
    th0 = -(0.5d0*ks/(gR*gR))*(tmp(:,2) +0.5d0*K*K*tmp(:,1)*(1d0+0.25d0*tmp(:,1)))
    IIC = (0d0, 0d0)
    do n=1,self%nHarm
      h=self%harm_tab(n)
      if(h==1) then
        !IIc  = (i1*twopi*self%IC(1,n) - self%IC(2,n))/ku
        IIC    = -self%IC(2,n)/ku
        II_th0 = self%IC(1,n)+0.5d0*self%IC(3,n)*xiR +0.5d0*IIC*th0
        exit
      endif
    enddo
    Keff =  K*(1d0+0.5d0*tmp(:,1))
    th1  =  (0.5d0*ks/(gR*gR))*(1d0+tmp(:,2)+0.5d0*Keff*Keff)
    xi1  =  xiR*(1d0+tmp(:,1))
    !===============================================================
    !==== Longitudinal RK4 =========================================
    !f1-------------------------------------------------------------
    delta = gamma/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))/ks*gamma
    f1(:,2) = 0d0
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -B*(((self%IC(1,n)*0.5d0*h)*(Keff*zetaR/((1d0+delta)**2)) &
                             +self%I0(1,n))*px/gamma)
        f1(:,1) = f1(:,1) + ((self%IC(1,n)*h*zetaR/gR)*px/(1d0+delta)**3)*real(B)
      else
        if(h==1) then
          Veff = (II_th0 + (IIc*(th1*(2d0+delta)*delta)-self%IC(3,n)*xi1)/(2d0*(1d0+delta)**2))*B
          f1(:,1) = f1(:,1) - real((th1*IIc+xi1*self%IC(3,n))*Keff/gR/(1d0+delta)**3*B)
          forall(i=1:Beam%npt/Beam%nBin)
            B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fldz(i,n)
          end forall
          Veff = -Keff/gamma*(Veff-0.5d0*i1*IIc*B)
        else
          Veff = -(self%IC(1,n)*Keff/gamma)*B
        endif
      endif
      f1(:,1) = f1(:,1)   - real(Veff)
      f1(:,2) = f1(:,2) +h*aimag(Veff)
    enddo
    f1(:,1) = f1(:,1)/gamma
    f1 = f1*ks
    !f2-------------------------------------------------------------
    tmp = Beam%pData(:,t_:g_) + f1*0.5d0*dz
    delta = tmp(:,2)/gR - 1d0
    f2(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))/ks*tmp(:,2)
    f2(:,2) = 0d0
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -B*(((self%IC(1,n)*0.5d0*h)*(Keff*zetaR/((1d0+delta)**2)) &
                             +self%I0(1,n))*px/tmp(:,2))
        f2(:,1) = f2(:,1) + ((self%IC(1,n)*h*zetaR/gR)*px/(1d0+delta)**3)*real(B)
      else
        if(h==1) then
          Veff = (II_th0 + (IIc*(th1*(2d0+delta)*delta)-self%IC(3,n)*xi1)/(2d0*(1d0+delta)**2))*B
          f2(:,1) = f2(:,1) - real((th1*IIc+xi1*self%IC(3,n))*Keff/gR/(1d0+delta)**3*B)
          forall(i=1:Beam%npt/Beam%nBin)
            B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fldz(i,n)
          end forall
          Veff = -Keff/tmp(:,2)*(Veff-0.5d0*i1*IIc*B)
        else
          Veff = -(self%IC(1,n)*Keff/tmp(:,2))*B
        endif
      endif
      f2(:,1) = f2(:,1)   - real(Veff)
      f2(:,2) = f2(:,2) +h*aimag(Veff)
    enddo
    f2(:,1) = f2(:,1)/tmp(:,2)
    f2 = f2*ks
    !k3-------------------------------------------------------------
    tmp = Beam%pData(:,t_:g_) + f2*0.5d0*dz
    f2 = f1+2*f2
    delta = tmp(:,2)/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))/ks*tmp(:,2)
    f1(:,2) = 0d0
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -B*(((self%IC(1,n)*0.5d0*h)*(Keff*zetaR/((1d0+delta)**2)) &
                             +self%I0(1,n))*px/tmp(:,2))
        f1(:,1) = f1(:,1) + ((self%IC(1,n)*h*zetaR/gR)*px/(1d0+delta)**3)*real(B)
      else
        if(h==1) then
          Veff = (II_th0 + (IIc*(th1*(2d0+delta)*delta)-self%IC(3,n)*xi1)/(2d0*(1d0+delta)**2))*B
          f1(:,1) = f1(:,1) - real((th1*IIc+xi1*self%IC(3,n))*Keff/gR/(1d0+delta)**3*B)
          forall(i=1:Beam%npt/Beam%nBin)
            B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fldz(i,n)
          end forall
          Veff = -Keff/tmp(:,2)*(Veff-0.5d0*i1*IIc*B)
        else
          Veff = -(self%IC(1,n)*Keff/tmp(:,2))*B
        endif
      endif
      f1(:,1) = f1(:,1)   - real(Veff)
      f1(:,2) = f1(:,2) +h*aimag(Veff)
    enddo
    f1(:,1) = f1(:,1)/tmp(:,2)
    f1 = f1*ks
    !k4-------------------------------------------------------------
    tmp = Beam%pData(:,t_:g_) + f1*dz
    f2 = f2+2*f1
    delta = tmp(:,2)/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))/ks*tmp(:,2)
    f1(:,2) = 0d0
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -B*(((self%IC(1,n)*0.5d0*h)*(Keff*zetaR/((1d0+delta)**2)) &
                             +self%I0(1,n))*px/tmp(:,2))
        f1(:,1) = f1(:,1) + ((self%IC(1,n)*h*zetaR/gR)*px/(1d0+delta)**3)*real(B)
      else
        if(h==1) then
          Veff = (II_th0 + (IIc*(th1*(2d0+delta)*delta)-self%IC(3,n)*xi1)/(2d0*(1d0+delta)**2))*B
          f1(:,1) = f1(:,1) - real((th1*IIc+xi1*self%IC(3,n))*Keff/gR/(1d0+delta)**3*B)
          forall(i=1:Beam%npt/Beam%nBin)
            B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fldz(i,n)
          end forall
          Veff = -Keff/tmp(:,2)*(Veff-0.5d0*i1*IIc*B)
        else
          Veff = -(self%IC(1,n)*Keff/tmp(:,2))*B
        endif
      endif
      f1(:,1) = f1(:,1)   - real(Veff)
      f1(:,2) = f1(:,2) +h*aimag(Veff)
    enddo
    f1(:,1) = f1(:,1)/tmp(:,2)
    f1 = f1*ks
    !sum(f1,f2,k3,k4)-----------------------------------------------
    Beam%pData(:,t_:g_) = Beam%pData(:,t_:g_) + dz*(f1 + f2)/6d0
    !===============================================================
  end associate
end subroutine undulator_tMap1_RK4

subroutine undulator_tMap1_RK2(self,Beam,fld,fldz,dz)
!=================================================================
! Longitudinal map under undulator and radiation field
!-----------------------------------------------------------------
  implicit none
  type(undulator),intent(in) :: self
  type(eBeam),    intent(inout) :: Beam
  real*8,         intent(in) :: dz
  complex(8),     intent(in),&
    dimension(Beam%npt/Beam%nBin,self%nHarm) :: fld, fldz

  integer :: h,n,i
  real*8  :: gR,xiR,zetaR
  real*8,    dimension(Beam%npt)  :: Keff,th0,th1,xi1,delta
  complex(8)                      :: IIc
  complex(8),dimension(Beam%npt)  :: Veff,B,II_th0
  real*8,    dimension(Beam%npt,2):: f1,tmp
  
  associate(x =>Beam%pData(:,x_), px =>Beam%pData(:,px_),&
            y =>Beam%pData(:,y_), py =>Beam%pData(:,py_),&
            gamma=>Beam%pData(:,g_),  ks =>Beam%ks ,&
            K=>self%K, kx2=>self%kx2, ky2=>self%ky2, ku=>self%ku)
  
    gR = sqrt(0.5d0*ks/ku*(1d0+0.5d0*K*K))
    xiR = 0.125d0*ks*K*K/(ku*gR*gR)
    zetaR = ks*K/(ku*gR*gR)
    tmp(:,1) = Kx2*x*x + Ky2*y*y
    tmp(:,2) = px*px+py*py
    th0 = -(0.5d0*ks/(gR*gR))*(tmp(:,2) +0.5d0*K*K*tmp(:,1)*(1d0+0.25d0*tmp(:,1)))
    IIC = (0d0, 0d0)
    do n=1,self%nHarm
      h=self%harm_tab(n)
      if(h==1) then
        IIC    = -self%IC(2,n)/ku
        II_th0 = self%IC(1,n)+0.5d0*self%IC(3,n)*xiR +0.5d0*IIC*th0
        exit
      endif
    enddo
    Keff =  K*(1d0+0.5d0*tmp(:,1))
    th1  =  (0.5d0*ks/(gR*gR))*(1d0+tmp(:,2)+0.5d0*Keff*Keff)
    xi1  =  xiR*(1d0+tmp(:,1))
    !===============================================================
    !==== Longitudinal RK4 =========================================
    !f1-------------------------------------------------------------
    delta = gamma/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))/ks*gamma
    f1(:,2) = 0d0
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -B*(((self%IC(1,n)*0.5d0*h)*(Keff*zetaR/((1d0+delta)**2)) &
                             +self%I0(1,n))*px/gamma)
        f1(:,1) = f1(:,1) + ((self%IC(1,n)*h*zetaR/gR)*px/(1d0+delta)**3)*real(B)
      else
        if(h==1) then
          Veff = (II_th0 + (IIc*(th1*(2d0+delta)*delta)-self%IC(3,n)*xi1)/(2d0*(1d0+delta)**2))*B
          f1(:,1) = f1(:,1) - real((th1*IIc+xi1*self%IC(3,n))*Keff/gR/(1d0+delta)**3*B)
          forall(i=1:Beam%npt/Beam%nBin)
            B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fldz(i,n)
          end forall
          Veff = -Keff/gamma*(Veff-0.5d0*i1*IIc*B)
        else
          Veff = -(self%IC(1,n)*Keff/gamma)*B
        endif
      endif
      f1(:,1) = f1(:,1)   - real(Veff)
      f1(:,2) = f1(:,2) +h*aimag(Veff)
    enddo
    f1(:,1) = f1(:,1)/gamma
    f1 = f1*ks
    !f2-------------------------------------------------------------
    tmp = Beam%pData(:,t_:g_) + f1*dz
    Beam%pData(:,t_:g_) = Beam%pData(:,t_:g_) + 0.5*f1*dz
    delta = tmp(:,2)/gR - 1d0
    f1(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))/ks*tmp(:,2)
    f1(:,2) = 0d0
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -B*(((self%IC(1,n)*0.5d0*h)*(Keff*zetaR/((1d0+delta)**2)) &
                             +self%I0(1,n))*px/tmp(:,2))
        f1(:,1) = f1(:,1) + ((self%IC(1,n)*h*zetaR/gR)*px/(1d0+delta)**3)*real(B)
      else
        if(h==1) then
          Veff = (II_th0 + (IIc*(th1*(2d0+delta)*delta)-self%IC(3,n)*xi1)/(2d0*(1d0+delta)**2))*B
          f1(:,1) = f1(:,1) - real((th1*IIc+xi1*self%IC(3,n))*Keff/gR/(1d0+delta)**3*B)
          forall(i=1:Beam%npt/Beam%nBin)
            B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*tmp((i-1)*Beam%nBin+1:i*Beam%nBin,1))*fldz(i,n)
          end forall
          Veff = -Keff/tmp(:,2)*(Veff-0.5d0*i1*IIc*B)
        else
          Veff = -(self%IC(1,n)*Keff/tmp(:,2))*B
        endif
      endif
      f1(:,1) = f1(:,1)   - real(Veff)
      f1(:,2) = f1(:,2) +h*aimag(Veff)
    enddo
    f1(:,1) = f1(:,1)/tmp(:,2)
    f1 = f1*ks
    !sum(f1,f2,k3,k4)-----------------------------------------------
    Beam%pData(:,t_:g_) = Beam%pData(:,t_:g_) + 0.5d0*dz*f1
    !===============================================================
  end associate
end subroutine undulator_tMap1_RK2

subroutine undulator_tMap1_RK1(self,Beam,fld,fldz,dz)
!=================================================================
! Longitudinal map under undulator and radiation field
!-----------------------------------------------------------------
  implicit none
  type(undulator),intent(in) :: self
  type(eBeam),    intent(inout) :: Beam
  real*8,         intent(in) :: dz
  complex(8),     intent(in),&
    dimension(Beam%npt/Beam%nBin,self%nHarm) :: fld, fldz

  integer :: h,n,i
  real*8  :: gR,xiR,zetaR
  real*8,    dimension(Beam%npt)  :: Keff,th0,th1,xi1,delta
  complex(8)                      :: IIc
  complex(8),dimension(Beam%npt)  :: Veff,B,II_th0
  real*8,    dimension(Beam%npt,2):: tmp
  
  associate(x =>Beam%pData(:,x_), px =>Beam%pData(:,px_),&
            y =>Beam%pData(:,y_), py =>Beam%pData(:,py_),&
            gamma=>Beam%pData(:,g_),  ks =>Beam%ks ,&
            K=>self%K, kx2=>self%kx2, ky2=>self%ky2, ku=>self%ku)
  
    gR = sqrt(0.5d0*ks/ku*(1d0+0.5d0*K*K))
    xiR = 0.125d0*ks*K*K/(ku*gR*gR)
    zetaR = ks*K/(ku*gR*gR)
    tmp(:,1) = Kx2*x*x + Ky2*y*y
    tmp(:,2) = px*px+py*py
    th0 = -(0.5d0*ks/(gR*gR))*(tmp(:,2) +0.5d0*K*K*tmp(:,1)*(1d0+0.25d0*tmp(:,1)))
    IIC = (0d0, 0d0)
    do n=1,self%nHarm
      h=self%harm_tab(n)
      if(h==1) then
        IIC    = -self%IC(2,n)/ku
        II_th0 = self%IC(1,n)+0.5d0*self%IC(3,n)*xiR +0.5d0*IIC*th0
        exit
      endif
    enddo
    Keff =  K*(1d0+0.5d0*tmp(:,1))
    th1  =  (0.5d0*ks/(gR*gR))*(1d0+tmp(:,2)+0.5d0*Keff*Keff)
    xi1  =  xiR*(1d0+tmp(:,1))
    !===============================================================
    !==== Longitudinal RK4 =========================================
    !f1-------------------------------------------------------------
    delta = gamma/gR - 1d0
    tmp(:,1) = (th0 + th1*delta*(2d0+delta)/((1d0+delta)**2))/ks*gamma
    tmp(:,2) = 0d0
    do n=1,self%nHarm
      h=self%harm_tab(n)
      forall(i=1:Beam%npt/Beam%nBin)
        B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fld(i,n)
      end forall
      if(mod(h,2)==0) then
        Veff = -B*(((self%IC(1,n)*0.5d0*h)*(Keff*zetaR/((1d0+delta)**2)) &
                             +self%I0(1,n))*px/gamma)
        tmp(:,1) = tmp(:,1) + ((self%IC(1,n)*h*zetaR/gR)*px/(1d0+delta)**3)*real(B)
      else
        if(h==1) then
          Veff = (II_th0 + (IIc*(th1*(2d0+delta)*delta)-self%IC(3,n)*xi1)/(2d0*(1d0+delta)**2))*B
          tmp(:,1) = tmp(:,1) - real((th1*IIc+xi1*self%IC(3,n))*Keff/gR/(1d0+delta)**3*B)
          forall(i=1:Beam%npt/Beam%nBin)
            B((i-1)*Beam%nBin+1:i*Beam%nBin) = exp(i1*h*Beam%pData((i-1)*Beam%nBin+1:i*Beam%nBin,t_))*fldz(i,n)
          end forall
          Veff = -Keff/gamma*(Veff-0.5d0*i1*IIc*B)
        else
          Veff = -(self%IC(1,n)*Keff/gamma)*B
        endif
      endif
      tmp(:,1) = tmp(:,1)   - real(Veff)
      tmp(:,2) = tmp(:,2) +h*aimag(Veff)
    enddo
    tmp(:,1) = tmp(:,1)/gamma
    !sum------------------------------------------------------------
    Beam%pData(:,t_:g_) = Beam%pData(:,t_:g_) + dz*(tmp*ks)
    !===============================================================
  end associate
end subroutine undulator_tMap1_RK1


subroutine track_undulator(self,beam,rad,dz)
  implicit none
  class(Undulator),intent(in)   :: self
  type(eBeam),intent(inout):: beam
  type(Radiation),intent(inout):: rad
  real*8, optional,intent(in)   :: dz
  integer :: i
  
  if(Rad%ks /= Beam%ks) then
    error stop 'Error <- track_undulator :: Rad%ks and Beam%ks does not agree. Program terminating...'
  endif
  if(self%ku /= Beam%ku) then
    error stop 'Error <- track_undulator :: undulator%ku and Beam%ku does not agree. Program terminating...'
  endif
  if(self%nHarm /= Rad%nHarm) then
    error stop 'Error <- track_undulator :: undulator%nHarm and Rad%nHarm does not agree. Program terminating...'
  endif
  if(all(self%harm_tab /= Rad%harm_tab)) then
    error stop 'Error <- track_undulator :: undulator%harm_tab and Rad%harm_tab does not agree. Program terminating...'
  endif

  if(present(dz)) then
    call self%track_slipsplit(beam,rad,dz)
  else
    do i=1,self%nStep
      call self%track_slipsplit(beam,rad,self%dz)
    enddo
  endif

end subroutine track_undulator

subroutine undulator_track_leapfrog(self,beam,rad,dz)
  implicit none
  class(Undulator),intent(in)   :: self
  type(eBeam),     intent(inout):: beam
  type(Radiation), intent(inout):: rad
  real*8,          intent(in)   :: dz
  call undulator_xyMap(self,Beam,0.25d0*dz)
  call Beam%Reorder()
  call undulator_tMap0(self,Beam,Rad,dz)
  call undulator_xyMap(self,Beam,0.25d0*dz)
  call Beam%Reorder()
  call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
  call Rad%FldSolve(dz)
  call Rad%Slip(dz)
end subroutine undulator_track_leapfrog

subroutine undulator_track_schmitt(self,beam,rad,dz)
  implicit none
  class(Undulator),intent(in)   :: self
  type(eBeam),     intent(inout):: beam
  type(Radiation), intent(inout):: rad
  real*8,          intent(in)   :: dz
  
  call undulator_xyMap(self,Beam,0.25d0*dz)
  call Beam%yReorder()
  call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
  call Rad%FldSolve(dz)
  call undulator_tMap0(self,Beam,Rad,dz)
  call Rad%Slip(dz)
  call Beam%tReorder()
  call undulator_xyMap(self,Beam,0.25d0*dz)
end subroutine undulator_track_schmitt

subroutine undulator_track_slipsplit(self,beam,rad,dz)
  implicit none
  class(Undulator),intent(in)   :: self
  type(eBeam),     intent(inout):: beam
  type(Radiation), intent(inout):: rad
  real*8,          intent(in)   :: dz
  if(self%HeffID==0) then 
    call undulator_track_slipsplit_Heff0(self,beam,rad,dz)
  else
    call undulator_track_slipsplit_Heff1(self,beam,rad,dz)
  endif
end subroutine

subroutine undulator_track_slipsplit_Heff0(self,beam,rad,dz)
  implicit none
  class(Undulator),intent(in)   :: self
  type(eBeam),     intent(inout):: beam
  type(Radiation), intent(inout):: rad
  real*8,          intent(in)   :: dz
  call undulator_xyMap(self,Beam,0.25d0*dz)
  call Beam%yReorder()
  call Rad%Slip(0.5d0*dz)
  call Beam%tReorder()
  call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
  call Rad%FldSolve(dz)
  call Rad%Slip(0.5d0*dz)
  call Beam%tReorder()
  call undulator_tMap0(self,Beam,Rad,dz)
  call undulator_xyMap(self,Beam,0.25d0*dz)
end subroutine undulator_track_slipsplit_Heff0

subroutine undulator_track_slipsplit_Heff1(self,beam,rad,dz)
  implicit none
  class(Undulator),intent(in)   :: self
  type(eBeam),     intent(inout):: beam
  type(Radiation), intent(inout):: rad
  real*8,          intent(in)   :: dz
  complex(8), allocatable :: FldTmp(:,:)
  
  call undulator_xyMap(self,Beam,0.25d0*dz)
  call Beam%yReorder()
  call Rad%Slip(0.5d0*dz)
  call Beam%tReorder()
  allocate(FldTmp(beam%npt/beam%nBin,rad%nHarm))
  call Beam%scatter(Rad,FldTmp)
  call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
  call Rad%FldSolve(dz)
  call Rad%Slip(0.5d0*dz)
  call Beam%tReorder(FldTmp)
  call undulator_tMap1(self,Beam,Rad,FldTmp,dz)
  call undulator_xyMap(self,Beam,0.25d0*dz)
end subroutine undulator_track_slipsplit_Heff1

subroutine dummy4comment()
! subroutine undulator_track_steady_leapfrog(self,beam,rad,dz)
  ! implicit none
  ! class(Undulator),intent(in)   :: self
  ! type(eBeam),     intent(inout):: beam
  ! type(Radiation), intent(inout):: rad
  ! real*8,          intent(in)   :: dz
  ! call undulator_xyMap(self,Beam,0.25d0*dz)
  ! call Beam%Reorder()
  ! call undulator_tMap0(self,Beam,Rad,dz)
  ! call undulator_xyMap(self,Beam,0.25d0*dz)
  ! call Beam%Reorder()
  ! call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
  ! call Rad%FldSolve(dz)
! end subroutine undulator_track_steady_leapfrog

! subroutine undulator_track_steady(self,beam,rad,dz)
  ! implicit none
  ! class(Undulator),intent(in)   :: self
  ! type(eBeam),     intent(inout):: beam
  ! type(Radiation), intent(inout):: rad
  ! real*8, optional,intent(in)   :: dz
    
  ! if(self%HeffID==0) then
    ! call undulator_track_steady_Heff0(self,beam,rad,dz)
  ! else
    ! call undulator_track_steady_Heff1(self,beam,rad,dz)
  ! endif
! end subroutine undulator_track_steady

! subroutine undulator_track_steady_Heff0(self,beam,rad,dz)
  ! implicit none
  ! class(Undulator),intent(in)   :: self
  ! type(eBeam),     intent(inout):: beam
  ! type(Radiation), intent(inout):: rad
  ! real*8, optional,intent(in)   :: dz
  ! integer :: i
  ! if(present(dz)) then
      ! call undulator_xyMap(self,Beam,0.25d0*dz)
      ! call Beam%yReorder()
      ! call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
      ! call Rad%FldSolve(dz)
      ! call undulator_tMap0(self,Beam,Rad,dz)
      ! call Beam%tReorder()
      ! call undulator_xyMap(self,Beam,0.25d0*dz)
  ! else
      ! call undulator_xyMap(self,Beam,0.25d0*self%dz)
    ! do i=1,self%nStep-1
      ! call Beam%yReorder()
      ! call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
      ! call Rad%FldSolve(self%dz)
      ! call undulator_tMap0(self,Beam,Rad,self%dz)
      ! call undulator_xyMap(self,Beam,0.5d0*self%dz)
    ! enddo
      ! call Beam%yReorder()
      ! call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
      ! call Rad%FldSolve(self%dz)
      ! call undulator_tMap0(self,Beam,Rad,self%dz)
      ! call undulator_xyMap(self,Beam,0.25d0*self%dz)
  ! endif
! end subroutine undulator_track_steady_Heff0

! subroutine undulator_track_steady_Heff1(self,beam,rad,dz)
  ! implicit none
  ! class(Undulator),intent(in)   :: self
  ! type(eBeam),     intent(inout):: beam
  ! type(Radiation), intent(inout):: rad
  ! real*8, optional,intent(in)   :: dz
  ! complex(8) :: fldz(Beam%npt/Beam%nBin,Rad%nHarm)
  ! integer :: i
  
  ! ! if(present(dz)) then
      ! ! call undulator_xyMap(self,Beam,0.25d0*dz)
      ! ! call Beam%yReorder()
      ! ! call Beam%scatter(Rad,fldz)
      ! ! call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
      ! ! call Rad%FldSolve(dz)
      ! ! call undulator_tMap1(self,Beam,Rad,fldz,dz)
      ! ! call undulator_xyMap(self,Beam,0.25d0*dz)
  ! ! else
      ! ! call undulator_xyMap(self,Beam,0.25d0*self%dz)
    ! ! do i=1,self%nStep-1
      ! ! call Beam%yReorder()
      ! ! call Beam%scatter(Rad,fldz)
      ! ! call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
      ! ! call Rad%FldSolve(self%dz)
      ! ! call undulator_tMap1(self,Beam,Rad,fldz,self%dz)
      ! ! call undulator_xyMap(self,Beam,0.5d0*self%dz)
    ! ! enddo
      ! ! call Beam%yReorder()
      ! ! call Beam%scatter(Rad,fldz)
      ! ! call Beam%Gather(Rad,[self%K,self%kx2,self%ky2],self%I0,self%IC,self%HeffID)
      ! ! call Rad%FldSolve(self%dz)
      ! ! call undulator_tMap1(self,Beam,Rad,fldz,self%dz)
      ! ! call undulator_xyMap(self,Beam,0.25d0*self%dz)
  ! ! endif
! end subroutine undulator_track_steady_Heff1

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMM   get_M   MMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! subroutine undulator_get_M(self,M,gamma)
! ! get transverse transfer matrix
  ! implicit none
  ! class(Undulator), intent(in) :: self
  ! real(8), intent(in) :: gamma
  ! real(8), intent(out):: M(4,4)
  ! real(8) :: bg, Gx, c,s,Mx(2,2), My(2,2)
  
  ! bg = sqrt(gamma**2-1d0)
  ! ! ==== Mx ====
  ! Gx = abs(self%K*self%kx/bg)
  ! if(Gx>tiny(Gx)) then
    ! c = cos(Gx*self%L)
    ! s = sin(Gx*self%L)
    ! Mx = reshape([c,s/Gx,-Gx*s,c],[2,2],order=[2,1])
  ! else
    ! Mx = reshape([1d0,self%L/bg,0d0,1d0],[2,2],order=[2,1])
  ! endif
  ! ! ==== My ====
  ! Gx = abs(self%K*self%ky/bg)
  ! if(Gx>tiny(Gx)) then
    ! c = cos(Gx*self%L)
    ! s = sin(Gx*self%L)
    ! My = reshape([c,s/Gx,-Gx*s,c],[2,2],order=[2,1])
  ! else
    ! My = reshape([1d0,self%L/bg,0d0,1d0],[2,2],order=[2,1])
  ! endif  
  ! M = 0d0
  ! M(1:2,1:2) = Mx
  ! M(3:4,3:4) = My
  ! return
! end subroutine undulator_get_M
end subroutine dummy4comment

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMM dummy deffered subroutines MMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine dummy_rad(self,rad,dz)
  implicit none
  class(Undulator),intent(in)   :: self
  type(Radiation),intent(inout) :: rad
  real*8, optional,intent(in)   :: dz
end subroutine dummy_rad

subroutine dummy_beam(self,beam,dz)
  implicit none
  class(Undulator),intent(in)  :: self
  type(eBeam),intent(inout)    :: beam
  real*8, optional,intent(in)  :: dz
end subroutine dummy_beam

end module undulatorModule
