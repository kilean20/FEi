module CompDomModule
  use ConstDataModule
  use MathModule
  use pGrid2dModule
  implicit none
  private
  type, public :: CompDom
    ! spatial range in each dimension.
    real*8, dimension(3,2) :: Rng    !(xMin,xMax,yMin,yMax,thetaMin,thetaMax)
    real*8, dimension(3,2) :: RngLoc ! local processor range
    ! mesh size i.e., distance between mesh points 
    real*8, dimension(3) :: MshSize  !(dx,dy,dtheta)
    real*8  :: dV
    integer, dimension(3) :: MshNum    !# of total mesh points
    integer, dimension(3) :: MshNumLoc !# of mesh points in each processor
    ! table of local meshes
    integer, allocatable, dimension(:) :: xMshTab  
    integer, allocatable, dimension(:) :: yMshTab
    integer, allocatable, dimension(:) :: tMshTab
    !========== Mesh Info ==========
    ! Mesh Table is, e.g., xMshTab(rank+1) = xMshTab(rank+1) +(MshNumLoc(1)-1)
    ! Mesh Range is, e.q., RngLoc(1,1) = Rng(1,1) + xMshTab(mycol)*MshSize(1)
    !===============================
    ! PML layer
    integer :: nxPML, nyPML
    ! weight / shape funtion ID
    integer :: wID    ! wID = 11  ! transversely and longitudinally Linear
                      ! wID = 00  ! transversely and longitudinally Uniform
                      ! wID = 10  ! transversely Linear and longitudinally Uniform
    ! MPI info
    contains
      procedure :: destruct => CompDom_destructor
      !procedure :: yLoadBalance, tLoadBalance
  end type CompDom
  
  interface CompDom
    procedure CompDom_constructor
  end interface 
  
contains


!=================================================================
!=====================3D CompDom Constructors=====================
!=================================================================
function CompDom_constructor(pGrd,Rng,MshNum,distID,Std,nslp,wID,nxPML,nyPML) result(cDom)
  !============================================================
  ! Construct Gaussian 2D Computation Domain over y and theta dimension and 1D uniform over x dimension
  ! input :
  !   pGrd = pGrid2D class, the grid of computational processors
  !   Rng = domain(xMin,yMin,thetaMin,xMax,yMax,,thetaMax) of interest excluding PML layer
  !   MshNum = number of mesh points such that, for example, computational resolution of dx is Rng(1,2)-Rng(1,1)/(MshNum-1)
  !   distID = 1 : uniform,  2: gaussian  for y-theta dimension
  !   Std = y[meter], theta[radian] standard deviation of e-Beam
  !   nslp = expected number of slippage
  !   nxPML,nyPML = number of mesh in PML layer, 7 by default. Need large Rng(~at least 9 times of e-beam xy rms size) when xyPML=0
  ! output :
  !   CompDom
  !============================================================
  implicit none
  class(CompDom),  pointer    :: cDom
  type(pGrid2D)               :: pGrd
  integer,         intent(in) :: MshNum(3), distID(2), nslp
  integer,optional,intent(in) :: wID, nxPML, nyPML
  real*8 ,         intent(in) :: Rng(3,2), Std(2)
  allocate(cDom)
  allocate(cDom%xMshTab(0:pGrd%npCol))
  allocate(cDom%yMshTab(0:pGrd%npCol))
  allocate(cDom%tMshTab(0:pGrd%npRow))
  if(present(nxPML)) then
    cDom%nxPML = nxPML
  else
    cDom%nxPML = 0
  endif
  if(present(nyPML)) then
    cDom%nyPML = nyPML
  else
    cDom%nyPML = 0
  endif
  if(present(wID)) then
    cDom%wID = wID
  else
    cDom%wID = 0
  endif
  
  call construct_xCompDom_uniform(cDom,pGrd,Rng(1,1),Rng(1,2),MshNum(1),cDom%nxPML)
  select case (distID(1))
  case (1)
  !if(Std(1)<=0) then
    call construct_yCompDom_uniform(cDom,pGrd,Rng(2,1),Rng(2,2),MshNum(2),cDom%nyPML)
  case (2) 
  !else
    call construct_yCompDom_gaussian(cDom,pGrd,Rng(2,1),Rng(2,2),Std(1),MshNum(2),cDom%nyPML)
  end select
  !endif
  select case (distID(2))
  case (2)
    call construct_tCompDom_gaussian(cDom,pGrd,Rng(3,1),Rng(3,2),Std(2),MshNum(3),nslp)
  case default
    call construct_tCompDom_uniform (cDom,pGrd,Rng(3,1),Rng(3,2),       MshNum(3),nslp)
  end select
  !endif
  cDom%dV = cDom%MshSize(1)*cDom%MshSize(2)*cDom%MshSize(3)
end function CompDom_constructor

subroutine CompDom_destructor(self)
  class(CompDom), intent(inout) :: self
  if(allocated(self%xMshTab))  deallocate(self%xMshTab)
  if(allocated(self%yMshTab))  deallocate(self%yMshTab)
  if(allocated(self%tMshTab))  deallocate(self%tMshTab)
end subroutine


!=================================================================
!=========Uniform CompDom Constructors 4 each dimension===========
!=================================================================
subroutine construct_xCompDom_uniform(self,pGrd,xmin,xmax,nx,nxPML)
  !============================================================
  ! Construct uniform computation domain over x domain
  ! input :
  !   pGrd = pGrid2D type, the grid of computational processors
  !   xmin,xmax[meter] = domain of interest excluding PML layer
  !   nx = number of mesh points such that computational resolution of dx is (xMax-xMin)/(nx-1)
  !   nxPML = number of mesh in PML layer, 7 by default. Need large (xMax-xMin)(~at least 9 times of e-beam x rms size) when xyPML=0
  ! in/output :
  !   self = CompDom type, the computational domain
  !============================================================
  implicit none
  type(pGrid2D), intent(in)   :: pGrd
  type(CompDom), intent(inout):: self
  integer,       intent(in)   :: nx
  integer,       intent(in)   :: nxPML
  real*8 ,       intent(in)   :: xmin,xmax
  integer                     :: i, nEach, remain
  self%nxPML = nxPML
  self%MshNum(1) = nx + 2*self%nxPML
  self%MshSize(1) = (xmax-xmin)/(nx-1)
  self%Rng(1,1) = xmin - self%nxPML*self%MshSize(1)
  self%Rng(1,2) = xmax + self%nxPML*self%MshSize(1)
  ! ==== Domain decomposition onto the grid of computational processors ====
  nEach = (self%MshNum(1)-1)/pGrd%npCol
  remain = (self%MshNum(1)-1) - pGrd%npCol*nEach
  !allocate(self%xMshTab(0:pGrd%npCol))
  self%xMshTab(0) = 0
  do i=1,(pGrd%npCol-remain)/2
    self%xMshTab(i) = self%xMshTab(i-1) + nEach
  enddo
  do i=(pGrd%npCol-remain)/2+1,(pGrd%npCol-remain)/2+remain
    self%xMshTab(i) = self%xMshTab(i-1) + nEach + 1
  enddo
  do i=(pGrd%npCol+remain)/2+1,pGrd%npCol
    self%xMshTab(i) = self%xMshTab(i-1) + nEach
  enddo
  self%RngLoc(1,1) = self%Rng(1,1) +self%xMshTab(pGrd%mycol  )*self%MshSize(1)
  self%RngLoc(1,2) = self%Rng(1,1) +self%xMshTab(pGrd%mycol+1)*self%MshSize(1)
  self%MshNumLoc(1) = self%xMshTab(pGrd%mycol+1) -self%xMshTab(pGrd%mycol) +1
end subroutine construct_xCompDom_uniform

subroutine construct_yCompDom_uniform(self,pGrd,yMin,yMax,ny,nyPML)
  !============================================================
  ! Construct uniform computation domain over y domain
  ! input :
  !   pGrd = pGrid2D type, the grid of computational processors
  !   yMin,yMax[meter] = domain of interest excluding PML layer
  !   ny = number of mesh points such that computational resolution of dy is (yMax-yMin)/(ny-1)
  !   nyPML = number of mesh in PML layer, 7 by default. Need large (yMax-yMin)(~at least 9 times of e-beam x rms size) when xyPML=0
  ! in/output :
  !   self = CompDom type, the computational domain
  !============================================================
  implicit none
  type(pGrid2D), intent(in)   :: pGrd
  type(CompDom), intent(inout):: self
  integer,       intent(in)   :: ny
  integer,       intent(in)   :: nyPML
  real*8 ,       intent(in)   :: yMin,yMax
  integer                     :: i, nEach, remain
  self%nyPML = nyPML
  self%MshNum(2) = ny + 2*self%nyPML
  self%MshSize(2) = (yMax-yMin)/(ny-1)
  self%Rng(2,1) = yMin - self%nyPML*self%MshSize(2)
  self%Rng(2,2) = yMax + self%nyPML*self%MshSize(2)
  ! ==== Domain decomposition onto the grid of computational processors ====
  nEach = (self%MshNum(2)-1)/pGrd%npCol
  remain = (self%MshNum(2)-1) - pGrd%npCol*nEach
  if(.not. allocated(self%yMshTab)) allocate(self%yMshTab(0:pGrd%npCol))
  self%yMshTab(0) = 0
  do i=1,(pGrd%npCol-remain)/2
    self%yMshTab(i) = self%yMshTab(i-1) + nEach
  enddo
  do i=(pGrd%npCol-remain)/2+1,(pGrd%npCol-remain)/2+remain
    self%yMshTab(i) = self%yMshTab(i-1) + nEach + 1
  enddo
  do i=(pGrd%npCol+remain)/2+1,pGrd%npCol
    self%yMshTab(i) = self%yMshTab(i-1) + nEach
  enddo
  self%RngLoc(2,1) = self%Rng(2,1) +self%yMshTab(pGrd%mycol  )*self%MshSize(2)
  self%RngLoc(2,2) = self%Rng(2,1) +self%yMshTab(pGrd%mycol+1)*self%MshSize(2)
  self%MshNumLoc(2) = self%yMshTab(pGrd%mycol+1) -self%yMshTab(pGrd%mycol) +1
end subroutine construct_yCompDom_uniform

subroutine construct_tCompDom_uniform(self,pGrd,tMin,tMax,nt,nslp)
  !============================================================
  ! Construct uniform computation domain over theta domain without slippage. 
  ! myCol = 0 have additional domain which includes future slippage
  ! input :
  !   pGrd = pGrid2D type, the grid of computational processors
  !   tMin,tMax[radian] = domain of interest
  !   nt = number of mesh points such that computational resolution is dt ~= (tMax-tMin)/(nt-1)
  !   nslp = (int) expected slippage length in unit of wavelength
  ! in/output :
  !   self = CompDom type, the computational domain
  !============================================================
  implicit none
  type(pGrid2D),intent(in)   :: pGrd
  type(CompDom),intent(inout):: self
  integer,      intent(in)   :: nt, nslp
  real*8,       intent(in)   :: tMin,tMax
  integer                    :: i, nMshSlip, nEach, remain
  ! ==== total longitudinal mesh ====
  self%Rng(3,2) = tMax
  self%MshSize(3) = (tMax -tMin)/(nt-1)
  ! correction on mesh number in consideration of slippage
  self%MshNum(3) = ceiling((tMax -tMin +twopi*nslp-1d-6)/self%MshSize(3)) + 1
  ! correction on computational domain in consideration of slippage
  self%Rng(3,1) = tMax - (self%MshNum(3)-1)*self%MshSize(3)
  ! ==== Domain decomposition onto the grid of computational processors ====
  if(.not. allocated(self%tMshTab)) allocate(self%tMshTab(0:pGrd%npRow))
  self%tMshTab(0) = 0
  if(pGrd%npRow > 1) then
    nMshSlip = int(twopi*nslp/self%MshSize(3)+1d-6)
    nEach = (self%MshNum(3)-1-nMshSlip)/pGrd%npRow
    remain = (self%MshNum(3)-1-nMshSlip) - pGrd%npRow*nEach
    do i=1,(pGrd%npRow-remain)/2
      self%tMshTab(i) = self%tMshTab(i-1) + nEach
    enddo
    do i=(pGrd%npRow-remain)/2+1,(pGrd%npRow-remain)/2+remain
      self%tMshTab(i) = self%tMshTab(i-1) + nEach + 1
    enddo
    do i=(pGrd%npRow+remain)/2+1,pGrd%npRow
      self%tMshTab(i) = self%tMshTab(i-1) + nEach
    enddo
    self%tMshTab(1:) = self%tMshTab(1:) + nMshSlip
  endif
  self%tMshTab(pGrd%npRow) = self%MshNum(3) -1
  self%RngLoc(3,1) = self%Rng(3,1) +self%tMshTab(pGrd%myrow  )*self%MshSize(3)
  self%RngLoc(3,2) = self%Rng(3,1) +self%tMshTab(pGrd%myrow+1)*self%MshSize(3)
  self%MshNumLoc(3) = self%tMshTab(pGrd%myrow+1) -self%tMshTab(pGrd%myrow) +1
end subroutine construct_tCompDom_uniform


!=================================================================
!=========Gaussian CompDom Constructors 4 each dimension==========
!=================================================================
subroutine construct_yCompDom_gaussian(self,pGrd,yMin,yMax,yStd,ny,nyPML)
  !============================================================
  ! Construct Gaussian computation domain over y domain
  ! input :
  !   pGrd = pGrid2D type, the grid of computational processors
  !   yMin,yMax[meter] = domain of interest excluding PML layer
  !   yStd[meter] = y-standard deviation of e-Beam
  !   ny = number of mesh points such that computational resolution is dy=(yMax-yMin)/(ny-1)
  !   nyPML = number of mesh in PML layer, 7 by default. Need large (yMax-yMin)(~at least 9 times of e-beam y rms size) when xyPML=0
  ! in/output :
  !   self = CompDom type, the computational domain
  !============================================================
  implicit none
  type(pGrid2D), intent(in)   :: pGrd
  type(CompDom), intent(inout):: self
  integer,       intent(in)   :: ny
  integer,       intent(in)   :: nyPML
  real*8,        intent(in)   :: yMin,yMax,yStd
  integer                     :: i, nEach
  real*8                      :: low, high, erfL, erfR, estimate
  self%nyPML = nyPML
  self%MshNum(2) = ny +2*self%nyPML
  self%MshSize(2) = (yMax-yMin)/(ny-1)
  self%Rng(2,1) = yMin - self%nyPML*self%MshSize(2)
  self%Rng(2,2) = yMax + self%nyPML*self%MshSize(2)  
  ! ==== Domain decomposition onto the grid of computational processors ====
  erfL = Erf(-self%Rng(2,1)/yStd/sqrt(2d0))
  erfR = Erf( self%Rng(2,2)/yStd/sqrt(2d0))
  nEach = (self%MshNum(2)-1)/pGrd%npCol
  !allocate(self%yMshTab(0:pGrd%npCol))
  self%yMshTab(0) = 0
  low = self%Rng(2,1)/yStd
  high = self%Rng(2,2)/yStd
  do i=1,pGrd%npCol-1
    estimate = sqrt(2d0)*math%iErf(dble(i)/pGrd%npCol*(erfR+erfL)-erfL)
    estimate = (estimate - low)/(high-low)*(self%MshNum(2)-1)
    self%yMshTab(i) = nint(estimate)
  enddo
  !**correction when zero mesh is assigned to a processor
  do i=pGrd%npCol/2,2,-1
    if (self%yMshTab(i)==self%yMshTab(i-1)) then
      self%yMshTab(1:i-1) = self%yMshTab(1:i-1)-1
    endif
  enddo
  do i=pGrd%npCol/2+1,pGrd%npCol-1
    if (self%yMshTab(i)==self%yMshTab(i-1)) then
      self%yMshTab(i:pGrd%npCol-1) = self%yMshTab(i:pGrd%npCol-1)+1
    endif
  enddo
  !** end of correction
  self%yMshTab(pGrd%npCol) = self%MshNum(2) -1
  self%RngLoc(2,1) = self%Rng(2,1) +self%yMshTab(pGrd%mycol  )*self%MshSize(2)
  self%RngLoc(2,2) = self%Rng(2,1) +self%yMshTab(pGrd%mycol+1)*self%MshSize(2)
  self%MshNumLoc(2) = self%yMshTab(pGrd%mycol+1) -self%yMshTab(pGrd%mycol) +1
end subroutine construct_yCompDom_gaussian 

subroutine construct_tCompDom_gaussian(self,pGrd,tMin,tMax,tStd,nt,nslp)
  !============================================================
  ! Construct Gaussian 2D Computation Domain over theta dimension
  ! input :
  !   pGrd = pGrid2D type, the grid of computational processors
  !   tMin,tMax[meter] = domain of interest
  !   tStd[radian] =  theta-standard deviation of e-Beam
  !   nt = number of mesh points such that computational resolution is dt ~= (tMax-tMin)/(nt-1)
  !   nslp = expected slippage number
  ! in/output :
  !   self = CompDom type, the computational domain
  !============================================================
  implicit none
  type(pGrid2D),intent(in)   :: pGrd
  type(CompDom),intent(inout):: self
  integer,      intent(in)   :: nt, nslp
  real*8,       intent(in)   :: tMin,tMax,tStd
  integer                    :: i, nEach
  real*8                     :: low, high, erfL, erfR, estimate
  ! ---- total longitudinal mesh ----
  self%Rng(3,2) = tMax
  ! correction on mesh size in unit of buckets
  self%MshSize(3) = twopi*nint((tMax -tMin)/twopi/(nt-1))
  if (self%MshSize(3)<twopi) self%MshSize(3) = twopi
  ! correction on mesh number in consideration of slippage
  self%MshNum(3) = ceiling((tMax -tMin +twopi*nslp)/self%MshSize(3)) + 1
  ! correction on computational domain in consideration of slippage
  self%Rng(3,1) = tMax - (self%MshNum(3)-1)*self%MshSize(3)
  ! ==== Domain decomposition onto the grid of computational processors ====
  ! ---- theta domain : Gaussian decomposition
  erfL = Erf(-self%Rng(3,1)/tStd/sqrt(2d0))
  erfR = Erf(self%Rng(3,2)/tStd/sqrt(2d0))
  nEach = (self%MshNum(3)-1)/pGrd%npRow
  if(.not. allocated(self%tMshTab)) allocate(self%tMshTab(0:pGrd%npRow))
  self%tMshTab(0) = 0
  low = self%Rng(3,1)/tStd
  high = self%Rng(3,2)/tStd
  do i=1,pGrd%npRow-1
    estimate = sqrt(2d0)*math%iErf(dble(i)/pGrd%npRow*(erfR+erfL)-erfL)
    estimate = (estimate - low)/(high-low)*(self%MshNum(3)-1)
    self%tMshTab(i) = nint(estimate)
  enddo
  !**correction when zero mesh is assigned to a processor
  do i=pGrd%npRow/2,2,-1
    if (self%tMshTab(i)==self%tMshTab(i-1)) then
      self%tMshTab(1:i-1) = self%tMshTab(1:i-1)-1
    endif
  enddo
  do i=pGrd%npRow/2+1,pGrd%npRow-1
    if (self%tMshTab(i)==self%tMshTab(i-1)) then
      self%tMshTab(i:pGrd%npRow-1) = self%tMshTab(i:pGrd%npRow-1)+1
    endif
  enddo
  !**end of correction
  self%tMshTab(pGrd%npRow) = self%MshNum(3) -1
  self%RngLoc(3,1) = self%Rng(3,1) +self%tMshTab(pGrd%myrow  )*self%MshSize(3)
  self%RngLoc(3,2) = self%Rng(3,1) +self%tMshTab(pGrd%myrow+1)*self%MshSize(3)
  self%MshNumLoc(3) = self%tMshTab(pGrd%myrow+1) -self%tMshTab(pGrd%myrow) +1
end subroutine construct_tCompDom_gaussian

!=================================================================
!======================= Load Balancing ==========================
!=================================================================
end module CompDomModule
