module eBeamModule
  use ConstDataModule  ! included in CompDomModule
  use MathModule
  use pGrid2DModule
  use CompDomModule
  use RadiationModule
  use fileIO_Module
 

  implicit none
  private
  
  type, public :: eBeam
    integer :: npt,nBin
    real*8  :: ks,ku,z
    real*8, allocatable :: pData(:,:)
    integer, allocatable :: beamletID(:)
    type(pGrid2D), pointer :: pGrd => null()
    type(CompDom), pointer :: cDom => null()
    contains
      procedure :: correct_centroid => correct_centroid_5D, &
                                       correct_centroid_slices_5D
      procedure :: destruct => eBeam_destructor
      procedure :: assign_beamletID, populate_beamlet, Load_Twiss, set=>eBeam_set
      procedure :: ReOrder, yReOrder, tReOrder, tReorder_middle, add_shotnoise
      procedure :: yLoadBalance, tLoadBalance
      procedure :: release_beamlet
      procedure :: gather, scatter
      procedure :: get_sample,Write_Slices,Read_Slices,Write_Loc,Read_Loc
      procedure :: get_mean, get_std
  end type eBeam

  interface eBeam 
    procedure eBeam_constructor, &
              eBeam_constructor_from_pData, &
              eBeam_constructor_distID, &
              eBeam_constructor_real_enum
  end interface

  ! === example use ===================
  ! type(eBeam),pointer :: beam
  ! beam -> eBeam(pGrd,cDom,n_beamlet,n_bin,ks,ku)
  ! beam%yLoadBalance(rad)
  ! ===================================
contains


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMM   beam loadingMM  MMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
function eBeam_constructor(pGrd,cDom,nLet,nBin,ks,ku,z) result(self)
  class(eBeam), pointer :: self
  type(pGrid2D),pointer :: pGrd
  type(CompDom),pointer :: cDom
  integer,intent(in) :: nLet,nBin
  real*8, intent(in) :: ks,ku
  real*8, intent(in), optional :: z
  
  allocate(self)
  self%pGrd => pGrd
  self%cDom => cDom
  if(present(z)) then
    self%z = z
  else
    self%z = 0
  endif
  self%ks  = ks
  self%ku  = ku
  self%npt = nLet
  self%nBin= nBin
  allocate(self%pData(nLet,7))
  allocate(self%beamletID(nLet))
  
end function eBeam_constructor

function eBeam_constructor_from_pData(pGrd,cDom,pData,beamletID,npt,nBin,ks,ku,z) result(self)
  class(eBeam), pointer :: self
  type(pGrid2D),pointer :: pGrd
  type(CompDom),pointer :: cDom
  integer,intent(in) :: npt,nBin
  real*8, intent(in) :: ks,ku,pData(npt,7)
  integer,optional,intent(in) :: beamletID(npt/nBin)
  real*8, intent(in), optional :: z
  
  allocate(self)
  self%pGrd => pGrd
  self%cDom => cDom
  call self%set(pData,beamletID,npt,nBin,ks,ku,z)
end function eBeam_constructor_from_pData

function eBeam_constructor_distID(pGrd,cDom,nLet,nBin,ks,ku,Ip,&
                                  distID,Mean,Std,Alpha,Emit,lowCL,upCL,z,shotID,iHamsl) result(self)
!=================================================================================================
!!!!!!
! input
!   nLet   [int    ] : number of beamlets  per MPI processor
!   nBin   [int    ] : number of particles per beamlet
!   ks,ku  [dble   ] : reference wave number of radiation and undulator respectively, 
!                       in order to define theta = ks*(z-ct)+ku*z
!   z      [dble   ] : longitudinal coordinate of the beam
!   Ip     [dble   ] : peak current (Ampere)
!   distID [int(6) ] : distribution id (1:uniform, 2:normal)
!   Mean   [dble(6)] : centroid offset    of x[meter],px,y[meter],py,theta[rad],gamma
!   Std    [dble(3)] : standard deviation of x[meter],y[meter],theta[rad]
!   alpha  [dble(3)] : twiss alpha
!   emit   [dble(3)] : emittance of x[meter],y[meter],theta[rad]
!   lowCL  [dble,optional] : lower bound confidence level
!   upCL   [dble,optional] : upper bound confidence level
!   shotID [int, optional] : (1: temporal coordinate perturbation[Fawley],
!                             2: Charge perturbation[M.Neil],
!                             3: weighted charge Perturbation[M.Neil] )
!   iHamsl [int, optional] : starting sequence number of Hammersley set
!-------------------------------------------------------------------------------------------------
  implicit none
  class(eBeam), pointer :: self
  type(pGrid2D),pointer :: pGrd
  type(CompDom),pointer :: cDom
  
  integer,intent(in) :: nLet,nBin,distID(6)
  integer,intent(in), optional :: shotID,iHamsl
  
  real*8, intent(in) :: ks,ku,Ip
  real*8, intent(in) :: Std(3),Alpha(3),Emit(3),Mean(6)
  real*8, intent(in), optional :: lowCL,upCL,z
  
  real*8  :: meanSlices(1,6)

  ! if(pGrd%myRank==0) print*, '==eBeam_constructor=='
  ! if(pGrd%myRank==0) print*, 'Std=',Std
  ! if(pGrd%myRank==0) print*, 'Alpha=',Alpha
  ! if(pGrd%myRank==0) print*, 'emit=',emit
  allocate(self)
  self%pGrd => pGrd
  self%cDom => cDom
  if(present(z)) then
    self%z = z
  else
    self%z = 0
  endif
  self%ks  = ks
  self%ku  = ku
  self%npt = nLet
  self%nBin= 1
  allocate(self%pData(nLet,7))
  
  ! load beamlet
  call Load_Twiss(self,distID,Mean,Std,Alpha,Emit,iHamsl,lowCL,upCL)
  self%pData(:,q_) = Get_Q(Ip,std(3),ks,distID(t_))/(nLet*pGrd%np)
  call Reorder(self)
  call tReorder_middle(self)
  
  
  ! correction numerical shotnoise on beamlet centroid 
  if(nBin>1) then
    meanSlices(1,1) = 0d0
    meanSlices(1,2:5) = Mean(1:4)
    meanSlices(1,6) = Mean(6)
    call correct_centroid_slices_5D(self,meanSlices,1)
  endif
  
  ! populate particles of beamlet
  self%nBin= nBin
  if(self%nBin>1) then
    call populate_beamlet(self)
  endif 
  if(present(shotID)) then
    call add_shotnoise(self,shotID,distID(t_),[Mean(t_),Std(3)])
  endif
end function eBeam_constructor_distID

function eBeam_constructor_real_enum(pGrd,cDom,ks,ku,Ip,&
                                     distID,Mean,Std,lowCL,upCL,z,seed) result(self)
  class(eBeam), pointer :: self
  type(pGrid2D),pointer :: pGrd
  type(CompDom),pointer :: cDom
  integer,intent(in) :: distID(6)
  real*8, intent(in) :: ks,ku,Ip
  real*8, intent(in) :: Mean(6),Std(6),lowCL(6),upCL(6)  
  integer,intent(in), optional ::seed
  real*8, intent(in), optional :: z
  
  integer:: npt,n
  integer,allocatable:: seedarr(:)
  
  if(present(seed)) then
    call random_seed(size = n)
    allocate(seedarr(n))
    call random_seed(get=seedarr)
    seedarr = seedarr*(seed+self%pGrd%myRank+1)
    call random_seed(put=seedarr)
  endif
  
  npt = int(Get_Q(Ip,std(t_),ks,distID(t_))/PhysConst%e/self%pGrd%np)
  self => eBeam(pGrd,cDom,npt,1,ks,ku,Ip,&
                distID,Mean,Std,lowCL,upCL,z)
                
end function eBeam_constructor_real_enum

subroutine eBeam_destructor(self)
  class(eBeam) :: self
  if(allocated(self%pData))  deallocate(self%pData)
  self%pGrd => null()
  self%cDom => null()
  self%npt=0
end subroutine eBeam_destructor

real*8 function Get_Q(Ip,std,ks,distID)
!=====================================================
! Get total charge from given peak current
! input :
!  Ip : peak current
!  std: theta standard deviation in unit of radian
!  ks : radiation wavenumber used to define theta
!  distID : longitudinal distribution ID
!           1 : uniform
!           2 : normal
!-----------------------------------------------------
  implicit none
  integer,intent(in) :: distID
  real*8, intent(in) :: Ip,ks,std
  select case(distID)
    case(2)
      Get_Q = Ip*sqrt(twopi)*std/(ks*PhysConst%c)
    case default
      Get_Q = Ip*sqrt(12d0) *std/(ks*PhysConst%c)
  end select
end function Get_Q

subroutine eBeam_set(self,pData,beamletID,npt,nBin,ks,ku,z)
  class(eBeam) :: self
  integer,intent(in) :: npt,nBin
  real*8, intent(in) :: pData(npt,7), ks, ku
  integer,intent(in), optional :: beamletID(npt)
  real*8, intent(in), optional :: z

  if(present(z)) then
    self%z = z
  else
    self%z = 0
  endif
  if(present(beamletID)) then
    self%beamletID = beamletID
  else
    call self%assign_beamletID()
  endif
  self%ks  = ks
  self%ku  = ku
  self%npt = npt
  self%nBin= nBin
  if(allocated(self%pData)) deallocate(self%pData)
  allocate(self%pData(npt,7))
  self%pData = pData
end subroutine

subroutine assign_beamletID(self)
  implicit none
  class(eBeam) :: self
end subroutine 

subroutine populate_beamlet(self)
!=================================================================================================
!  populate beamlet from exiting macro-particles in beam
!-------------------------------------------------------------------------------------------------
  implicit none
  class(eBeam) :: self
  real*8 :: pDataTmp(self%npt,7)
  integer :: i,j
  
  associate(nBin=>self%nBin, npt=>self%npt)
    if(nBin==1) then
      return
    endif
    if(.not. size(self%pData) == npt*7) then
      print*, 'Error <- beam%populate_beamlet :',&
              'Mismatch b/w size of pData and npt. Enforcing npt=len(pData)'
      npt=size(self%pData)/7
    endif
    pDataTmp = self%pData(1:npt,1:7)
    deallocate(self%pData)
    allocate(self%pData(npt*nBin,7))
    do i=1,npt
      do j=1,4
        self%pData((i-1)*nBin+1:i*nBin,j) = pDataTmp(i,j)
      enddo
      self%pData((i-1)*nBin+1:i*nBin,g_) = pDataTmp(i,g_)
      self%pData((i-1)*nBin+1:i*nBin,q_) = pDataTmp(i,q_)/nBin
      do j=1,nBin
        self%pData((i-1)*nBin+j,t_) = pDataTmp(i,t_) -pi +twopi*(j-0.5d0)/dble(nBin)
      enddo
    enddo
    npt = npt*nBin
  end associate
end subroutine populate_beamlet

subroutine Load_Twiss(self,distID,Mean,Std,Alpha,Emit,iHamsl,lowCL,upCL)
!=================================================================================================
!             load phase-space using twiss parameters
! input
!   distID[int(6)]  : distribution id for x,px,y,py,theta,gamma
!   Mean  [dble(6)] : offset of x,px,y,py,theta,gamma
!   Std   [dble(3)] : standard deviation of x[meter], y[meter],theta[radian]
!   Alpha [dble(3)] : twiss paramter alpha
!   Emit  [dble(3)] : (normalzied)emittance of x-px[meter-rad], y-py[meter-rad], theta-gamma[rad-1])
!   lowCL    [dble,optional] : lower bound confidence level
!   upCL     [dble,optional] : upper bound confidence level
!   iHamsl   [int,    optional] : optional. starting sequence number of Hammersley set.                               
!-------------------------------------------------------------------------------------------------
  implicit none
  class(eBeam) :: self
  integer,intent(in) :: distID(6)
  real*8, intent(in) :: Std(3),Alpha(3),Emit(3),Mean(6)
  integer,intent(in),optional :: iHamsl
  real*8, intent(in),optional :: lowCL,upCL
  integer:: i,j,iHamslLoc
  real*8 :: ibet_x, std_p, std_xp
  
  ! if(self%pGrd%myRank==0) print*, '==Load_Twiss=='
  ! if(self%pGrd%myRank==0) print*, 'Std=',Std
  ! if(self%pGrd%myRank==0) print*, 'Alpha=',Alpha
  ! if(self%pGrd%myRank==0) print*, 'emit=',emit
  ! if(self%pGrd%myRank==0) print*, 'Mean=',Mean
  ! if(self%pGrd%myRank==0) print*, 'distID=',distID
  
  iHamslLoc = 0
  if(present(iHamsl)) then
    if(iHamsl>0) iHamslLoc = iHamsl + self%pGrd%myRank * self%npt
  endif
  !print*, 'iHamslLoc',iHamslLoc
  do j=1,6
    do i=1,self%npt
      self%pData(i,j) = math%randID(distID(j),j,iHamslLoc,lowCL,upCL)
    enddo
  enddo
  do i=1,3
    ibet_x = Emit(i)/(Std(i)**2)
    std_p = Std(i)*ibet_x
    std_xp = -Alpha(i)*ibet_x
    self%pData(:,2*i-1) = Std(i)*self%pData(:,2*i-1)
    self%pData(:,2*i  ) = std_xp*self%pData(:,2*i-1) +std_p*self%pData(:,2*i)
  enddo
 
  forall(j=1:6)
    self%pData(:,j) = self%pData(:,j) + Mean(j)
  end forall
end subroutine Load_Twiss


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMM   correct_centroid  MMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine correct_centroid_5D(self,Mean)
!===============================================================
!  correction on numerical centroid shot error over whole bunch
!---------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  real*8, intent(in) :: Mean(5)
  real*8 :: sBuff(6),rBuff(6)
  integer :: j,ierr
  
  associate(pData=>self%pData, npt=>self%npt, comm_2d=>self%pGrd%comm_2d)
    do j=1,4
      sBuff(j) = sum(pData(:,j )*pData(:,q_),dim=1)
    enddo
      sBuff(5) = sum(pData(:,g_)*pData(:,q_),dim=1)
      sBuff(6) = sum(pData(:,q_))
    call MPI_ALLreduce(sBuff,rBuff,6,MPI_double,MPI_sum,comm_2d,ierr)
    if(npt==0) return
    rBuff(1:5) = rBuff(1:5)/rBuff(6)
    do j=1,4
      pData(:,j ) = pData(:,j ) -rBuff(j) +Mean(j)
    enddo
      pData(:,g_) = pData(:,g_) -rBuff(5) +Mean(5)
  end associate
  
end subroutine correct_centroid_5D

subroutine correct_centroid_slices_5D(self,Mean,nPoints)
!===============================================================
! correction on numerical centroid shot error slice by slice
!---------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  integer,intent(in)   :: nPoints
  real*8, intent(in)   :: Mean(nPoints,6)
  real*8, dimension(self%cDom%MshNumLoc(3),6) :: sBuff,rBuff,MeanSlice
  integer, allocatable :: gt(:)
  integer :: ierr,i,j
  ! type(logIO),    pointer :: fLog

  call tReorder_middle(self)
  associate(pData=>self%pData, npt=>self%npt, &
            cDom=>self%cDom,   pGrd=>self%pGrd)
    ! interpolate Mean data to Grid points
    do i=1,cDom%MshNumLoc(3)
      MeanSlice(i,1) = cDom%RngLoc(3,1) + (i-1)*cDom%MshSize(3)
    enddo
    call math%linear_interpolation(Mean,nPoints,MeanSlice,cDom%MshNumLoc(3),6)
    allocate(gt(npt))
    gt = int((pData(:,t_)-cDom%RngLoc(3,1))/cDom%MshSize(3)+0.5d0)+1
    sBuff = 0d0
    do i=1,npt
      do j=1,4
        sBuff(gt(i),j) = sBuff(gt(i),j) +pData(i,j )*pData(i,q_)
      enddo
        sBuff(gt(i),5) = sBuff(gt(i),5) +pData(i,g_)*pData(i,q_)
        sBuff(gt(i),6) = sBuff(gt(i),6) +pData(i,q_)
    enddo
    sBuff = sBuff/physConst%e
    if(pGrd%npCol>1) then
      call MPI_ALLreduce(sBuff,rBuff,cDom%MshNumLoc(3)*6,MPI_double,MPI_sum,pGrd%comm_col,ierr)
    else
      rBuff = sBuff
    endif
    forall(j=1:5, i=1:cDom%MshNumLoc(3), rBuff(i,6) .ne. 0d0)
      rBuff(i,j) = rBuff(i,j)/rBuff(i,6)
    end forall
    
    ! fLog => logIO('log.',pGrd%myRank)
    ! call fLog%reopen()
    ! write(fLog%unit,*), 'cDom%RngLoc(3,:)=',cDom%RngLoc(3,:)
    ! write(fLog%unit,*), 'minval(MeanSlice(:,1))=',minval(MeanSlice(:,1))
    ! write(fLog%unit,*), 'maxval(MeanSlice(:,1))=',maxval(MeanSlice(:,1))
    ! write(fLog%unit,*), 'minval(MeanSlice(:,2))=',minval(MeanSlice(:,2))
    ! write(fLog%unit,*), 'maxval(MeanSlice(:,2))=',maxval(MeanSlice(:,2))
    ! write(fLog%unit,*), 'maxval(abs(rBuff(:,x_)))=',maxval(abs(rBuff(:,x_)))
    ! write(fLog%unit,*), 'maxval(abs(rBuff(:,y_)))=',maxval(abs(rBuff(:,y_)))
    ! write(fLog%unit,*), 'minval(rBuff(:,5))=',minval(rBuff(:,5))
    ! write(fLog%unit,*), 'maxval(rBuff(:,5))=',maxval(rBuff(:,5))
    ! write(fLog%unit,*), 'cDom%MshNumLoc(3)=',cDom%MshNumLoc(3)
    ! write(fLog%unit,*), 'minval(gt)=',minval(gt)
    ! write(fLog%unit,*), 'maxval(gt)=',maxval(gt)
    ! call fLog%close()
    
    ! correct centroid
    forall(i=1:npt, gt(i) <= cDom%MshNumLoc(3) .and. gt(i) > 0)
      forall(j=1:4)
        pData(i,j ) = pData(i,j ) -rBuff(gt(i),j) +MeanSlice(gt(i),j+1)
      end forall
        pData(i,g_) = pData(i,g_) -rBuff(gt(i),5) +MeanSlice(gt(i),6  )
    end forall
  end associate
  call tReorder(self)
end subroutine correct_centroid_slices_5D


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMM   add_shotnoise  MMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine add_shotnoise(self,shotID,distID,args)
!=================================================================================================
!!!!!!  add_shotnoise
! input
!   shotID: choice of shot-model method
!           1 : perturbation on theta (Fawley's method)
!           2 : perturbation on charge(uniform average over a macro-particle)
!           3 : perturbation on charge(average e-number is weighted among 
!               the micro-particles cosisting a macro-particle)
! optional input
!    distID : distribution ID
!    args   : arguments of distribution parameters    
!-------------------------------------------------------------------------------------------------
  implicit none
  class(eBeam) :: self
  integer, intent(in) :: shotID
  integer, optional, intent(in) :: distID
  real*8,  optional, intent(in) :: args(:)
  select case (shotID)
    case (1)
      call add_shotnoise_position(self)
    case (3)
      if(present(distID))  then
        call add_shotnoise_charge_weighted(self,distID,args)
      else
        print*, 'Error <- eBeam%add_shotnoise :: distID is required for shotID=3. We change it to shotID=2 and continue to run'
        call add_shotnoise_charge(self)
      endif
    case default
      call add_shotnoise_charge(self)
  end select
end subroutine add_shotnoise

subroutine add_shotnoise_position(self)
!=================================================================================================
!  shot-model by adding random temporal perturbation
!-------------------------------------------------------------------------------------------------
  implicit none
  type(eBeam) :: self
  integer :: i,j,k
  real*8 :: ne, dtheta(self%nBin)
  real*8, dimension(ceiling(self%nBin/2.0)) :: a,phi
  
  associate(pData=>self%pData, nBin=>self%nBin, npt=>self%npt)
    do i=1,npt/nBin
      ne = pData(i*nBin,q_)*nBin/PhysConst%e
      do j=1,ceiling(nBin/2.0)
        call math%rande(a(j),sqrt(2d0/ne)/j)
      enddo
      call random_number(phi)
      phi = (phi-0.5d0)*twopi
      dtheta = 0d0
      do j=1,nBin
        do k=1,ceiling(nBin/2.0)
          dtheta(j) = dtheta(j) + a(k)*cos(k*pData((i-1)*nBin+j,t_)+phi(k))
        enddo
      enddo
      pData((i-1)*nBin+1:i*nBin,t_) = pData((i-1)*nBin+1:i*nBin,t_) + dtheta
    enddo
  end associate
end subroutine add_shotnoise_position

subroutine add_shotnoise_charge(self)
!=================================================================================================
!  shot-model by adding random charge perturbation
!-------------------------------------------------------------------------------------------------
  implicit none
  type(eBeam) :: self
  integer :: i,j
  real*8 :: q
  associate(pData=>self%pData, nBin=>self%nBin, npt=>self%npt)
    do i=1,npt/nBin
      do j=1,nBin
        q = pData((i-1)*nBin+j,q_)
        pData((i-1)*nBin+j,q_) = q  + math%randn()*sqrt(PhysConst%e*q)
      enddo
    enddo
  end associate
end subroutine add_shotnoise_charge

subroutine add_shotnoise_charge_weighted(self,distID,args)
!=================================================================================================
!  shot-model by adding weighted random charge perturbation
!  the weight is defined by particle distribution function -> may include CSR effect 
!-------------------------------------------------------------------------------------------------
  implicit none
  type(eBeam) :: self
  integer, intent(in) :: distID
  real*8,  intent(in) :: args(:)
  integer :: i,j
  real*8 :: q, w(self%nBin)
  
  associate(pData=>self%pData, nBin=>self%nBin, npt=>self%npt)
    do i=1,npt/nBin
      call math%get_density_at(w,pData((i-1)*nBin+1:i*nBin,t_),nBin,distID,args)
      w = nBin*w/sum(w)
      do j=1,nBin
        q = pData((i-1)*nBin+j,q_)*w(j)
        pData((i-1)*nBin+j,q_) = q  + math%randn()*sqrt(PhysConst%e*q)
      enddo
    enddo
  end associate
  
end subroutine add_shotnoise_charge_weighted


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMM   release particles bound to beamlet  MMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine release_beamlet(Beam1,Beam2,limit)
!==============================================================
! release micro-particles from the bound to the corresponding macro-particle.
!   Beam1 : Beam composed of beamlets. i.e., Beam1%nBin > 1
!   Beam2 : Beam composed of micro-particles. i.e., Beam1%nBin = 1
!   ** Beam1%npt+Beam2%npt is conserved
!--------------------------------------------------------------
  implicit none
  class(eBeam):: Beam1
  type(eBeam) :: Beam2
  real*8, intent(in) :: limit
  real*8, allocatable :: pDataTmp(:,:)
  integer :: i,j,nMove
  logical :: mask(Beam1%npt/Beam1%nBin)
  
  if(Beam1%npt==0) return
  
  do i=1,Beam1%npt/Beam1%nBin
    mask(i) = limit <= math%std(Beam1%pData((i-1)*Beam1%nBin+1:i*Beam1%nBin,g_),Beam1%nBin)
  enddo  
  nMove = count(mask)*Beam1%nBin
  if(nMove==0) return
  if(Beam2%npt>0) then
    allocate(pDataTmp(Beam2%npt,7))
    pDataTmp = Beam2%pData
  endif
  if(allocated(Beam2%pData)) deallocate(Beam2%pData)
  allocate(Beam2%pData(Beam2%npt+nMove,7))
  if(Beam2%npt>0) Beam2%pData(1:Beam2%npt,1:7) = pDataTmp
  j=1
  do i=1,Beam1%npt/Beam1%nBin
    if(mask(i)) then
      Beam2%pData(Beam2%npt+(j-1)*Beam1%nBin+1:Beam2%npt+j*Beam1%nBin,:) = Beam1%pData((i-1)*Beam1%nBin+1:i*Beam1%nBin,:)
      j=j+1
    endif
  enddo
  if(allocated(pDataTmp)) deallocate(pDataTmp)
  allocate(pDataTmp(Beam1%npt,7))
  pDataTmp = Beam1%pData
  deallocate(Beam1%pData)
  allocate(Beam1%pData(Beam1%npt-nMove,7))
  j=0
  do i=1,Beam1%npt/Beam1%nBin
    if(.not. mask(i)) then
      j=j+1
      Beam1%pData((j-1)*Beam1%nBin+1:j*Beam1%nBin,:) = pDataTmp((i-1)*Beam1%nBin+1:i*Beam1%nBin,:)
    endif
  enddo
  Beam1%npt = Beam1%npt -nMove
  Beam2%npt = Beam2%npt +nMove
  
end subroutine release_beamlet


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMM   Reorder i.e., domain decomposiion  MMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine yReorder2adj(self,nSend)
! =============================================================================
! move outside-of-domain-particles to adjacent column processors 
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(eBeam) :: self
  integer, intent(out) :: nSend
  real*8, allocatable, dimension(:) :: sUpBuff,sDnBuff,rUpBuff,rDnBuff
  real*8, allocatable, dimension(:,:) :: pDataTemp
  
  
  integer :: i, myCol, iTemp, ierr, req(4)
  integer :: iSendUp,iSendDn
  integer :: nSendUp,nSendDn,nRecvUp,nRecvDn,nRecv
  logical, dimension(self%npt) :: isUp, isDn
  
  myCol = self%pGrd%myCol
  isUp = self%pData(:,y_) > self%cDom%RngLoc(2,2)
  isDn = self%pData(:,y_) < self%cDom%RngLoc(2,1)
  if(myCol == 0)            isDn = .false.
  if(myCol == self%pGrd%npcol-1) isUp = .false.

  nSendUp = count(isUp)
  nSendDn = count(isDn)
  nRecvUp=0
  nRecvDn=0
  ! -- send/recv size of buff --
  call MPI_IRECV(nRecvUp,1,MPI_integer,self%pGrd%myU,0,self%pGrd%comm_col,req(1),ierr)
  call MPI_IRECV(nRecvDn,1,MPI_integer,self%pGrd%myD,0,self%pGrd%comm_col,req(2),ierr)
  call MPI_ISEND(nSendUp,1,MPI_integer,self%pGrd%myU,0,self%pGrd%comm_col,req(3),ierr)
  call MPI_ISEND(nSendDn,1,MPI_integer,self%pGrd%myD,0,self%pGrd%comm_col,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  allocate(sUpBuff(7*nSendUp))
  allocate(sDnBuff(7*nSendDn))
  allocate(rUpBuff(7*nRecvUp))
  allocate(rDnBuff(7*nRecvDn))
  nSend = nSendUp+nSendDn
  nRecv = nRecvUp+nRecvDn
  
  ! -- prepare sBuff and save self%pData --
  sUpBuff = 0d0
  sDnBuff = 0d0
  rUpBuff = 0d0
  rDnBuff = 0d0
  allocate(pDataTemp(self%npt-nSend,7))
  pDataTemp = 0d0
  
  iSendUp = 0
  iSendDn = 0
  iTemp = 1
  do i=1,self%npt
    if(isUp(i)) then
      sUpBuff(iSendUp+1:iSendUp+7) = self%pData(i,1:7)
      iSendUp = iSendUp + 7
    elseif(isDn(i)) then
      sDnBuff(iSendDn+1:iSendDn+7) = self%pData(i,1:7)
      iSendDn = iSendDn + 7
    else
      pDataTemp(iTemp,1:7) = self%pData(i,1:7)
      iTemp = iTemp +1
    endif
  enddo

  call MPI_IRECV(rUpBuff,nRecvUp*7,MPI_DOUBLE_PRECISION,self%pGrd%myU,0,self%pGrd%comm_col,req(1),ierr)
  call MPI_ISEND(sDnBuff,nSendDn*7,MPI_DOUBLE_PRECISION,self%pGrd%myD,0,self%pGrd%comm_col,req(2),ierr)
  call MPI_IRECV(rDnBuff,nRecvDn*7,MPI_DOUBLE_PRECISION,self%pGrd%myD,0,self%pGrd%comm_col,req(3),ierr)
  call MPI_ISEND(sUpBuff,nSendUp*7,MPI_DOUBLE_PRECISION,self%pGrd%myU,0,self%pGrd%comm_col,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  if(nSend==0 .and. nRecv==0) return
  deallocate(self%pData)
  allocate(self%pData(self%npt-nSend+nRecv,7))
  if(self%npt-nSend/=0) self%pData(1:self%npt-nSend,1:7) = pDataTemp(1:self%npt-nSend,1:7)
  do i=1,nRecvUp
    self%pData(self%npt-nSend+i,1:7) = rUpBuff((i-1)*7+1:i*7)
  enddo
  do i=1,nRecvDn
    self%pData(self%npt-nSend+nRecvUp+i,1:7) = rDnBuff((i-1)*7+1:i*7)
  enddo
  self%npt = self%npt-nSend+nRecv
end subroutine yReorder2adj

subroutine yReorder2adjBin(self,nSend)
! =============================================================================
! move outside-of-domain-beamlet to adjacent column processors 
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(eBeam) :: self
  integer, intent(out) :: nSend
  real*8,allocatable, dimension(:) :: sUpBuff,sDnBuff,rUpBuff,rDnBuff
  real*8,allocatable, dimension(:,:) :: pDataTemp
  
  integer :: i, j, myCol, iTemp, ierr, req(4)
  integer :: iSendUp,iSendDn
  integer :: nSendUp,nSendDn,nRecvUp,nRecvDn,nRecv
  real*8,  dimension(self%npt/self%nBin) :: yBeamLet
  logical, dimension(self%npt/self%nBin) :: isUp, isDn
  
  myCol = self%pGrd%myCol
  do i=1,self%npt/self%nBin
    yBeamLet(i) = sum(self%pData((i-1)*self%nBin+1:i*self%nBin,y_))/self%nBin
  end do
  isUp = yBeamLet > self%cDom%RngLoc(2,2)
  isDn = yBeamLet < self%cDom%RngLoc(2,1)
  if(myCol == 0)            isDn = .false.
  if(myCol == self%pGrd%npcol-1) isUp = .false.

  nSendUp = count(isUp)*self%nBin
  nSendDn = count(isDn)*self%nBin
  nRecvUp=0
  nRecvDn=0
  call MPI_IRECV(nRecvUp,1,MPI_integer,self%pGrd%myU,0,self%pGrd%comm_col,req(1),ierr)
  call MPI_IRECV(nRecvDn,1,MPI_integer,self%pGrd%myD,0,self%pGrd%comm_col,req(2),ierr)
  call MPI_ISEND(nSendUp,1,MPI_integer,self%pGrd%myU,0,self%pGrd%comm_col,req(3),ierr)
  call MPI_ISEND(nSendDn,1,MPI_integer,self%pGrd%myD,0,self%pGrd%comm_col,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  allocate(sUpBuff(7*nSendUp))
  allocate(sDnBuff(7*nSendDn))
  allocate(rUpBuff(7*nRecvUp))
  allocate(rDnBuff(7*nRecvDn))
  nSend = nSendUp+nSendDn
  nRecv = nRecvUp+nRecvDn
  sUpBuff = 0d0
  sDnBuff = 0d0
  rUpBuff = 0d0
  rDnBuff = 0d0
  allocate(pDataTemp(self%npt-nSend,7))
  pDataTemp = 0d0
  
  iSendUp = 0
  iSendDn = 0
  iTemp = 0
  do i=1,self%npt/self%nBin
    if(isUp(i)) then
      do j=1,self%nBin
        sUpBuff(iSendUp+1:iSendUp+7) = self%pData((i-1)*self%nBin+j,1:7)
        iSendUp = iSendUp + 7
      enddo
    elseif(isDn(i)) then
      do j=1,self%nBin
        sDnBuff(iSendDn+1:iSendDn+7) = self%pData((i-1)*self%nBin+j,1:7)
        iSendDn = iSendDn + 7
      enddo
    else
      pDataTemp(iTemp+1:iTemp+self%nBin,1:7) = self%pData((i-1)*self%nBin+1:i*self%nBin,1:7)
      iTemp = iTemp +self%nBin
    endif
  enddo
  call MPI_IRECV(rUpBuff,nRecvUp*7,MPI_DOUBLE_PRECISION,self%pGrd%myU,0,self%pGrd%comm_col,req(1),ierr)
  call MPI_ISEND(sDnBuff,nSendDn*7,MPI_DOUBLE_PRECISION,self%pGrd%myD,0,self%pGrd%comm_col,req(2),ierr)
  call MPI_IRECV(rDnBuff,nRecvDn*7,MPI_DOUBLE_PRECISION,self%pGrd%myD,0,self%pGrd%comm_col,req(3),ierr)
  call MPI_ISEND(sUpBuff,nSendUp*7,MPI_DOUBLE_PRECISION,self%pGrd%myU,0,self%pGrd%comm_col,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  if(nSend==0 .and. nRecv==0) return
  deallocate(self%pData)
  allocate(self%pData(self%npt-nSend+nRecv,7))
  if(self%npt-nSend/=0) self%pData(1:self%npt-nSend,1:7) = pDataTemp(1:self%npt-nSend,1:7)
  do i=1,nRecvUp
    self%pData(self%npt-nSend+i,1:7) = rUpBuff((i-1)*7+1:i*7)
  enddo
  do i=1,nRecvDn
    self%pData(self%npt-nSend+nRecvUp+i,1:7) = rDnBuff((i-1)*7+1:i*7)
  enddo
  self%npt = self%npt-nSend+nRecv
end subroutine yReorder2adjBin

subroutine yReorder(self)
! =============================================================================
! move outside-of-domain-particles to corresponding column processors
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  integer :: nSend, nSendTot,ierr

  call MPI_BARRIER(self%pGrd%comm_row,ierr)
  if(self%nBin>1) then
    do
      call yReorder2adjBin(self,nSend)
      call MPI_ALLREDUCE(nSend,nSendTot,1,MPI_INT,MPI_SUM,self%pGrd%comm_col,ierr)
      if(nSendTot==0) exit
    enddo
  else
    do
      call yReorder2adj(self,nSend)
      call MPI_ALLREDUCE(nSend,nSendTot,1,MPI_INT,MPI_SUM,self%pGrd%comm_col,ierr)
      if(nSendTot==0) exit
    enddo
  endif
end subroutine yReorder

subroutine tReorder2adj(self,nSend)
! =============================================================================
! move outside-of-domain-particles to adjacent column processors 
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(eBeam) :: self
  integer, intent(out) :: nSend
  real*8, allocatable, dimension(:) :: sRigtBuff,sLeftBuff,rRigtBuff,rLeftBuff
  real*8, allocatable, dimension(:,:) :: pDataTemp
  
  
  integer :: i, myRow, iTemp, ierr, req(4)
  integer :: iSendRigt,iSendLeft
  integer :: nSendRigt,nSendLeft,nRecvRigt,nRecvLeft,nRecv
  logical, dimension(self%npt) :: isRigt, isLeft
  
  myRow = self%pGrd%myRow
  isRigt = self%pData(:,t_) > self%cDom%RngLoc(3,2)
  isLeft = self%pData(:,t_) < self%cDom%RngLoc(3,1)
  if(myRow == 0)            isLeft = .false.
  if(myRow == self%pGrd%npRow-1) isRigt = .false.

  nSendRigt = count(isRigt)
  nSendLeft = count(isLeft)
  nRecvRigt=0
  nRecvLeft=0
  ! -- send/recv size of buff --
  call MPI_IRECV(nRecvRigt,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_IRECV(nRecvLeft,1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_ISEND(nSendRigt,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(nSendLeft,1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  allocate(sRigtBuff(7*nSendRigt))
  allocate(sLeftBuff(7*nSendLeft))
  allocate(rRigtBuff(7*nRecvRigt))
  allocate(rLeftBuff(7*nRecvLeft))
  nSend = nSendRigt+nSendLeft
  nRecv = nRecvRigt+nRecvLeft
  
  ! -- prepare sBuff and save self%pData --
  sRigtBuff = 0d0
  sLeftBuff = 0d0
  rRigtBuff = 0d0
  rLeftBuff = 0d0
  allocate(pDataTemp(self%npt-nSend,7))
  pDataTemp = 0d0
  
  iSendRigt = 0
  iSendLeft = 0
  iTemp = 1
  do i=1,self%npt
    if(isRigt(i)) then
      sRigtBuff(iSendRigt+1:iSendRigt+7) = self%pData(i,1:7)
      iSendRigt = iSendRigt + 7
    elseif(isLeft(i)) then
      sLeftBuff(iSendLeft+1:iSendLeft+7) = self%pData(i,1:7)
      iSendLeft = iSendLeft + 7
    else
      pDataTemp(iTemp,1:7) = self%pData(i,1:7)
      iTemp = iTemp +1
    endif
  enddo
  call MPI_IRECV(rRigtBuff,nRecvRigt*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_ISEND(sLeftBuff,nSendLeft*7,MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_IRECV(rLeftBuff,nRecvLeft*7,MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(sRigtBuff,nSendRigt*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  if(nSend==0 .and. nRecv==0) return
  deallocate(self%pData)
  allocate(self%pData(self%npt-nSend+nRecv,7))
  if(self%npt-nSend/=0) self%pData(1:self%npt-nSend,1:7) = pDataTemp(1:self%npt-nSend,1:7)
  do i=1,nRecvRigt
    self%pData(self%npt-nSend+i,1:7) = rRigtBuff((i-1)*7+1:i*7)
  enddo
  do i=1,nRecvLeft
    self%pData(self%npt-nSend+nRecvRigt+i,1:7) = rLeftBuff((i-1)*7+1:i*7)
  enddo
  self%npt = self%npt-nSend+nRecv
end subroutine tReorder2adj

subroutine tReorder2adjBin(self,nSend)
! =============================================================================
! move outside-of-domain-beamlet to adjacent column processors 
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(eBeam) :: self
  integer, intent(out) :: nSend
  real*8, allocatable, dimension(:) :: sRightBuff,sLeftBuff,rRightBuff,rLeftBuff
  real*8, allocatable, dimension(:,:) :: pDataTemp
  
  
  integer :: i, j, myRow, iTemp, ierr, req(4)
  integer :: iSendRight,iSendLeft
  integer :: nSendRight,nSendLeft,nRecvRight,nRecvLeft,nRecv
  real*8,  dimension(self%npt/self%nBin) :: tBeamLet
  logical, dimension(self%npt/self%nBin) :: isRight, isLeft
  
  myRow = self%pGrd%myRow
  do i=1,self%npt/self%nBin
    tBeamLet(i) = sum(self%pData((i-1)*self%nBin+1:i*self%nBin,t_))/self%nBin
  end do
  isRight = tBeamLet > self%cDom%RngLoc(3,2)
  isLeft = tBeamLet < self%cDom%RngLoc(3,1)
  if(myRow == 0)            isLeft  = .false.
  if(myRow == self%pGrd%npRow-1) isRight = .false.

  nSendRight = count(isRight)*self%nBin
  nSendLeft  = count(isLeft)*self%nBin
  nRecvRight = 0
  nRecvLeft  = 0
  !call write_log('comm pNum')
  call MPI_IRECV(nRecvRight,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_IRECV(nRecvLeft, 1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_ISEND(nSendRight,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(nSendLeft, 1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  allocate(sRightBuff(7*nSendRight))
  allocate(sLeftBuff (7*nSendLeft ))
  allocate(rRightBuff(7*nRecvRight))
  allocate(rLeftBuff (7*nRecvLeft ))
  nSend = nSendRight+nSendLeft
  nRecv = nRecvRight+nRecvLeft
  sRightBuff = 0d0
  sLeftBuff  = 0d0
  rRightBuff = 0d0
  rLeftBuff  = 0d0
  allocate(pDataTemp(self%npt-nSend,7))
  pDataTemp = 0d0
  
  iSendRight= 0
  iSendLeft = 0
  iTemp = 0

  do i=1,self%npt/self%nBin
    if(isRight(i)) then
      do j=1,self%nBin
        sRightBuff(iSendRight+1:iSendRight+7) = self%pData((i-1)*self%nBin+j,1:7)
        iSendRight = iSendRight +7
      enddo
    elseif(isLeft(i)) then
      do j=1,self%nBin
        sLeftBuff(iSendLeft+1:iSendLeft+7) = self%pData((i-1)*self%nBin+j,1:7)
        iSendLeft = iSendLeft +7
      enddo
    else
      pDataTemp(iTemp+1:iTemp+self%nBin,1:7) = self%pData((i-1)*self%nBin+1:i*self%nBin,1:7)
      iTemp = iTemp +self%nBin
    endif
  enddo
  
  call MPI_IRECV(rRightBuff,nRecvRight*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_ISEND(sLeftBuff, nSendLeft*7, MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_IRECV(rLeftBuff, nRecvLeft*7, MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(sRightBuff,nSendRight*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  if(nSend==0 .and. nRecv==0) return
  deallocate(self%pData)
  allocate(self%pData(self%npt-nSend+nRecv,7))
  if(self%npt-nSend/=0) self%pData(1:self%npt-nSend,1:7) = pDataTemp(1:self%npt-nSend,1:7)
  do i=1,nRecvRight
    self%pData(self%npt-nSend+i,1:7) = rRightBuff((i-1)*7+1:i*7)
  enddo
  do i=1,nRecvLeft
    self%pData(self%npt-nSend+nRecvRight+i,1:7) = rLeftBuff((i-1)*7+1:i*7)
  enddo
  self%npt = self%npt-nSend+nRecv
end subroutine tReorder2adjBin

subroutine tReorder2adjBinArr(self,nSend,Arr)
! =============================================================================
! move outside-of-domain-beamlet to adjacent column processors 
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(eBeam) :: self
  integer, intent(out) :: nSend
  complex(8), allocatable, intent(inout) :: Arr(:,:)
  complex(8), allocatable, dimension(:)  :: sRightBuffArr,sLeftBuffArr,rRightBuffArr,rLeftBuffArr
  complex(8), allocatable, dimension(:,:):: ArrTmp
  real*8, allocatable, dimension(:) :: sRightBuff,sLeftBuff,rRightBuff,rLeftBuff
  real*8, allocatable, dimension(:,:) :: pDataTemp
  
  
  integer :: i, j, myRow, iTemp, ierr, req(4), N
  integer :: iSendRight,iSendLeft
  integer :: nSendRight,nSendLeft,nRecvRight,nRecvLeft,nRecv
  real*8,  dimension(self%npt/self%nBin) :: tBeamLet
  logical, dimension(self%npt/self%nBin) :: isRight, isLeft
  
  if(size(Arr,1) /= self%npt/self%nBin) then
    print*, 'Error :: eBeam%tReorder2adjBinArr -> wrong Arr shape. Program exit'
    call exit(1)
  endif
  N = size(Arr,2)  
  
  myRow = self%pGrd%myRow
  do i=1,self%npt/self%nBin
    tBeamLet(i) = sum(self%pData((i-1)*self%nBin+1:i*self%nBin,t_))/self%nBin
  end do
  isRight = tBeamLet > self%cDom%RngLoc(3,2)
  isLeft = tBeamLet < self%cDom%RngLoc(3,1)
  if(myRow == 0)                 isLeft  = .false.
  if(myRow == self%pGrd%npRow-1) isRight = .false.

  nSendRight = count(isRight)*self%nBin
  nSendLeft  = count(isLeft)*self%nBin
  nRecvRight = 0
  nRecvLeft  = 0
  !call write_log('comm pNum')
  call MPI_IRECV(nRecvRight,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_IRECV(nRecvLeft, 1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_ISEND(nSendRight,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(nSendLeft, 1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  allocate(sRightBuff(7*nSendRight))
  allocate(sLeftBuff (7*nSendLeft ))
  allocate(rRightBuff(7*nRecvRight))
  allocate(rLeftBuff (7*nRecvLeft ))
  allocate(sRightBuffArr(N*nSendRight/self%nBin))
  allocate(sLeftBuffArr (N*nSendLeft /self%nBin))
  allocate(rRightBuffArr(N*nRecvRight/self%nBin))
  allocate(rLeftBuffArr (N*nRecvLeft /self%nBin))
   
  nSend = nSendRight+nSendLeft
  nRecv = nRecvRight+nRecvLeft
  sRightBuff = 0d0
  sLeftBuff  = 0d0
  rRightBuff = 0d0
  rLeftBuff  = 0d0
  sRightBuffArr = (0d0,0d0)
  sLeftBuffArr  = (0d0,0d0)
  rRightBuffArr = (0d0,0d0)
  rLeftBuffArr  = (0d0,0d0)
  allocate(pDataTemp(self%npt-nSend,7))
  pDataTemp = 0d0
  allocate(ArrTmp((self%npt-nSend)/self%nBin,N))
  ArrTmp = (0d0,0d0)
  
  iSendRight= 0
  iSendLeft = 0
  iTemp = 0
  do i=1,self%npt/self%nBin
    if(isRight(i)) then
      sRightBuffArr(iSendRight/(7*self%nBin)*N+1:iSendRight/(7*self%nBin)*N+N) = Arr(i,1:N)
      do j=1,self%nBin
        sRightBuff(iSendRight+1:iSendRight+7) = self%pData((i-1)*self%nBin+j,1:7)
        iSendRight = iSendRight +7
      enddo
    elseif(isLeft(i)) then
      sLeftBuffArr(iSendLeft/(7*self%nBin)*N+1:iSendLeft/(7*self%nBin)*N+N) = Arr(i,1:N)
      do j=1,self%nBin
        sLeftBuff(iSendLeft+1:iSendLeft+7) = self%pData((i-1)*self%nBin+j,1:7)
        iSendLeft = iSendLeft +7
      enddo
    else
      ArrTmp(iTemp/self%nBin+1,1:N) = Arr(i,1:N)
      pDataTemp(iTemp+1:iTemp+self%nBin,1:7) = self%pData((i-1)*self%nBin+1:i*self%nBin,1:7)
      iTemp = iTemp +self%nBin
    endif
  enddo
  
  call MPI_IRECV(rRightBuff,nRecvRight*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_ISEND(sLeftBuff, nSendLeft*7, MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_IRECV(rLeftBuff, nRecvLeft*7, MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(sRightBuff,nSendRight*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  
  call MPI_IRECV(rRightBuffArr,nRecvRight/self%nBin*N,MPI_DOUBLE_COMPLEX,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_ISEND(sLeftBuffArr, nSendLeft/self%nBin *N,MPI_DOUBLE_COMPLEX,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_IRECV(rLeftBuffArr, nRecvLeft/self%nBin *N,MPI_DOUBLE_COMPLEX,self%pGrd%myL,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(sRightBuffArr,nSendRight/self%nBin*N,MPI_DOUBLE_COMPLEX,self%pGrd%myR,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  
  if(nSend==0 .and. nRecv==0) return
  deallocate(self%pData, Arr)
  allocate(self%pData(self%npt-nSend+nRecv,7))
  allocate(Arr((self%npt-nSend+nRecv)/self%nBin,N))
  if(self%npt-nSend/=0) then
    self%pData(1:self%npt-nSend,1:7) = pDataTemp(1:self%npt-nSend,1:7)
    Arr(1:(self%npt-nSend)/self%nBin,1:N) = ArrTmp(1:(self%npt-nSend)/self%nBin,1:N)
  endif
  do i=1,nRecvRight
    self%pData(self%npt-nSend+i,1:7) = rRightBuff((i-1)*7+1:i*7)
  enddo
  do i=1,nRecvRight/self%nBin
    Arr((self%npt-nSend)/self%nBin+i,1:N) = rRightBuffArr((i-1)*N+1:i*N)
  enddo
  do i=1,nRecvLeft
    self%pData(self%npt-nSend+nRecvRight+i,1:7) = rLeftBuff((i-1)*7+1:i*7)
  enddo
  do i=1,nRecvLeft/self%nBin
    Arr((self%npt-nSend+nRecvRight)/self%nBin+i,1:N) = rRightBuffArr((i-1)*N+1:i*N)
  enddo
  self%npt = self%npt-nSend+nRecv
end subroutine tReorder2adjBinArr



subroutine tReorder(self,Arr)
! =============================================================================
! move outside-of-domain-particles to corresponding row processors
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  complex(8), allocatable, optional, intent(inout) :: Arr(:,:)
  integer :: nSend, nSendTot,ierr

  call MPI_BARRIER(self%pGrd%comm_col,ierr)
  if(present(Arr)) then
    do
      call tReorder2adjBinArr(self,nSend,Arr)
      call MPI_ALLREDUCE(nSend,nSendTot,1,MPI_INT,MPI_SUM,self%pGrd%comm_row,ierr)
      if(nSendTot==0) exit
    enddo
  elseif(self%nBin>1) then
    do
      call tReorder2adjBin(self,nSend)
      call MPI_ALLREDUCE(nSend,nSendTot,1,MPI_INT,MPI_SUM,self%pGrd%comm_row,ierr)
      if(nSendTot==0) exit
    enddo
  else
    do
      call tReorder2adj(self,nSend)
      call MPI_ALLREDUCE(nSend,nSendTot,1,MPI_INT,MPI_SUM,self%pGrd%comm_row,ierr)
      if(nSendTot==0) exit
    enddo
  endif
end subroutine tReorder

subroutine Reorder(self)
! =============================================================================
! move outside-of-domain-particles to corresponding row processors
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  integer :: ierr
  call self%yReorder()
  call self%tReorder()
  call MPI_BARRIER(self%pGrd%comm_2d,ierr)
end subroutine Reorder

subroutine tReorder2adj_middle(self,nSend)
! =============================================================================
! move outside-of-domain-particles to adjacent column(theta domain) processors 
! the domain used here is re-defined such that the mesh points are 
! not on domain boundary but at the middle of boundary
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(eBeam) :: self
  integer, intent(out) :: nSend
  real*8, allocatable, dimension(:) :: sRigtBuff,sLeftBuff,rRigtBuff,rLeftBuff
  real*8, allocatable, dimension(:,:) :: pDataTemp
  
  
  integer :: i, myRow, iTemp, ierr, req(4)
  integer :: iSendRigt,iSendLeft
  integer :: nSendRigt,nSendLeft,nRecvRigt,nRecvLeft,nRecv
  logical, dimension(self%npt) :: isRigt, isLeft
  
  myRow = self%pGrd%myRow
  isRigt = self%pData(:,t_) > self%cDom%RngLoc(3,2) - 0.5d0*self%cDom%MshSize(3)
  isLeft = self%pData(:,t_) < self%cDom%RngLoc(3,1) - 0.5d0*self%cDom%MshSize(3)
  if(myRow == 0)            isLeft = .false.
  if(myRow == self%pGrd%npRow-1) isRigt = .false.

  nSendRigt = count(isRigt)
  nSendLeft = count(isLeft)
  nRecvRigt=0
  nRecvLeft=0
  ! -- send/recv size of buff --
  call MPI_IRECV(nRecvRigt,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_IRECV(nRecvLeft,1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_ISEND(nSendRigt,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(nSendLeft,1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  allocate(sRigtBuff(7*nSendRigt))
  allocate(sLeftBuff(7*nSendLeft))
  allocate(rRigtBuff(7*nRecvRigt))
  allocate(rLeftBuff(7*nRecvLeft))
  nSend = nSendRigt+nSendLeft
  nRecv = nRecvRigt+nRecvLeft
  
  ! -- prepare sBuff and save self%pData --
  sRigtBuff = 0d0
  sLeftBuff = 0d0
  rRigtBuff = 0d0
  rLeftBuff = 0d0
  allocate(pDataTemp(self%npt-nSend,7))
  pDataTemp = 0d0
  
  iSendRigt = 0
  iSendLeft = 0
  iTemp = 1
  do i=1,self%npt
    if(isRigt(i)) then
      sRigtBuff(iSendRigt+1:iSendRigt+7) = self%pData(i,1:7)
      iSendRigt = iSendRigt + 7
    elseif(isLeft(i)) then
      sLeftBuff(iSendLeft+1:iSendLeft+7) = self%pData(i,1:7)
      iSendLeft = iSendLeft + 7
    else
      pDataTemp(iTemp,1:7) = self%pData(i,1:7)
      iTemp = iTemp +1
    endif
  enddo
  call MPI_IRECV(rRigtBuff,nRecvRigt*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_ISEND(sLeftBuff,nSendLeft*7,MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_IRECV(rLeftBuff,nRecvLeft*7,MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(sRigtBuff,nSendRigt*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  if(nSend==0 .and. nRecv==0) return
  deallocate(self%pData)
  allocate(self%pData(self%npt-nSend+nRecv,7))
  if(self%npt-nSend/=0) self%pData(1:self%npt-nSend,1:7) = pDataTemp(1:self%npt-nSend,1:7)
  do i=1,nRecvRigt
    self%pData(self%npt-nSend+i,1:7) = rRigtBuff((i-1)*7+1:i*7)
  enddo
  do i=1,nRecvLeft
    self%pData(self%npt-nSend+nRecvRigt+i,1:7) = rLeftBuff((i-1)*7+1:i*7)
  enddo
  self%npt = self%npt-nSend+nRecv
end subroutine tReorder2adj_middle

subroutine tReorder2adjBin_middle(self,nSend)
! =============================================================================
! move outside-of-domain-beamlet to adjacent column processors 
! the domain used here is re-defined such that domain boundary is 
! not passing the mesh points but between the mesh point.
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(eBeam) :: self
  integer, intent(out) :: nSend
  real*8, allocatable, dimension(:) :: sRightBuff,sLeftBuff,rRightBuff,rLeftBuff
  real*8, allocatable, dimension(:,:) :: pDataTemp
  
  integer :: i, j, myRow, iTemp, ierr, req(4)
  integer :: iSendRight,iSendLeft
  integer :: nSendRight,nSendLeft,nRecvRight,nRecvLeft,nRecv
  real*8,  dimension(self%npt/self%nBin) :: tBeamLet
  logical, dimension(self%npt/self%nBin) :: isRight, isLeft
  
  !call write_log('tReorder2adjBin start!')
  myRow = self%pGrd%myRow
  do i=1,self%npt/self%nBin
    tBeamLet(i) = sum(self%pData((i-1)*self%nBin+1:i*self%nBin,t_))/self%nBin
  end do
  isRight = tBeamLet > self%cDom%RngLoc(3,2) - 0.5d0*self%cDom%MshSize(3)
  isLeft  = tBeamLet < self%cDom%RngLoc(3,1) - 0.5d0*self%cDom%MshSize(3)
  if(myRow == 0)            isLeft  = .false.
  if(myRow == self%pGrd%npRow-1) isRight = .false.

  !call write_log('self%npt,self%nBin,self%npt/self%nBin'//num2str(self%npt)//num2str(self%nBin)//num2str(self%npt/self%nBin))
  nSendRight = count(isRight)*self%nBin
  nSendLeft  = count(isLeft)*self%nBin
  nRecvRight = 0
  nRecvLeft  = 0
  !call write_log('comm pNum')
  call MPI_IRECV(nRecvRight,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_IRECV(nRecvLeft, 1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_ISEND(nSendRight,1,MPI_integer,self%pGrd%myR,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(nSendLeft, 1,MPI_integer,self%pGrd%myL,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  allocate(sRightBuff(7*nSendRight))
  allocate(sLeftBuff (7*nSendLeft ))
  allocate(rRightBuff(7*nRecvRight))
  allocate(rLeftBuff (7*nRecvLeft ))
  nSend = nSendRight+nSendLeft
  nRecv = nRecvRight+nRecvLeft
  sRightBuff = 0d0
  sLeftBuff  = 0d0
  rRightBuff = 0d0
  rLeftBuff  = 0d0
  allocate(pDataTemp(self%npt-nSend,7))
  pDataTemp = 0d0
  
  iSendRight= 0
  iSendLeft = 0
  iTemp = 0
  do i=1,self%npt/self%nBin
    if(isRight(i)) then
      do j=1,self%nBin
        sRightBuff(iSendRight+1:iSendRight+7) = self%pData((i-1)*self%nBin+j,1:7)
        iSendRight = iSendRight +7
      enddo
    elseif(isLeft(i)) then
      do j=1,self%nBin
        sLeftBuff(iSendLeft+1:iSendLeft+7) = self%pData((i-1)*self%nBin+j,1:7)
        iSendLeft = iSendLeft +7
      enddo
    else
      pDataTemp(iTemp+1:iTemp+self%nBin,1:7) = self%pData((i-1)*self%nBin+1:i*self%nBin,1:7)
      iTemp = iTemp +self%nBin
    endif
  enddo
  
  call MPI_IRECV(rRightBuff,nRecvRight*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(1),ierr)
  call MPI_ISEND(sLeftBuff, nSendLeft*7, MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(2),ierr)
  call MPI_IRECV(rLeftBuff, nRecvLeft*7, MPI_DOUBLE_PRECISION,self%pGrd%myL,0,self%pGrd%comm_row,req(3),ierr)
  call MPI_ISEND(sRightBuff,nSendRight*7,MPI_DOUBLE_PRECISION,self%pGrd%myR,0,self%pGrd%comm_row,req(4),ierr)
  call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
  if(nSend==0 .and. nRecv==0) return
  deallocate(self%pData)
  allocate(self%pData(self%npt-nSend+nRecv,7))
  if(self%npt-nSend/=0) self%pData(1:self%npt-nSend,1:7) = pDataTemp(1:self%npt-nSend,1:7)
  do i=1,nRecvRight
    self%pData(self%npt-nSend+i,1:7) = rRightBuff((i-1)*7+1:i*7)
  enddo
  do i=1,nRecvLeft
    self%pData(self%npt-nSend+nRecvRight+i,1:7) = rLeftBuff((i-1)*7+1:i*7)
  enddo
  self%npt = self%npt-nSend+nRecv
end subroutine tReorder2adjBin_middle

subroutine tReorder_middle(self)
! =============================================================================
! move outside-of-domain-particles to corresponding row processors
! -----------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  integer :: nSend, nSendTot,ierr
  
  call MPI_BARRIER(self%pGrd%comm_col,ierr)
  if(self%nBin>1) then
    do
      call tReorder2adjBin_middle(self,nSend)
      call MPI_ALLREDUCE(nSend,nSendTot,1,MPI_INT,MPI_SUM,self%pGrd%comm_row,ierr)
      if(nSendTot==0) exit
    enddo
  else
    do
      call tReorder2adj_middle(self,nSend)
      call MPI_ALLREDUCE(nSend,nSendTot,1,MPI_INT,MPI_SUM,self%pGrd%comm_row,ierr)
      if(nSendTot==0) exit
    enddo
  endif
end subroutine tReorder_middle


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMM   LoadBalance  MMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine yLoadBalance(self,Rad)
  implicit none
  include 'mpif.h'
  class(eBeam)    :: self
  type(Radiation) :: rad
  integer :: i,j,k,n,nptMean
  integer, dimension(self%cDom%MshNum(2)-1) :: nptTab,nptTabLoc
  integer :: yMshTabOld(0:self%pGrd%npCol)
  
  associate(pGrd =>self%pGrd,  cDom=>self%cDom,&
            pData=>self%pData, npt=>self%npt,   nBin=>self%nBin)
    yMshTabOld = cDom%yMshTab
    nptTabLoc = 0 
    nptTab = 0
    do i=1,npt/nBin
      j = int((sum(pData((i-1)*nBin+1:i*nBin,y_))-cDom%RngLoc(2,1))/cDom%MshSize(2))
      nptTabLoc(j+cDom%yMshTab(pGrd%myCol)+1) = nptTabLoc(j+cDom%yMshTab(pGrd%myCol)+1) + nBin
    enddo
    call MPI_AllReduce(nptTabLoc,nptTab,cDom%MshNum(2)-1,MPI_INT,MPI_SUM,pGrd%comm_2d,i)
    nptMean = sum(nptTab)/pGrd%npCol
    n=0
    k=0
    do i=1,pGrd%npCol-2
      do j=k+1,cDom%MshNum(2)-1
        n = n+nptTab(j)
        if(n>=i*nptMean) then
          if (n-i*nptMean < i*nptMean-n+nptTab(j)) then
            k = j
          else 
            k = j-1
            n = n-nptTab(j)
          endif
          if(k == cDom%yMshTab(i-1)) then
            k = k+1
            n = n+nptTab(k)
          endif
          cDom%yMshTab(i) = k
          exit
        endif
      enddo
    enddo
    cDom%RngLoc(2,1) = cDom%Rng(2,1) +cDom%yMshTab(pGrd%mycol  )*cDom%MshSize(2)
    cDom%RngLoc(2,2) = cDom%Rng(2,1) +cDom%yMshTab(pGrd%mycol+1)*cDom%MshSize(2)
    cDom%MshNumLoc(2)= cDom%yMshTab(pGrd%mycol+1) -cDom%yMshTab(pGrd%mycol) +1
  end associate
  call Rad%update_Radiation_yMsh(yMshTabOld)
end subroutine yLoadBalance

subroutine tLoadBalance(self,Rad)
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  type(Radiation) :: Rad
  integer :: i,j,k,n,nptMean
  integer, dimension(self%cDom%MshNum(3)-1) :: nptTab,nptTabLoc
  integer :: tMshTabOld(0:self%pGrd%npRow)

  associate(pGrd =>self%pGrd,  cDom=>self%cDom,&
            pData=>self%pData, npt=>self%npt,   nBin=>self%nBin)
    tMshTabOld = cDom%tMshTab
    nptTabLoc = 0 
    nptTab = 0
    do i=1,npt/nBin
      j = int((sum(pData((i-1)*nBin+1:i*nBin,y_))-cDom%RngLoc(3,1))/cDom%MshSize(3))
      nptTabLoc(j+cDom%tMshTab(pGrd%myRow)+1) = nptTabLoc(j+cDom%tMshTab(pGrd%myRow)+1) + nBin
    enddo
    call MPI_AllReduce(nptTabLoc,nptTab,cDom%MshNum(3)-1,MPI_INT,MPI_SUM,pGrd%comm_2d,i)
    nptMean = sum(nptTab)/pGrd%npRow
    n=0
    k=0
    do i=1,pGrd%npCol-2
      do j=k+1,cDom%MshNum(3)-1
        n = n+nptTab(j)
        if(n>=i*nptMean) then
          if (n-i*nptMean < i*nptMean-n+nptTab(j)) then
            k = j
          else 
            k = j-1
            n = n-nptTab(j)
          endif
          if(k == cDom%tMshTab(i-1)) then
            k = k+1
            n = n+nptTab(k)
          endif
          cDom%tMshTab(i) = k
          exit
        endif
      enddo
    enddo
    cDom%RngLoc(3,1) = cDom%Rng(3,1) +cDom%tMshTab(pGrd%myRow  )*cDom%MshSize(3)
    cDom%RngLoc(3,2) = cDom%Rng(3,1) +cDom%tMshTab(pGrd%myRow+1)*cDom%MshSize(3)
    cDom%MshNumLoc(3)= cDom%tMshTab(pGrd%myRow+1) -cDom%tMshTab(pGrd%myRow) +1
  end associate
  call Rad%update_Radiation_tMsh(tMshTabOld)
end subroutine tLoadBalance


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMM   Gather Source  MMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine cal_monopole0(self,W0,h,Kprm,I0,IC,dV)!,qWeight)
!================================================================
! calculate monopole of effective source. 
! Kprm are
! K  = Kprm(1)
! kx2= Kprm(2)
! ky2= Kprm(3)
!----------------------------------------------------------------
  implicit none
  class(eBeam) :: self
  integer,    intent(in) :: h
  real*8,     intent(in) :: Kprm(3),I0(1),IC(1),dV
  complex(8), intent(out):: W0(self%npt)
!  real*8,optional,intent(in) :: qWeight
    
  associate(npt => self%npt, &
            K    =>Kprm(1),     ku  =>self%ku,    ks=>self%ks,&
            px   =>self%pData(:,px_),                         &
            theta=>self%pData(:,t_),                          &
            gamma=>self%pData(:,g_),                          &
            q    =>self%pData(:,q_),                          &
            Keff => Kprm(1)*(1d0+0.5d0*(Kprm(2)*self%pData(:,x_)*self%pData(:,x_)&
                                       +Kprm(3)*self%pData(:,y_)*self%pData(:,y_))))
    if(mod(h,2)==0) then
      W0 = -(2d0*PhysConst%Zeff/dV*(ku+ks))*(q*px/gamma)&
           *(I0(1)+(0.5d0*h*K*ks/ku*IC(1))*Keff/(gamma*gamma))&
           *exp(-i1*h*theta)
    else
      W0 = -(2d0*PhysConst%Zeff*(ku+ks)*IC(1)/dV)&
           *(q*Keff/gamma)*exp(-i1*h*theta)
    endif
!    if(present(qWeight)) W0=W0*qWeight
  end associate
end subroutine cal_monopole0

subroutine cal_monopole1(self,W0,h,Kprm,I0,IC,dV)!,qWeight)
!================================================================
! calculate monopole of effective source. 
! Kprm are
! K  = Kprm(1)
! kx2= Kprm(2)
! ky2= Kprm(3)
!----------------------------------------------------------------
  implicit none  
  class(eBeam) :: self  
  integer,    intent(in) :: h
  real*8,     intent(in) :: Kprm(3),I0(1),IC(3),dV
  complex(8), intent(out):: W0(self%npt)
!  real*8,optional,intent(in) :: qWeight
  real*8     :: gR
  complex(8) :: IIC
  
  associate(npt  =>self%npt,&
            ku   =>self%ku,           ks   =>self%ks,         &
            K    =>Kprm(1),           kx2  =>Kprm(2),           ky2  =>Kprm(3),&
            px   =>self%pData(:,px_),                         &
            theta=>self%pData(:,t_),  gamma=>self%pData(:,g_),&
            q    =>self%pData(:,q_),                          &
            Keff =>Kprm(1)*(1d0+0.5d0*(Kprm(2)*self%pData(:,x_)*self%pData(:,x_)&
                                      +Kprm(3)*self%pData(:,y_)*self%pData(:,y_))))
    
    if(mod(h,2)==0) then
      W0 = -(2d0*PhysConst%Zeff/dV*(ks+ku))*(q*px/gamma)     &
           *(I0(1) +(0.5d0*h*K*ks/ku*IC(1))*Keff/(gamma*gamma))&
           * exp(-i1*h*theta)
    else
      if(h==1) then
        gR = sqrt(0.5d0*ks/ku*(1d0+0.5d0*K*K))
        associate(delta=> gamma/gR -1d0,&
                  xiR  => 0.25d0*K*K/(1d0+0.5d0*K*K),               &
                  x2   => kx2*self%pData(:,x_)*self%pData(:,x_)+ky2*self%pData(:,y_)*self%pData(:,y_),&
                  p2   => px*px+self%pData(:,py_)*self%pData(:,py_))
          IIc  = (i1*twopi*IC(1) - IC(2))*(0.25d0*ks/(gR*gR)/ku)
          ! W0 = int_C
          W0 = IC(1)+0.5d0*IC(3)*xiR&
              -IIc*(p2 +0.5d0*K*K*x2*(1d0+0.25d0*x2)) &! th0
             +(IIc*((2d0+delta)*delta*(1d0 +p2 +0.5d0*Keff*Keff))   &!th1
               -(IC(3)*xiR)*(1d0+x2))/(2d0*(1d0+delta)**2)
          ! W0 = mono pole
          W0 = -(2d0*PhysConst%Zeff/dV*(ks+ku))&
               *(q*Keff/gamma)*W0*exp(-i1*h*theta)
        end associate
      else
        W0 = -(2d0*PhysConst%Zeff/dV*(ks+ku)*IC(1))&
             *(q*Keff/gamma)*exp(-i1*h*theta)
      endif
    endif
!  if(present(qWeight)) W0=W0*qWeight
  end associate
end subroutine cal_monopole1

subroutine Gather(self,Rad,Kprm,I0,IC,HeffID,src_init)!,qWeight)
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  type(Radiation) :: Rad
  integer,    intent(in) :: HeffID
  real*8,     intent(in) :: Kprm(3),I0(:,:),IC(:,:)
  logical,optional,intent(in) :: src_init
  !real*8,optional,intent(in) :: qWeight
  
  if(.not. present(src_init)) then
    Rad%Src = (0d0, 0d0) ! need to be set elsewhere in order to gather multiple beam
  elseif(src_init) then
    Rad%Src = (0d0, 0d0)
  endif
  if(Kprm(1)==0) return
  select case(self%cDom%wID)
    case(11)  ! transversely and longitudinally Linear
      call Gather_Linear (self,Rad,Kprm,I0,IC,HeffID)!,qWeight)
    case(00)  ! transversely and longitudinally Uniform
      call Gather_Uniform(self,Rad,Kprm,I0,IC,HeffID)!,qWeight)
    case(10)  ! transversely Linear and longitudinally Uniform
      call Gather_xyL_tU (self,Rad,Kprm,I0,IC,HeffID)!,qWeight)
  end select
end subroutine Gather

subroutine Gather_Linear (self,Rad,Kprm,I0,IC,HeffID)!,qWeight)
!================================================================
!     get source deposition with linear weight
!     K  = Kprm(1)
!     kx2= Kprm(2)
!     ky2= Kprm(3)
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  type(Radiation) :: Rad
  integer,    intent(in) :: HeffID
  real*8,     intent(in) :: Kprm(3),I0(:,:),IC(:,:)
  !real*8,optional,intent(in) :: qWeight
   
  integer  :: nx,ny,nt,n,h,i,j,ierr,req(2)
  integer,   dimension(self%npt/self%nBin)  :: gx,gy,gt
  real*8,    dimension(self%npt/self%nBin)  :: fx,fy,ft
  complex(8),dimension(self%npt)            :: W0
  complex(8),allocatable                    :: sBuff(:,:,:), rBuff(:,:,:)
  

  associate(pData=>self%pData, npt =>self%npt,  nBin=>self%nBin,&
            pGrd =>self%pGrd,  cDom=>self%cDom, &
            Src  =>Rad%Src,    nH  =>Rad%nHarm, harm_tab=>Rad%harm_tab)
    nx = cDom%MshNum   (1)
    ny = cDom%MshNumLoc(2)
    nt = cDom%MshNumLoc(3)
    do i=1,npt/nBin
      fx(i) = sum(pData((i-1)*nBin+1:i*nBin,x_))
      fy(i) = sum(pData((i-1)*nBin+1:i*nBin,y_))
      ft(i) = sum(pData((i-1)*nBin+1:i*nBin,t_))
    enddo
    fx = (fx/nBin -cDom%Rng   (1,1))/cDom%MshSize(1)
    fy = (fy/nBin -cDom%RngLoc(2,1))/cDom%MshSize(2)
    ft = (ft/nBin -cDom%RngLoc(3,1))/cDom%MshSize(3)
    gx = int(fx)
    gy = int(fy)
    gt = int(ft)
    fx = fx - dble(gx)
    fy = fy - dble(gy)
    ft = ft - dble(gt)
    gx = gx + 1
    gy = gy + 1
    gt = gt + 1

    do n=1,nH
      h = harm_tab(n)
      if(HeffID==0) then
        call cal_monopole0(self,W0,h,Kprm,I0(:,n),IC(:,n),cDom%dV)!,qWeight)
      else
        call cal_monopole1(self,W0,h,Kprm,I0(:,n),IC(:,n),cDom%dV)!,qWeight)
      endif
      do i=1,npt
        j = (i-1+nBin)/nBin
        if(0<gx(j).and.gx(j)<nx.and.0<gy(j).and.gy(j)<ny.and.0<gt(j).and.gt(j)<nt) then
          Src(gx(j)  ,gy(j)  ,gt(j)  ,n) = Src(gx(j)  ,gy(j)  ,gt(j)  ,n) +W0(i)*(1d0-fx(j))*(1d0-fy(j))*(1d0-ft(j))
          Src(gx(j)+1,gy(j)  ,gt(j)  ,n) = Src(gx(j)+1,gy(j)  ,gt(j)  ,n) +W0(i)*(    fx(j))*(1d0-fy(j))*(1d0-ft(j))
          Src(gx(j)+1,gy(j)+1,gt(j)  ,n) = Src(gx(j)+1,gy(j)+1,gt(j)  ,n) +W0(i)*(    fx(j))*(    fy(j))*(1d0-ft(j))
          Src(gx(j)  ,gy(j)+1,gt(j)  ,n) = Src(gx(j)  ,gy(j)+1,gt(j)  ,n) +W0(i)*(1d0-fx(j))*(    fy(j))*(1d0-ft(j))
          Src(gx(j)  ,gy(j)  ,gt(j)+1,n) = Src(gx(j)  ,gy(j)  ,gt(j)+1,n) +W0(i)*(1d0-fx(j))*(1d0-fy(j))*(    ft(j))
          Src(gx(j)+1,gy(j)  ,gt(j)+1,n) = Src(gx(j)+1,gy(j)  ,gt(j)+1,n) +W0(i)*(    fx(j))*(1d0-fy(j))*(    ft(j))
          Src(gx(j)+1,gy(j)+1,gt(j)+1,n) = Src(gx(j)+1,gy(j)+1,gt(j)+1,n) +W0(i)*(    fx(j))*(    fy(j))*(    ft(j))
          Src(gx(j)  ,gy(j)+1,gt(j)+1,n) = Src(gx(j)  ,gy(j)+1,gt(j)+1,n) +W0(i)*(1d0-fx(j))*(    fy(j))*(    ft(j))
        endif
      enddo
    enddo  
    ! === fill the src(nx,1:ny-1,1:nt-1,nH) ==== 
    ! -- y direction communication 
    allocate(sBuff(nx,nt,nH))
    allocate(rBuff(nx,nt,nH))
    sBuff(:,:,:) = Src(:,ny,:,:)
    n = nx*nt*nH 
    call MPI_IRECV(rBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myD,0,pGrd%comm_col,req(1),ierr)
    call MPI_ISEND(sBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myU,0,pGrd%comm_col,req(2),ierr)
    call MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,ierr)
    if(pGrd%myCol/=0) Src(:,1,:,:) = Src(:,1,:,:) + rBuff(:,:,:)  
    deallocate(sBuff,rBuff)
    
    ! -- theta direction communication 
    allocate(sBuff(nx,ny,nH))
    allocate(rBuff(nx,ny,nH))
    sBuff(:,:,:) = Src(:,:,nt,:)
    n = nx*ny*nH 
    call MPI_IRECV(rBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myL,0,pGrd%comm_row,req(1),ierr)
    call MPI_ISEND(sBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myR,0,pGrd%comm_row,req(2),ierr)
    call MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,ierr)
    if(pGrd%myRow/=0) Src(:,:,1,:) = Src(:,:,1,:) + rBuff(:,:,:)
  end associate
end subroutine Gather_Linear

subroutine Gather_Uniform(self,Rad,Kprm,I0,IC,HeffID)!,qWeight)
!================================================================
!     get source deposition with uniform weight
!     K  = Kprm(1)
!     kx2= Kprm(2)
!     ky2= Kprm(3)
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  type(Radiation) :: Rad
  integer,    intent(in) :: HeffID
  real*8,     intent(in) :: Kprm(3),I0(:,:),IC(:,:)
  !real*8,optional,intent(in) :: qWeight
  integer  :: nx,ny,nt,n,h,i,j,ierr,req(2)
  integer,   dimension(self%npt/self%nBin)  :: gx,gy,gt
  real*8,    dimension(self%npt/self%nBin)  :: fx,fy,ft
  complex(8),dimension(self%npt)            :: W0
  complex(8),allocatable                    :: sBuff(:,:,:), rBuff(:,:,:)
  
  associate(pData=>self%pData, npt =>self%npt,  nBin=>self%nBin,&
            pGrd =>self%pGrd,  cDom=>self%cDom, &
            Src  =>Rad%Src,    nH  =>Rad%nHarm, harm_tab=>Rad%harm_tab)
    nx = cDom%MshNum   (1)
    ny = cDom%MshNumLoc(2)
    nt = cDom%MshNumLoc(3)
    do i=1,npt/nBin
      fx(i) = sum(pData((i-1)*nBin+1:i*nBin,x_))
      fy(i) = sum(pData((i-1)*nBin+1:i*nBin,y_))
      ft(i) = sum(pData((i-1)*nBin+1:i*nBin,t_))
    enddo
    fx = (fx/nBin -cDom%Rng   (1,1))/cDom%MshSize(1) + 0.5d0
    fy = (fy/nBin -cDom%RngLoc(2,1))/cDom%MshSize(2) + 0.5d0
    ft = (ft/nBin -cDom%RngLoc(3,1))/cDom%MshSize(3) + 0.5d0
    gx = int(fx)
    gy = int(fy)
    gt = int(ft)
    gx = gx + 1
    gy = gy + 1
    gt = gt + 1

    do n=1,nH
      h = harm_tab(n)
      if(HeffID==0) then
        call cal_monopole0(self,W0,h,Kprm,I0(:,n),IC(:,n),cDom%dV)!,qWeight)
      else
        call cal_monopole1(self,W0,h,Kprm,I0(:,n),IC(:,n),cDom%dV)!,qWeight)
      endif
      do i=1,npt
        j = (i-1+nBin)/nBin
        if(0<gx(j).and.gx(j)<=nx.and.0<gy(j).and.gy(j)<=ny.and.0<gt(j).and.gt(j)<=nt) then
          Src(gx(j),gy(j),gt(j),n) = Src(gx(j),gy(j),gt(j),n) +W0(i)
        endif
      enddo
    enddo
    
    
    ! === fill the src(nx,1:ny-1,1:nt-1,nH) ==== 
    ! -- y direction communication 
    allocate(sBuff(nx,nt,nH))
    allocate(rBuff(nx,nt,nH))
    sBuff(:,:,:) = Src(:,ny,:,:)
    n = nx*nt*nH 
    
    call MPI_IRECV(rBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myD,0,pGrd%comm_col,req(1),ierr)
    call MPI_ISEND(sBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myU,0,pGrd%comm_col,req(2),ierr)
    call MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,ierr)
    if(pGrd%myCol/=0) Src(:,1,:,:) = Src(:,1,:,:) + rBuff(:,:,:)  
    deallocate(sBuff,rBuff)
    
    ! -- theta direction communication 
    allocate(sBuff(nx,ny,nH))
    allocate(rBuff(nx,ny,nH))
    sBuff(:,:,:) = Src(:,:,nt,:)
    n = nx*ny*nH 
    call MPI_IRECV(rBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myL,0,pGrd%comm_row,req(1),ierr)
    call MPI_ISEND(sBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myR,0,pGrd%comm_row,req(2),ierr)
    call MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,ierr)
    if(pGrd%myRow/=0) Src(:,:,1,:) = Src(:,:,1,:) + rBuff(:,:,:)
  end associate
end subroutine Gather_Uniform

subroutine Gather_xyL_tU (self,Rad,Kprm,I0,IC,HeffID)!,qWeight)
!===============================================================================
! get source deposition with transverse linear and longitudinal uniform weight 
!     K  = Kprm(1)
!     kx2= Kprm(2)
!     ky2= Kprm(3)
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  type(Radiation) :: Rad
  integer,    intent(in) :: HeffID
  real*8,     intent(in) :: Kprm(3),I0(:,:),IC(:,:)
  !real*8,optional,intent(in) :: qWeight  
  integer  :: nx,ny,nt,n,h,i,j,ierr,req(2)
  integer,   dimension(self%npt/self%nBin)  :: gx,gy,gt
  real*8,    dimension(self%npt/self%nBin)  :: fx,fy,ft
  complex(8),dimension(self%npt)            :: W0
  complex(8),allocatable                    :: sBuff(:,:,:), rBuff(:,:,:)
  
  associate(pData=>self%pData, npt =>self%npt,  nBin=>self%nBin,&
            pGrd =>self%pGrd,  cDom=>self%cDom, &
            Src  =>Rad%Src,    nH  =>Rad%nHarm, harm_tab=>Rad%harm_tab)
              
    nx = cDom%MshNum   (1)
    ny = cDom%MshNumLoc(2)
    nt = cDom%MshNumLoc(3)

    do i=1,npt/nBin
      fx(i) = sum(pData((i-1)*nBin+1:i*nBin,x_))
      fy(i) = sum(pData((i-1)*nBin+1:i*nBin,y_))
      ft(i) = sum(pData((i-1)*nBin+1:i*nBin,t_))
    enddo
    fx = (fx/nBin -cDom%Rng   (1,1))/cDom%MshSize(1)
    fy = (fy/nBin -cDom%RngLoc(2,1))/cDom%MshSize(2)
    ft = (ft/nBin -cDom%RngLoc(3,1))/cDom%MshSize(3) + 0.5d0
    ! print*, 'cDom%RngLoc,cDom%MshSize',cDom%RngLoc,cDom%MshSize
    ! print*, 'npt=',npt
    gx = int(fx)
    gy = int(fy)
    gt = int(ft)
    fx = fx - dble(gx)
    fy = fy - dble(gy)
    gx = gx + 1
    gy = gy + 1
    gt = gt + 1
    
    do n=1,nH
      h = harm_tab(n)
      if(HeffID==0) then
        call cal_monopole0(self,W0,h,Kprm,I0(:,n),IC(:,n),cDom%dV)!,qWeight)
      else
        call cal_monopole1(self,W0,h,Kprm,I0(:,n),IC(:,n),cDom%dV)!,qWeight)
      endif
      do i=1,npt
        j = (i-1+nBin)/nBin
        if(0<gx(j).and.gx(j)<nx.and.0<gy(j).and.gy(j)<ny.and.0<gt(j).and.gt(j)<=nt) then
          Src(gx(j)  ,gy(j)  ,gt(j)  ,n) = Src(gx(j)  ,gy(j)  ,gt(j),n) +W0(i)*(1d0-fx(j))*(1d0-fy(j))
          Src(gx(j)+1,gy(j)  ,gt(j)  ,n) = Src(gx(j)+1,gy(j)  ,gt(j),n) +W0(i)*(    fx(j))*(1d0-fy(j))
          Src(gx(j)+1,gy(j)+1,gt(j)  ,n) = Src(gx(j)+1,gy(j)+1,gt(j),n) +W0(i)*(    fx(j))*(    fy(j))
          Src(gx(j)  ,gy(j)+1,gt(j)  ,n) = Src(gx(j)  ,gy(j)+1,gt(j),n) +W0(i)*(1d0-fx(j))*(    fy(j))
  !      else
  !        !call write_log('gx,gy,gt='//num2str(gx(j))//','//num2str(gy(j))//','//num2str(gt(j)))
        endif
      enddo
    enddo 
    !call write_log('Src yBuff ready!!')
    ! === fill the src(nx,1:ny-1,1:nt-1,nH) ==== 
    ! -- y direction communication 
    allocate(sBuff(nx,nt,nH))
    allocate(rBuff(nx,nt,nH))
    sBuff(:,:,:) = Src(:,ny,:,:)
    
    n = nx*nt*nH 
    call MPI_IRECV(rBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myD,0,pGrd%comm_col,req(1),ierr)
    call MPI_ISEND(sBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myU,0,pGrd%comm_col,req(2),ierr)
    call MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,ierr)
    if(pGrd%myCol/=0) Src(:,1,:,:) = Src(:,1,:,:) + rBuff(:,:,:)  
    deallocate(sBuff,rBuff)
    ! -- theta direction communication 
    !call write_log('Src tBuff ready!!')
    allocate(sBuff(nx,ny,nH))
    allocate(rBuff(nx,ny,nH))
    sBuff(:,:,:) = Src(:,:,nt,:)
    n = nx*ny*nH 
    call MPI_IRECV(rBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myL,0,pGrd%comm_row,req(1),ierr)
    call MPI_ISEND(sBuff,n,MPI_DOUBLE_COMPLEX,pGrd%myR,0,pGrd%comm_row,req(2),ierr)
    call MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,ierr)
    if(pGrd%myRow/=0) Src(:,:,1,:) = Src(:,:,1,:) + rBuff(:,:,:)
  end associate
end subroutine Gather_xyL_tU


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMM Scatter Field  MMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine Scatter(self,Rad,fld4ptcl)
!==================================================================
!     linear field interpolation
!------------------------------------------------------------------
  implicit none
  class(eBeam) :: self
  type(Radiation) :: Rad
  complex(8),   intent(out):: fld4ptcl(self%npt/self%nBin,Rad%nHarm)
  
  if(self%npt==0) return
  select case(self%cDom%wID)
    case(11)  ! transversely and longitudinally Linear
      call Scatter_Linear  (self,Rad,fld4ptcl)
    case(00)  ! transversely and longitudinally Uniform
      call Scatter_Uniform (self,Rad,fld4ptcl)
    case(10)  ! transversely Linear and longitudinally Uniform
      call Scatter_xyL_tU  (self,Rad,fld4ptcl)
  end select
end subroutine Scatter

subroutine Scatter_Linear(self,Rad,fld4ptcl)
!==================================================================
!     linear field interpolation
!------------------------------------------------------------------
  implicit none
  class(eBeam) :: self
  type(Radiation) :: Rad
  complex(8),   intent(out):: fld4ptcl(self%npt/self%nBin,Rad%nHarm)
  
  integer,dimension(self%npt/self%nBin) :: gx,gy,gt
  real*8, dimension(self%npt/self%nBin) :: fx,fy,ft
  integer :: i
 
  associate(pData=>self%pData, npt =>self%npt,  nBin=>self%nBin,&
            cDom=>self%cDom,   Fld  =>Rad%Fld)
            
    do i=1,npt/nBin
      fx(i) = sum(pData((i-1)*nBin+1:i*nBin,x_))
      fy(i) = sum(pData((i-1)*nBin+1:i*nBin,y_))
      ft(i) = sum(pData((i-1)*nBin+1:i*nBin,t_))
    enddo
    fx = (fx/nBin -cDom%Rng   (1,1))/cDom%MshSize(1)
    fy = (fy/nBin -cDom%RngLoc(2,1))/cDom%MshSize(2)
    ft = (ft/nBin -cDom%RngLoc(3,1))/cDom%MshSize(3)
    gx = int(fx)
    gy = int(fy)
    gt = int(ft)
    fx = fx - dble(gx)
    fy = fy - dble(gy)
    ft = ft - dble(gt)
    gx = gx + 1
    gy = gy + 1
    gt = gt + 1
    
    forall(i=1:npt/nBin, 0<gx(i).and.gx(i)<cDom%MshNum   (1).and.&
                         0<gy(i).and.gy(i)<cDom%MshNumLoc(2).and.&
                         0<gt(i).and.gt(i)<cDom%MshNumLoc(3))
      fld4ptcl(i,:) =(1d0-fx(i))*(1d0-fy(i))*(1d0-ft(i))*Fld(gx(i)  ,gy(i)  ,gt(i)  ,:) +&
                     (    fx(i))*(1d0-fy(i))*(1d0-ft(i))*Fld(gx(i)+1,gy(i)  ,gt(i)  ,:) +&
                     (    fx(i))*(    fy(i))*(1d0-ft(i))*Fld(gx(i)+1,gy(i)+1,gt(i)  ,:) +&
                     (1d0-fx(i))*(    fy(i))*(1d0-ft(i))*Fld(gx(i)  ,gy(i)+1,gt(i)  ,:) +&
                     (1d0-fx(i))*(1d0-fy(i))*(    ft(i))*Fld(gx(i)  ,gy(i)  ,gt(i)+1,:) +&
                     (    fx(i))*(1d0-fy(i))*(    ft(i))*Fld(gx(i)+1,gy(i)  ,gt(i)+1,:) +&
                     (    fx(i))*(    fy(i))*(    ft(i))*Fld(gx(i)+1,gy(i)+1,gt(i)+1,:) +&
                     (1d0-fx(i))*(    fy(i))*(    ft(i))*Fld(gx(i)  ,gy(i)+1,gt(i)+1,:)
    end forall
  end associate
end subroutine Scatter_Linear

subroutine Scatter_Uniform(self,Rad,fld4ptcl)
!==================================================================
!     uniform field interpolation
!------------------------------------------------------------------
  implicit none
  class(eBeam) :: self
  type(Radiation) :: Rad
  complex(8),   intent(out):: fld4ptcl(self%npt/self%nBin,Rad%nHarm)
  
  integer,dimension(self%npt/self%nBin) :: gx,gy,gt
  real*8, dimension(self%npt/self%nBin) :: fx,fy,ft
  integer :: i
 
  associate(pData=>self%pData, npt =>self%npt,  nBin=>self%nBin,&
            cDom=>self%cDom,   Fld  =>Rad%Fld)
 
   do i=1,npt/nBin
      fx(i) = sum(pData((i-1)*nBin+1:i*nBin,x_))
      fy(i) = sum(pData((i-1)*nBin+1:i*nBin,y_))
      ft(i) = sum(pData((i-1)*nBin+1:i*nBin,t_))
    enddo
    fx = (fx/nBin -cDom%Rng   (1,1))/cDom%MshSize(1) + 0.5d0
    fy = (fy/nBin -cDom%RngLoc(2,1))/cDom%MshSize(2) + 0.5d0
    ft = (ft/nBin -cDom%RngLoc(3,1))/cDom%MshSize(3) + 0.5d0
    gx = int(fx)
    gy = int(fy)
    gt = int(ft)
    gx = gx + 1
    gy = gy + 1
    gt = gt + 1
    
    forall(i=1:npt/nBin, 0<gx(i).and.gx(i)<=cDom%MshNum   (1).and.&
                         0<gy(i).and.gy(i)<=cDom%MshNumLoc(2).and.&
                         0<gt(i).and.gt(i)<=cDom%MshNumLoc(3))
      fld4ptcl(i,:) = Fld(gx(i),gy(i),gt(i),:)
    end forall
  end associate
end subroutine Scatter_Uniform

subroutine Scatter_xyL_tU(self,Rad,fld4ptcl)
!==================================================================
!     linear field interpolation
!------------------------------------------------------------------
  implicit none
  class(eBeam) :: self
  type(Radiation) :: Rad
  complex(8),   intent(out):: fld4ptcl(self%npt/self%nBin,Rad%nHarm)
  
  integer,dimension(self%npt/self%nBin) :: gx,gy,gt
  real*8, dimension(self%npt/self%nBin) :: fx,fy,ft
  integer :: i
 
  associate(pData=>self%pData, npt =>self%npt,  nBin=>self%nBin,&
            cDom=>self%cDom,   Fld =>Rad%Fld)
 
   do i=1,npt/nBin
      fx(i) = sum(pData((i-1)*nBin+1:i*nBin,x_))
      fy(i) = sum(pData((i-1)*nBin+1:i*nBin,y_))
      ft(i) = sum(pData((i-1)*nBin+1:i*nBin,t_))
    enddo
    fx = (fx/nBin -cDom%Rng   (1,1))/cDom%MshSize(1)
    fy = (fy/nBin -cDom%RngLoc(2,1))/cDom%MshSize(2)
    ft = (ft/nBin -cDom%RngLoc(3,1))/cDom%MshSize(3) + 0.5d0
    gx = int(fx)
    gy = int(fy)
    gt = int(ft)
    fx = fx - dble(gx)
    fy = fy - dble(gy)
    gx = gx + 1
    gy = gy + 1
    gt = gt + 1

    
    forall(i=1:npt/nBin, 0<gx(i).and.gx(i)<cDom%MshNum   (1).and.&
                         0<gy(i).and.gy(i)<cDom%MshNumLoc(2).and.&
                         0<gt(i).and.gt(i)<=cDom%MshNumLoc(3))
      fld4ptcl(i,:) =(1d0-fx(i))*(1d0-fy(i))*Fld(gx(i)  ,gy(i)  ,gt(i),:) +&
                     (    fx(i))*(1d0-fy(i))*Fld(gx(i)+1,gy(i)  ,gt(i),:) +&
                     (    fx(i))*(    fy(i))*Fld(gx(i)+1,gy(i)+1,gt(i),:) +&
                     (1d0-fx(i))*(    fy(i))*Fld(gx(i)  ,gy(i)+1,gt(i),:) 
    end forall
  end associate
end subroutine Scatter_xyL_tU

subroutine Get_sample(self,pData,nSample)
!========================================================================
! collect particle data from all processors to root
!------------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  integer, intent(in) :: nSample
  real*8, intent(out) :: pData(nSample,7)
  integer :: ierr,i,j,nSampleLoc,npt_tot
  integer, dimension(0:self%pGrd%np-1) :: npt_tab, npt_cumm
  real*8 :: sRate
  integer,allocatable :: randIndex(:)
  real*8, allocatable :: pDataTmp(:)

  call MPI_ALLreduce(self%npt,npt_tot,1,MPI_INT,MPI_sum,self%pGrd%comm_2d,ierr)
  sRate = nSample/dble(npt_tot)
  call MPI_AllGather(ceiling(self%npt*sRate),1,MPI_INT,npt_tab,1,MPI_INT,self%pGrd%comm_2d,ierr)
  ierr = sum(npt_tab)-nSample
  j=0
  do i=0,self%pGrd%np-1
    if(j==ierr) exit
    if(npt_tab(i)>1) then
      npt_tab(i) = npt_tab(i) -1
      j=j+1
    else
      cycle
    endif
    
  enddo
  npt_cumm(0)=0
  do i=0,self%pGrd%np-2
    npt_cumm(i+1) = npt_cumm(i)+npt_tab(i)
  enddo
  nSampleLoc = npt_tab(self%pGrd%myRank)
  allocate(pDataTmp (nSampleLoc))
  allocate(randIndex(self%npt))
  call math%shuffle(randIndex,self%npt,nSampleLoc)
  do i=1,7
    do j=1,nSampleLoc
      pDataTmp(j) = self%pData(randIndex(j),i)
    enddo
    call MPI_GatherV(pDataTmp,nSampleLoc,&
                     MPI_DOUBLE,pData(:,i),npt_tab,npt_cumm,MPI_DOUBLE,0,self%pGrd%comm_2d,ierr)
  enddo
end subroutine Get_sample

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMM  I/O  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
subroutine Write_Loc(self,fname)
  implicit none
  class(eBeam),     intent(in) :: self
  character(len=*), intent(in) :: fname
  integer :: i,u
  
  u=fileIO%get_free_unit()
  open(u,form="unformatted",action="write",&
        file=fname//trim(fileIO%num2str(self%pGrd%myRank)))
  write(u) self%npt,self%nBin
  do i=1,self%npt
    write(u) self%pData(i,:)
  enddo
  close(u)
end subroutine Write_Loc

subroutine Read_Loc(self,fname)
!================================================================
!     Read particle data written slice by slice 
!----------------------------------------------------------------
  implicit none
  class(eBeam),     intent(inout):: self
  character(len=*), intent(in)   :: fname
  integer :: i,u

  if(allocated(self%pData)) deallocate(self%pData)
  u=fileIO%get_free_unit()
  open(u,form="unformatted",action="read",&
      file=fname//trim(fileIO%num2str(self%pGrd%myRank)))
  read(u) self%npt,self%nBin
  allocate(self%pData(self%npt,7))
  do i=1,self%npt
    read(u) self%pData(i,x_:q_)
  enddo
  close(u)
end subroutine Read_Loc


subroutine Write_Slices(self,fname)
!==================================================================
!     write eBeam data slice by slice
!  output [file] : fname+'islice'  where  islice = 1,2,3,...'nSlice'
!                                         nSlice is determined based on grid number and longotudinal particle data
!------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  character(len=*), intent(in) :: fname
  integer, dimension(0:self%pGrd%npCol-1) :: npt_tab, npt_cumm
  integer, allocatable :: gt(:),u(:)
  real*8,  allocatable :: pDataTmp(:,:),ft(:)
  integer :: nptTot,nt,ierr,i,j,min_gt
  
  call self%tReorder_middle()
  call MPI_AllGather(self%npt,1,MPI_INT,npt_tab,1,MPI_INT,self%pGrd%comm_col,ierr)
  nptTot = sum(npt_tab)
  npt_cumm(0)=0
  do i=0,self%pGrd%npCol-2
    npt_cumm(i+1) = npt_cumm(i)+npt_tab(i)
  enddo 
  allocate(pDataTmp(nptTot,7))
  do i=1,7
    call MPI_GatherV(self%pData(:,i),self%npt,MPI_DOUBLE,pDataTmp(:,i),npt_tab,npt_cumm,MPI_DOUBLE,0,self%pGrd%comm_col,ierr)
  enddo 
  if(self%pGrd%myCol==0) then
    allocate(ft(nptTot/self%nBin))
    allocate(gt(nptTot/self%nBin))
    do i=1,nptTot/self%nBin
      ft(i) = sum(pDataTmp((i-1)*self%nBin+1:i*self%nBin,t_))
    enddo
    ft = (ft/self%nBin -self%cDom%RngLoc(3,1))/self%cDom%MshSize(3) + 0.5d0
    gt = int(ft)
    if(nptTot>0) then
      j = minval(gt+self%cDom%tMshTab(self%pGrd%myRow))
    else
      j = huge(j)
    endif
    call MPI_AllReduce(j,min_gt,1,MPI_INT,MPI_MIN,self%pGrd%comm_row,ierr)
    nt = self%cDom%MshNumLoc(3)-1
    if(self%pGrd%myRow==self%pGrd%npRow-1) nt = nt+1
    if(self%cDom%tMshTab(self%pGrd%myRow+1) <= min_gt) then
      nt = 0
    elseif(self%cDom%tMshTab(self%pGrd%myRow) <= min_gt .and. &
           min_gt < self%cDom%tMshTab(self%pGrd%myRow+1)) then
      nt = nt - (min_gt-self%cDom%tMshTab(self%pGrd%myRow))
      gt = gt - (min_gt-self%cDom%tMshTab(self%pGrd%myRow))
    endif
    gt = gt+1
    allocate(u(nt))
    do i=1,nt
      u(i)=fileIO%get_free_unit()
      if(self%cDom%tMshTab(self%pGrd%myRow) <= min_gt .and. &
           min_gt < self%cDom%tMshTab(self%pGrd%myRow+1)) then
        open(u(i),form="unformatted",action="write",&
                  file=fname//trim(fileIO%num2str(i)))
      else
        open(u(i),form="unformatted",action="write",&
                  file=fname//trim(fileIO%num2str(i+self%cDom%tMshTab(self%pGrd%myRow)-min_gt)))
      endif
    enddo
    do i=1,nptTot/self%nBin
      do j=1,self%nBin
        write(u(gt(i))) pDataTmp((i-1)*self%nBin+j,:)
      enddo
    enddo
    do i=1,nt
      close(u(i))
    enddo
  endif
  call self%tReorder()
end subroutine Write_Slices

subroutine Read_Slices(self,fname,nSlice)
!================================================================
!     Read particle data written slice by slice 
! input
!   nSlice [int] : total number of particle data file named as fname+'iSlice'
!                  where iSlice=1,2,3,...nSlince
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  integer, intent(in) :: nSlice
  character(len=*), intent(in) :: fname
  integer, allocatable :: npt(:)
  integer :: u,nS,nS1,iS,i,ip,IO,nptSum
  character(len=5) :: str_iS
  character(len=*), parameter  :: fmt_ = "(I5)"   

  nS = nSlice/self%pGrd%np
  nS1= nSlice-nS*self%pGrd%np
  iS = self%pGrd%myrank*nS
  if(self%pGrd%myrank >= self%pGrd%np-nS1) then
    nS = nS +1
    iS = iS +self%pGrd%myrank-(self%pGrd%np-nS1)
  endif
  allocate(npt(nS))
  npt=0
  do i=1,nS
    u=fileIO%get_free_unit()
    write(str_iS,fmt_) iS+i
    open(u,file=fname//trim(ADJUSTL(str_iS)),&
           form="unformatted",action="read")
    do
      read(u,iostat=IO)
      if (IO/=0) exit
      npt(i) = npt(i) +1
    end do
    close(u)
  enddo
  
  if(allocated(self%pData)) deallocate(self%pData)
  nptSum = sum(npt)
  allocate(self%pData(nptSum,7))
  self%npt=nptSum
  nptSum = 0
  do i=1,nS
    u=fileIO%get_free_unit()
    write(str_iS,fmt_) iS+i
    open(u,file=fname//trim(ADJUSTL(str_iS)),&
           form="unformatted",action="read")
    if(i>1) nptSum = nptSum + npt(i-1)
    do ip=nptSum+1,nptSum+npt(i)
      read(u) self%pData(ip,x_:q_)
    enddo
    close(u)
  enddo
  call self%Reorder()
end subroutine Read_Slices

real*8 function get_mean(self,i_)
  implicit none
  include 'mpif.h'
  class(eBeam) :: self
  integer, intent(in) :: i_
  integer:: npt,ierr 
  call MPI_ALLreduce(sum(self%pData(:,i_)),get_mean,1,MPI_double,MPI_sum,self%pGrd%comm_2d,ierr)
  call MPI_ALLreduce(self%npt,npt,1,MPI_integer,MPI_sum,self%pGrd%comm_2d,ierr)
  get_mean = get_mean/npt 
  return
end function get_mean

real*8 function get_std(self,i_)
implicit none
  include 'mpif.h'
  class(eBeam) :: self
  integer, intent(in) :: i_
  real*8 :: mean
  integer:: npt,ierr 
  mean = self%get_mean(i_)
  call MPI_ALLreduce(self%npt,npt,1,MPI_integer,MPI_sum,self%pGrd%comm_2d,ierr)
  call MPI_ALLreduce(sum((self%pData(:,i_)-mean)**2),get_std,1,MPI_double,MPI_sum,self%pGrd%comm_2d,ierr)
  get_std = sqrt(get_std/(npt-1))
  return
end function get_std

end module eBeamModule
