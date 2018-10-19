module RadiationModule
  use ConstDataModule  ! included in CompDomModule
  use MathModule
  use pGrid2DModule
  use CompDomModule
  implicit none
  private
  
  type, public :: Radiation
    integer :: nx,ny,nt
    integer :: nxT,nyT,ntT
    integer :: nxT2,ny2
    ! integer :: wID  ! shape function 
               ! !wID = 00 transversly, longitudinally uniform
               ! !wID = 10 transversly linear, longitudinally uniform
               ! !wID = 11 transversly, longitudinally linear
    integer :: nHarm
    integer,allocatable :: harm_tab(:)
    real*8 :: ks,ku,dx,dy,dt  ! dt : Msh size of theta domain not time
    complex(8), allocatable :: Fld(:,:,:,:),FldT(:,:,:,:)
    complex(8), allocatable :: Src(:,:,:,:),SrcT(:,:,:,:)
    complex(8), allocatable :: sBuff_y2x(:),rBuff_y2x(:)
    ! save pre-calculated table for MPI computation of xy transpose
    integer, allocatable, dimension(:) :: sCnt_y2x,rCnt_y2x,sDsp_y2x,rDsp_y2x
    ! MPI and domain info
    type(pGrid2D), pointer :: pGrd => null()
    type(CompDom), pointer :: cDom => null()
    
    contains
      procedure :: destruct => Radiation_destructor
      procedure :: FldSolve, Slip, Get_Pwr, Get_Energy
      ! need to be used together with loadbalance
      procedure :: update_Radiation_tMsh, update_Radiation_yMsh
  end type Radiation
  
  interface Radiation 
    procedure Radiation_constructor
  end interface
  
  ! === example use ============================
  ! type(Radiation),pointer :: rad
  ! rad -> Radiation(pGrd,cDom,harm_tab,nHarm,ks,ku)  ! constructor
  ! call rad%FldSolve(dz)  ! steady-state fld solver
  ! call rad%Slip(dz)      ! slippage
  ! pwr = rad%get_pwr()    ! get slice-by-slice power for each harmonic
  ! ===========================================
  
contains

function Radiation_constructor(pGrd,cDom,harm_tab,nHarm,ks,ku) result(self)
!================================================================
!    Initialize field class.
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(Radiation),pointer :: self
  type(CompDom),   pointer :: cDom
  type(pGrid2D),   pointer :: pGrd
  integer,intent(in) :: nHarm, harm_tab(nHarm)
  real*8, intent(in) :: ks,ku
  integer :: n
  
  allocate(self)
  self%pGrd => pGrd
  self%cDom => cDom
  self%ks  = ks
  self%ku  = ku
  self%nx  = cDom%MshNum   (1)
  self%ny  = cDom%MshNumLoc(2)
  self%nt  = cDom%MshNumLoc(3)
  self%nxT = cDom%MshNumLoc(1)
  self%nyT = cDom%MshNum   (2)

  if(pGrd%myRow==pGrd%npRow-1) then
    self%ntT = self%nt
  else
    self%ntT = self%nt -1
  endif

  if(pGrd%myCol==pGrd%npCol-1) then
    self%nxT2 = self%nxT
    self%ny2  = self%ny
  else
    self%nxT2 = self%nxT-1
    self%ny2  = self%ny-1
  endif
  self%nHarm = nHarm
  allocate(self%harm_tab(nHarm))
  self%harm_tab = harm_tab
  allocate(self%Fld (  self%nx, 0:self%ny, self%nt, nHarm))
  allocate(self%Src (  self%nx,   self%ny, self%nt, nHarm))
  allocate(self%FldT(0:self%nxT,  self%nyT,self%ntT,nHarm))
  allocate(self%SrcT(  self%nxT,  self%nyT,self%ntT,nHarm))
  self%Fld  = (0d0,0d0)
  self%FldT = (0d0,0d0)
  self%Src  = (0d0,0d0)
  self%SrcT = (0d0,0d0)
  
!  do i=1,self%nx
!    do j=0,self%ny
!      do k=1,self%nt
!        self%Fld(i,j,k,:) = i*(1d0+i1) &
!                          +(10d0+0.1d0*i1)*(j+cDom%yMshTab(pGrd%myCol))&
!                          +(100d0+0.01d0*i1)*(k+cDom%tMshTab(pGrd%myRow))
!      enddo
!    enddo
!  enddo
!  if(pGrd%myCol==0) self%Fld(:,0,:,:) = (0d0,0d0)
  
  self%dx = cDom%MshSize(1)
  self%dy = cDom%MshSize(2)
  self%dt = cDom%MshSize(3)

  !--------------------------------------------------------------
  !-------- pre-calculate y2x transformation buffer info --------
  allocate(self%sCnt_y2x(0:pGrd%npCol-1),self%rCnt_y2x(0:pGrd%npCol-1))
  do n=0,pGrd%npCol-2
    self%sCnt_y2x(n) = (cDom%xMshTab(n+1)-cDom%xMshTab(n))*self%ny2
    self%rCnt_y2x(n) = self%nxT2*(cDom%yMshTab(n+1)-cDom%yMshTab(n))
  enddo
  n = pGrd%npCol-1
  self%sCnt_y2x(n) = (cDom%xMshTab(n+1)-cDom%xMshTab(n)+1)*self%ny2
  self%rCnt_y2x(n) = self%nxT2*(cDom%yMshTab(n+1)-cDom%yMshTab(n)+1)
  self%sCnt_y2x = self%sCnt_y2x*self%ntT*nHarm
  self%rCnt_y2x = self%rCnt_y2x*self%ntT*nHarm
  
  allocate(self%sDsp_y2x(0:pGrd%npCol-1),self%rDsp_y2x(0:pGrd%npCol-1))
  self%sDsp_y2x(0)=0
  self%rDsp_y2x(0)=0
  do n=1,pGrd%npCol-1
    self%sDsp_y2x(n) = self%sDsp_y2x(n-1)+self%sCnt_y2x(n-1)
    self%rDsp_y2x(n) = self%rDsp_y2x(n-1)+self%rCnt_y2x(n-1)
  enddo
  allocate(self%sBuff_y2x(self%ntT*nHarm*self%ny2*self%nx  ))
  allocate(self%rBuff_y2x(self%ntT*nHarm*self%nyT*self%nxT2))
end function Radiation_constructor

subroutine update_Radiation_tMsh(self,tMshTabOld)
!================================================================
!    re-construct radiation data structure 
!    after theta domain CompDom change (by load balancing or radiation slip)
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(Radiation)    :: self
  integer, intent(in) :: tMshTabOld(0:self%pGrd%npRow)
  
  integer, dimension(0:self%pGrd%npRow-1) :: sCnt,sDsp,rCnt,rDsp
  complex(8) :: sBuff(self%nx,0:self%ny,self%ntT,self%nHarm)
  integer :: i,k,myStart,myEnd
  
  associate(cDom => self%cDom, pGrd=>self%pGrd)
    !call write_log('update tRad')
    if(all(cDom%tMshTab(1:pGrd%npRow-1)==tMshTabOld(1:pGrd%npRow-1))) then
      !call write_log('tMsh == old_tMsh ????')
      return
    endif
    ! ==== Prepare Buff ==== 
    sBuff(:,:,:,:) = self%Fld(:,:,1:self%ntT,:)
    ! ---- send counts/displacements -----
    myStart = tMshTabOld(pGrd%myRow)
    myEnd   = tMshTabOld(pGrd%myRow+1)
    sCnt=0
    do k=pGrd%npRow,0,-1
      if(cDom%tMshTab(k) <= myStart) exit
    enddo
    if(cDom%tMshTab(k+1) >= myEnd) then
      sCnt(k) = myEnd - myStart
    else
      sCnt(k) = cDom%tMshTab(k+1) - myStart
      do i=k+1,pGrd%npRow-1
        if(cDom%tMshTab(i+1) >= myEnd) then
          sCnt(i) = myEnd - cDom%tMshTab(i)
          exit
        else
          sCnt(i) = cDom%tMshTab(i+1) - cDom%tMshTab(i)
        endif
      enddo
    endif
    if(pGrd%myRow == pGrd%npRow-1) sCnt(pGrd%npRow-1)=sCnt(pGrd%npRow-1)+1
    sCnt = sCnt*self%nx*(self%ny+1)
    sDsp=0
    do i=1,pGrd%npRow-1
      sDsp(i) = sDsp(i-1) + sCnt(i-1)
    enddo
    ! ----- recv counts/displacements ----
    myStart = cDom%tMshTab(pGrd%myRow)
    myEnd   = cDom%tMshTab(pGrd%myRow+1)
    rCnt = 0
    do k=pGrd%npRow,0,-1
      if(tMshTabOld(k) <= myStart) exit
    enddo
    if(tMshTabOld(k+1) >= myEnd) then
      rCnt(k) = myEnd-myStart
    else
      rCnt(k) = tMshTabOld(k+1)-myStart
      do i=k+1,pGrd%npRow
        if(tMshTabOld(i+1) >= myEnd) then
          rCnt(i) = myEnd-tMshTabOld(i)
          exit
        else
          rCnt(i) = tMshTabOld(i+1)-tMshTabOld(i)
        endif 
      enddo
    endif
    if(pGrd%myRow == pGrd%npRow-1) rCnt(pGrd%npRow-1)=rCnt(pGrd%npRow-1)+1
    rCnt = rCnt*self%nx*(self%ny+1)
    rDsp=0
    do i=1,pGrd%npRow-1
      rDsp(i) = rDsp(i-1) + rCnt(i-1)
    enddo
    ! ==== re-struct radiation data =====
    self%sCnt_y2x = self%sCnt_y2x/self%ntT
    self%rCnt_y2x = self%rCnt_y2x/self%ntT
    self%sDsp_y2x = self%sDsp_y2x/self%ntT
    self%rDsp_y2x = self%rDsp_y2x/self%ntT
    
    self%nt  = cDom%MshNumLoc(3)
    if(pGrd%myRow==pGrd%npRow-1) then
      self%ntT = self%nt
    else
      self%ntT = self%nt -1
    endif
    
    self%sCnt_y2x = self%sCnt_y2x*self%ntT
    self%rCnt_y2x = self%rCnt_y2x*self%ntT
    self%sDsp_y2x = self%sDsp_y2x*self%ntT
    self%rDsp_y2x = self%rDsp_y2x*self%ntT
    
    deallocate(self%sBuff_y2x,self%rBuff_y2x)
    allocate(self%sBuff_y2x(self%ntT*self%nHarm*self%ny2*self%nx  ))
    allocate(self%rBuff_y2x(self%ntT*self%nHarm*self%nyT*self%nxT2))
    
    deallocate(self%Fld,self%Src,self%FldT,self%SrcT)
    allocate(self%Fld (  self%nx, 0:self%ny, self%nt, self%nHarm))
    allocate(self%Src (  self%nx,   self%ny, self%nt, self%nHarm))
    allocate(self%FldT(0:self%nxT,  self%nyT,self%ntT,self%nHarm))
    allocate(self%SrcT(  self%nxT,  self%nyT,self%ntT,self%nHarm))
    self%Fld  = (0d0,0d0)
    self%FldT = (0d0,0d0)
    self%Src  = (0d0,0d0)
    self%SrcT = (0d0,0d0)
    ! ==== communicate buffs ====
    do i=1,self%nHarm
      call MPI_AlltoAllv(sBuff  (:,:,:,i), sCnt,sDsp,MPI_DOUBLE_COMPLEX,&
                         self%Fld(:,:,:,i),rCnt,rDsp,MPI_DOUBLE_COMPLEX, pGrd%comm_row, k)
    enddo
    call fill_t_guard_cell(self)
  end associate
end subroutine update_Radiation_tMsh

subroutine update_Radiation_yMsh(self,yMshTabOld)
!================================================================
!    re-construct radiation data structure 
!    after theta domain CompDom change (by load balancing or radiation slip)
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  class(Radiation)    :: self
  integer, intent(in) :: yMshTabOld(0:self%pGrd%npCol)
  
  integer, dimension(0:self%pGrd%npCol-1) :: sCnt,sDsp,rCnt,rDsp
  complex(8) :: sBuff(self%nx,self%ntT,self%nHarm,self%ny2)
  integer :: i,k,myStart,myEnd
  
  associate(cDom => self%cDom, pGrd=>self%pGrd)
  
    if(any(cDom%yMshTab/=yMshTabOld)) return
    ! ==== Prepare Buff ==== 
    do i=1,self%ny2
      sBuff(:,:,:,i) = self%Fld(:,i,:,:)
    enddo
    ! ---- send counts/displacements -----
    myStart = yMshTabOld(pGrd%myCol)
    myEnd   = yMshTabOld(pGrd%myCol+1)
    sCnt=0
    do k=pGrd%npCol,0,-1
      if(cDom%yMshTab(k) <= myStart) exit
    enddo
    if(cDom%yMshTab(k+1) >= myEnd) then
      sCnt(k) = myEnd - myStart
    else
      sCnt(k) = cDom%yMshTab(k+1) - myStart
      do i=k+1,pGrd%npCol-1
        if(cDom%yMshTab(i+1) >= myEnd) then
          sCnt(i) = myEnd - cDom%yMshTab(i)
          exit
        else
          sCnt(i) = cDom%yMshTab(i+1) - cDom%yMshTab(i)
        endif
      enddo
    endif
    if(pGrd%myCol == pGrd%npCol-1) sCnt(pGrd%npCol-1)=sCnt(pGrd%npCol-1)+1
    sCnt = sCnt*self%nx*self%nt*self%nHarm
    sDsp=0
    do i=1,pGrd%npCol-1
      sDsp(i) = sDsp(i-1) + sCnt(i-1)
    enddo
    ! ----- recv counts/displacements ----
    myStart = cDom%yMshTab(pGrd%myCol)
    myEnd   = cDom%yMshTab(pGrd%myCol+1)
    rCnt = 0
    do k=pGrd%npCol,0,-1
      if(yMshTabOld(k) <= myStart) exit
    enddo
    if(yMshTabOld(k+1) >= myEnd) then
      rCnt(k) = myEnd-myStart
    else
      rCnt(k) = yMshTabOld(k+1)-myStart
      do i=k+1,pGrd%npCol
        if(yMshTabOld(i+1) >= myEnd) then
          rCnt(i) = myEnd-yMshTabOld(i)
          exit
        else
          rCnt(i) = yMshTabOld(i+1)-yMshTabOld(i)
        endif 
      enddo
    endif
    if(pGrd%myCol == pGrd%npCol-1) rCnt(pGrd%npCol-1)=rCnt(pGrd%npCol-1)+1
    rCnt = rCnt*self%nx*self%nt*self%nHarm
    rDsp=0
    do i=1,pGrd%npCol-1
      rDsp(i) = rDsp(i-1) + rCnt(i-1)
    enddo
    ! ==== re-struct radiation data =====
    self%ny  = cDom%MshNumLoc(2)
    if(pGrd%myCol==pGrd%npCol-1) then
      self%ny2  = self%ny
    else
      self%ny2  = self%ny-1
    endif
    do i=0,pGrd%npCol-2
      self%sCnt_y2x(i) = (cDom%xMshTab(i+1)-cDom%xMshTab(i))*self%ny2
      self%rCnt_y2x(i) = self%nxT2*(cDom%yMshTab(i+1)-cDom%yMshTab(i))
    enddo
    i = pGrd%npCol-1
    self%sCnt_y2x(i) = (cDom%xMshTab(i+1)-cDom%xMshTab(i)+1)*self%ny2
    self%rCnt_y2x(i) = self%nxT2*(cDom%yMshTab(i+1)-cDom%yMshTab(i)+1)
    self%sCnt_y2x = self%sCnt_y2x*self%ntT*self%nHarm
    self%rCnt_y2x = self%rCnt_y2x*self%ntT*self%nHarm

    self%sDsp_y2x(0)=0
    self%rDsp_y2x(0)=0
    do i=1,pGrd%npCol-1
      self%sDsp_y2x(i) = self%sDsp_y2x(i-1)+self%sCnt_y2x(i-1)
      self%rDsp_y2x(i) = self%rDsp_y2x(i-1)+self%rCnt_y2x(i-1)
    enddo
    
    deallocate(self%sBuff_y2x,self%rBuff_y2x)
    allocate(self%sBuff_y2x(self%ntT*self%nHarm*self%ny2*self%nx  ))
    allocate(self%rBuff_y2x(self%ntT*self%nHarm*self%nyT*self%nxT2))
    
    deallocate(self%Fld,self%Src,self%FldT,self%SrcT)
    allocate(self%Fld (  self%nx, 0:self%ny, self%nt, self%nHarm))
    allocate(self%Src (  self%nx,   self%ny, self%nt, self%nHarm))
    allocate(self%FldT(0:self%nxT,  self%nyT,self%ntT,self%nHarm))
    allocate(self%SrcT(  self%nxT,  self%nyT,self%ntT,self%nHarm))
    self%Fld  = (0d0,0d0)
    self%FldT = (0d0,0d0)
    self%Src  = (0d0,0d0)
    self%SrcT = (0d0,0d0)
    ! ==== communicate buffs ====
    do i=1,self%nHarm
      call MPI_AlltoAllv(sBuff  (:,:,:,i),sCnt,sDsp,MPI_DOUBLE_COMPLEX,&
                         self%Fld(:,:,:,i),rCnt,rDsp,MPI_DOUBLE_COMPLEX, pGrd%comm_col, k)
    enddo
    call fill_t_guard_cell(self)
  end associate
end subroutine update_Radiation_yMsh

subroutine Radiation_destructor(self)
  implicit none
  class(Radiation) :: self
  if(allocated(self%harm_tab))  deallocate(self%harm_tab)
  if(allocated(self%Fld))       deallocate(self%Fld)
  if(allocated(self%FldT))      deallocate(self%FldT)
  if(allocated(self%Src))       deallocate(self%Src)
  if(allocated(self%SrcT))      deallocate(self%SrcT)
  if(allocated(self%sBuff_y2x)) deallocate(self%sBuff_y2x)
  if(allocated(self%rBuff_y2x)) deallocate(self%rBuff_y2x)
  self%pGrd => null()
  self%cDom => null()
end subroutine Radiation_destructor


subroutine xADI(self,dz)
  implicit none
  type(Radiation) :: self
  real*8, intent(in) :: dz
  integer :: n,h,iy,it
  real*8 :: dx2,dy2
  complex(8) :: dz2,off_diag
  complex(8), dimension(self%nx,self%ny2,self%ntT,self%nHarm) :: TempSrc
  

  TempSrc=(0d0,0d0)
  dx2 = self%dx*self%dx
  dy2 = self%dy*self%dy
  do n=1,self%nHarm
    h=self%harm_tab(n)
    dz2 = -i1*0.25d0*dz/h/self%ks
    do iy=1,self%ny-1 !----------------------------------------------------------------
      TempSrc(:,iy,:,n) = self%Fld(:,iy,1:self%ntT,n)+dz2*self%Src(:,iy,1:self%ntT,n) -dz2/dy2&
                        *(self%Fld(:,iy-1,1:self%ntT,n)+self%Fld(:,iy+1,1:self%ntT,n)&
                          -2d0*self%Fld(:,iy,1:self%ntT,n))
    enddo
    if(self%pGrd%myCol==self%pGrd%npCol-1) then
      iy=self%ny !------------------------------------------------------------------------
      TempSrc(:,iy,:,n) = self%Fld(:,iy,1:self%ntT,n)+dz2*self%Src(:,iy,1:self%ntT,n) -dz2/dy2&
                        *(self%Fld(:,iy-1,1:self%ntT,n) &
                          -2d0*self%Fld(:,iy,1:self%ntT,n))
    endif
  enddo
  do n=1,self%nHarm
    h=self%harm_tab(n)
    dz2 = -i1*0.25d0*dz/h/self%ks
    off_diag = dz2/dx2
    do iy=1,self%ny2
      do it=1,self%ntT
        call math%TriDiag(self%Fld(:,iy,it,n),TempSrc(:,iy,it,n),off_diag,self%nx)
      enddo
    enddo
  enddo
end subroutine xADI

subroutine yADI(self,dz)
  implicit none
  include 'mpif.h'
  type(Radiation) :: self
  real*8, intent(in) :: dz
  
  integer :: n,h,ix,it
  real*8 :: dx2,dy2
  complex(8) :: dz2,off_diag
  complex(8), dimension(self%nxT2,self%nyT,self%ntT,self%nHarm) :: TempSrc

  dx2 = self%dx*self%dx
  dy2 = self%dy*self%dy
  do n=1,self%nHarm
    h=self%harm_tab(n)
    dz2 = -i1*0.25d0*dz/h/self%ks
    do ix=1,self%nxT-1
      TempSrc(ix,:,:,n) = self%FldT(ix,:,:,n)+dz2*self%SrcT(ix,:,:,n)&
                         -dz2/dx2*(self%FldT(ix-1,:,:,n)-2d0*self%FldT(ix,:,:,n)+self%FldT(ix+1,:,:,n))
    enddo
    if(self%pGrd%myCol==self%pGrd%npCol-1) then
      ix=self%nxT
      TempSrc(ix,:,:,n) = self%FldT(ix,:,:,n)+dz2*self%SrcT(ix,:,:,n)&
                          -dz2/dx2*(self%FldT(ix-1,:,:,n)-2d0*self%FldT(ix,:,:,n))
    endif
  enddo
  
  do n=1,self%nHarm
    h=self%harm_tab(n)
    dz2 = -i1*0.25d0*dz/h/self%ks
    off_diag = dz2/dy2
    do ix=1,self%nxT2
      do it=1,self%ntT
        call math%TriDiag(self%FldT(ix,:,it,n),TempSrc(ix,:,it,n),off_diag,self%nyT)
      enddo
    enddo
  enddo
end subroutine yADI


subroutine y2x_transpose(self)
!================================================================
! transpose y-chunked radiation matrix to x-chunked matirx
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(Radiation) :: self
  integer :: i,n,ierr
  associate(pGrd=>self%pGrd, cDom=>self%cDom)
    !----------------------- Fld Transpose --------------------------
    n=self%nHarm*self%ntT*self%ny2
    do i=0,pGrd%npCol-2
      self%sBuff_y2x(n*cDom%xMshTab(i)+1:n*cDom%xMshTab(i+1)) &
        = reshape(self%Fld(cDom%xMshTab(i)+1:cDom%xMshTab(i+1),1:self%ny2,1:self%ntT,:),&
                  [n*(cDom%xMshTab(i+1)-cDom%xMshTab(i))])
    enddo
    i=pGrd%npCol-1
    self%sBuff_y2x(n*cDom%xMshTab(i)+1:) &
      = reshape(self%Fld(cDom%xMshTab(i)+1:,1:self%ny2,1:self%ntT,:),&
                [n*(self%nx-cDom%xMshTab(i))])
    call MPI_AlltoAllv(self%sBuff_y2x,self%sCnt_y2x,self%sDsp_y2x,MPI_DOUBLE_COMPLEX,&
                       self%rBuff_y2x,self%rCnt_y2x,self%rDsp_y2x,MPI_DOUBLE_COMPLEX,Pgrd%comm_col,ierr)
    do i=0,pGrd%npCol-2
      self%FldT(1:self%nxT2,cDom%yMshTab(i)+1:cDom%yMshTab(i+1),1:self%ntT,1:self%nHarm) &
        = reshape(self%rBuff_y2x(self%rDsp_y2x(i)+1:self%rDsp_y2x(i+1)) &
                  ,[self%nxT2,cDom%yMshTab(i+1)-cDom%yMshTab(i),self%ntT,self%nHarm])
    enddo
    i=pGrd%npCol-1
    self%FldT(1:self%nxT2,cDom%yMshTab(i)+1:self%nyT,1:self%ntT,1:self%nHarm) &
        = reshape(self%rBuff_y2x(self%rDsp_y2x(i)+1:self%rDsp_y2x(i)+self%rCnt_y2x(i)) &
                  ,[self%nxT2,self%nyT-cDom%yMshTab(i),self%ntT,self%nHarm])
    !----------------------- Src Transpose --------------------------
    do i=0,pGrd%npCol-2
      self%sBuff_y2x(n*cDom%xMshTab(i)+1:n*cDom%xMshTab(i+1)) &
        = reshape(self%Src(cDom%xMshTab(i)+1:cDom%xMshTab(i+1),1:self%ny2,1:self%ntT,:),&
                  [n*(cDom%xMshTab(i+1)-cDom%xMshTab(i))])
    enddo
    i=pGrd%npCol-1
    self%sBuff_y2x(n*cDom%xMshTab(i)+1:) &
      = reshape(self%Src(cDom%xMshTab(i)+1:,1:self%ny2,1:self%ntT,:),&
                [n*(self%nx-cDom%xMshTab(i))])

    call MPI_AlltoAllv(self%sBuff_y2x,self%sCnt_y2x,self%sDsp_y2x,MPI_DOUBLE_COMPLEX,&
                       self%rBuff_y2x,self%rCnt_y2x,self%rDsp_y2x,MPI_DOUBLE_COMPLEX,Pgrd%comm_col,ierr)
                      
    do i=0,pGrd%npCol-2
      self%SrcT(1:self%nxT2,cDom%yMshTab(i)+1:cDom%yMshTab(i+1),1:self%ntT,1:self%nHarm) &
        = reshape(self%rBuff_y2x(self%rDsp_y2x(i)+1:self%rDsp_y2x(i+1)) &
                  ,[self%nxT2,cDom%yMshTab(i+1)-cDom%yMshTab(i),self%ntT,self%nHarm])
    enddo
    i=pGrd%npCol-1
    self%SrcT(1:self%nxT2,cDom%yMshTab(i)+1:self%nyT,1:self%ntT,1:self%nHarm) &
        = reshape(self%rBuff_y2x(self%rDsp_y2x(i)+1:self%rDsp_y2x(i)+self%rCnt_y2x(i)) &
                  ,[self%nxT2,self%nyT-cDom%yMshTab(i),self%ntT,self%nHarm])
                  
    !------------ fill x guard cell ------------------
    call fill_x_guard_cell(self)
  end associate
end subroutine y2x_transpose

subroutine x2y_transpose(self)
!================================================================
! transpose y-chunked radiation matrix to x-chunked matirx
!----------------------------------------------------------------
  implicit none
  include 'mpif.h'
  type(Radiation) :: self
  !complex(8) :: yBuff(self%nx,self%ntT,self%nHarm)
  integer :: i,n,ierr
  !----------------------- Fld Transpose --------------------------
  associate(pGrd=>self%pGrd, cDom=>self%cDom)
    n=self%nHarm*self%ntT*self%nxT2
    do i=0,pGrd%npCol-2
      self%rBuff_y2x(n*cDom%yMshTab(i)+1:n*cDom%yMshTab(i+1)) &
        = reshape(self%FldT(1:self%nxT2,cDom%yMshTab(i)+1:cDom%yMshTab(i+1),:,:),&
                  [n*(cDom%yMshTab(i+1)-cDom%yMshTab(i))])
    enddo
    i=pGrd%npCol-1
    self%rBuff_y2x(n*cDom%yMshTab(i)+1:) &
      = reshape(self%FldT(1:self%nxT2,cDom%yMshTab(i)+1:,:,:),&
                [n*(self%nyT-cDom%yMshTab(i))])
    call MPI_AlltoAllv(self%rBuff_y2x,self%rCnt_y2x,self%rDsp_y2x,MPI_DOUBLE_COMPLEX,&
                       self%sBuff_y2x,self%sCnt_y2x,self%sDsp_y2x,MPI_DOUBLE_COMPLEX,Pgrd%comm_col,ierr)
    do i=0,pGrd%npCol-2
      self%Fld(cDom%xMshTab(i)+1:cDom%xMshTab(i+1),1:self%ny2,1:self%ntT,:) &
        = reshape(self%sBuff_y2x(self%sDsp_y2x(i)+1:self%sDsp_y2x(i+1)) &
                  ,[cDom%xMshTab(i+1)-cDom%xMshTab(i),self%ny2,self%ntT,self%nHarm])
    enddo
    i=pGrd%npCol-1
    self%Fld(cDom%xMshTab(i)+1:self%nx,1:self%ny2,1:self%ntT,:) &
        = reshape(self%sBuff_y2x(self%sDsp_y2x(i)+1:self%sDsp_y2x(i)+self%sCnt_y2x(i)) &
                  ,[self%nx-cDom%xMshTab(i),self%ny2,self%ntT,self%nHarm])
  end associate
end subroutine x2y_transpose

subroutine fill_t_guard_cell(self)
!------------ fill t guard cell ------------------
  implicit none
  include 'mpif.h'
  type(Radiation) :: self
  complex(8) :: sBuff(self%nx,0:self%ny,self%nHarm), rBuff(self%nx,0:self%ny,self%nHarm)
  integer :: n,ierr,req(2)
  
  associate(pGrd=>self%pGrd)
    sBuff = (0d0,0d0)
    rBuff = (0d0,0d0)
    n = self%nx*(self%ny+1)*self%nHarm
    !if(pGrd%myRow/=pGrd%npRow-1)
    sBuff(:,:,:) = self%Fld(:,:,1,:)
    call MPI_ISEND(sBuff,n,MPI_DOUBLE_COMPLEX,Pgrd%myL,0,Pgrd%comm_row,req(1),ierr)
    call MPI_IRECV(rBuff,n,MPI_DOUBLE_COMPLEX,Pgrd%myR,0,Pgrd%comm_row,req(2),ierr)
    call MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,ierr)
    if(pGrd%myRow/=pGrd%npRow-1)  self%Fld(:,:,self%nt,:) = rBuff
  end associate
end subroutine fill_t_guard_cell

subroutine fill_y_guard_cell(self)
!------------ fill y guard cell ------------------
  implicit none
  include 'mpif.h'
  type(Radiation) :: self
  integer :: n,ierr,req(4)
  complex(8) :: sBuff(self%nx,self%nt,self%nHarm,2),rBuff(self%nx,self%nt,self%nHarm,2)
  
  associate(pGrd=>self%pGrd)
    sBuff = (0d0,0d0)
    rBuff = (0d0,0d0)
    n = self%nx*self%nt*self%nHarm
    !if(pGrd%myCol/=pGrd%npCol-1)
    sBuff(:,:,:,1) = self%Fld(:,self%ny-1,:,:)
    !if(pGrd%myCol/=0)
    sBuff(:,:,:,2) = self%Fld(:,1,:,:)
    call MPI_ISEND(sBuff(:,:,:,2),n,MPI_DOUBLE_COMPLEX,pGrd%myD,0,pGrd%comm_col,req(1),ierr)
    call MPI_IRECV(rBuff(:,:,:,1),n,MPI_DOUBLE_COMPLEX,pGrd%myU,0,pGrd%comm_col,req(2),ierr)
    call MPI_ISEND(sBuff(:,:,:,1),n,MPI_DOUBLE_COMPLEX,pGrd%myU,0,pGrd%comm_col,req(3),ierr)
    call MPI_IRECV(rBuff(:,:,:,2),n,MPI_DOUBLE_COMPLEX,pGrd%myD,0,pGrd%comm_col,req(4),ierr)
    call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
    if(pGrd%myCol/=0)              self%Fld(:,0,:,:)      = rBuff(:,:,:,2)
    if(pGrd%myCol/=pGrd%npCol-1)   self%Fld(:,self%ny,:,:) = rBuff(:,:,:,1)
  end associate
end subroutine fill_y_guard_cell

subroutine fill_x_guard_cell(self)
!------------ fill y guard cell ------------------
  implicit none
  include 'mpif.h'
  type(Radiation) :: self
  integer :: n,ierr,req(4)
  complex(8) :: sBuff(self%nyT,self%ntT,self%nHarm,2),rBuff(self%nyT,self%ntT,self%nHarm,2)
  associate(pGrd=>self%pGrd)
    !call write_log('filling_x_guard')
    sBuff = (0d0,0d0)
    rBuff = (0d0,0d0)
    n = self%nyT*self%ntT*self%nHarm
    !if(pGrd%myCol/=pGrd%npCol-1)
    sBuff(:,:,:,1) = self%FldT(self%nxT-1,:,:,:)
    !if(pGrd%myCol/=0)
    sBuff(:,:,:,2) = self%FldT(1,:,:,:)
    !call write_log('fill x_guard comm')
    call MPI_ISEND(sBuff(:,:,:,2),n,MPI_DOUBLE_COMPLEX,pGrd%myD,0,pGrd%comm_col,req(1),ierr)
    call MPI_IRECV(rBuff(:,:,:,1),n,MPI_DOUBLE_COMPLEX,pGrd%myU,0,pGrd%comm_col,req(2),ierr)
    call MPI_ISEND(sBuff(:,:,:,1),n,MPI_DOUBLE_COMPLEX,pGrd%myU,0,pGrd%comm_col,req(3),ierr)
    call MPI_IRECV(rBuff(:,:,:,2),n,MPI_DOUBLE_COMPLEX,pGrd%myD,0,pGrd%comm_col,req(4),ierr)
    call MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,ierr)
    !call write_log('comm done')
    if(pGrd%myCol/=0)              self%FldT(0,:,:,:)      = rBuff(:,:,:,2)
    if(pGrd%myCol/=pGrd%npCol-1)   self%FldT(self%nxT,:,:,:)= rBuff(:,:,:,1)
  end associate
end subroutine fill_x_guard_cell

subroutine FldSolve(self,dz)
! ADI method. self%Src need to be prepared priori
  implicit none
  include 'mpif.h'
  class(Radiation)   :: self
  real*8, intent(in) :: dz
  logical, save :: flag_xy = .true.
  integer :: ierr
  associate(comm_2d => self%pGrd%comm_2d)
    if (flag_xy) then
      !call write_log('xADI')
      call mpi_barrier(comm_2d,ierr)
      call xADI(self,dz)
      !call write_log('y2x_transpose')
      call mpi_barrier(comm_2d,ierr)
      call y2x_transpose(self)
      !call write_log('yADI')
      call mpi_barrier(comm_2d,ierr)
      call yADI(self,dz)
      !call write_log('x2y_transpose')
      call mpi_barrier(comm_2d,ierr)
      call x2y_transpose(self)
      flag_xy = .false.
    else
      !call write_log('y2x_transpose')
      call mpi_barrier(comm_2d,ierr)
      call y2x_transpose(self)
      !call write_log('yADI')
      call mpi_barrier(comm_2d,ierr)
      call yADI(self,dz)
      !call write_log('x2y_transpose')
      call mpi_barrier(comm_2d,ierr)
      call x2y_transpose(self)
      call mpi_barrier(comm_2d,ierr)
      call fill_y_guard_cell(self)
      !call write_log('xADI')
      call mpi_barrier(comm_2d,ierr)
      call xADI(self,dz)
      flag_xy = .true.
    endif
    !call write_log('fill_y_guard_cell')
    call mpi_barrier(comm_2d,ierr)
    call fill_y_guard_cell(self)
    !call write_log('fill_t_guard_cell')
    call mpi_barrier(comm_2d,ierr)
    call fill_t_guard_cell(self)
  end associate
  self%Src = (0d0, 0d0)
end subroutine FldSolve

subroutine Slip(self,dz,reset_round)
  !=======================================================================
  ! Slippage through moving window
  !----------------------------------------------------------------------
  implicit none
  class(Radiation)   :: self
  real*8, intent(in) :: dz
  logical,intent(in), optional :: reset_round
  integer :: nt, tMshTabOld(0:self%pGrd%npRow)
  real*8 :: dt
  real*8, save :: round_dt = 0d0
  if(present(reset_round)) round_dt = 0d0
  
  associate(cDom => self%cDom, pGrd=>self%pGrd)
    tMshTabOld = cDom%tMshTab
    dt = self%ku*dz
    nt = nint((dt+round_dt)/cDom%MshSize(3))
    round_dt = dt - nt*cDom%MshSize(3)
    
    if(cDom%tMshTab(1) - nt > 0 .and. nt/=0 ) then
      cDom%tMshTab(1:pGrd%npRow-1) = cDom%tMshTab(1:pGrd%npRow-1) - nt
      cDom%RngLoc(3,1) = cDom%Rng(3,1) +cDom%tMshTab(pGrd%myRow  )*cDom%MshSize(3)
      cDom%RngLoc(3,2) = cDom%Rng(3,1) +cDom%tMshTab(pGrd%myRow+1)*cDom%MshSize(3)
      cDom%MshNumLoc(3)= cDom%tMshTab(pGrd%myRow+1) -cDom%tMshTab(pGrd%myRow) +1
      call update_Radiation_tMsh(self,tMshTabOld)
    endif
    cDom%Rng   (3,:) = cDom%Rng   (3,:) + dt
    cDom%RngLoc(3,:) = cDom%RngLoc(3,:) + dt
  end associate
end subroutine Slip


function Get_Pwr_tLoc(self) result(Pwr)
  implicit none
  include 'mpif.h'
  type(Radiation)  :: self
  real*8           :: Pwr(self%ntT,self%nHarm)
  integer :: ih,it,h
  real*8  :: dummy(self%ntT,self%nHarm)
  
  associate(cDom => self%cDom, pGrd=>self%pGrd)
    Pwr = 0d0
    do ih=1,self%nHarm
      h = self%harm_tab(ih)
      do it=1,self%ntT
        dummy(it,ih) = 0.5d0*h*h*sum(abs(self%Fld(1:self%nx,1:self%ny2,it,ih))**2)
      enddo
    enddo
    dummy = dummy*cDom%MshSize(1)*cDom%MshSize(2)*(physConst%Me*self%ks)**2/physConst%Z
    call MPI_BARRIER(pGrd%comm_row,h)
    call MPI_AllReduce(dummy,Pwr,self%nHarm*self%ntT,MPI_DOUBLE,MPI_SUM,pGrd%comm_col,h)
  end associate
  
end function Get_Pwr_tLoc

function Get_Pwr(self) result(Pwr)
  implicit none
  include 'mpif.h'
  class(Radiation) :: self
  real*8           :: Pwr(self%cDom%MshNum(3),self%nHarm)
  integer :: irow,it,ih,nCount,ierr
  integer :: nSend(0:self%pGrd%npRow-1),displs(0:self%pGrd%npRow-1)
  real*8  :: PwrChunk(self%ntT,self%nHarm),RecvBuff(self%cDom%MshNum(3)*self%nHarm)
  
  associate(cDom => self%cDom, pGrd=>self%pGrd)
    PwrChunk = Get_Pwr_tLoc(self)
    displs(0)=0
    do irow=0,pGrd%npRow-2
      nSend(irow)=cDom%tMshTab(irow+1)-cDom%tMshTab(irow)
      displs(irow+1) = displs(irow) + nSend(irow)
    enddo
    irow=pGrd%npRow-1
    nSend(irow)=cDom%tMshTab(irow+1)-cDom%tMshTab(irow)+1
     
    call MPI_AllGatherV(PwrChunk,nSend(pGrd%myRow) *self%nHarm,MPI_DOUBLE,&
                        RecvBuff,nSend*self%nHarm,displs*self%nHarm,MPI_DOUBLE,pGrd%comm_row,ierr)
    nCount = 1
    do irow=0,pGrd%npRow-1
      do ih=1,self%nHarm
        do it=displs(irow)+1,displs(irow)+nSend(irow)
          Pwr(it,ih) = RecvBuff(nCount)
          nCount = nCount+1
        enddo
      enddo
    enddo
  end associate
end function Get_Pwr

function Get_energy(self) result(energy)
  implicit none
  include 'mpif.h'
  class(Radiation) :: self
  real*8           :: energy(self%nHarm)
  real*8           :: Pwr(self%cDom%MshNum(3),self%nHarm)
  pwr = self%get_pwr()
  energy = sum(pwr,1)
  energy = energy*self%dt/(self%ks*physConst%c)
end function

end module RadiationModule
