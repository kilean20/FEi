module DiagnoElemModule
  use ElementModule
  use ConstDataModule
  use eBeamModule
  use RadiationModule
  implicit none
  private
  type, public, extends(Element) :: write_field_csv
    integer :: fileID
    contains
      procedure :: act_rad => write_field_csv_func
      !procedure :: get_M => drift_get_M
      !======== dummy overloadings =========
      procedure :: act_Beam => dummy_Beam_write_field_csv
      procedure :: act_BeamRad => dummy_BeamRad_write_field_csv
      !====== end of dummy overloading =====
  end type write_field_csv
  
  type, public, extends(Element) :: write_ptcl_csv
    integer :: fileID
    contains
      procedure :: act_beam => write_ptcl_csv_func
      !======== dummy overloadings =========
      procedure :: act_rad => dummy_Rad_write_ptcl_csv
      procedure :: act_BeamRad => dummy_BeamRad_write_ptcl_csv
      !====== end of dummy overloading =====
  end type write_ptcl_csv
  
  type, public, extends(Element) :: write_csv
    integer :: fileID
    contains
      procedure :: act_BeamRad => write_csv_func
      !======== dummy overloadings =========
      procedure :: act_rad => dummy_Rad_write_csv
      procedure :: act_beam => dummy_Beam_write_csv
      !====== end of dummy overloading =====
  end type write_csv
  
  interface write_field_csv
    procedure write_field_csv_constructor
  end interface
  
  interface write_ptcl_csv
    procedure write_ptcl_csv_constructor
  end interface
  
  interface write_csv
    procedure write_csv_constructor
  end interface

contains

function write_field_csv_constructor(fID)
  class(write_field_csv),pointer :: write_field_csv_constructor
  integer,intent(in) :: fID
  allocate(write_field_csv_constructor)
  write_field_csv_constructor%L  = 0.0
  write_field_csv_constructor%nStep = 1
  write_field_csv_constructor%dz = 0.0
  write_field_csv_constructor%flagBeam = .false.
  write_field_csv_constructor%flagRad = .true.
  write_field_csv_constructor%fileID = fID
end function write_field_csv_constructor

function write_ptcl_csv_constructor(fID)
  class(write_ptcl_csv),pointer :: write_ptcl_csv_constructor
  integer,intent(in) :: fID
  allocate(write_ptcl_csv_constructor)
  write_ptcl_csv_constructor%L  = 0.0
  write_ptcl_csv_constructor%nStep = 1
  write_ptcl_csv_constructor%dz = 0.0
  write_ptcl_csv_constructor%flagBeam = .true.
  write_ptcl_csv_constructor%flagRad = .false.
  write_ptcl_csv_constructor%fileID = fID
end function write_ptcl_csv_constructor

function write_csv_constructor(fID)
  class(write_csv),pointer :: write_csv_constructor
  integer,intent(in) :: fID
  allocate(write_csv_constructor)
  write_csv_constructor%L  = 0.0
  write_csv_constructor%nStep = 1
  write_csv_constructor%dz = 0.0
  write_csv_constructor%flagBeam = .true.
  write_csv_constructor%flagRad = .true.
  write_csv_constructor%fileID = fID
end function write_csv_constructor

subroutine write_field_csv_Func(self,Rad,dz)
  use fileIO_module
  implicit none
  include 'mpif.h'
  class(write_field_csv),intent(in) :: self
  type(Radiation),intent(inout):: Rad
  real*8, optional,intent(in) :: dz
  
  integer :: i,j,k,iH,ip,nx,ny,nt,nH,nMax,myIndex
  integer :: myCol,myRow,root,ierr,iUnit
  integer,dimension(0:Rad%Pgrd%np-1) :: nyLocList,ntLocList
  real*8 :: x,y,ct,amp,phs,eps
  complex(8), allocatable :: RecvBuff(:)
  character(7) :: cfID
  print*, 'calling write_field_csv_Func',self%fileID
  root = Rad%pGrd%np-1
  myCol = Rad%pGrd%myCol
  myRow = Rad%pGrd%myRow
  nH = Rad%nHarm
  call MPI_GATHER(Rad%ny2,1,MPI_INTEGER,nyLocList,1,&
                  MPI_INTEGER,root,Rad%pGrd%comm_2d,ierr)
  call MPI_GATHER(Rad%ntT,1,MPI_INTEGER,ntLocList,1,&
                  MPI_INTEGER,root,Rad%pGrd%comm_2d,ierr)
  nx=Rad%nx  

  if(Rad%pGrd%myRank==root) then
    iUnit = fileIO%get_free_unit()
    write(cfID,'(I0)') self%fileID
    open(iUnit,file='field.'//trim(adjustl(cfID))//'.csv',action='write')
    write(iUnit,'("x, y, ct")',advance="no") 
    do iH=1,nH-1
      write(iUnit,'(", amplitude",I1,", phase",I1,", log_amplitude",I1)',advance="no")&
           Rad%harm_tab(iH),Rad%harm_tab(iH),Rad%harm_tab(iH)
    enddo
    write(iUnit,'(", amplitude",I1,", phase",I1,", log_amplitude",I1)')&
         Rad%harm_tab(nH),Rad%harm_tab(nH),Rad%harm_tab(nH)
    eps = 1.0d-30
    do i=1,nx
      x = Rad%cDom%RngLoc(1,1) +(i-1)*Rad%dx
      do j=1,Rad%ny2
        y = Rad%cDom%RngLoc(2,1) + (j-1)*Rad%dy
        do k=1,Rad%ntT
          ct = Rad%cDom%RngLoc(3,1) + (k-1)*Rad%dt
          ct = ct/Rad%ks
          write(iUnit,'(3(ES14.7,","))',advance="no") x,y,ct
          do iH=1,nH-1
            amp = abs(Rad%Fld(i,j,k,iH))
            phs = atan2(aimag(Rad%Fld(i,j,k,iH)),real(Rad%Fld(i,j,k,iH)))
            write(iUnit,'(3(ES14.7,","))',advance="no") amp,phs,log(amp+eps)
          enddo
          amp = abs(Rad%Fld(i,j,k,nH))
          phs = atan2(aimag(Rad%Fld(i,j,k,nH)),real(Rad%Fld(i,j,k,nH)))
          write(iUnit,'(2(ES14.7,","),ES14.7)') amp,phs,log(amp+eps)
        enddo
      enddo
    enddo
    nMax = maxval(nyLocList*ntLocList)
    allocate(RecvBuff(nx*nMax*nH))
    do ip=0,root-1 !recall that root = np-1
      ny=nyLocList(ip)
      nt=ntLocList(ip)
      call MPI_RECV(RecvBuff,nx*ny*nt*nH,MPI_DOUBLE_COMPLEX,&
                    ip,1,Rad%pGrd%comm_2d,MPI_STATUS_IGNORE,ierr)
      do i=1,nx
        x = Rad%cDom%Rng(1,1) +(i-1)*Rad%dx
        do j=1,ny
          y = Rad%cDom%Rng(2,1) +(Rad%cDom%yMshTab(myCol)+j-1)*Rad%dy
          do k=1,nt
            ct = Rad%cDom%Rng(3,1) +(Rad%cDom%tMshTab(myRow)+k-1)*Rad%dt
            ct = ct/Rad%ks
            write(iUnit,'(3(ES14.7,","))',advance="no") x,y,ct
            do iH=1,nH-1
              myIndex = i+(j-1)*nx+(k-1)*nx*ny+(iH-1)*nx*ny*nt
              amp = abs(RecvBuff(myIndex))
              phs = atan2(aimag(RecvBuff(myIndex)),real(RecvBuff(myIndex)))
              write(iUnit,'(3(ES14.7,","))',advance="no") amp,phs,log(amp+eps)
            enddo
            myIndex = i+(j-1)*nx+(k-1)*nx*ny+(nH-1)*nx*ny*nt
            amp = abs(RecvBuff(myIndex))
            phs = atan2(aimag(RecvBuff(myIndex)),real(RecvBuff(myIndex)))
            write(iUnit,'(2(ES14.7,","),ES14.7)') amp,phs,log(amp+eps)
          enddo
        enddo
      enddo
    enddo
  else
    call MPI_SEND(Rad%Fld(:,1:Rad%ny2,1:Rad%ntT,:),&
                  nx*Rad%ny2*Rad%ntT*nH,MPI_DOUBLE_COMPLEX,&
                  root,1,Rad%pGrd%comm_2d,ierr)
  endif
end subroutine write_field_csv_Func

subroutine write_ptcl_csv_Func(self,beam,dz)
  use fileIO_module
  implicit none
  include 'mpif.h'
  class(write_ptcl_csv),intent(in) :: self
  type(eBeam),intent(inout) :: beam
  real*8,optional,intent(in) :: dz

  integer :: i,nMax,ip,root,ierr,iUnit
  integer,dimension(0:beam%Pgrd%np-1) :: nptList
  real*8, allocatable :: RecvBuff(:)
  character(7) :: cfID
  
  root = beam%pGrd%np-1
  call MPI_GATHER(beam%npt,1,MPI_INTEGER,nptList,1,&
                  MPI_INTEGER,root,beam%pGrd%comm_2d,ierr)

  if(beam%pGrd%myRank==root) then
    iUnit = fileIO%get_free_unit()
    write(cfID,'(I0)') self%fileID
    open(iUnit,file='ptcl.'//trim(adjustl(cfID))//'.csv',action='write')
    write(iUnit,'("x, y, ct, charge")') 

    do i=1,beam%npt
      write(iUnit,'(3(ES14.7,","),ES14.7)') beam%pData(i,[x_,y_]),&
                                            beam%pData(i,t_)/beam%ks,&
                                            beam%pData(i,q_)
    enddo
    nMax = maxval(nptList)
    allocate(RecvBuff(nMax*4))
    do ip=0,root-1 !recall that root = np-1
      call MPI_RECV(RecvBuff,nptList(ip)*4,MPI_DOUBLE_PRECISION,&
                    ip,1,beam%pGrd%comm_2d,MPI_STATUS_IGNORE,ierr)
      do i=1,beam%npt
        write(iUnit,'(3(ES14.7,","),ES14.7)') RecvBuff((i-1)*4+1:(i-1)*4+1), &
                                              RecvBuff((i-1)*4+3)/beam%ks,&
                                              RecvBuff(i*4)
      enddo
    enddo
  else
    call MPI_SEND(beam%pData(:,[x_,y_,t_,q_]),beam%npt*4,&
                  MPI_DOUBLE_PRECISION,root,1,beam%pGrd%comm_2d,ierr)
  endif
end subroutine write_ptcl_csv_Func

subroutine write_csv_Func(self,beam,rad,dz)
  use fileIO_module
  implicit none
  class(write_csv),intent(in) :: self
  type(eBeam),intent(inout) :: beam
  type(Radiation),intent(inout) :: rad
  real*8,optional,intent(in) :: dz
  
  type(write_ptcl_csv), pointer :: ptclCSV
  type(write_field_csv), pointer :: fieldCSV
  
  ptclCSV => write_ptcl_csv(self%fileID)
  fieldCSV=> write_field_csv(self%fileID)
  call fieldCSV%act(Beam,Rad,dz)
  call ptclCSV %act(Beam,Rad,dz)
  
end subroutine write_csv_Func
  
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
subroutine dummy_Beam_write_field_csv(self,beam,dz)
  implicit none
  class(write_field_csv),intent(in) :: self
  type(eBeam),intent(inout):: beam
  real*8,optional,intent(in)   :: dz
end subroutine

subroutine dummy_BeamRad_write_field_csv(self,beam,rad,dz)
  implicit none
  class(write_field_csv),intent(in) :: self
  type(eBeam),intent(inout)    :: beam
  type(Radiation),intent(inout):: rad
  real*8,optional,intent(in)   :: dz
end subroutine

subroutine dummy_Rad_write_ptcl_csv(self,rad,dz)
  implicit none
  class(write_ptcl_csv),intent(in) :: self
  type(Radiation),intent(inout):: rad
  real*8,optional,intent(in)   :: dz
end subroutine

subroutine dummy_BeamRad_write_ptcl_csv(self,beam,rad,dz)
  implicit none
  class(write_ptcl_csv),intent(in) :: self
  type(eBeam),intent(inout)    :: beam
  type(Radiation),intent(inout):: rad
  real*8,optional,intent(in)   :: dz
end subroutine

subroutine dummy_Beam_write_csv(self,beam,dz)
  implicit none
  class(write_csv),intent(in) :: self
  type(eBeam),intent(inout):: beam
  real*8,optional,intent(in)   :: dz
end subroutine

subroutine dummy_Rad_write_csv(self,rad,dz)
  implicit none
  class(write_csv),intent(in) :: self
  type(Radiation),intent(inout):: rad
  real*8,optional,intent(in)   :: dz
end subroutine

end module DiagnoElemModule

