module pGrid2dModule
  implicit none
  private
  type, public :: pGrid2D
    integer :: comm_2d,&          ! comunicator for entire grid
               comm_col,&         ! comunicator for my col
               comm_row,&         ! comunicator for my row
               np,&               ! total number of processors
               npCol,&            ! number of processors in my col
               npRow,&            ! number of processors in my row
               myRank,&           ! my rank in the grid comm
               myCol,&            ! my column number
               myRow,&            ! my row number
               myU, myD, myR, myL ! adjacent processors
    integer, allocatable :: rank2cart(:,:)  ! table of row,col rank from comm_2d rank
    contains
      procedure :: destruct => pGrid2D_destructor
  end type pGrid2D
  interface pGrid2D
    procedure pGrid2D_constructor
  end interface
contains
! set up 2-D Cartisian grid.
function pGrid2D_constructor(comm,nrow,ncol) result(pGrd)
  implicit none
  include 'mpif.h'
  class(pGrid2D), pointer :: pGrd
  integer, intent(in) :: comm,nrow,ncol
  integer, dimension(2) :: dims,local
  integer :: irank,ierr
  integer, allocatable :: seed(:)
  logical, dimension(2) :: period,remaindims
  allocate(pGrd)
  ! set up global grid information.
  call MPI_COMM_SIZE(comm,pGrd%np,ierr)

  ! set up the size of 2-D logical processor grid.
  pGrd%npCol = ncol
  pGrd%npRow = nrow

  ! create 2-D grid communicator.
  dims(1) = ncol
  dims(2) = nrow
  period(1) = .false.
  period(2) = .false.

  call MPI_CART_CREATE(comm,2,dims,period,.true., &
                       pGrd%comm_2d,ierr)

  ! get my rank, my column number and my row number in grid.
  call MPI_COMM_RANK(pGrd%comm_2d, pGrd%myRank, ierr)
  call MPI_CART_COORDS(pGrd%comm_2d,pGrd%myRank,2,local,ierr)
  pGrd%myCol = local(1)
  pGrd%myRow = local(2)

  ! set up column communicators.
  remaindims(1) = .true.
  remaindims(2) = .false.
  call MPI_CART_SUB(pGrd%comm_2d,remaindims,pGrd%comm_col,ierr)

  ! set up row communicators.
  remaindims(1) = .false.
  remaindims(2) = .true.
  call MPI_CART_SUB(pGrd%comm_2d,remaindims,pGrd%comm_row,ierr)
  
  if(pGrd%myCol == 0) then
    pGrd%myD = MPI_PROC_NULL
  else
    pGrd%myD = pGrd%myCol-1
  endif
  if(pGrd%myCol == pGrd%npCol-1) then
    pGrd%myU = MPI_PROC_NULL
  else
    pGrd%myU = pGrd%myCol+1
  endif
  
  if(pGrd%myRow == 0) then
    pGrd%myL = MPI_PROC_NULL
  else
    pGrd%myL = pGrd%myRow-1
  endif
  if(pGrd%myRow == pGrd%npRow-1) then
    pGrd%myR = MPI_PROC_NULL
  else
    pGrd%myR = pGrd%myRow+1
  endif
  
  
  allocate(pGrd%rank2cart(2,0:pGrd%np-1))
  do irank=0, pGrd%np-1
    call MPI_Cart_Coords(pGrd%comm_2d,irank,2,pGrd%rank2cart(:,irank),ierr)
  enddo
  
  call random_seed(size = irank)
  allocate(seed(irank))
  call random_seed(get=seed)
  seed = seed*(pGrd%myRank+1)
  call random_seed(put=seed)
end function pGrid2D_constructor


subroutine pGrid2D_destructor(self)
  class(pGrid2D), intent(inout) :: self
  if(allocated(self%rank2cart))  deallocate(self%rank2cart)
end subroutine pGrid2D_destructor


end module pGrid2dModule
