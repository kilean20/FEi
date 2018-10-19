module fileIO_Module
  implicit none
  private 
  interface logIO
    procedure logIO_constructor
  end interface 
  
  type,public :: logIO  ! file IO in append mode for log.
                 ! close file at each append.  Useful for debugging.
                 ! (when file is not closed before program exit due to error, file is not guranteed to be wrote)
    integer :: unit
    character(len=20) :: fname
    contains
      procedure :: append => write_dble, write_int, write_char
      procedure :: reopen
      procedure :: close => logIO_close
  end type

  type :: fileIO_Template
    contains
      procedure :: get_free_unit, save, load
      procedure :: num2str => num2str_int, num2str_dble, num2str_darr, num2str_iarr
  end type
  
  type(fileIO_Template), public :: fileIO
  
contains

function logIO_constructor(fname,fID)
  class(logIO), pointer :: logIO_constructor
  character(len=*),intent(in) :: fname
  integer,optional,intent(in) :: fID
  allocate(logIO_constructor)
  logIO_constructor%unit = fileIO%get_free_unit()
  logIO_constructor%fname = fname
  if(present(fID))  logIO_constructor%fname=fname//trim(fileIO%num2str(fID))
  open(logIO_constructor%unit,file=logIO_constructor%fname)
  close(logIO_constructor%unit)
end function logIO_constructor

subroutine  logIO_close(self)
  class(logIO) :: self
  close(self%unit)
end subroutine  logIO_close

integer function get_free_unit(self,initUnit)
  implicit none 
  class(fileIO_Template) :: self
  integer, optional, intent(in) :: initUnit
  integer :: lun 
  logical :: file_open 
  if(present(initUnit)) then
    lun = initUnit
  else
    lun = 10
  endif
  file_open = .true. 
  do while ( file_open ) 
    lun = lun + 1 
    inquire( lun, opened = file_open ) 
  end do 
  get_free_unit = lun 
  return 
end function get_free_unit

subroutine reopen(self)
  implicit none
  class(logIO),intent(in) :: self
  integer :: ifail
  open(self%unit, file=self%fname, access='SEQUENTIAL', status='OLD', position='APPEND', iostat=ifail) 
  if(ifail /= 0)  STOP '--- Error in opening log file in APPEND mode ---'
end subroutine

subroutine write_dble(self,d)
  implicit none
  class(logIO),intent(in) :: self
  real*8, intent(in) :: d(:)
  call reopen(self)
  write(self%unit,*) d
  close(self%unit)
end subroutine

subroutine write_char(self,c)
  implicit none
  class(logIO),intent(in) :: self
  character(len=*), intent(in) :: c
  call reopen(self)
  write(self%unit,*) c
  close(self%unit)
end subroutine

subroutine write_int(self,c)
  implicit none
  class(logIO),intent(in) :: self
  integer, intent(in) :: c
  call reopen(self)
  write(self%unit,*) c
  close(self%unit)
end subroutine

function num2str_dble(self,num)  
  implicit none
  class(fileIO_Template) :: self
  real*8, intent(in) :: num
  character(len=10) :: num2str_dble
  character(len=*), parameter  :: fmt_ = "(ES10.3)"
  write(num2str_dble,fmt_) num
end function

function num2str_darr(self,num,n)
  implicit none
  class(fileIO_Template) :: self
  integer, intent(in) :: n
  real*8, dimension(n), intent(in) :: num
  character(len=20*n) :: num2str_darr
  character(len=12),dimension(n) :: tempStr
  integer :: i
  do i=1,n
    tempStr(i) = num2str_dble(self,num(i))
  enddo
  write(num2str_darr,*) tempStr
  num2str_darr = adjustl(num2str_darr)
end function

function num2str_int(self,num)
  implicit none
  class(fileIO_Template) :: self
  integer, intent(in) :: num
  character(len=10) :: num2str_int
  character(len=*), parameter  :: fmt_ = "(I8)"
  write(num2str_int,fmt_) num
  num2str_int = adjustl(num2str_int)
end function

function num2str_iarr(self,num,n)
  implicit none
  class(fileIO_Template) :: self
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: num
  character(len=20*n) :: num2str_iarr
  character(len=10),dimension(n) :: tempStr
  integer :: i
  do i=1,n
    tempStr(i) = num2str_int(self,num(i))
  enddo
  write(num2str_iarr,*) tempStr
  num2str_iarr = adjustl(num2str_iarr)
end function
  

  
subroutine save(self,fname,arr,arrsize)
  ! save arbitrary array to file in txt format
  implicit none
  class(fileIO_Template) :: self
  character(len=*), intent(in) :: fname
  integer, intent(in) :: arrsize
  real*8,  intent(in) :: arr(*)
  integer :: u
  u = fileIO%get_free_unit() 
  open(u,file=fname)
  write(u,*) arr(1:arrsize)
  close(u)
end subroutine save

subroutine load(self,fname,arr,arrsize)
  ! save arbitrary array to file in txt format
  implicit none
  class(fileIO_Template) :: self
  character(len=*), intent(in) :: fname
  integer, intent(in) :: arrsize
  real*8,  intent(out) :: arr(*)
  integer :: u
  u = fileIO%get_free_unit() 
  open(u,file=fname)
  read(u,*) arr(1:arrsize)
  close(u)
end subroutine load

end module fileIO_Module
