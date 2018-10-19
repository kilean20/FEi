program test_run
  use ImpactFEL
  
  implicit none
  include 'mpif.h'
  integer :: i, j, n, m, wID, shotID
  integer :: npRow, npCol, distID(6), nLet, nBin, nBucket, dz_lu
  integer :: nPeriod, nStep, nTail
  integer :: nHarm, harm_tab(1)
  integer :: MshNum(3), distID2D(2),iUnit
  real*8  :: gamma0, K, kx, ky, ks, ku, ls, lu
  real*8  :: sigmaX,sigmaY,sigmaG,alphax,alphay,rangeT,sigmaZ,sigmaT
  real*8  :: currPk, Std2D(2),Std(3), alpha(3), Emit(3)
  real*8  :: dgrid, dz,z, meanx, meany, stdx, stdy
  real*8,dimension(6) :: Mean
  real*8,dimension(3,2) :: Rng
  real*8,allocatable :: pwr(:,:)
  integer, parameter :: nSample = 1024
  real*8 :: pData_sample(nSample,7)
  
  type(pGrid2D), pointer :: pGrd
  type(CompDom), pointer :: cDom
  type(Radiation), pointer :: Rad
  type(eBeam),pointer :: Beam
  
  integer,parameter :: nFODO = 2
  type(arrElem) :: lattice(99*2*nFODO)
  type(undulator), pointer :: wigl
  type(propagator),pointer :: Pro
  type(drift),     pointer :: d0,d1
  type(quad),      pointer :: QF,QD
  
!! control parameters
  npRow  = 1
  npCol  = 1
  shotID = 1
  wID    = 10
  nBin   = 4
  nLet   = 10000/(npRow*npCol*nBin)
  nBucket= 5
  dz_lu  = 5
  nPeriod= 165*2*nFODO + 20
  nTail  = 50
  nStep  = nPeriod/dz_lu
  gamma0 = 2.4d9/physConst%Me
  nHarm  = 1
  harm_tab = [1]
  dgrid  = 5d-4
  MshNum = [51,51,2*nTail]
  Rng(:2,1)=-dgrid
  Rng(:2,2)= dgrid
  Rng(3,1) =-twopi*(nTail-0.5)*nBucket
  Rng(3,2) = twopi*(nTail-0.5)*nBucket

!! undulator parameters
  lu = 2d-2
  ku = twopi/lu
  ls = 1d-9
  ks = twopi/ls
  K  = sqrt(4d0*gamma0**2*ku/ks-2d0)
  kx = 0d0
  ky = ku
  dz = nPeriod*lu/dble(nStep) 


!! electron beam parameters
  distID = [2,2,2,2,1,2]
  
  sigmaX = 4.07895d-5 
  sigmaY = 3.05758d-5 
  alphax = 1.3223
  alphay =-0.809419
  
  sigmaG = 0.2935421
  RangeT = twopi*(nBucket*2*nTail)
  sigmaZ = rangeT/sqrt(12d0)/ks
  sigmaT = sigmaZ*ks

  Mean  = [0d0,0d0,0d0,0d0,0d0,gamma0]
  Std   = [sigmaX,sigmaY,sigmaT]
  Alpha = [alphax,alphay,0d0]
  Emit  = [6.0d-7,6.0d-7,sigmaG*sigmaT]
  
  Std2D = [Std(2),Std(3)]
  distID2D = [distID(y_),distID(t_)]
    
  currPk= 500d0

!! initilaize beam
  call MPI_INIT(i)
  pGrd => pGrid2D(MPI_COMM_WORLD,npRow,npCol)
  cDom => CompDom(pGrd,Rng,MshNum,distID2D,Std2D,nPeriod)
  Rad  => Radiation(pGrd,cDom,harm_tab,nHarm,ks,ku)
  Beam => eBeam(pGrd,cDom,nLet,nBin,ks,ku,currPk,&
                distID,Mean,Std,Alpha,Emit,shotID=shotID,iHamsl=2)

!! check inital pData
  call Beam%get_sample(pData_sample,nSample)
  if(pGrd%myRank==0) then
    iUnit = fileIO%get_free_unit()
    open(unit=iUnit,file='pDataSample.txt')
    do i=1,nSample
      write(iUnit,*) pData_sample(i,:)
    enddo
  close(iUnit)
  endif
 

!! initilaize lattice (FODO cell)
  Wigl => undulator(5,1,K,kx,ky,ku,harm_tab,nHarm) ! *33
  QF => quad(L=0.1d0,nStep=1,B1= 7.86435d0)        ! *2
  QD => quad(L=0.1d0,nStep=1,B1=-7.19682d0)        ! *2
  D0 => drift(L=0.1d0,nStep=1)
  D1 => drift(L=0.1d0,nStep=1)                     ! *8
  Pro=> propagator(L = 0.1d0, nStep = 1)           ! *11
  n=1
  m=1
  do i=1,nFODO
    do j=1,33
      lattice(n)%op => Wigl
      n=n+1
      lattice(n)%op => write_csv(m)
      n=n+1
      m=m+1
    enddo
    lattice(n)%op => D0
    n=n+1
    lattice(n)%op => Pro
    n=n+1
    lattice(n)%op => write_csv(m)
    n=n+1
    m=m+1
    do j=1,2
      lattice(n)%op => QD
      n=n+1
      lattice(n)%op => Pro
      n=n+1
      lattice(n)%op => write_csv(m)
      n=n+1
      m=m+1
    enddo
    do j=1,8
      lattice(n)%op => D1
      n=n+1
      lattice(n)%op => Pro
      n=n+1
      lattice(n)%op => write_csv(m)
      n=n+1
      m=m+1
    enddo
    do j=1,33
      lattice(n)%op => Wigl
      n=n+1
      lattice(n)%op => write_csv(m)
      n=n+1
      m=m+1
    enddo
    lattice(n)%op => D0
    n=n+1
    lattice(n)%op => Pro
    n=n+1
    lattice(n)%op => write_csv(m)
    n=n+1
    m=m+1
    do j=1,2
      lattice(n)%op => QF
      n=n+1
      lattice(n)%op => Pro
      n=n+1
      lattice(n)%op => write_csv(m)
      n=n+1
      m=m+1
    enddo
    do j=1,8
      lattice(n)%op => D1
      n=n+1
      lattice(n)%op => Pro
      n=n+1
      lattice(n)%op => write_csv(m)
      n=n+1
      m=m+1
    enddo
  enddo


!! prepare pwr output
  allocate(pwr(Rad%cDom%MshNum(3),Rad%nHarm))
  z = 0d0
  do i=1,99*2*nFODO
    call lattice(i)%op%act(Beam,Rad,dz)
    pwr(:,:) = Rad%get_pwr()
    if(pGrd%myRank==0) print*, sum(pwr(:,1))
  enddo

  call MPI_finalize(i)
end program test_run
