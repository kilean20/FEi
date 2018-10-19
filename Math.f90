module MathModule
  implicit none
  private
  real*8, parameter :: pi    = 3.141592653589793d0
  real*8, parameter :: twopi = 6.283185307179586d0
  type :: mathTemplate
    contains
      procedure :: hamsl,iErf,cdfn,randn,randp
      procedure :: randID,cdfID,get_density_at
      procedure :: shuffle,sort
      procedure :: TriDiag,aTriDiag
      procedure :: std
      procedure :: rande => rande_scalar, rande_vector
      procedure :: linear_interpolation
  end type
  type(mathTemplate), public :: math
  ! === example use ==================================
  ! x = math%randn()             ! normal distribution
  ! x_std = math%std(x,size(x))  ! get s.t.d. of x_std
  ! ==================================================
contains

real*8 function hamsl(self,j,at)
!==================================================================
! same with hammersley except that the sequence now start from "at"
!------------------------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer, optional, intent(in) :: j,at
  integer, parameter            :: jmax=10
  real*8                        :: xs(jmax),xsi(jmax)
  integer                       :: i(jmax),nbase(jmax),i1(jmax),i2(jmax)
  data nbase/2,3,5,7,11,13,17,19,23,29/
  data i/jmax*1/
  if((.not. present(j)) .or. (.not. present(at))) then
    call random_number(hamsl) 
    return
  endif
  if(j<1) then
    call random_number(hamsl) 
    return
  endif
  if(at<=0) then
    call random_number(hamsl) 
    return
  endif
  if(at>i(j)) i(j) = at  ! starting index
  xs (j)=0.d0
  xsi(j)=1.0d0
  i2(j)= i(j)
  do
    xsi(j)=xsi(j)/float(nbase(j))
    i1 (j)= i2(j)/nbase(j)
    xs (j)= xs(j)+(i2(j)-nbase(j)*i1(j))*xsi(j)
    i2 (j)= i1(j)
    if(i2(j)<=0) exit
  enddo
  hamsl=xs(j)
  i(j)= i(j)+1
  return
end function 

real*8 function iErf(self,x)
! ==========================================================
! inverted error function
! original author:  Mark Law
! reference: https://github.com/markthelaw/GoStatHelper
! ------------------------------------------------------------      
  implicit none
  class(mathTemplate) :: self
  real*8, intent(in) :: x
  real*8 :: y, w, p
  
  y = x
  if(x>1d0-1d-16) y=1d0-1d-16    ! divergence correction 
  if(x<-1d0+1d-16) y=-1d0+1d-16  ! divergence correction 
  w = -log((1d0-y)*(1d0+y))      ! 1.0 - x * x would induce rounding errors near the boundaries +/-1
  if(w<6.25d0) then
    w = w - 3.125
    p =  -3.6444120640178196996d-21
    p =   -1.685059138182016589d-19 + p * w
    p =   1.2858480715256400167d-18 + p * w
    p =    1.115787767802518096d-17 + p * w
    p =   -1.333171662854620906d-16 + p * w
    p =   2.0972767875968561637d-17 + p * w
    p =   6.6376381343583238325d-15 + p * w
    p =  -4.0545662729752068639d-14 + p * w
    p =  -8.1519341976054721522d-14 + p * w
    p =   2.6335093153082322977d-12 + p * w
    p =  -1.2975133253453532498d-11 + p * w
    p =  -5.4154120542946279317d-11 + p * w
    p =    1.051212273321532285d-09 + p * w
    p =  -4.1126339803469836976d-09 + p * w
    p =  -2.9070369957882005086d-08 + p * w
    p =   4.2347877827932403518d-07 + p * w
    p =  -1.3654692000834678645d-06 + p * w
    p =  -1.3882523362786468719d-05 + p * w
    p =    0.0001867342080340571352 + p * w
    p =  -0.00074070253416626697512 + p * w
    p =   -0.0060336708714301490533 + p * w
    p =      0.24015818242558961693 + p * w
    p =       1.6536545626831027356 + p * w 
  elseif(w<16d0) then
    w = sqrt(w) - 3.25d0
    p =   2.2137376921775787049d-09
    p =   9.0756561938885390979d-08 + p * w
    p =  -2.7517406297064545428d-07 + p * w
    p =   1.8239629214389227755d-08 + p * w
    p =   1.5027403968909827627d-06 + p * w
    p =   -4.013867526981545969d-06 + p * w
    p =   2.9234449089955446044d-06 + p * w
    p =   1.2475304481671778723d-05 + p * w
    p =  -4.7318229009055733981d-05 + p * w
    p =   6.8284851459573175448d-05 + p * w
    p =   2.4031110387097893999d-05 + p * w
    p =   -0.0003550375203628474796 + p * w
    p =   0.00095328937973738049703 + p * w
    p =   -0.0016882755560235047313 + p * w
    p =    0.0024914420961078508066 + p * w
    p =   -0.0037512085075692412107 + p * w
    p =     0.005370914553590063617 + p * w
    p =       1.0052589676941592334 + p * w
    p =       3.0838856104922207635 + p * w      
  else
    w =  sqrt(w) - 5d0
    p =  -2.7109920616438573243d-11
    p =  -2.5556418169965252055d-10 + p * w
    p =   1.5076572693500548083d-09 + p * w
    p =  -3.7894654401267369937d-09 + p * w
    p =   7.6157012080783393804d-09 + p * w
    p =  -1.4960026627149240478d-08 + p * w
    p =   2.9147953450901080826d-08 + p * w
    p =  -6.7711997758452339498d-08 + p * w
    p =   2.2900482228026654717d-07 + p * w
    p =  -9.9298272942317002539d-07 + p * w
    p =   4.5260625972231537039d-06 + p * w
    p =  -1.9681778105531670567d-05 + p * w
    p =   7.5995277030017761139d-05 + p * w
    p =  -0.00021503011930044477347 + p * w
    p =  -0.00013871931833623122026 + p * w
    p =       1.0103004648645343977 + p * w
    p =       4.8499064014085844221 + p * w
  endif
  iErf = p*x
end function      

real*8 function cdfn(self,low,up,lowCut,upCut)
  implicit none
  class(mathTemplate) :: self
  real*8, intent(in) :: low,up
  real*8, optional, intent(in) :: lowCut,upCut
  real*8 :: Erf_lc,Erf_hc,ErfL,ErfH
  if(present(lowCut)) then
    Erf_lc = Erf(lowCut/sqrt(2d0))
  else
    Erf_lc = -1d0
  endif
  if(present(upCut)) then
    Erf_hc = Erf(upCut/sqrt(2d0))
  else
    Erf_hc =  1d0
  endif
  ErfL = Erf(low/sqrt(2d0))
  ErfH = Erf(up/sqrt(2d0))
  cdfn = (ErfH-ErfL)/(Erf_hc-Erf_lc)
end function cdfn

real*8 function randn(self, jHamm, at, low, up)
!==========================================
! hammersley sampling of local gaussian distribution 
! starting at "at" using inverse error function
!------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer,optional, intent(in) :: jHamm,at
  real*8, optional, intent(in) :: low,up
  real*8                       :: erfL,erfH
  if(present(low)) then
    erfL = Erf(low/sqrt(2d0))
  else
    erfL = -1d0
  endif
  if(present(up)) then
    erfH = Erf(up/sqrt(2d0))
  else
    erfH = 1d0
  endif  
  randn = sqrt(2d0)*self%iErf((erfH-erfL)*self%hamsl(jHamm,at)+erfL)
end function randn

subroutine rande_scalar(self,x,mean)
!===================================================
!   random generator of exponential distribution. scalar version
!---------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  real*8, intent(in) :: mean
  real*8, intent(out):: x
  call random_number(x)
  x = -log(x)*mean
end subroutine

subroutine rande_vector(self,x,mean)
!===================================================
!   random generator of exponential distribution. vector version
!---------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  real*8, intent(in) :: mean
  real*8, intent(out):: x(:)
  call random_number(x)
  x = -log(x)*mean
end subroutine

integer function randP(self,lambda)
  !==================================================================
  ! Poisson random generator
  ! input :
  !   lambda = average
  !==================================================================
  implicit none
  class(mathTemplate) :: self
  real*8, intent(in) :: lambda
  !integer :: PoissonSmall, PoissonLarge
  if (lambda<30d0) then
    randP = PoissonSmall(lambda)
  else
    randP = PoissonLarge(lambda)
  endif
end function randP

integer function PoissonSmall(lambda)
  implicit none
  real*8, intent(in) :: lambda
  integer :: k
  real*8 :: p,L,dummy
  
  p = 1d0
  L = exp(-lambda)
  k=0
  do
    k = k+1
    call random_number(dummy)
    p = p*dummy
    if(p <= L) exit
  enddo
  PoissonSmall = k-1
end function PoissonSmall

integer function PoissonLarge(lambda)
  !"Rejection method PA" from "The Computer Generation of Poisson Random Variables" by A. C. Atkinson
  !Journal of the Royal Statistical Society Series C (Applied Statistics) Vol. 28, No. 1. (1979)
  !The article is on pages 29-35. The algorithm given here is on page 32.
  implicit none
  real*8, intent(in) :: lambda
  integer :: n
  real*8 :: c,beta,alpha,k,u,x,v,y,temp,lhs,rhs
  
  c = 0.767d0 - 3.36d0/lambda
  beta = pi/sqrt(3d0*lambda)
  alpha = beta*lambda
  k = log(c) - lambda - log(beta)

  do
    call random_number(u)
    x = (alpha - log((1d0 - u)/u))/beta
    n = floor(x + 0.5d0)
    if (n < 0)  cycle
    call random_number(v)
    y = alpha - beta*x
    temp = 1d0 + exp(y)
    lhs = y + log(v/(temp*temp))
    rhs = k + n*log(lambda) - LogFactorial(n)
    if (lhs <= rhs) then
      PoissonLarge = n
      return
    endif
  enddo
end function PoissonLarge

real*8 function randID(self,distID,jHamm,at,low,up)
!==================================================
! random number generator of distribution specified by distID
!--------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer,         intent(in) :: distID
  integer,optional,intent(in) :: jHamm,at
  real*8, optional,intent(in) :: low,up
  select case (distID)
    case (1) !uniform
      randID = sqrt(12d0)*(self%hamsl(jHamm,at)-0.5d0)
    case (2) !normal
      randID = self%randn(jHamm,at,low,up)
    case default
      randID = 0d0
      print*, 'Error :: math%randID -> unknown distID '
  end select
end function

real*8 function cdfID(self,distID,low,up,lowCut,upCut)
  implicit none
  class(mathTemplate) :: self
  integer, intent(in) :: distID  
  real*8,  intent(in) :: low,up,lowCut,upCut
  select case (distID)
    case (1)
      cdfID = (up-low)/(upCut-lowCut)
    case (2)
      cdfID = self%cdfn(low,up,lowCut,upCut)
    case default
      cdfID = 0d0
      print*, 'Error :: math%cdfID -> unknown distID '
  end select
end function cdfID

subroutine get_density_at(self,density,theta,n,distID,args)
!=============================================================================
! get density of given location and density function ID
! input
!   : theta(n) = locations
!     distID   = ID of density function
!     args     = argument needed for corresponding density function
!-----------------------------------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer, intent(in) :: distID,n
  real*8,  intent(in) :: theta(n)
  real*8,  intent(in) :: args(:)
  real*8,  intent(out):: density(n)
  logical,save :: first = .true.
  real*8 :: norm
  integer :: i
  density = 1d0
  select case(distID)
    case (1)
      if(size(args)<2) then
        print*, 'Error <- math%get_density_at :: args for distID=1 are not enough'
        density = 1d0
      endif
      forall(i=1:n, theta(i) < args(1)-sqrt(3d0)*args(2) &
               .or. theta(i) > args(1)+sqrt(3d0)*args(2))
        density(i) = 0d0
      end forall
      if(first) then
        first = .false.
        print*, 'args(1), sqrt(3d0)*args(2)'
      endif
      if(sum(density)<n) print*, 'density=0'
    case (2)
      if(size(args)<2) then
        print*, 'Error <- math%get_density_at :: args for distID=2 are not enough'
        density = 1d0
      endif
      norm = exp(-0.5d0*(sum(theta)/n-args(1))**2/(args(2)*args(2)))
      density = exp(-0.5d0*(theta-args(1))**2/(args(2)*args(2)))/norm
    case default
      print*, 'Error <- math%get_density_at :: unknown density ID'
  end select
end subroutine get_density_at


real*8 function LogFactorial(n)
  implicit none
  integer, intent(in) :: n
  real*8 :: x
  real*8, parameter :: PI = 3.14159265359
  real*8, dimension(255), parameter :: &
    lf=(/ 0.000000000000000d0,&
          0.000000000000000d0,&
          0.693147180559945d0,&
          1.791759469228055d0,&
          3.178053830347946d0,&
          4.787491742782046d0,&
          6.579251212010101d0,&
          8.525161361065415d0,&
          10.604602902745251d0,&
          12.801827480081469d0,&
          15.104412573075516d0,&
          17.502307845873887d0,&
          19.987214495661885d0,&
          22.552163853123421d0,&
          25.191221182738683d0,&
          27.899271383840894d0,&
          30.671860106080675d0,&
          33.505073450136891d0,&
          36.395445208033053d0,&
          39.339884187199495d0,&
          42.335616460753485d0,&
          45.380138898476908d0,&
          48.471181351835227d0,&
          51.606675567764377d0,&
          54.784729398112319d0,&
          58.003605222980518d0,&
          61.261701761002001d0,&
          64.557538627006323d0,&
          67.889743137181526d0,&
          71.257038967168000d0,&
          74.658236348830158d0,&
          78.092223553315307d0,&
          81.557959456115029d0,&
          85.054467017581516d0,&
          88.580827542197682d0,&
          92.136175603687079d0,&
          95.719694542143202d0,&
          99.330612454787428d0,&
          102.968198614513810d0,&
          106.631760260643450d0,&
          110.320639714757390d0,&
          114.034211781461690d0,&
          117.771881399745060d0,&
          121.533081515438640d0,&
          125.317271149356880d0,&
          129.123933639127240d0,&
          132.952575035616290d0,&
          136.802722637326350d0,&
          140.673923648234250d0,&
          144.565743946344900d0,&
          148.477766951773020d0,&
          152.409592584497350d0,&
          156.360836303078800d0,&
          160.331128216630930d0,&
          164.320112263195170d0,&
          168.327445448427650d0,&
          172.352797139162820d0,&
          176.395848406997370d0,&
          180.456291417543780d0,&
          184.533828861449510d0,&
          188.628173423671600d0,&
          192.739047287844900d0,&
          196.866181672889980d0,&
          201.009316399281570d0,&
          205.168199482641200d0,&
          209.342586752536820d0,&
          213.532241494563270d0,&
          217.736934113954250d0,&
          221.956441819130360d0,&
          226.190548323727570d0,&
          230.439043565776930d0,&
          234.701723442818260d0,&
          238.978389561834350d0,&
          243.268849002982730d0,&
          247.572914096186910d0,&
          251.890402209723190d0,&
          256.221135550009480d0,&
          260.564940971863220d0,&
          264.921649798552780d0,&
          269.291097651019810d0,&
          273.673124285693690d0,&
          278.067573440366120d0,&
          282.474292687630400d0,&
          286.893133295426990d0,&
          291.323950094270290d0,&
          295.766601350760600d0,&
          300.220948647014100d0,&
          304.686856765668720d0,&
          309.164193580146900d0,&
          313.652829949878990d0,&
          318.152639620209300d0,&
          322.663499126726210d0,&
          327.185287703775200d0,&
          331.717887196928470d0,&
          336.261181979198450d0,&
          340.815058870798960d0,&
          345.379407062266860d0,&
          349.954118040770250d0,&
          354.539085519440790d0,&
          359.134205369575340d0,&
          363.739375555563470d0,&
          368.354496072404690d0,&
          372.979468885689020d0,&
          377.614197873918670d0,&
          382.258588773060010d0,&
          386.912549123217560d0,&
          391.575988217329610d0,&
          396.248817051791490d0,&
          400.930948278915760d0,&
          405.622296161144900d0,&
          410.322776526937280d0,&
          415.032306728249580d0,&
          419.750805599544780d0,&
          424.478193418257090d0,&
          429.214391866651570d0,&
          433.959323995014870d0,&
          438.712914186121170d0,&
          443.475088120918940d0,&
          448.245772745384610d0,&
          453.024896238496130d0,&
          457.812387981278110d0,&
          462.608178526874890d0,&
          467.412199571608080d0,&
          472.224383926980520d0,&
          477.044665492585580d0,&
          481.872979229887900d0,&
          486.709261136839360d0,&
          491.553448223298010d0,&
          496.405478487217580d0,&
          501.265290891579240d0,&
          506.132825342034830d0,&
          511.008022665236070d0,&
          515.890824587822520d0,&
          520.781173716044240d0,&
          525.679013515995050d0,&
          530.584288294433580d0,&
          535.496943180169520d0,&
          540.416924105997740d0,&
          545.344177791154950d0,&
          550.278651724285620d0,&
          555.220294146894960d0,&
          560.169054037273100d0,&
          565.124881094874350d0,&
          570.087725725134190d0,&
          575.057539024710200d0,&
          580.034272767130800d0,&
          585.017879388839220d0,&
          590.008311975617860d0,&
          595.005524249382010d0,&
          600.009470555327430d0,&
          605.020105849423770d0,&
          610.037385686238740d0,&
          615.061266207084940d0,&
          620.091704128477430d0,&
          625.128656730891070d0,&
          630.172081847810200d0,&
          635.221937855059760d0,&
          640.278183660408100d0,&
          645.340778693435030d0,&
          650.409682895655240d0,&
          655.484856710889060d0,&
          660.566261075873510d0,&
          665.653857411105950d0,&
          670.747607611912710d0,&
          675.847474039736880d0,&
          680.953419513637530d0,&
          686.065407301994010d0,&
          691.183401114410800d0,&
          696.307365093814040d0,&
          701.437263808737160d0,&
          706.573062245787470d0,&
          711.714725802289990d0,&
          716.862220279103440d0,&
          722.015511873601330d0,&
          727.174567172815840d0,&
          732.339353146739310d0,&
          737.509837141777440d0,&
          742.685986874351220d0,&
          747.867770424643370d0,&
          753.055156230484160d0,&
          758.248113081374300d0,&
          763.446610112640200d0,&
          768.650616799717000d0,&
          773.860102952558460d0,&
          779.075038710167410d0,&
          784.295394535245690d0,&
          789.521141208958970d0,&
          794.752249825813460d0,&
          799.988691788643450d0,&
          805.230438803703120d0,&
          810.477462875863580d0,&
          815.729736303910160d0,&
          820.987231675937890d0,&
          826.249921864842800d0,&
          831.517780023906310d0,&
          836.790779582469900d0,&
          842.068894241700490d0,&
          847.352097970438420d0,&
          852.640365001133090d0,&
          857.933669825857460d0,&
          863.231987192405430d0,&
          868.535292100464630d0,&
          873.843559797865740d0,&
          879.156765776907600d0,&
          884.474885770751830d0,&
          889.797895749890240d0,&
          895.125771918679900d0,&
          900.458490711945270d0,&
          905.796028791646340d0,&
          911.138363043611210d0,&
          916.485470574328820d0,&
          921.837328707804890d0,&
          927.193914982476710d0,&
          932.555207148186240d0,&
          937.921183163208070d0,&
          943.291821191335660d0,&
          948.667099599019820d0,&
          954.046996952560450d0,&
          959.431492015349480d0,&
          964.820563745165940d0,&
          970.214191291518320d0,&
          975.612353993036210d0,&
          981.015031374908400d0,&
          986.422203146368590d0,&
          991.833849198223450d0,&
          997.249949600427840d0,&
          1002.670484599700300d0,&
          1008.095434617181700d0,&
          1013.524780246136200d0,&
          1018.958502249690200d0,&
          1024.396581558613400d0,&
          1029.838999269135500d0,&
          1035.285736640801600d0,&
          1040.736775094367400d0,&
          1046.192096209724900d0,&
          1051.651681723869200d0,&
          1057.115513528895000d0,&
          1062.583573670030100d0,&
          1068.055844343701400d0,&
          1073.532307895632800d0,&
          1079.012946818975000d0,&
          1084.497743752465600d0,&
          1089.986681478622400d0,&
          1095.479742921962700d0,&
          1100.976911147256000d0,&
          1106.478169357800900d0,&
          1111.983500893733000d0,&
          1117.492889230361000d0,&
          1123.006317976526100d0,&
          1128.523770872990800d0,&
          1134.045231790853000d0,&
          1139.570684729984800d0,&
          1145.100113817496100d0,&
          1150.633503306223700d0,&
          1156.170837573242400d0 /)

  if(n > 254) then
    x = n + 1
    LogFactorial = (x - 0.5)*log(x) - x + 0.5*log(2*PI) + 1.0/(12.0*x)
  else
    LogFactorial = lf(n+1)
  endif
end function LogFactorial

subroutine shuffle(self,arr,n,m)
!=============================================================================
! shuffle m integers out of n integers
!-----------------------------------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer, intent(in) :: n,m
  integer, intent(out):: arr(n)
  integer :: i,j,k
  real*8  :: r
  
  do i=1,n
    arr(i)=i
  enddo
  do i=1,m
    call random_number(r)
    k=int((n-i+1)*r)+i
    j=arr(k)
    arr(k)=arr(i)
    arr(i)=j
  enddo
    
end subroutine shuffle

recursive subroutine sort(self,data2sort, iCol4Sort, nRow, nCol, iRowStart, iRowEnd)
!==============================================================================
! quicksort algorithm
!------------------------------------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  real*8 data2sort(nRow, nCol), mid, temp(nCol)
  integer iCol4Sort, nCol, nRow, iRowStart, iRowEnd
  integer i, j
 
  mid = data2sort((iRowStart+iRowEnd)/2, iCol4Sort)
  i = iRowStart
  j = iRowEnd
  do
    do while (data2sort(i, iCol4Sort) < mid)
      i=i+1
    end do
    do while (mid < data2sort(j, iCol4Sort))
      j=j-1
    end do
    if (i >= j) exit
    temp = data2sort(i,:)
    data2sort(i,:) = data2sort(j,:)
    data2sort(j,:) = temp
    i=i+1
    j=j-1
  end do
  if (iRowStart < i-1) call self%sort(data2sort, iCol4Sort, nRow, nCol, iRowStart, i-1)
  if (j+1 < iRowEnd)   call self%sort(data2sort, iCol4Sort, nRow, nCol, j+1, iRowEnd)
end subroutine sort

real*8 function std(self,arr,n)
  !==================================================================
  ! get standard deviation of an array
  !==================================================================
  implicit none
  class(mathTemplate) :: self
  integer, intent(in) :: n
  real*8,intent(in) :: arr(n)
  real*8 :: ave
  ave = sum(arr)/n
  std = sqrt(sum((arr-ave)*(arr-ave))/(n-1))
end function

subroutine TriDiag(self,solution,Src,off_diag,ndim)
!====================================================================
!  (symmetric) Tri-diagonal matrix algorithm 
!  https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
!  a=c=off_diag, b=1-2*off_diag
!  d=Src
!--------------------------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer, intent(in) :: ndim
  complex(8), intent(in) :: off_diag
  complex(8), dimension(ndim), intent(in) :: Src
  complex(8), dimension(ndim), intent(out) :: solution
 
  integer :: n
  complex(8), dimension(ndim) :: cp, dp !see wikipedia for definition
 
  cp(1) = off_diag/(1d0-2d0*off_diag)
  dp(1) = Src(1)/(1d0-2d0*off_diag)
  do n=2,ndim
    cp(n) = 1d0/(1d0-2d0*off_diag-off_diag*cp(n-1))
    dp(n) = (Src(n)-off_diag*dp(n-1))*cp(n)
    cp(n) = off_diag*cp(n)
  enddo
  solution(ndim) = dp(ndim)
  do n=ndim-1,1,-1
    solution(n) = dp(n)-cp(n)*solution(n+1)
  enddo
end subroutine TriDiag 

subroutine aTriDiag(self,solution,Src,off_diag,ndim)
!====================================================================
!  (antisymmetric) Tri-diagonal matrix algorithm 
!  https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
!  a=-c=-off_diag, b=1
!  d=Src
!--------------------------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer, intent(in) :: ndim
  real*8, intent(in) :: off_diag
  complex(8), dimension(ndim), intent(in) :: Src
  complex(8), dimension(ndim), intent(out) :: solution
 
  integer :: n
  real*8, dimension(ndim) :: cp !see wikipedia for definition
  complex(8), dimension(ndim) :: dp !see wikipedia for definition
 
  cp(1) = off_diag
  dp(1) = Src(1)
  do n=2,ndim
    cp(n) = 1d0/(1d0+off_diag*cp(n-1))
    dp(n) = (Src(n)+off_diag*dp(n-1))*cp(n)
    cp(n) = off_diag*cp(n)
  enddo
  solution(ndim) = dp(ndim)
  do n=ndim-1,1,-1
    solution(n) = dp(n)-cp(n)*solution(n+1)
  enddo
end subroutine aTriDiag

subroutine linear_interpolation(self,arr1,n1,arr2,n2,narg)
!==================================================
!  linear interpolation from arr1 to arr2
!  arr1(:,1 )  : locations(from)
!  arr1(:,2:)  : values to be interpolated
!  arr2(:,1 )  : locations(to)
!  arr2(:,2:)  : values interpoalted
!--------------------------------------------------
  class(mathTemplate) :: self
  integer,intent(in)  :: n1,n2,narg
  real*8, intent(in)  :: arr1(n1,narg)
  real*8, intent(out) :: arr2(n2,narg)
  integer :: i1,i2,k
  real*8  :: R,L
  
  k=1
  do i2=1,n2
    do i1=k,n1
      if(arr1(i1,1)>arr2(i2,1)) then
        k = i1
        R = arr1(i1,1)-arr2(i2,1)
        if(i1>1) then
          L = arr2(i2,1)-arr1(i1-1,1)
          arr2(i2,2:) = (arr1(i1,2:)*L + arr1(i1-1,2:)*R)/(R+L)
        else
          arr2(i2,2:) = arr1(i1,2:)
        endif
        exit
      endif
      if(i1==n1) arr2(i2,2:) = arr1(i1,2:)
    enddo
  enddo
  
end subroutine

end module MathModule
