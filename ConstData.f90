module ConstDataModule
  implicit none
!    real*8, parameter :: cLight= 299792458.0d0    !speed of light in vaccum [m/s]
!    real*8, parameter :: VacImp= 376.73031346177d0  !vaccum impedance [V/A]
!    real*8, parameter :: eMass = 510998.9461d0      !electron mass [eV]
!    real*8, parameter :: ImpEff = VacImp/eMass*Clight
!    real*8, parameter :: e1    = 1.60217662d-19   !charge of 1 electron [Coulomb]
  
  type, private :: PhysConstTemplate
    real*8  :: c   = 299792458.0d0    !speed of light in vaccum [m/s]
    real*8  :: Z   = 376.73031346177d0!vaccum impedance [V/A]
    real*8  :: Me  = 510998.9461d0    !electron mass [eV]
    real*8  :: Zeff= 221019.842678d0  !effective Impedance  = VacImp/eMass*Clight
    real*8  :: e   = 1.60217662d-19   !charge of 1 electron [Coulomb]
  end type
  type(PhysConstTemplate), public :: PhysConst
  
  complex(8), parameter :: i1=(0d0, 1d0)
  integer, parameter :: x_=1,px_=2,y_=3,py_=4,t_=5,g_=6,q_=7  ! particle phase-space index
  integer, parameter :: uniform_=1,gaussian_=2                ! distribution index
  integer, parameter :: Undulator_=1, Drift_=2, Quad_=3       ! Element type
  real*8,  parameter :: pi    = 3.141592653589793d0
  real*8,  parameter :: twopi = 6.283185307179586d0

end module ConstDataModule
