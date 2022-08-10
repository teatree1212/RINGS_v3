!======================================================================!
MODULE PARAMS
!----------------------------------------------------------------------!
USE DOUBLE
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
! Maximum allowed ZBRENT iterations.
!----------------------------------------------------------------------!
INTEGER :: ITMAX = 100
!----------------------------------------------------------------------!
! Machine floating-point precision.
! Not sure if this is be best value to use as now DP.
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: EPS = 3.0D-8
!----------------------------------------------------------------------!
! ZBRENT accuracy (used for equilibrium carbohydrate concentrations).
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: TOL = 0.001_DP
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: one     = 1.0_DP
REAL(DP), PARAMETER :: tf      = 273.15_DP ! Freezing point of water (K)
!----------------------------------------------------------------------!
! Boltzmann constant (eV K-1)
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: kb = 8.617D-5
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: mu_b    = 8.8_DP
REAL(DP), PARAMETER :: rd      = 3.14150_DP / 180.0_DP ! Deg -> radians
REAL(DP), PARAMETER :: pz_a    = -50.86372_DP
REAL(DP), PARAMETER :: pz_b    = 7.466117_DP
REAL(DP), PARAMETER :: ez_a    = -743.5649_DP
REAL(DP), PARAMETER :: ez_b    = 51.72611_DP
REAL(DP), PARAMETER :: tz_a    = -427.3073_DP
REAL(DP), PARAMETER :: tz_b    = 68.81009_DP
REAL(DP), PARAMETER :: sigma_a = 0.105_DP
REAL(DP), PARAMETER :: sigma   = 0.227_DP    ! Noise on size at division
REAL(DP), PARAMETER :: gasym   = 0.43_DP     ! Asymmetry of division
REAL(DP), PARAMETER :: fd      = 0.48_DP     ! Mode of division control
REAL(DP), PARAMETER :: tal      = 2680.0_DP  ! Cell axial length           (µm)
REAL(DP), PARAMETER :: ttl      = 50.1_DP    ! Cell tangential length      (µm)
REAL(DP), PARAMETER :: t_pw     = 0.8_DP     ! Primary cell wall thickness (µm)
REAL(DP), PARAMETER :: f_phloem = 0.315_DP
REAL(DP), PARAMETER :: Tfc      = 0.0_DP
!----------------------------------------------------------------------!
! Radial growth constant (µm µm-1 day-1)
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: r0 = 0.05705_DP
!----------------------------------------------------------------------!
! Wall growth constant   (mg ml-1 day-1)
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: m0 = 19.4_DP
!----------------------------------------------------------------------!
! Effective activation energy for radial growth (eV)
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: Ea = 0.3736_DP
!----------------------------------------------------------------------!
! Effective activation energy for wall building (eV)
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: Eaw = 1.425_DP
!----------------------------------------------------------------------!
! Primary cell-wall mass density (mg[DM] ml-1)
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: rhopw = 1.540D3
!----------------------------------------------------------------------!
! Secondary cell-wall mass density (mg[DM] ml-1).
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: rho_w = 1.43D3
!----------------------------------------------------------------------!
! Phloem carbohydrate concentration (mg ml-1).
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: Cbase = 95.05_DP
!----------------------------------------------------------------------!
! Effective Michaelis constant for wall growth (mg ml-1)
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: Km = 5.1_DP
!----------------------------------------------------------------------!
! Resistance to carbohydrate flux between cells (day ml-1)
!----------------------------------------------------------------------!
REAL(DP), PARAMETER :: res = 4359.0_DP
!----------------------------------------------------------------------!
END MODULE PARAMS
!======================================================================!
