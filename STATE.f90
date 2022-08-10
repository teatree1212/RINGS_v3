!======================================================================!
MODULE STATE
!----------------------------------------------------------------------!
USE DOUBLE
USE CONTROL
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
REAL(DP), DIMENSION (nfi, ncells_max) :: L     ! Radial cell length (µm)
REAL(DP), DIMENSION (nfi, ncells_max) :: M     ! Cell DM        (mg[DM])
REAL(DP), DIMENSION (nfi, ncells_max) :: rasym ! Asymmetry on enlarg.
REAL(DP), DIMENSION (nfi, ncells_max) :: I     ! Div. inhib. vol  (µm^3)
REAL(DP), DIMENSION (nfi, ncells_max) :: D ! Cell distance from phloem (µm)
!----------------------------------------------------------------------!
! Carbohydrate concentration in cell lumen (mg ml-1).
!----------------------------------------------------------------------!
REAL(DP), DIMENSION (nfi, ncells_max) :: SUC
!----------------------------------------------------------------------!
REAL(DP) :: rlat              ! Latitude of site               (radians)
REAL(DP) :: pz_min            ! Minimum width of proliferation zone (µm)
REAL(DP) :: pz                ! Width of proliferation zone         (µm)
REAL(DP) :: ez                ! Width of enlargement zone           (µm)
REAL(DP) :: tz                ! Width of wall thickening zone       (µm)
REAL(DP) :: T                 ! Mean daily temperature               (K)
REAL(DP) :: TC                ! Mean daily temperature              (oC)
REAL(DP) :: dd                ! Degree-days                         (oC)
!----------------------------------------------------------------------!
INTEGER, DIMENSION (nfi) :: ncells       ! No. cells/file            (n)
INTEGER :: ncells_alive                  ! No. cells in cambial zone (n)
!----------------------------------------------------------------------!
INTEGER :: kyr                                  ! Year              (CE)
INTEGER :: kday                                 ! Day of year      (day)
INTEGER :: cd                                   ! Chilling days   (days)
!----------------------------------------------------------------------!
LOGICAL :: dorm                                 ! Domancy state      (-)
!----------------------------------------------------------------------!
END MODULE STATE
!======================================================================!
