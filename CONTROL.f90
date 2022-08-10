!======================================================================!
MODULE CONTROL
!----------------------------------------------------------------------!
USE DOUBLE
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER :: fi
INTEGER :: ic
INTEGER :: kt
!----------------------------------------------------------------------!
REAL(DP) :: dT ! For temperature experiments (K)
!****adf was 100
INTEGER, PARAMETER :: nfi        =  100 ! No. radial files           (n)
INTEGER, PARAMETER :: ncells_max = 1000 ! Max. no. cells/file        (n)
INTEGER, PARAMETER :: ndays      =  365 ! No. days in year        (days)
!----------------------------------------------------------------------!
INTEGER :: syr                          ! Start year of simulation  (CE)
INTEGER :: eyr                          ! End year of simulation    (CE)
!----------------------------------------------------------------------!
END MODULE CONTROL
!======================================================================!
