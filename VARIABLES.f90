!======================================================================!
MODULE VARIABLES
!----------------------------------------------------------------------!
USE DOUBLE
USE CONTROL
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
! Wall mass growth with unlimiting carbohydrates (mg cell-1 day-1).
!----------------------------------------------------------------------!
REAL(DP), DIMENSION (nfi, ncells_max) :: dM_max
!----------------------------------------------------------------------!
! Wall mass growth (mg cell-1 day-1).
!----------------------------------------------------------------------!
REAL(DP), DIMENSION (nfi, ncells_max) :: dM
!----------------------------------------------------------------------!
REAL(DP) :: Sp             ! Phloem carbohydrate concentration (mg ml-1)
REAL(DP) :: SG ! Cell lumen carbohydrate concentration to test (mg ml-1)
REAL(DP) :: zbrent_out   ! Carbohydrate concentration solution (mg ml-1)
REAL(DP) :: FS                    ! Error on solution using SG (mg ml-1)
!----------------------------------------------------------------------!
END MODULE VARIABLES
!======================================================================!
