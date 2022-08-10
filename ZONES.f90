!======================================================================!
SUBROUTINE ZONES
!----------------------------------------------------------------------!
! Compute zone edge distances from phloem based on year-day and
! daylength (µm).
!----------------------------------------------------------------------!
USE DOUBLE
USE PARAMS ! MODULE declaring model parameters.
USE STATE  ! MODULE declaring model variables.
USE ENVIRON
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
REAL(DP) :: a_ez
REAL(DP) :: b_ez
!----------------------------------------------------------------------!
IF (kday < 185) THEN
 IF (dorm) THEN
  pz = pz_a + pz_b * dlength (kday)
  ez = pz
 ELSE
  pz = pz_a + pz_b * dlength (185)
  !-----------------------------------------------------------------!
  ! Slope.
  !-----------------------------------------------------------------!
  b_ez = ((ez_a + ez_b * dlength(210)) - &
          (ez_a + ez_b * dlength(185))) / &
          (210.0_DP - 185.0_DP + 1.0_DP)
  !-----------------------------------------------------------------!
  ! Intercept.
  !-----------------------------------------------------------------!
  a_ez = (ez_a + ez_b * dlength(185)) - b_ez * FLOAT(185)
  ez = a_ez + b_ez * FLOAT(kday)
  !-----------------------------------------------------------------!
 END IF
ELSE
 pz = pz_a + pz_b * dlength (kday)
 ez = ez_a + ez_b * dlength (kday)
END IF
pz = MAX (pz_min, pz)
ez = MAX (pz    , ez)
tz = MAX (ez    , tz_a + tz_b * dlength (kday))
!----------------------------------------------------------------------!
END SUBROUTINE ZONES
!======================================================================!
