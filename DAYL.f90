!======================================================================!
SUBROUTINE DAYL
!----------------------------------------------------------------------!
USE DOUBLE
USE CONTROL
USE PARAMS
USE ENVIRON
USE STATE
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
REAL(DP) :: delta ! Solar declination                          (radians)
REAL(DP) :: x     ! Intermediate variable         (cosine of hour angle)
REAL(DP) :: ha    ! Hour angle of sunrise (sunset) from noon   (radians)
!----------------------------------------------------------------------!
DO kday = 1, ndays
 !---------------------------------------------------------------------!
 ! http://fer3.com/arc/img/111019.life-raft-declination.pdf
 !---------------------------------------------------------------------!
 delta = rd * 23.45_DP * SIN ((rd * 360.0_DP / 365.25_DP) * &
         FLOAT (kday - 80))
 !---------------------------------------------------------------------!
 ! https://physics.stackexchange.com/questions/28563/
 ! hours-of-light-per-day-based-on-latitude-longitude-formula/320092
 !---------------------------------------------------------------------!
 x = COS (rd * 90.833_DP) / (COS (rlat) * COS (delta)) - &
     TAN (rlat) * TAN (delta)
 x = MAX (-1.0_DP, x)
 x = MIN ( 1.0_DP, x)
 !---------------------------------------------------------------------!
 ha = ACOS (x) ! Hour angle of sunrise (sunset) from noon      (radians)
 !---------------------------------------------------------------------!
 dlength (kday) = 2.0_DP * (ha / rd) / 15.0_DP ! Daylength       (hours)
 !---------------------------------------------------------------------!
END DO ! kday = 1, ndays
!----------------------------------------------------------------------!
END SUBROUTINE DAYL
!======================================================================!
