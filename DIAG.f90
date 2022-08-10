!======================================================================!
SUBROUTINE DIAG (TRW, MXD)
!----------------------------------------------------------------------!
! End-of-year diagnostics.
!----------------------------------------------------------------------!
USE DOUBLE
USE CONTROL
USE PARAMS
USE STATE  ! MODULE declaring model variables.
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
REAL(DP) :: fMXD, den
! Max. ring of 10 cm
INTEGER, PARAMETER :: lmax = 10 * 10 * 1000
REAL(DP), DIMENSION (lmax) :: profile
REAL(DP), DIMENSION (2022-1901+1) :: TRW ! Annual growth (mm)
REAL(DP), DIMENSION (2022-1901+1) :: MXD ! Annual max. den (g[DM] cm-3)
REAL(DP) :: tM_mean, den_mean, N_mean
INTEGER :: lt, ia, ib, iyr
!----------------------------------------------------------------------!
iyr = kyr - syr + 1
TRW (iyr) = 0.0_DP
MXD (iyr) = 0.0_DP
DO fi = 1, nfi
 TRW (iyr) = TRW (iyr) + SUM (L(fi,1:ncells(fi)))
END DO
! Work up profile as sum across files of 1 um lengths.
profile = 0.0_DP
tM_mean = 0.0_DP
N_mean = 0.0_DP
DO fi = 1, nfi
 N_mean = N_mean + FLOAT (ncells (fi))
 DO ic = 1, ncells (fi)
  ia = NINT (D (fi,ic) - L (fi,ic) / 2.0_DP)
  ib = NINT (D (fi,ic) + L (fi,ic) / 2.0_DP) - 1
  den = 1.0D-3 * M (fi,ic) / (tal * ttl * L (fi,ic) / 1.0D12)
  profile (ia:ib) = profile (ia:ib) + den
  tM_mean = tM_mean + 1.0D3 * M (fi,ic)
 END DO ! ic
END DO ! fi
profile = profile / FLOAT (nfi)
tM_mean = tM_mean / FLOAT (nfi)
N_mean = N_mean / FLOAT (nfi)
MXD (iyr) = MAXVAL (profile(1:lmax))
TRW (iyr) = TRW (iyr) / FLOAT (nfi)
WRITE (32,*) kyr, TRW (iyr), MXD (iyr)
den_mean = 1.0D-6 * tM_mean / (tal * ttl * TRW (iyr) / 1.0D12)
IF (kyr == 1995) THEN
 WRITE (93,'(6f12.4)') dT, TRW (iyr), tM_mean, den_mean, MXD (iyr), N_mean
 WRITE (*,'(6f12.4)') dT, TRW (iyr), tM_mean, den_mean, MXD (iyr), N_mean
END IF
!----------------------------------------------------------------------!
END SUBROUTINE DIAG
!======================================================================!
