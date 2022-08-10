!======================================================================!
SUBROUTINE DIAG_CELLS
!----------------------------------------------------------------------!
! Diagnostics for mean cells on particular days.
!----------------------------------------------------------------------!
USE DOUBLE
USE CONTROL
USE PARAMS
USE ENVIRON
USE STATE  ! MODULE declaring model variables.
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
! Max. ring of 10 cm
INTEGER, PARAMETER :: lmax = 10 * 10 * 1000
REAL(DP) :: rncells_mean
REAL(DP), DIMENSION (ncells_max) :: D_mean, L_mean, M_mean
REAL(DP) :: z_phloem, f_Vl, Vc, Vl, len
INTEGER :: nincm, ncells_mean
!----------------------------------------------------------------------!
! For Figure 1, density profiles.
!----------------------------------------------------------------------!
IF ((kyr > 1975) .AND. (kday == 365)) THEN
 DO fi = 1, nfi
  len = D (fi,ncells(fi)) + L (fi,ncells(fi)) / 2.0_DP - pz
  DO ic = 1, ncells (fi)
   ! Subtract pz to get final ring.
   IF (D(fi,ic) > pz) THEN
    WRITE (97,*) (D(fi,ic) - pz) / len, L(fi,ic), 1.0D3*M(fi,ic)
   END IF
   !-------------------------------------------------------------------!
   ! For Figure 2 density profiles.
   !-------------------------------------------------------------------!
   IF (kyr == 1995) THEN
    WRITE (94,*) (D(fi,ic) - pz), L(fi,ic), 1.0D3*M(fi,ic)
   END IF
   !-------------------------------------------------------------------!
  END DO
 END DO
END IF
!----------------------------------------------------------------------!
IF (kyr > syr) WRITE (30,*) kyr, kday, cd, dd
IF (kyr == 1995) WRITE (31,*) kday, pz, ez, tz
IF (kyr == 1995) THEN
 ! For zones figure.
 IF (dorm) THEN
  WRITE (20,'(2I5,6F12.4)') kyr,kday,dlength(kday),pz,ez,tz,1.0_DP,TC
   ELSE
  WRITE (20,'(2I5,6F12.4)') kyr,kday,dlength(kday),pz,ez,tz,0.0_DP,TC
 END IF
END IF
rncells_mean = 0.0_DP
DO fi = 1, nfi
 rncells_mean = rncells_mean + FLOAT (ncells (fi))
END DO
ncells_mean = NINT (rncells_mean / FLOAT (nfi))
L_mean = 0.0_DP
M_mean = 0.0_DP
nincm = 0
DO fi = 1, nfi
 IF (ncells (fi) == ncells_mean) THEN
  L_mean (:) = L_mean (:) + L (fi,:)
  M_mean (:) = M_mean (:) + M (fi,:)
  nincm = nincm + 1
 END IF
END DO
L_mean = L_mean / FLOAT (nincm)
M_mean = M_mean / FLOAT (nincm)
D_mean = 0.0_DP
z_phloem = f_phloem * pz
D_mean (1) = z_phloem + L_mean (1) / 2.0_DP
DO ic = 2, ncells_mean
 D_mean (ic) = D_mean (ic-1) + L_mean (ic-1) / 2.0_DP + L_mean (ic) / 2.0_DP
END DO ! ic
IF ((kday == 185) .OR. (kday == 210) .OR. (kday == 231) .OR. &
 (kday == 365)) THEN
 DO ic = 1, ncells_mean
  Vc = tal * ttl * L_mean (ic) / 1.D12 ! Cell volume (ml).
  ! Lumen volume (ml).
  IF (D_mean (ic) <= pz) THEN
   Vl = Vc - M_mean (ic) / rhopw
  ELSE
   Vl = Vc - M_mean (ic) / rho_w
  END IF
  IF (Vc > 0.0_DP)THEN
   f_Vl = Vl / Vc
  ELSE
   f_Vl = 0.0_DP
  END IF
  ! For zones figure.
  WRITE (21,'(2I5,4F12.4)') kday, ic, D_mean(ic), L_mean(ic), &
                            1.0D3*M_mean(ic), f_Vl
 END DO
END IF
!----------------------------------------------------------------------!
END SUBROUTINE DIAG_CELLS
!======================================================================!
