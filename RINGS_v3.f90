!======================================================================!
PROGRAM RINGS_v3
!----------------------------------------------------------------------!
! Update of RINGS for resubmission of Nature Communications paper.
!----------------------------------------------------------------------!
USE DOUBLE
USE PARAMS
USE ENVIRON
USE STATE
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
REAL(DP) :: dlat                            ! Latitude of site (degrees)
REAL(DP) :: iT, T_save
!****adf
!----------------------------------------------------------------------!
! Following to compute slope of relationship to calibrate Ea from
! Torn forcings and obs.
!----------------------------------------------------------------------!
REAL(DP), DIMENSION (2022-1901+1) :: jja ! Mean annual JJA temperature (degC)
REAL(DP), DIMENSION (2022-1901+1) :: TRW ! Annual growth (mm)
REAL(DP), DIMENSION (2022-1901+1) :: MXD ! Annual max. den (g[DM] cm-3)
real(dp) :: sum_x
real(dp) :: sum_y
real(dp) :: sum_xs
real(dp) :: sum_ys
real(dp) :: sum_xy
real(dp) :: a, b
real(dp) :: err, step, err_old
real(dp), dimension (ncells_max) :: Dmean, Smean, Mmean, WTHmean
real(dp) :: Vc, Vl, Vw, WTH
integer, dimension (ncells_max) :: nDmean
integer :: ii, jj, n, eyr_calib, i4
!****adf
!----------------------------------------------------------------------!
CHARACTER (LEN=200) :: Site ! Site                                (name)
!----------------------------------------------------------------------!
INTEGER :: iyr, ikt
!----------------------------------------------------------------------!
WRITE (*,*) "Running RINGS_v2..."
!----------------------------------------------------------------------!
OPEN (20,FILE='zones_fig.out',STATUS='UNKNOWN')
OPEN (21,FILE='cells_fig.out',STATUS='UNKNOWN')
OPEN (30,FILE='test.out',STATUS='UNKNOWN')
OPEN (31,FILE='zones.out',STATUS='UNKNOWN')
OPEN (32,FILE='TRW.dat',STATUS='UNKNOWN')
!----------------------------------------------------------------------!
OPEN (10,FILE='driver.txt',STATUS='OLD')
READ (10,*) Site
CLOSE (10)
!----------------------------------------------------------------------!
SELECT CASE (TRIM(Site))
 CASE ('Uggla')
  OPEN (11,FILE='uggla_clm_crujrav2.2.txt',STATUS='OLD')
  !OPEN (11,FILE='CLIM_Uggla_1951-1995.dat',STATUS='OLD')
  !read (11,*)
  !read (11,*)
  dlat = 64.0_DP + 21.0_DP / 60.0_DP
  syr = 1951 ! 1901 in file, so will need to skip if later
  eyr = 1995
  eyr_calib = 1995 !****adf
 CASE ('Tornetrask')
  OPEN (11,FILE='torn_clm_crujrav2.2.txt',STATUS='OLD')
  dlat = 68.26_DP
  syr = 1901
  eyr = 2020
  eyr_calib = 2004 !****adf
END SELECT
WRITE (*,*) 'Running for ', TRIM(Site), syr, eyr
!----------------------------------------------------------------------!
! Convert site latitude to radians                             (radians)
!----------------------------------------------------------------------!
rlat = rd * dlat
!----------------------------------------------------------------------!
! Compute daylength for each day (hr)
!----------------------------------------------------------------------!
CALL DAYL
!----------------------------------------------------------------------!
! For Figure 4.
!----------------------------------------------------------------------!
dT = -5.5_DP
!DO i4 = 1, 25
!dT = dT + 0.5_DP
!----------------------------------------------------------------------!
! Initialise radial files, dormancy state, and phenology signals.
!----------------------------------------------------------------------!
CALL INIT
!----------------------------------------------------------------------!
! Loop through timepoints and compute development of radial file.
!----------------------------------------------------------------------!
DO kyr = syr, eyr
 DO kday = 1, ndays
  !--------------------------------------------------------------------!
  T = 0.0_DP
  !--------------------------------------------------------------------!
  ! Loop over timepoints in day.
  !--------------------------------------------------------------------!
  DO kt = 1, 4
   READ (11,*) iyr, ikt, iT
   T = T + iT
  END DO ! kt = 1, 4
  !--------------------------------------------------------------------!
  T = T / 4.0_DP
  !--------------------------------------------------------------------!
  TC = T - tf
  !--------------------------------------------------------------------!
  ! Following adjusted from original.
  !--------------------------------------------------------------------!
  IF (dorm) THEN
    IF (((kday >= 306) .OR. (kday < 182)) .AND. (TC <= tfc)) cd = cd + 1
    IF ((kday > 32) .AND. (kday < 182)) &
      dd = dd + MAX (0.0_DP, (TC - Tfc))
    IF (dd >= (15.0_DP + 4401.8_DP * EXP (-0.042_DP * FLOAT (cd)))) THEN
     dorm = .FALSE.
     dd = 0.0_DP
     cd = 0
    END IF
  END IF
  IF (kday == 231) dorm = .TRUE.
  !--------------------------------------------------------------------!
  ! Compute zone widths. Fixing DOY at 171 for fort.97_20_fzw
  ! No autumn domancy for fort.97_20_fzw_nd
  !--------------------------------------------------------------------!
  !IF (kyr <= 1975) THEN
   CALL ZONES
  !ELSE
  ! pz =  104.18_DP
  ! ez =  400.73_DP
  ! tz = 1040.70_DP
  ! ! Following for _nd
  ! IF (kday > 182) dorm = .FALSE.
  !END IF
  !--------------------------------------------------------------------!
  ! Grow volume and mass of cells. T_save used for fort.97_20_nT
  !--------------------------------------------------------------------!
  !IF (kyr > 1975) T = 10.0_DP + tf
  !T_save = T
  !--------------------------------------------------------------------!
  ! For Figure 4.
  !--------------------------------------------------------------------!
  !T = T + dT
  !--------------------------------------------------------------------!
  CALL GROW
  !--------------------------------------------------------------------!
  !T = T_save
  !--------------------------------------------------------------------!
  ! Recalculate cell positions.
  !--------------------------------------------------------------------!
  CALL DIST
  !--------------------------------------------------------------------!
  ! Proliferate cells.
  !--------------------------------------------------------------------!
  IF (.NOT. dorm) CALL DIVIDE
  !--------------------------------------------------------------------!
  ! Recalculate cell positions.
  !--------------------------------------------------------------------!
  CALL DIST
  !--------------------------------------------------------------------!
  ! Daily cell diagnostics.
  !--------------------------------------------------------------------!
  CALL DIAG_CELLS
  !--------------------------------------------------------------------!
 END DO ! kday = 1, ndays
 !---------------------------------------------------------------------!
 ! Calculate and output end-of-year diagnostics.
 !---------------------------------------------------------------------!
 CALL DIAG (TRW, MXD)
 !---------------------------------------------------------------------!
 ! Before beginning of following year, remove all non-pz cells.
 !---------------------------------------------------------------------!
 DO fi = 1, nfi
  ic = 1
  DO WHILE (D (fi,ic) <= pz)
   ic = ic + 1
  END DO
  ncells (fi) = ic - 1
 END DO ! fi
 !---------------------------------------------------------------------!
END DO ! kyr = syr, eyr
!----------------------------------------------------------------------!
! For Figure 4.
!----------------------------------------------------------------------!
!REWIND (11)
!END DO ! i4
!----------------------------------------------------------------------!
CLOSE (11) ! Climate input file.
CLOSE (20) ! zones_fig.out
CLOSE (21) ! cells_fig.out
CLOSE (30) ! Testing output file.
CLOSe (31) ! Zone widths.
CLOSE (32) ! Mean annual final radial file length (µm).
!----------------------------------------------------------------------!
WRITE (*,*) "RINGS_v2 finished cleanly..."
!----------------------------------------------------------------------!
END PROGRAM RINGS_v3
!======================================================================!


