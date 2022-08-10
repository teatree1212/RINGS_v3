!======================================================================!
SUBROUTINE CARB
!----------------------------------------------------------------------!
! Solve for equilibrium distributions of carbohydrates (SUC) across the
! radial files and associated cell wall mass growth rates (dM).
!----------------------------------------------------------------------!
USE DOUBLE
USE CONTROL
USE PARAMS
USE STATE
USE VARIABLES
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
! Initialise wall mass growth (mg cell-1 day-1).
!----------------------------------------------------------------------!
dM (:,:) = 0.0_DP
!----------------------------------------------------------------------!
! Loop over radial files.
!----------------------------------------------------------------------!
DO fi = 1, nfi
 !---------------------------------------------------------------------!
 ! Initialise number of living cells.
 !---------------------------------------------------------------------!
 ncells_alive = 0
 !---------------------------------------------------------------------!
 ! Count number of living cells.
 !---------------------------------------------------------------------!
 DO ic = 1, ncells (fi)
  IF (D (fi,ic) <= tz) ncells_alive = ncells_alive + 1
 END DO ! ic
 !---------------------------------------------------------------------!
 ! Set phloem carbohydrate concentration (mg ml-1).
 !---------------------------------------------------------------------!
 Sp = Cbase
 ! For Figure 2 and 4.
 !Sp = 1.0D6 * Sp
 !---------------------------------------------------------------------!
 ! Solve equilibrium carbohydrate concentration in last living cell.
 !---------------------------------------------------------------------!
 CALL ZBRENT
 !---------------------------------------------------------------------!
 ! Carbohydrate concentration in last living cell (mg ml-1).
 !---------------------------------------------------------------------!
 SG = zbrent_out
 !---------------------------------------------------------------------!
 ! Make sure dM are the final equilibrium values at SG.
 !---------------------------------------------------------------------!
 CALL FUNC
 !---------------------------------------------------------------------!
END DO ! fi = 1, nfi
!----------------------------------------------------------------------!
END SUBROUTINE CARB
!======================================================================!

!======================================================================!
SUBROUTINE ZBRENT
!----------------------------------------------------------------------!
USE DOUBLE
USE PARAMS
USE STATE
USE VARIABLES
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER :: iter
REAL(DP) :: x1,x2,a,b,fa,fb,c,fc,dz,e,tol1,xm,s,q,p,rz
!----------------------------------------------------------------------!
x1 = 0.0_DP
!x1 = EPS !****
x2 = Sp
a = x1
b = x2
SG = a
CALL FUNC
fa = FS
SG = b
CALL FUNC
fb = FS
if(fb * fa > 0.0_DP)then
 write(*,*) 'Root must be bracketed for ZBRENT'
 write(*,*)'a b fa fb',a,b,fa,fb
 write(*,*)'fi ncells ncells_alive',fi,ncells(fi),ncells_alive
 write(*,*)'SUC     ',SUC(fi,1:ncells_alive)
 write(*,*)'pz ez tz',pz,ez,tz
 write(*,*)'L       ',L(fi,1:ncells(fi))
 write(*,*)'D       ',D(fi,1:ncells(fi))
 write(*,*)'dM      ',dM    (fi,1:ncells_alive)
 write(*,*)'dM_max  ',dM_max(fi,1:ncells_alive)
 stop
endif
c = b !****adf seems to be have been forgotten in _NC! (14/12/21)
fc = fb
do iter=1, ITMAX
 if(fb*fc>0.0_DP)then
  c=a
  fc=fa
  dz=b-a
  e=dz
 endif
 if(abs(fc)<abs(fb))then
  a=b
  b=c
  c=a
  fa=fb
  fb=fc
  fc=fa
 endif
 tol1 = 2.0_DP * EPS * ABS (b) + 0.5_DP * TOL
 xm=0.5_DP*(c-b)
 if((abs(xm)<=tol1).or.(fb==0.0_DP))then
  zbrent_out=b
  goto 8000
 endif
 if((abs(e)>=tol1).and.(abs(fa)>abs(fb)))then
  s=fb/fa
  if(a==c)then
   p=2.0_DP*xm*s
   q=1.0_DP-s
  else
   q=fa/fb
   rz=fb/fc
   p=s*(2.0_DP*xm*q*(q-rz)-(b-a)*(rz-1.0_DP))
   q=(q-1.0_DP)*(rz-1.0_DP)*(s-1.0_DP)
  endif
  if(p>0.0_DP)q=-q
  p=abs(p)
  if((2.0_DP*p)<(min(3.0_DP*xm*q-abs(tol1*q),abs(e*q))))then
   e=dz
   dz=p/q
  else
   dz=xm
   e=dz
  endif
 else
  dz=xm
  e=dz
 endif
 a=b
 fa=fb
 if(abs(dz)>tol1)then
  b=b+dz
 else
  b=b+sign(tol1,xm)
 endif
 SG=b
 call func
 fb=FS
enddo !iter
write(*,*)'ZBRENT exceeded maximum iterations'
stop
8000 continue
zbrent_out=b
!----------------------------------------------------------------------!
END SUBROUTINE ZBRENT
!======================================================================!

!======================================================================!
SUBROUTINE FUNC
!----------------------------------------------------------------------!
! Determine FS, the difference between the prescribed sugar
! concentration in the last cell (S1=SG) and that inferred from the wall
! growth imposed on the concentration gradient from the phloem to the
! final living cell (S2).
!----------------------------------------------------------------------!
USE DOUBLE
USE CONTROL
USE PARAMS
USE STATE
USE VARIABLES
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER :: jc
!----------------------------------------------------------------------!
REAL(DP) :: S1
REAL(DP) :: S2
REAL(DP) :: temp
!----------------------------------------------------------------------!
! Cell number of last living cell.
!----------------------------------------------------------------------!
ic = ncells_alive
!----------------------------------------------------------------------!
! Prescribed carbohydrate concentration in last cell to test (mg ml-1).
!----------------------------------------------------------------------!
S1 = SG
SUC (fi,ic) = SG
!----------------------------------------------------------------------!
! Growth of cell walls based on sugar content of cell (mg cell-1 day-1).
!----------------------------------------------------------------------!
dM (fi,ic) = dM_max (fi,ic) * S1 / (S1 + Km)
!----------------------------------------------------------------------!
temp = dM (fi,ic)
!----------------------------------------------------------------------!
! Loop through cell, last in file first.
!----------------------------------------------------------------------!
DO ic = ncells_alive-1, 1, -1
 !---------------------------------------------------------------------!
 ! SUC for this cell based on required flux to next cell to balance
 ! growth given SUC in last cell = SG (mg ml-1).
 !---------------------------------------------------------------------!
 SUC (fi,ic) = SUC (fi,ic+1) + res * temp
 !---------------------------------------------------------------------!
 ! Wall mass growth based on SUC (mg cell-1 day-1).
 !---------------------------------------------------------------------!
 dM (fi,ic) = dM_max (fi,ic) * SUC (fi,ic) / (SUC (fi,ic) + Km)
 !---------------------------------------------------------------------!
 ! Accumulate wall growth along file (mg day-1).
 !---------------------------------------------------------------------!
 temp = temp + dM (fi,ic)
 !---------------------------------------------------------------------!
END DO ! ic
!----------------------------------------------------------------------!
! Inferred sugar concentration in last cell from growth and phloem
! concentration (mg ml-1).
!----------------------------------------------------------------------!
S2 = .0
DO ic = 1, ncells_alive
 !****adf looks like following is wrong, so changed!
 !DO jc = 1, ncells_alive
 DO jc = ic, ncells_alive
 !****adf
  S2 = S2 + res * dM (fi,jc)
 END DO ! jc
END DO ! ic
S2 = Sp - S2
!----------------------------------------------------------------------!
! Error in last cell's concentration between supply and demand
! calculations (mg ml-1).
!----------------------------------------------------------------------!
FS = S1 - S2
!----------------------------------------------------------------------!
END SUBROUTINE FUNC
!======================================================================!
