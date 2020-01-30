!-----------------------------------------------------------------------
!   SUBROUTINE: FINDLAI
!
!   Description: Find current LAI and previous LAI from LAIS46 and IDATE
!
!   Call: None
!
!   Require: None
!
!   Input:
!            1) IDATE : current julian date
!
!   Output:  1) LAIp_I
!            2) LAIc_I
!
!   Created by Tan 07/28/11
!   Updated by Ling Huang 02/18/17 for MEGAN3: LAI data is saved as
!   LAI1, LAI2, ... LAIS92, instead of one variable with multiple time
!   step.
!             
!-----------------------------------------------------------------------
      SUBROUTINE FINDLAI( IDATE, MXLAI, LAIp_I, LAIc_I)

      IMPLICIT NONE

! input
      INTEGER,INTENT(IN) :: IDATE  ! YYYYJJJ
      INTEGER,INTENT(IN) :: MXLAI
! output
      INTEGER,INTENT(OUT) :: LAIp_I, LAIc_I
! Local
      INTEGER :: JJJ
      REAL    :: XXX

! Calculation

      JJJ = MOD(IDATE,1000)
      XXX = JJJ/8.0
      LAIc_I = CEILING(XXX)

      IF (LAIc_I .EQ. 1) THEN
        LAIp_I = MXLAI
      ELSE
        LAIp_I = LAIc_I - 1
      ENDIF

      RETURN
      END SUBROUTINE FINDLAI
!-----------------------------------------------------------------------
