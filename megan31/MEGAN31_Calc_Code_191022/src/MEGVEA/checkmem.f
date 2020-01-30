
        SUBROUTINE CHECKMEM( MSTATUS, ONVAR, CALLER )
 
!***********************************************************************
!  subroutine body starts at line  105
!
!  DESCRIPTION:
!       Reports an error and exits if memory status flag is non-zero.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Adapted 10/98 by M Houyoux
!
!***********************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id: checkmem.f,v 1.2 1999/11/23 12:46:12 mhouyoux Exp $
!
! COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
! All Rights Reserved
!
! See file COPYRIGHT for conditions of use.
!
! Environmental Programs Group
! MCNC--North Carolina Supercomputing Center
! P.O. Box 12889
! Research Triangle Park, NC  27709-2889
!
! env_progs@mcnc.org
!
! Pathname: $Source: /afs/isis/depts/cep/emc/apps/archive/edss_tools/edss_tools/src/lib/checkmem.f,v $
! Last updated: $Date: 1999/11/23 12:46:12 $ 
!
!***************************************************************************
 
      IMPLICIT NONE

!...........   ARGUMENTS and their descriptions:

       INTEGER       MSTATUS !  ALLOCATE function exit status
       CHARACTER*(*) ONVAR   !  Variable name of previous ALLOCATE statement
       CHARACTER*(*) CALLER  !  Name of calling program

!...........   ARGUMENTS and their descriptions:
       INTEGER      TRIMLEN
       EXTERNAL     TRIMLEN

!...........   Local variables

       INTEGER         L1
       INTEGER         L2
       CHARACTER*256   MESG

       CHARACTER*16 :: PROGNAME = 'CHECKMEM' ! program name

!***********************************************************************
!   begin body of function CHECKMEM

!.........  Get lengths of input character strings
        L1 = TRIMLEN( ONVAR )
        L2 = TRIMLEN( CALLER )

!.........  Abort if memory status is non-zero

        IF( MSTATUS .GT. 0 ) THEN            
            MESG = 'Failure allocating memory for "' // ONVAR( 1:L1 ) // 
     &             '" variable'
            CALL M3EXIT( CALLER( 1:L2 ), 0, 0, MESG, 2 )
        ENDIF

        RETURN

        END

