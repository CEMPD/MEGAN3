
      SUBROUTINE CNVT_EFS(GDNAM,EFSTXTF,EFMAPS)
!***********************************************************************
!  This subroutine converts EF txt to EF IOAPI
!
!  Called by TXT2IOAPI
!
!  Call:
!     PARSLINE
!     CHECKMEM
!
!  Created by Tan  07/07/2011 for MEGAN v2.10
!
!  History:
!***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
      INCLUDE 'PARMS3.EXT'    !  I/O API parameters
      INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

      INCLUDE 'SPC_MGN.EXT'   !  SPC

C...........   EXTERNAL FUNCTIONS and their descriptions:
      INTEGER      PROMPTFFILE, STR2INT
      REAL         STR2REAL
      LOGICAL      DSCGRID
      INTEGER      INDEX1
      EXTERNAL     PROMPTFFILE, DSCGRID
     &             STR2INT, STR2REAL, INDEX1

C.........  LOCAL VARIABLES and their descriptions:
      INTEGER       IOS            ! i/o status
      INTEGER       ILINE          ! current line
      INTEGER       ADEV           ! unit no. for input data file
      INTEGER       I, J, S, N     ! Counters
      INTEGER       ITIME          ! Counters
      INTEGER       IDATE          ! Counters
      CHARACTER*16    CNAME        ! Coord name
      CHARACTER(LEN=1000) LINE     ! input ine buffer
      CHARACTER*256    MESG        ! message buffer
      INTEGER       CID, INX, INY  ! Input grid x and y

      CHARACTER*16  :: EFSTXTF
      CHARACTER*16  :: EFMAPS
      CHARACTER*16  :: GDNAM

      ! CID, 20 EF, 1 invalid
      INTEGER, PARAMETER :: MXTCOL = 22      ! Columns in an input line
      CHARACTER*30     SEGMENT( MXTCOL )     ! Input line fields

      INTEGER, PARAMETER :: NFLD  = MXTCOL-2 ! Number of output variables
      INTEGER,         XREF(NFLD)
      INTEGER       :: NSTEP = NFLD
      INTEGER, PARAMETER :: NVARS = NFLD         ! Number of output variables
      INTEGER, PARAMETER :: NEF = NFLD      ! Number of EF

      REAL, ALLOCATABLE :: IBUFF( :, :, : )  ! input, output buffer

      CHARACTER*16  :: PROGNAME = 'CNVT_EFS'

!***********************************************************************
!... Begin program

!... Get output file grid parameters
      GDNAM3D = GDNAM
      IF( .NOT. DSCGRID( GDNAM3D, CNAME, GDTYP3D,
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!... Allocate memory
      ALLOCATE ( IBUFF( NCOLS3D, NROWS3D, NSTEP ),
     &           STAT = IOS )
      CALL CHECKMEM ( IOS, 'IBUFF', PROGNAME )
      IBUFF = 0.

!... Get input file unit
      MESG = 'Enter logical name of EF data file.'
      ADEV = PROMPTFFILE( MESG, .TRUE., .TRUE., EFSTXTF, PROGNAME )
      IF( ADEV < 0 ) THEN
         MESG = 'Error opening input EF file.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!... Read TXT input file
      ILINE = 0
      DO
         READ( ADEV, 1000, IOSTAT = IOS ) LINE
         ILINE = ILINE + 1
         IF( IOS < 0 ) EXIT       ! End of file, exit loop
         IF( LINE == ' ' ) CYCLE  ! Skip blank line
         IF( IOS > 0 ) THEN
            WRITE( MESG, 1010 ) 'I/O error', IOS,
     &            'reading input file at line', ILINE
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
         ENDIF

         CALL PARSLINE( LINE, MXTCOL, SEGMENT )  ! break lines into segments
         IF ( ILINE .EQ. 1 ) THEN

            DO S = 1 , NEF
               XREF( S ) = S
               WRITE(MESG,'(I0.2)') S
               PRINT*,'READ EF'//TRIM(MESG)
            ENDDO

            CYCLE  ! go back to read 2nd line
         ENDIF

         CID = STR2INT( SEGMENT(1) )   ! convert character to integer
         INY = (CID - 1)/NCOLS3D + 1
         INX = CID - (INY-1)*NCOLS3D
         DO S = 1, NFLD
            IBUFF(INX,INY,XREF(S)) = STR2REAL( SEGMENT( 1+S ) )  ! convert char to real
         ENDDO
      ENDDO
      WHERE (IBUFF(:,:,:) .LT.0 ) 
           IBUFF = 0.0
      ENDWHERE

      WRITE( MESG, 1010 ) 'Total number of input lines read:',ILINE
      CALL M3MESG( MESG )

!... Setup additional outfile parameters and open outfile
      FTYPE3D = GRDDED3
      NVARS3D = NVARS
!... Set VNAME3D and UNITS3D
      DO N = 1,NEF
        VNAME3D( N ) = 'EF_'//TRIM(MGN_SPC(N))
        UNITS3D( N ) = 'nanomol/m^2/s'
      ENDDO
      VTYPE3D(1:NVARS) = M3REAL
      TSTEP3D = 10000
      MXREC3D = 1
      NLAYS3D = 1
      CALL NAMEVAL ( EFSTXTF, MESG )  ! get input file name
      FDESC3D( 1 ) = 'Converted from MEGAN input text file: '
      FDESC3D( 2 ) = TRIM(MESG)
      IF ( .NOT. OPEN3( EFMAPS, FSCREA3, PROGNAME ) ) THEN
           MESG = 'Could not open file "'//TRIM(EFMAPS)//'" for output'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!... Write output
      IDATE = 0
      ITIME = 0
      DO S = 1, NEF
         IF ( .NOT. WRITE3(EFMAPS,VNAME3D(S),IDATE,ITIME,
     &                                              IBUFF(:,:,S))) THEN
            MESG = 'Error writing to  file "'//TRIM(EFMAPS)//'"'
            CALL M3EXIT(PROGNAME,0,0,MESG,2)
         ENDIF
      ENDDO

      IF ( CLOSE3(EFMAPS) ) THEN
        CALL M3MESG( 'Closing EF IOAPI file')
      ENDIF
      CLOSE(ADEV)

1000  FORMAT( 40A )
1010  FORMAT( 40( A, :, I8, :, 1X ) )

      END SUBROUTINE
