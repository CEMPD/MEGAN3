	PROGRAM DAYMET

!   New for MEGANv3: generate daily average temperature/solar
!                    radiation, maximum/minimum temperature/
!                    wind speed and write to an output file
!   Created by Ling Huang 
!*****************************************************************
!   Input varibles
!
!   Day                  Julian day
!   TEMP                 Temperature [K]
!   PPFD           Incoming photosynthetic active radiation [umol/m2.s]
!   Wind                 Wind speed [m/s]
!
!*****************************************************************
! OUTPUT
!
! For each day and each location
!   MaxT                 Daily maximum temperature [K]
!   MinT                 Daily minimum temperature [K]
!   MaxWS                Daily maximum wind speed [m/s]
!   D_TEMP               Daily average temperature [K]
!   D_PPFD               Daily average solar radiation [umol/m2.s]

!*****************************************************************

      IMPLICIT NONE

!...  INCLUDES:
      INCLUDE 'PARMS3.EXT'          !  I/O API parameters
      INCLUDE 'IODECL3.EXT'         !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'          !  I/O API file description data structures

!...  EXTERNAL FUNCTIONS and their descriptions:
      INTEGER, EXTERNAL       ::   ENVINT
      LOGICAL      DSCGRID
      EXTERNAL     DSCGRID

! input

!...  Program I/O files: From run script
! Program name
      CHARACTER*16  :: PROGNAME = 'DAYMET'
! Met files
      CHARACTER*16 :: MGNMET

! Output file (daily meteorology)
      CHARACTER*16  :: OUTMET = 'DailyMET'     ! Output file logical name

!...  Parameters for file units
      INTEGER  LOGDEV                      ! Logfile unit number

! ... External parameters
! From run script
      INTEGER       SDATE          ! Start date YYYYDDD
      INTEGER       STIME          ! Start time HHMMSS
      INTEGER       RLENG          ! Run length HHMMSS
      INTEGER       NDAYS          ! Episode length, i.e. number of days

! I/O API file parameters
      INTEGER :: JDATE        ! Date YYYYDDD from inpname
      INTEGER :: JTIME        ! Time HHMMSS from inpname
      INTEGER :: NCOLS        ! Number of columns
      INTEGER :: NROWS        ! Number of rows
      INTEGER :: MXREC        ! Total number of timesteps
      INTEGER :: TSTEP        ! Time step

!... Internal parameters
! Internal parameters (status and buffer)
      INTEGER       IOS                    ! i/o status
      CHARACTER*256 MESG                   ! message buffer

! Parameters for output species
      INTEGER, PARAMETER :: NOUTPUT = 5
                          ! number of output variables: D_TEMP, D_PPFD,
                          ! MaxT, MinT, MaxWS

      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

! variables from input file

      REAL, ALLOCATABLE :: TEMP( :,:,: )    ! Temperature (K)
      REAL, ALLOCATABLE :: PPFD( :,:,: )    ! Calculated PAR (umol/m2.s)
      REAL, ALLOCATABLE :: WIND( :,:,: )    ! Wind speed (m/s)

      ! loop indices
      INTEGER :: IDATE             ! Looping
      INTEGER :: ITIME             ! Looping
      INTEGER :: N, T, I

! output
      REAL, ALLOCATABLE :: D_TEMP  ( :,: )   ! Daily average temperature [K]
      REAL, ALLOCATABLE :: D_PPFD  ( :,: )   ! Daily average solar radiation [umol/m2.s]
      REAL, ALLOCATABLE :: MaxT    ( :,: )   ! Daily maximum temperature [K]
      REAL, ALLOCATABLE :: MinT    ( :,: )   ! Daily minimum temperature [K]
      REAL, ALLOCATABLE :: MaxWS   ( :,: )   ! Daily maximum wind speed [m/s]

! local variables and their descriptions:

      INTEGER       AVEBY           ! Divider for daily average

!**************************************************************************

!--========================================================================
!... Begin program
!--========================================================================

!--------------------------------------------------------------------------
!.....1) File set up and assign I/O parameters
!--------------------------------------------------------------------------
!...  Initialize log file unit
      LOGDEV = INIT3()
           !  Now I/O API is set up, and LOGUNIT is the unit number
           !  for the log file (or it 6 for st'd output).

!...  Get input parameters from run script
      CALL ENVSTR( 'GDNAM3D', MESG, 'ASACA36km', GDNAM, IOS )
      IF( .NOT. DSCGRID( GDNAM, CNAME, GDTYP3D,
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      MESG = 'Model start date (YYYYDDD)'
      SDATE = ENVINT( 'EPISDATE', MESG, JDATE, IOS )

      MESG = 'Model start time (HHMMSS)'
      STIME = ENVINT( 'STIME', MESG, JTIME, IOS )

      MESG = 'Model run length (HHMMSS)'
      RLENG = ENVINT( 'RLENG', MESG, MXREC*10000, IOS )

      MESG = 'Number of days'
      NDAYS = ENVINT( 'NDAYS', MESG, 1, IOS )

      MGNMET = 'MGNMET001'
      IF ( .NOT. OPEN3( MGNMET, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      ! Check grid
      IF ( .NOT. FILCHK3 ( MGNMET,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'MGNMET has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      IF ( .NOT. DESC3( MGNMET ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Check start date, start time, end date, end time in MGNMET
      IDATE = SDATE; ITIME = STIME
      IF ( .NOT. CHECK3( MGNMET, 'TEMP2', IDATE, ITIME ) ) THEN
         MESG = 'Starting time not on met file'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      CALL NEXTIME ( IDATE, ITIME, RLENG-10000 )
      IF ( .NOT. CHECK3( MGNMET, 'TEMP2', IDATE, ITIME ) ) THEN
         MESG = 'Ending time not on met file'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Set output parameters that are different from met file and open file
      SDATE3D = SDATE                ! From run-script
      STIME3D = STIME                ! From run-script
      MXREC3D = RLENG/10000                ! From run-script
      MXREC = MXREC3D
      NLAYS3D = 1
      NVARS3D = NOUTPUT
      NCOLS = NCOLS3D
      NROWS = NROWS3D
      TSTEP = TSTEP3D
      TSTEP3D = 240000
      FDESC3D(:) = ' '
      FDESC3D(1) = 'Daily meteorology with daily average temperature/PPFD, daily
     &                                 max/min temperature/wind speed'

! Define output variables

      VNAME3D(1) = 'D_TEMP'
      UNITS3D(1) = 'K'
      VTYPE3D(1) = M3REAL
      VDESC3D(1) = 'Daily average temperature (K)'

      VNAME3D(2) = 'D_PPFD'
      UNITS3D(2) = 'umol/m2.s'
      VTYPE3D(2) = M3REAL
      VDESC3D(2) = 'Daily average PPFD (umol/m2.s)'

      VNAME3D(3) = 'MaxT'
      UNITS3D(3) = 'K'
      VTYPE3D(3) = M3REAL
      VDESC3D(3) = 'Daily maximum temperature (K)'

      VNAME3D(4) = 'MinT'
      UNITS3D(4) = 'K'
      VTYPE3D(4) = M3REAL
      VDESC3D(4) = 'Daily minimum temperature (K)'

      VNAME3D(5) = 'MaxWS'
      UNITS3D(5) = 'K'
      VTYPE3D(5) = M3REAL
      VDESC3D(5) = 'Daily maximum wind speed (m/s)'

      IF ( .NOT. OPEN3( OUTMET, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (OUTMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

! Allocate memory

      ALLOCATE ( TEMP(24,NCOLS,NROWS), STAT = IOS)
      CALL CHECKMEM    ( IOS, 'TEMP',     PROGNAME )
      ALLOCATE ( PPFD(24,NCOLS,NROWS), STAT = IOS)
      CALL CHECKMEM    ( IOS, 'PPFD',     PROGNAME )
      ALLOCATE ( WIND(24,NCOLS,NROWS), STAT = IOS)
      CALL CHECKMEM    ( IOS, 'WIND',     PROGNAME )
      ALLOCATE ( D_TEMP  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'D_TEMP',     PROGNAME )
      ALLOCATE ( D_PPFD ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'D_PPFD',     PROGNAME )
      ALLOCATE ( MaxT  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MaxT',     PROGNAME )
      ALLOCATE ( MinT  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MinT',     PROGNAME )
      ALLOCATE ( MaxWS   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MaxWS',     PROGNAME )

!--------------------------------------------------------------------------
!.....2) Read in data and calculate max/min/avg
!--------------------------------------------------------------------------

! Loop over episode days

      DO I = 1, NDAYS

        IDATE = SDATE; ITIME = STIME
        CALL NEXTIME ( IDATE, ITIME, (I-1)*240000 )
        
        WRITE(MESG,"I0.3") I
        MGNMET = 'MGNMET'//TRIM(MESG)

!...    Open files
        WRITE(MESG,'(A,10X, A)') 'Opening file', MGNMET
        CALL M3MESG( MESG )

! Met file
        IF ( .NOT. OPEN3( MGNMET, FSREAD3, PROGNAME ) ) THEN
           CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
           MESG = 'Could not open file '//TRIM(MESG)
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        ! Check grid
        IF ( .NOT. FILCHK3 ( MGNMET,
     &                GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
           MESG = 'MGNMET has differenet grid definition'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF
        IF ( .NOT. DESC3( MGNMET ) ) THEN
           CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
           MESG = 'Could not get description of '//TRIM(MESG)
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

!...  Check start date, start time, end date, end time in MGNMET
        IF ( .NOT. CHECK3( MGNMET, 'TEMP2', IDATE, ITIME ) ) THEN
          MESG = 'Starting time not on met file'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF
        CALL NEXTIME ( IDATE, ITIME, RLENG-10000 )
        IF ( .NOT. CHECK3( MGNMET, 'TEMP2', IDATE, ITIME ) ) THEN
          MESG = 'Ending time not on met file'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

!...  Read in hourly data
        AVEBY = MIN(24,MXREC)
        IDATE = SDATE;ITIME = STIME
        CALL NEXTIME ( IDATE, ITIME, (I-1)*240000 )
        TEMP = 0.0
        WIND = 0.0
        PPFD = 0.0
        ! Start the loop over the hours
        DO T = 1, AVEBY

          IF ( .NOT. READ3(MGNMET,'TEMP2',1,IDATE,ITIME,TEMP(T,:,:))) THEN
            MESG = 'Error reading TEMP2'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

          IF ( .NOT. READ3(MGNMET,'PAR',1,IDATE,ITIME,PPFD(T,:,:))) THEN
             MESG = 'Error reading PAR'
             CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

          IF ( .NOT. READ3(MGNMET,'WINDSPD',1,IDATE,ITIME,WIND(T,:,:))) THEN
             MESG = 'Error reading WIND'
             CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

          CALL NEXTIME( IDATE, ITIME, TSTEP )

        ENDDO ! End loop over hours

        MaxT = MAXVAL(TEMP,1)
        MinT = MINVAL(TEMP,1)
        MaxWS = MAXVAL(WIND,1)
        D_TEMP = (SUM(TEMP,1))/AVEBY
        ! Convert incoming PAR in W/m2 to umol/m2.s by multiplying 4.5
        D_PPFD = (SUM(PPFD,1))*4.5/AVEBY

!-----------------------------------------------------------------------
!.....3) Write out data
!-----------------------------------------------------------------------
!... Write met data to file
        IDATE = SDATE;ITIME = STIME
        CALL NEXTIME ( IDATE, ITIME, (I-1)*240000 )

        WRITE(MESG,1000) 'Writing met data at ',IDATE,ITIME
        CALL M3MESG( MESG )

! #1
        IF ( .NOT. WRITE3(OUTMET,'D_TEMP',IDATE,ITIME,
     &            D_TEMP)) THEN
          CALL NAMEVAL (OUTMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #2
        IF ( .NOT. WRITE3(OUTMET,'D_PPFD',IDATE,ITIME,
     &            D_PPFD)) THEN
          CALL NAMEVAL (OUTMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #3
        IF ( .NOT. WRITE3(OUTMET,'MaxT',IDATE,ITIME,
     &            MaxT)) THEN
          CALL NAMEVAL (OUTMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #4
        IF ( .NOT. WRITE3(OUTMET,'MinT',IDATE,ITIME,
     &            MinT)) THEN
          CALL NAMEVAL (OUTMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #5
        IF ( .NOT. WRITE3(OUTMET,'MaxWS',IDATE,ITIME,
     &            MaxWS)) THEN
          CALL NAMEVAL (OUTMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

      ENDDO ! End loop over days

      DEALLOCATE ( D_TEMP )
      DEALLOCATE ( D_PPFD )
      DEALLOCATE ( MaxT )
      DEALLOCATE ( MinT )
      DEALLOCATE ( MaxWS )
      DEALLOCATE ( TEMP )
      DEALLOCATE ( WIND )
      DEALLOCATE ( PPFD )
! ... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

!--=====================================================================
!...  FORMAT
!--=====================================================================
1000  FORMAT (A20,I8,X,I8)

!--=====================================================================
!...  End program
!--=====================================================================

       RETURN
       END PROGRAM DAYMET

