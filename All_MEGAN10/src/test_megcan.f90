	PROGRAM test_MEGCAN
	
! FSB This is the basic MEGCAN program except the all of the code processing
!    the canopy has been placed inside the subroutine MEGCANOPY
! FSB November 26,2019 UNC Institute of the Environment..



!   Based on code initiated by Alex Guenther in 1990s
!   Coded in FORTRAN by Xuemei Wang--Nov. 2007
!   Revised by Alex Guenther and Ling Huang in Feb 2017
!   to correct, modify, and update the code and make it
!   a stand-alone program
!
!*****************************************************************
!
!   Select Input/output files before running the program
!*****************************************************************
!   Input varibles
!
!   Day                  Julian day
!   Lat                  Latitude
!   Long                 Longitude
!   Hour                 Hour of the day
!   TEMP                 Temperature [K]
!   PPFD           Incoming photosynthetic active radiation [umol/m2/s1]
!   Wind                 Wind speed [m s-1]
!   Humidity             Relative humidity [%]
!   Cantype             Defines set of canopy characteristics
!   LAI                  Leaf area index [m2 per m2 ground area]
!   Pres                 Pressure [Pa]
!
!*****************************************************************
! Variables used
!
!   PPFDfrac             Fraction of solar radiation that is PPFD
!   Solar                Solar radiation [W/m2]
!   Maxsolar             Maximum of solar radiation
!   Sinbeta              Sin of solar angle above horizon
!   Beta                 Solar angle above horizon
!   TairK0               Above canopy air temperature [K]
!   TairK                Array of canopy air temperature [K]
!   Ws0                  Above canopy wind speed [m/s]
!   Ws                   Array of canopy wind speed [m/s]
!   HumidairPA0          Above canopy ambient humidity [Pa]
!   HumidairPa           Array of canopy ambient humidity in [Pa]
!   Transmis             Transmission of PPFD that is diffuse
!   Difffrac             Fraction of PPFD that is diffuse
!   PPFDfrac             Fraction of solar rad that is PPFD
!   Trate			   temperature vertical profile
!   QbAbsV, QbAbsN       Absorbed direct beam visible/near IR
!   QdAbsV, QdAbsN       Absorbed diffuse visible/near IR
!   QsAbsV, QsAbsN       Absorbed scattered visible//near IR
!   QBeamV, QBeamN       Above canopy direct beam visible/near IR
!   QDiffV, QDiffN       Above canopy diffuse visible/near IR
!
! Arrays with values for each canopy layer (vertical profile)
!   SunleafSH            sensible heat flux for sun leaves [W/m2]
!   SunleafLH            latent heat flux for sun leaves [W/m2]
!   SunleafIR            infrared flux for sun leaves [W/m2]
!   ShadeleafSH          sensible heat for shade leaves [W/m2]
!   ShadeleafLH          latent heat flux for shade leaves [W/m2]
!   ShadeleafIR          infrared flux for shade leaves [W/m2]
!   VPgausDis            gaussian weighting factors for distance
!   SunQv                visible radiation on sun leaves
!   ShadeQv              visible radiation on shade leaves
!   SunQn                near IR radiation on sun leaves
!   ShadeQn              near IR radiation on shade leaves
!   sun_ppfd             Array of incoming (NOT absorbed) PPFD on a sun leaf [umol/m2/s]
!   shade_ppfd           Array of incoming (NOT absorbed) PPFD on a shade leaf [umol/m2/s]
!   sun_tk               Array of leaf temperature for sun leaves [K]
!   shade_tk             Array of leaf temperature for shade leaves [K]
!   sun_frac             Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.

!*****************************************************************
! OUTPUT
! For each time step and location
! Each variable is an array with a value for each canopy layer
!			       (vertical profile)
! i = 1 is the top canopy layer, 2 is the next layer, etc.
!   ShadeleafTK          leaf temperature for shade leaves [K] (weighted by canopy type)
!   SunleafTK            leaf temperature for sun leaves [K] (weighted by canopy type)
!   SunFrac              fraction of sun leaves (weighted by canopy type)
!   SunPPFD              PPFD on a sun leaf [umol/m2/s] (weighted by canopy type)
!   ShadePPFD            PPFD on a shade leaf [umol/m2/s] (weighted by canopy type)
!
!*****************************************************************
! FUNCTIONS
!   Calcbeta             Calculation of solar zenith angle
!   WaterVapPres         Convert water mixing ratio (kg/kg) to water vapor
!   pressure
!   Stability            Temperature lapse rate in canopy
!   CalcEccentricity     Eccentricity of earth's orbit
!
           use all_megan
           USE M3UTILIO

      IMPLICIT NONE

!...  INCLUDES:
!      INCLUDE 'PARMS3.EXT'          !  I/O API parameters
!      INCLUDE 'IODECL3.EXT'         !  I/O API function declarations
!      INCLUDE 'FDESC3.EXT'          !  I/O API file description data structures
       
       
!...  EXTERNAL FUNCTIONS and their descriptions:
!      INTEGER, EXTERNAL       ::   ENVINT
!      INTEGER, EXTERNAL       ::   ENVYN
!      LOGICAL      DSCGRID
!      EXTERNAL     DSCGRID

! input

!...  Program I/O files: From run script
! Program name
      CHARACTER*16  :: PROGNAME = 'TEST_MEGCAN'
! Netcdf file
      CHARACTER*16  :: CANTYP = 'CANTYP'     ! canopy type file logical name
      CHARACTER*16  :: LAIS46 = 'LAIS46'     ! LAI file logical name
! Met files
      CHARACTER*16  :: MGNMET = 'MGNMET'     ! Met file logical name

! Output file (canopy meteorology)
      CHARACTER*16  :: CANMET = 'CANMET'     ! Output file logical name

!...  Parameters for file units
      INTEGER  LOGDEV                      ! Logfile unit number

! ... External parameters
! From run script
      INTEGER       SDATE          ! Start date YYYYDDD
      INTEGER       STIME          ! Start time HHMMSS
      INTEGER       RLENG          ! Run length HHMMSS

! I/O API file parameters
      INTEGER :: JDATE        ! Date YYYYDDD from inpname
      INTEGER :: JTIME        ! Time HHMMSS from inpname
      INTEGER :: NCOLS        ! Number of columns
      INTEGER :: NROWS        ! Number of rows
!     INTEGER :: LAYERS       ! Number of layers in canopy
      INTEGER :: MXREC        ! Total number of timesteps
      INTEGER :: TSTEP        ! Time step

!... Internal parameters
! Internal parameters (status and buffer)
      INTEGER       IOS                    ! i/o status
      CHARACTER*256 MESG                   ! message buffer

! Parameters for output species
      INTEGER, PARAMETER :: NOUTPUT = 5
                          ! number of output variables: sunleaftk, shadeleaftk,
                          ! sunPPFD, shadePPFD, sunfrac, minT, maxT, maxWS
!
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

! variables from input file

      REAL, ALLOCATABLE :: LAT( :,: )     ! Latitude of grid cell
      REAL, ALLOCATABLE :: LONG( :,: )    ! Longitude of grid cell

      REAL, ALLOCATABLE :: LAIc( :,: )    ! Current step LAI
      REAL, ALLOCATABLE :: TEMP( :,: )    ! Temperature (K)
      REAL, ALLOCATABLE :: PPFD( :,: )    ! Calculated PAR (umol/m2.s)
!     REAL, ALLOCATABLE :: PAR ( :,: )    ! Calculated PAR (W/m2)

      REAL, ALLOCATABLE :: WIND( :,: )
      REAL, ALLOCATABLE :: PRES( :,: )
      REAL, ALLOCATABLE :: QV( :,: )

      REAL, ALLOCATABLE :: CTF( :, :, : ) ! Canopy type factor array

      REAL :: TotalCT

      ! loop indices
      INTEGER :: IDATE             ! Looping
      INTEGER :: ITIME             ! Looping
      INTEGER :: I_CT, N, T, I, J
      INTEGER :: LAIp_I, LAIc_I
      INTEGER :: MXCT, MXLAI
      INTEGER :: YR
! output
      REAL, ALLOCATABLE :: SunleafTK  ( :,:,: )  ! Sun leaf temperature [K]
      REAL, ALLOCATABLE :: ShadeleafTK( :,:,: )  ! Shade leaf temperature [K]
      REAL, ALLOCATABLE :: SunPPFD    ( :,:,: )  ! PPFD on a sun leaf [umol/m2.s]
      REAL, ALLOCATABLE :: ShadePPFD  ( :,:,: )  ! PPFD on a shade leaf [(umol/m2.s]
      REAL, ALLOCATABLE :: SunFrac    ( :,:,: )   ! fraction of sun leaves

! local variables and their descriptions:

      INTEGER :: Day
      REAL    :: DOY,  ZTIME, YEAR, Hour           

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
      IF( .NOT. DSCGRID( GDNAM, CNAME, GDTYP3D,                     &
                   P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,     &
                   XORIG3D, YORIG3D, XCELL3D, YCELL3D,              &      
                   NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Open files
      WRITE(MESG,1030) 'Checking up files',0,0,0
      CALL M3MESG( MESG )
! Canopy type file
      IF ( .NOT. OPEN3( CANTYP, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (CANTYP, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( CANTYP,                                 &
                   GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'CANTYP has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( CANTYP ) ) THEN
         CALL NAMEVAL (CANTYP, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      MXCT = MXREC3D

! LAI file
      IF ( .NOT. OPEN3( LAIS46, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (LAIS46, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( LAIS46,                                 &
                   GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'LAIS46 has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( LAIS46 ) ) THEN
         CALL NAMEVAL (LAIS46, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      MXLAI = MXREC3D

! Met file
      IF ( .NOT. OPEN3( MGNMET, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( MGNMET,                                 &
                   GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'MGNMET has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( MGNMET ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      NCOLS = NCOLS3D
      NROWS = NROWS3D
      TSTEP = TSTEP3D

!...  Get input parameters from run script
      MESG = 'Model start date (YYYYDDD)'
      SDATE = ENVINT( 'SDATE', MESG, JDATE, IOS )

      MESG = 'Model start time (HHMMSS)'
      STIME = ENVINT( 'STIME', MESG, JTIME, IOS )

      MESG = 'Model run length (HHMMSS)'
      RLENG = ENVINT( 'RLENG', MESG, MXREC*10000, IOS )

!...  Check start date, start time, end date, end time in MGNMET
      WRITE(MESG,1030) 'Checking up MGNMET',0,0,0
      CALL M3MESG( MESG )
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
      IDATE = SDATE; ITIME = STIME

!...  Set output parameters that are different from met file and open file
      SDATE3D = SDATE                ! From run-script
      STIME3D = STIME                ! From run-script
      MXREC3D = RLENG / 10000        ! From run-script
      MXREC   = MXREC3D
      NVARS3D = NOUTPUT
      NLAYS3D = Layers
      write(*,*) SDATE3D,STIME3D,MXREC3D,MXREC,NVARS3D,NLAYS3D

      VNAME3D(1) = 'SunleafTK'
      UNITS3D(1) = 'K'
      VTYPE3D(1) = M3REAL
      VDESC3D(1) = 'Sun leaf temperature (K)'

      VNAME3D(2) = 'ShadeleafTK'
      UNITS3D(2) = 'K'
      VTYPE3D(2) = M3REAL
      VDESC3D(2) = 'Shade leaf temperature (K)'

      VNAME3D(3) = 'SunPPFD'
      UNITS3D(3) = 'umol/m2.s'
      VTYPE3D(3) = M3REAL
      VDESC3D(3) = 'Sun leaf PPFD (umol/m2.s)'

      VNAME3D(4) = 'ShadePPFD'
      UNITS3D(4) = 'umol/m2.s'
      VTYPE3D(4) = M3REAL
      VDESC3D(4) = 'Shade leaf PPFD (umol/m2.s)'

      VNAME3D(5) = 'SunFrac'
      UNITS3D(5) = 'fraction'
      VTYPE3D(5) = M3REAL
      VDESC3D(5) = 'Sun leaf fraction'

      CALL NAMEVAL (CANMET, MESG)  ! get output file name and path
      FDESC3D(:) = ' '
      FDESC3D(1) = 'Output CANMET file: '//TRIM(MESG)

      IF ( .NOT. OPEN3( CANMET, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (CANMET, MESG)  ! get output file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!-----------------------------------------------------------------------
!.....2) Process canopy model
!-----------------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( SunleafTK ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunleafTK',    PROGNAME )
      ALLOCATE ( ShadeleafTK ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'ShadeleafTK',    PROGNAME )
      ALLOCATE ( SunPPFD ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunPPFD',    PROGNAME )
      ALLOCATE ( ShadePPFD ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'ShadePPFD',    PROGNAME )
      ALLOCATE ( SunFrac ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunFrac',    PROGNAME )
      ALLOCATE ( LAT   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAT',    PROGNAME )
      ALLOCATE ( LONG  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LONG',   PROGNAME )
      ALLOCATE ( LAIc  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAIc',   PROGNAME )
      ALLOCATE ( PPFD  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PPFD',   PROGNAME )
      ALLOCATE ( TEMP  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'TEMP',   PROGNAME )
      ALLOCATE ( WIND  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'WIND',   PROGNAME )
      ALLOCATE ( PRES  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PRES',   PROGNAME )
      ALLOCATE ( QV    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'QV',   PROGNAME )
      ALLOCATE ( CTF( NRTYP, NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'CTF', PROGNAME )

!...  Read LAIS46

      IF ( .NOT. READ3(LAIS46,'LAT',1,0,0,LAT)) THEN
         MESG = 'Error reading LAT'
         CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. READ3(LAIS46,'LONG',1,0,0,LONG)) THEN
         MESG = 'Error reading LONG'
         CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
      ENDIF

!...  Read CANTYP
      DO N = 1, MXCT
         IF ( .NOT. READ3(CANTYP,'CTS',1,0,(N-1)*10000,CTF(N,:,:))) THEN
            MESG = 'Error reading CTS'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
         ENDIF
      ENDDO

!...  Start the loop over the time period
      IDATE = SDATE
      ITIME = STIME
      DO T = 1, MXREC
        WRITE(MESG,1030) 'Processing: ',T,IDATE,ITIME
        CALL M3MESG( MESG )
!...  Initialize hourly variables
        TEMP = 0.
        PPFD = 0.
        LAIc = 0.

        IF ( .NOT. READ3(MGNMET,'TEMP2', ALLAYS3,IDATE,ITIME,TEMP)) THEN
          MESG = 'Error reading temperature'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PAR', ALLAYS3,IDATE,ITIME,PPFD)) THEN
          MESG = 'Error reading PAR'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
        !PPFD = PPFD * 4.766
        PPFD = PPFD * 4.5

        IF( .NOT. READ3(MGNMET,'WINDSPD', ALLAYS3,IDATE,ITIME,WIND)) THEN
          MESG = 'Error reading wind speed'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PRES',  ALLAYS3,IDATE,ITIME,PRES)) THEN
          MESG = 'Error reading pressure'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'QV',    ALLAYS3,IDATE,ITIME,QV)) THEN
          MESG = 'Error reading QV'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        ! Find LAIc from date
        CALL FINDLAI(IDATE,MXLAI,LAIp_I,LAIc_I)
        write(*,*) IDATE,MXLAI,LAIp_I,LAIc_I
        WRITE(MESG,1020) 'Found LAI current period for YYYYJJJ : ',  &
                          IDATE,LAIc_I
        CALL M3MESG( MESG )
        WRITE(MESG,'(I0.2)') LAIc_I
        IF ( .NOT. READ3(LAIS46,'LAI'//TRIM(MESG),ALLAYS3,0,0,LAIc)) THEN
          MESG = 'Error reading LAI'//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

          WRITE(MESG,1030) 'Entering CANOPY: ',IDATE,ITIME
          CALL M3MESG( MESG )

! FSB The loops over NCOLS, and NROWS are now inside the subroutine
! FSB capture the year and day of the year from IDATE

                Day = MOD(IDATE,1000)
                YR  = (IDATE - Day) / 1000
                DOY = float(Day)
                YEAR = float(YR)
                
!            Convert from XXXXXX format to XX.XX (solar hour)
!            HOUR = 0 -> 23.xx
!            Solar hour
                Hour  = float(ITIME)/10000.0 + LONG(I,J) / 15.0
                IF ( Hour  .LT. 0.0 ) THEN
                  Hour  = Hour  + 24.0
                  Day  = Day  - 1
                ELSEIF ( Hour  .GT. 24.0 ) THEN
                  print*,'Invalid hour: HOUR  is ', Hour
                ENDIF

                ZTIME = Hour
                
        CALL MEGCANOPY( NCOLS, NROWS, LAYERS, YEAR, DOY, ZTIME,         &
                 LAT, LONG, LAIc, TEMP, PPFD, WIND, PRES, QV, CTF,       &
                 ShadeleafTK, SunleafTK, SunFrac, SunPPFD, ShadePPFD )
           
        WRITE(MESG,1030) 'Exited CANOPY: ',IDATE,ITIME
        CALL M3MESG( MESG )

!-----------------------------------------------------------------------
!.....3) Write out the calculated met data
!-----------------------------------------------------------------------
!... Write met data to file
        WRITE(MESG,1030) 'Writing met data at ',T,IDATE,ITIME
        CALL M3MESG( MESG )

! #1
        IF ( .NOT. WRITE3(CANMET,VNAME3D(1),IDATE,ITIME,         &
                 SunleafTK)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #2
        IF ( .NOT. WRITE3(CANMET,VNAME3D(2),IDATE,ITIME,        &
                 ShadeleafTK)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #3
        IF ( .NOT. WRITE3(CANMET,VNAME3D(3),IDATE,ITIME,        &
                 SunPPFD)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

!#4
        IF ( .NOT. WRITE3(CANMET,VNAME3D(4),IDATE,ITIME,         &
                 ShadePPFD)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #5
        IF ( .NOT. WRITE3(CANMET,VNAME3D(5),IDATE,ITIME,          &
                    SunFrac)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

       CALL NEXTIME( IDATE, ITIME, TSTEP )
      ENDDO ! End loop for time step (T)

!
! ... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

      DEALLOCATE( SunleafTK     )
      DEALLOCATE( ShadeleafTK   )
      DEALLOCATE( SunPPFD       )
      DEALLOCATE( ShadePPFD     )
      DEALLOCATE( SunFrac       )

      DEALLOCATE ( LAT     )   ! input latitude of grid cell
      DEALLOCATE ( LONG    )   ! input longitude of grid cell

      DEALLOCATE ( LAIc    )   ! current monthly LAI

      DEALLOCATE ( TEMP    )   ! input hourly temperature (K)
      DEALLOCATE ( PPFD    )   ! calculated PAR (umol/m2.s)

      DEALLOCATE ( WIND    )
      DEALLOCATE ( PRES    )
      DEALLOCATE ( QV      )
      DEALLOCATE ( CTF    )

!--=====================================================================
!...  FORMAT
!--=====================================================================
1000  FORMAT( A )
1010  FORMAT( 43( A, :, I8, :, 1X ) )
1020  FORMAT (A40,I8,X,I8,X,I8)
1030  FORMAT (A20,I8,X,I8,X,I8)

!--=====================================================================
!...  End program
!--=====================================================================      
   
       END 	PROGRAM test_MEGCAN

