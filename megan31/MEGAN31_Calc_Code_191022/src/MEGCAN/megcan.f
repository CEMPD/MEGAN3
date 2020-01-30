	PROGRAM MEGCAN

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

      IMPLICIT NONE

!...  INCLUDES:
      INCLUDE 'MEGCAN.EXT'
      INCLUDE 'PARMS3.EXT'          !  I/O API parameters
      INCLUDE 'IODECL3.EXT'         !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'          !  I/O API file description data structures

!...  EXTERNAL FUNCTIONS and their descriptions:
      INTEGER, EXTERNAL       ::   ENVINT
      INTEGER, EXTERNAL       ::   ENVYN
      LOGICAL      DSCGRID
      EXTERNAL     DSCGRID

! input

!...  Program I/O files: From run script
! Program name
      CHARACTER*16  :: PROGNAME = 'MEGCAN'
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

      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

! variables from input file

      REAL, ALLOCATABLE :: LAT( :,: )     ! Latitude of grid cell
      REAL, ALLOCATABLE :: LONG( :,: )    ! Longitude of grid cell

      REAL, ALLOCATABLE :: LAIc( :,: )    ! Current step LAI
      REAL, ALLOCATABLE :: TEMP( :,: )    ! Temperature (K)
      REAL, ALLOCATABLE :: PPFD( :,: )    ! Calculated PAR (umol/m2.s)

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

! output
      REAL, ALLOCATABLE :: SunleafTK  ( :,:,: )   ! Sun leaf temperature [K]
      REAL, ALLOCATABLE :: ShadeleafTK( :,:,: ) ! Shade leaf temperature [K]
      REAL, ALLOCATABLE :: SunPPFD    ( :,:,: )     ! PPFD on a sun leaf [umol/m2.s]
      REAL, ALLOCATABLE :: ShadePPFD  ( :,:,: )   ! PPFD on a shade leaf [(umol/m2.s]
      REAL, ALLOCATABLE :: SunFrac    ( :,:,: )     ! fraction of sun leaves

! local variables and their descriptions:

      INTEGER :: Day
      REAL :: Sinbeta,Beta
      REAL,DIMENSION(Layers) ::  VPgausWt, VPgausDis2, 
     &  VPgausDis, VPslwWT,  QdAbsV, QsAbsV, QdAbsn,     
     &  QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn,
     &  TairK, HumidairPa, Ws, SunleafSH, sun_ppfd,shade_ppfd,            
     &  SunleafLH,SunleafIR, ShadeleafSH, sun_tk,shade_tk,sun_frac,
     &  ShadeleafLH,ShadeleafIR, sun_ppfd_total, shade_ppfd_total,
     &  sun_tk_total, shade_tk_total, sun_frac_total

      REAL :: Hour, Solar, Maxsolar,  
     &         Difffrac, PPFDfrac, QbAbsn,                  
     &         Trate, Qbeamv,Qdiffv, Qbeamn, Qdiffn,  
     &         QbAbsV,Ea1tCanopy, Ea1pCanopy,                    
     &         TairK0, HumidairPa0,Ws0, SH
                       
      REAL ::  CalcEccentricity,WaterVapPres,        
     &          Stability, Calcbeta

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
      IF ( .NOT. FILCHK3 ( CANTYP,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
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
      IF ( .NOT. FILCHK3 ( LAIS46,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
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

        IF ( .NOT. READ3(MGNMET,'TEMP2',  ALLAYS3,IDATE,ITIME,TEMP)) THEN
          MESG = 'Error reading temperature'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PAR',   ALLAYS3,IDATE,ITIME,PPFD)) THEN
          MESG = 'Error reading PAR'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
        !PPFD = PPFD * 4.766
        PPFD = PPFD * 4.5

        IF( .NOT. READ3(MGNMET,'WINDSPD',ALLAYS3,IDATE,ITIME,WIND)) THEN
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
        WRITE(MESG,1020) 'Found LAI current period for YYYYJJJ : ',
     &                    IDATE,LAIc_I
        CALL M3MESG( MESG )
        WRITE(MESG,'(I0.2)') LAIc_I
        IF ( .NOT. READ3(LAIS46,'LAI'//TRIM(MESG),ALLAYS3,0,0,LAIc)) THEN
          MESG = 'Error reading LAI'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

          WRITE(MESG,1030) 'Entering CANOPY: ',IDATE,ITIME
          CALL M3MESG( MESG )

          DO I=1, NCOLS
           DO J=1, NROWS
            SunleafTK(I,J,:)   = TEMP(I,J)
            ShadeleafTK(I,J,:) = TEMP(I,J)
            SunPPFD(I,J,:)     = PPFD(I,J)
            ShadePPFD(I,J,:)   = PPFD(I,J)
            SunFrac(I,J,:)     = 1.0
            TotalCT           = 0.0
            DO I_CT = 1,NRTYP   !canopy types
              TotalCT = TotalCT + CTF(I_CT,I,J) * 0.01
            ENDDO   ! ENDDO I_CT

            IF (TotalCT .GT. 0.0 .AND. LAIc(I,J) .GT. 0.0) THEN
!           only invoke canopy model when both CT and LAI are valid

              sun_ppfd_total     = 0.0
              shade_ppfd_total   = 0.0
              sun_tk_total       = 0.0
              shade_tk_total     = 0.0
              sun_frac_total     = 0.0

              DO I_CT = 1,NRTYP   !canopy types
                IF (CTF(I_CT,I,J) .NE. 0.0) THEN
                sun_ppfd           = 0.0
                shade_ppfd         = 0.0
                sun_tk             = 0.0
                shade_tk           = 0.0
                sun_frac           = 0.0
                Day = MOD(IDATE,1000)
!            Convert from XXXXXX format to XX.XX (solar hour)
!            HOUR = 0 -> 23.xx
!            Solar hour
                Hour  = ITIME/10000 + LONG(I,J) /15
                IF ( Hour  .LT. 0.0 ) THEN
                  Hour  = Hour  + 24.0
                  Day  = Day  - 1
                ELSEIF ( Hour  .GT. 24.0 ) THEN
                  print*,'Invalid hour: HOUR  is ', Hour
                ENDIF
!            Solar angle
                Beta   = Calcbeta(Day , Lat(I,J) , Hour )
                Sinbeta    = SIN(Beta  / 57.29578)
                TairK0     = TEMP(I,J)
                Ws0        = WIND(I,J)
!               Solar      = PPFD(I,J)/ConvertWm2toUmolm2s*2
                Solar      = PPFD(I,J)/2.25
                Maxsolar   = Sinbeta  * SolarConstant * CalcEccentricity(Day )
                Call GaussianDist(VPgausDis, Layers)
                Call SolarFractions(Solar, Maxsolar, Qdiffv,Qbeamv,Qdiffn,Qbeamn)
                Call CanopyRad(VPgausDis, Layers, LAIc(I,J), Sinbeta, Qbeamv,
     &               Qdiffv, Qbeamn, Qdiffn,I_CT ,Canopychar, sun_frac,
     &               QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv,
     &               ShadeQv, SunQn, ShadeQn, sun_ppfd, shade_ppfd,
     &               NrCha,NrTyp)

                 HumidairPa0  =  WaterVapPres(QV(I,J), PRES(I,J), WaterAirRatio)
                 Trate    =  Stability(Canopychar, I_CT, Solar , NrCha, NrTyp)
                 Call CanopyEB(Trate, Layers, VPgausDis, Canopychar, I_CT,
     &                         TairK, HumidairPa, Ws, sun_ppfd,
     &                         shade_ppfd, SunQv, ShadeQv, SunQn, ShadeQn,
     &                         sun_tk, SunleafSH, SunleafLH, SunleafIR,
     &                         shade_tk,ShadeleafSH,ShadeleafLH,ShadeleafIR,
     &                         NrCha, NrTyp, Ws0, TairK0, HumidairPa0)
                     sun_ppfd_total(:)   = sun_ppfd_total(:) + 
     &                                0.01*CTF(I_CT,I,J)*sun_ppfd(:)
                     shade_ppfd_total(:) = shade_ppfd_total(:) +
     &                                0.01*CTF(I_CT,I,J)*shade_ppfd(:)
                     sun_tk_total(:)     = sun_tk_total(:) +
     &                                0.01*CTF(I_CT,I,J)*sun_tk(:)
                     shade_tk_total(:)   = shade_tk_total(:) + 
     &                                0.01*CTF(I_CT,I,J)*shade_tk(:)
                     sun_frac_total(:)   = sun_frac_total(:) + 
     &                                0.01*CTF(I_CT,I,J)*sun_frac(:)
                ENDIF
              ENDDO  ! ENDDO I_CT
              SunleafTK(I,J,:)   = sun_tk_total(:)/TotalCT
              ShadeleafTK(I,J,:) = shade_tk_total(:)/TotalCT
              SunPPFD(I,J,:)     = sun_ppfd_total(:)/TotalCT
              ShadePPFD(I,J,:)   = shade_ppfd_total(:)/TotalCT
              SunFrac(I,J,:)     = sun_frac_total(:)/TotalCT

            ELSEIF( TotalCT .LT. 0.0) THEN
              CALL M3ERR(PROGNAME,IDATE,ITIME,
     &              'TotalCT is less than 0.0',.TRUE.)
            ELSE
            ! total CT is zero
            SunleafTK(I,J,:)   = TEMP(I,J)
            ShadeleafTK(I,J,:) = TEMP(I,J)
            SunPPFD(I,J,:)     = PPFD(I,J)
            ShadePPFD(I,J,:)   = PPFD(I,J)
            SunFrac(I,J,:)     = 1

            ENDIF

           ENDDO   ! ENDDO J
          ENDDO   ! ENDDO I
          WRITE(MESG,1030) 'Exited CANOPY: ',IDATE,ITIME
          CALL M3MESG( MESG )

!-----------------------------------------------------------------------
!.....3) Write out the calculated met data
!-----------------------------------------------------------------------
!... Write met data to file
        WRITE(MESG,1030) 'Writing met data at ',T,IDATE,ITIME
        CALL M3MESG( MESG )

! #1
        IF ( .NOT. WRITE3(CANMET,VNAME3D(1),IDATE,ITIME,
     &            SunleafTK)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #2
        IF ( .NOT. WRITE3(CANMET,VNAME3D(2),IDATE,ITIME,
     &            ShadeleafTK)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #3
        IF ( .NOT. WRITE3(CANMET,VNAME3D(3),IDATE,ITIME,
     &            SunPPFD)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

!#4
        IF ( .NOT. WRITE3(CANMET,VNAME3D(4),IDATE,ITIME,
     &            ShadePPFD)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #5
        IF ( .NOT. WRITE3(CANMET,VNAME3D(5),IDATE,ITIME,
     &               SunFrac)) THEN
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

       RETURN
       END PROGRAM MEGCAN

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE GaussianDist
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE GaussianDist
     &           (Distgauss, Layers)
 
      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      INTEGER,INTENT(IN) ::  Layers

      REAL,DIMENSION(Layers),INTENT(OUT) :: 
     &                         Distgauss

! local variables
      INTEGER ::  i
!----------------------------------------------------------------

      IF (Layers .EQ. 1) THEN
        Distgauss(1)   = 0.5
      ELSEIF (Layers .EQ. 3) THEN
        Distgauss(1)   = 0.112702
        Distgauss(2)   = 0.5
        Distgauss(3)   = 0.887298
      ELSEIF (Layers .EQ. 5) THEN
        Distgauss(1)   = 0.0469101
        Distgauss(2)   = 0.2307534
        Distgauss(3)   = 0.5
        Distgauss(4)   = 0.7692465
        Distgauss(5)   = 0.9530899
      ELSE
        DO i = 1, Layers
          Distgauss(i)   = (i - 0.5) / Layers
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE GaussianDist


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE SolarFractions
!   Based on actual and potential max solar radiation:
!   Determine the fraction of solar radiation that is 
!   diffuse PPFD, direct PPFD, diffuse near IR, direct near IR 
!
!   Originally developed by Alex Guenther in 1990s
!   Modified in 2010
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE SolarFractions( Solar, Maxsolar,  
     &                          Qdiffv,Qbeamv, Qdiffn, Qbeamn)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
 
      REAL,INTENT(IN) :: Solar, Maxsolar
      REAL,INTENT(OUT) ::  Qdiffv,Qbeamv, Qdiffn, Qbeamn
      REAL :: FracDiff, PPFDfrac,PPFDdifFrac,Qv, Qn
! internal variables
      INTEGER :: I,J
      REAL ::  Transmis 
!-----------------------------------------------------
      IF (Maxsolar  <= 0) THEN
        Transmis  = 0.5
       ELSEIF (Maxsolar < Solar) THEN
      Transmis  = 1.0
        ELSE
       Transmis  = Solar  / Maxsolar 
      ENDIF

!FracDiff is based on Lizaso 2005
        FracDiff = 0.156 + 0.86/(1 + EXP(11.1*(Transmis -0.53)))

!PPFDfrac is based on Goudrian and Laar 1994
        PPFDfrac  = 0.55 -Transmis*0.12

!PPFDdifFrac is based on data in Jacovides 2007
        PPFDdifFrac = FracDiff *(1.06 + Transmis*0.4)  

! Calculate  Qdiffv,Qbeamv, Qdiffn, Qbeamn in the subroutine

        IF (PPFDdifFrac > 1.0) THEN
        PPFDdifFrac  = 1.0
        ENDIF
 
        Qv  = PPFDfrac * Solar
        Qdiffv = Qv * PPFDdifFrac
        Qbeamv = Qv - Qdiffv
        Qn = Solar - Qv
        Qdiffn =  Qn * FracDiff
        Qbeamn =  Qn - Qdiffn
            
      RETURN
      END SUBROUTINE SolarFractions
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!   Subroutine CanopyRad
!
!   Canopy light environment model
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!   based on Spitters et al. (1986), 
!   Goudrian and van Laar (1994), Leuning (1997)
!   Initial code 8-99, 
!   modified 7-2000, 12-2001, 1-2017
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE CanopyRad(Distgauss, Layers, LAI, Sinbeta,   
     &            Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype,     
     &            Canopychar, Sunfrac, QbAbsV, QdAbsV, QsAbsV, 
     &            QbAbsn, QdAbsn, QsAbsn, SunQv,               
     &            ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, 
     &            NrCha, NrTyp)

      IMPLICIT NONE

! input
      INTEGER,INTENT(IN) :: Layers, NrCha, NrTyp, Cantype
      REAL,INTENT(IN) :: Qbeamv,Qdiffv,Sinbeta,LAI,Qbeamn,Qdiffn
      REAL,DIMENSION(Layers),INTENT(IN) :: Distgauss

! output
      REAL,INTENT(OUT) :: QbAbsV, QbAbsn

      REAL,DIMENSION(Layers),INTENT(OUT) :: ShadePPFD, SunPPFD,  
     &                QdAbsv, QsAbsv, QsAbsn, ShadeQv,  SunQn,    
     &                QdAbsn, SunQv, ShadeQn, Sunfrac

      REAL,DIMENSION(NrCha,NrTyp),INTENT(OUT) :: Canopychar

! internal variables
      INTEGER :: i
      REAL :: ScatV, ScatN, RefldV, RefldN, ReflbV, ReflbN,      
     &         Kb, Kd, KbpV, KbpN, KdpV, KdpN, LAIdepth, Cluster, 
     &         QdAbsVL, QsAbsVL, QdAbsNL, QsAbsNL, CANTRAN, LAIadj

! Stefan-boltzman constant  W m-2 K-4
      REAL,PARAMETER :: Sb = 0.0000000567    
      REAL,PARAMETER :: ConvertShadePPFD = 4.6
      REAL,PARAMETER :: ConvertSunPPFD = 4.0
!---------------------------------------------------------------------
        
! adjust LAI for canopy transparency
      CANTRAN = Canopychar(17,Cantype)
      LAIadj = LAI / ( 1 - CANTRAN )

      IF (((Qbeamv  + Qdiffv ) > 0.001) .AND.    
     &     (Sinbeta  > 0.002) .AND.             
     &     (LAIadj  > 0.001)) THEN       ! Daytime

! Scattering coefficients (scatV,scatN), diffuse and beam reflection 
! coefficients (ref..) for visible or near IR
        ScatV   = Canopychar(5,Cantype)
        ScatN   = Canopychar(6,Cantype)
        RefldV  = Canopychar(7,Cantype)
        RefldN  = Canopychar(8,Cantype)
        Cluster = Canopychar(9,Cantype)
!        print*,'cluster',  Cluster
! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
        Kb = Cluster * 0.5 / Sinbeta 
! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
        Kd = 0.8 * Cluster              
! (0.8 assumes a spherical leaf angle distribution)

        Call CalcExtCoeff(Qbeamv,ScatV,Kb,Kd,ReflbV,KbpV,KdpV,QbAbsV)
        Call CalcExtCoeff(Qbeamn,ScatN,Kb,Kd,ReflbN,KbpN,KdpN,QbAbsn)

        DO i = 1,layers
! LAI depth at this layer
          LAIdepth   = LAIadj  * Distgauss(i)  
!fraction of leaves that are sunlit
          Sunfrac(i) = EXP(-Kb * LAIdepth)  

          Call CalcRadComponents(Qdiffv , Qbeamv , kdpV,  
     &                          kbpV, kb, scatV, refldV, 
     &                          reflbV, LAIdepth, QdAbsVL, QsAbsVL)

          Call CalcRadComponents(Qdiffn , Qbeamn , kdpN,  
     &                          kbpN, kb, scatN, refldN, 
     &                          reflbN, LAIdepth, QdAbsNL, QsAbsNL)

       ShadePPFD(i) = (QdAbsVL + QsAbsVL) * ConvertShadePPFD/(1 - scatV)
       SunPPFD(i) = ShadePPFD(i) + (QbAbsV* ConvertSunPPFD/(1 - scatV))
          QdAbsV(i) = QdAbsVL
          QsAbsV(i) = QsAbsVL
          QdAbsn(i) = QdAbsNL
          QsAbsn(i) = QsAbsNL
          ShadeQv(i) = QdAbsVL + QsAbsVL
          SunQv(i)   = ShadeQv(i) + QbAbsV 
          ShadeQn(i) = QdAbsNL + QsAbsNL
          SunQn(i)   = ShadeQn(i) + QbAbsn 
        ENDDO

      ELSE                           ! Night time

        QbAbsV  = 0
        QbAbsn   = 0

        Sunfrac(:)   = 0.2
        SunQn(:)     = 0
        ShadeQn(:)   = 0
        SunQv(:)     = 0
        ShadeQv(:)   = 0
        SunPPFD(:)   = 0
        ShadePPFD(:) = 0
        QdAbsV(:)    = 0
        QsAbsV(:)    = 0
        QdAbsn(:)    = 0
        QsAbsn(:)    = 0

      ENDIF
         
      RETURN
      END SUBROUTINE CanopyRad

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CalcExtCoeff
!   Calculate light extinction coefficients
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE CalcExtCoeff(Qbeam,scat,kb,kd,reflb,
     &                        kbp,kdp,QbeamAbsorb)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      REAL,INTENT(IN) :: Qbeam, scat, Kb, Kd
      REAL,INTENT(OUT) :: Reflb, Kbp, Kdp, QbeamAbsorb

! local variables
      REAL :: P
!-------------------------------------------------------------------
      
      P     = (1 - scat)**0.5
      Reflb = 1 - Exp((-2 * ((1 - P) / (1 + P)) * kb) / (1 + kb))

! Extinction coefficients
      Kbp   = Kb * P
      Kdp   = Kd * P
! Absorbed beam radiation
      QbeamAbsorb = kb * Qbeam * (1 - scat) 

      RETURN
      END SUBROUTINE CalcExtCoeff

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CalcRadComponents
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE CalcRadComponents(Qdiff, Qbeam, kdp, kbp, kb,    
     &        scat, refld, reflb, LAIdepth, QdAbs, QsAbs)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      REAL,INTENT(IN) :: Qdiff,Qbeam,kdp,kbp,kb,scat,
     &                   refld,reflb,LAIdepth
      REAL,INTENT(OUT) :: QdAbs, QsAbs
!-------------------------------------------------------------------

      QdAbs = Qdiff *   Kdp * (1 - Refld) * Exp(-Kdp * LAIdepth)
      QsAbs = Qbeam * ((Kbp * (1 - Reflb) * Exp(-Kbp * LAIdepth)) - 
     &                  (Kb * (1 - Scat) * Exp(-Kb * LAIdepth)))

      RETURN
      END SUBROUTINE CalcRadComponents

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CanopyEB
!
!   Canopy energy balance model for estimating leaf temperature
!   Coded into FORTRAN by Xuemei Wang
!   Code developed by Alex Guenther in 1990s
!   based on Goudrian and Laar (1994) and Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!   Modified in 1-2017 by Alex Guenther and Ling Huang
!   to correct IR balance and atmos. emissivity
!   Note: i denotes an array containing a vertical profile 
!         through the canopy with 0 (above canopy conditions) 
!         plus 1 to number of canopy layers
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE CanopyEB(Trate, Layers, Distgauss, Canopychar,    
     &             Cantype, TairK, HumidairPa, Ws,       
     &             SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn,  
     &             Sunleaftk, SunleafSH, SunleafLH,                 
     &             SunleafIR, Shadeleaftk, ShadeleafSH,             
     &             ShadeleafLH, ShadeleafIR, NrCha, NrTyp, Ws0,     
     &             TairK0, HumidairPa0)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

! inputs
      INTEGER,INTENT(IN) :: NrCha, NrTyp, Layers, Cantype
      REAL,INTENT(IN) :: Trate, TairK0, HumidairPa0, Ws0
      REAL,DIMENSION(Layers),INTENT(IN) ::  Distgauss, SunQv,ShadeQv, 
     &            SunQn, ShadeQn, SunPPFD, ShadePPFD
      REAL,DIMENSION(NrCha, NrTyp),INTENT(IN)  :: Canopychar

! outputs
      REAL,DIMENSION(Layers),INTENT(OUT) :: HumidairPa,            
     &       Ws, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, TairK, 
     &       Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR

! local variables
      INTEGER :: i
      REAL :: Cdepth, Lwidth, Llength, Cheight, Eps, TranspireType,   
!     &         Deltah, UnexposedLeafIRin, ExposedLeafIRin, IRin,IRout
     &         Deltah, EmissAtm, LeafIR, IRin,IRout
      REAL,DIMENSION(Layers) :: Ldepth, Wsh
!-----------------------------------------------------------------------
         
      Cdepth        = Canopychar(1, Cantype)
      Lwidth        = Canopychar(2, Cantype)
      Llength       = Canopychar(3, Cantype)
      Cheight       = Canopychar(4, Cantype)
      Eps           = Canopychar(10,Cantype)
      TranspireType = Canopychar(11,Cantype)

      IF (TairK0  > 288) THEN
! Pa m-1  (humidity profile for T < 288)
        Deltah =  Canopychar(14, Cantype) / Cheight
      ELSEIF (TairK0  > 278) THEN
        Deltah =(Canopychar(14,Cantype)-((288-TairK0)/10) *             
     &          (Canopychar(14,Cantype)-Canopychar(15,Cantype)))/Cheight
      ELSE
! Pa m-1  (humidity profile for T <278)
        Deltah = Canopychar(15, Cantype) / Cheight
      ENDIF

      Ldepth(:)     = Cdepth * Distgauss(:)
      TairK(:)      = TairK0  + (Trate  * Ldepth(:))      ! check this
      HumidairPa(:) = HumidairPa0  + (Deltah * Ldepth(:))

      Wsh(:) = (Cheight-Ldepth(:)) - (Canopychar(16,Cantype) * Cheight)
      Ws(:)  = (Ws0*LOG(Wsh(:))/LOG(Cheight-Canopychar(16,Cantype)
     &           *Cheight))
      WHERE (Wsh(:) < 0.001) Ws(:) = 0.05

      DO i=1,Layers

! REVISE - Replace UnexposedLeafIR with LeafIR

!        IRin     = UnexposedLeafIRin(TairK(i), Eps)
!        ShadeleafIR(i) = 2 * IRin
!        SunleafIR(i) = 0.5*ExposedLeafIRin(HumidairPa0,TairK0)+1.5*IRin

! Apparent atmospheric emissivity for clear skies: 
! function of water vapor pressure (Pa) 
! and ambient Temperature (K) based on Brutsaert(1975) 
! referenced in Leuning (1997)
         EmissAtm        = 0.642 * (HumidairPa(i) / TairK(i))**(1./7.)   
         IRin            = LeafIR (TairK(i), EmissAtm)
         ShadeleafIR(i)  = IRin
         SunleafIR(i)    = IRin 

      ! Sun
        CALL LeafEB(SunPPFD(i), SunQv(i) + SunQn(i),                    
     &               SunleafIR(i), Eps, TranspireType, Lwidth, Llength, 
     &               TairK(i), HumidairPa(i), Ws(i),                   
     &               Sunleaftk(i), SunleafSH(i),SunleafLH(i),          
     &               IRout )

         SunleafIR(i) = SunleafIR(i) - IRout

      ! Shade
        CALL LeafEB(ShadePPFD(i), ShadeQv(i)+ShadeQn(i),               
     &               ShadeleafIR(i),Eps,TranspireType, Lwidth,Llength,
     &               TairK(i), HumidairPa(i), Ws(i),        
     &               Shadeleaftk(i), ShadeleafSH(i),ShadeleafLH(i),
     &               IRout)

         ShadeleafIR(i) = ShadeleafIR(i) - IRout
      ENDDO

      RETURN
      END SUBROUTINE CanopyEB

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine LeafEB
!
!   Leaf energy balance
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE LeafEB(PPFD, Q, IRin, Eps, TranspireType,             
     &         Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf,           
     &         SH, LH, IRout)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      REAL,INTENT(IN) :: Eps, TranspireType, Lwidth, Llength,PPFD, Q,  
     &                    IRin, TairK, HumidairPa, Ws
      REAL,INTENT(OUT) :: IRout, Tleaf, SH, LH

! local variables

      INTEGER :: i
      REAL :: HumidAirKgm3,GHforced,StomRes,IRoutairT,LHV,LatHv,Ws1, 
!     &        LHairT,Tdelt,Balance,LeafBLC,LeafH,LeafLE,LeafIRout,   
     &        LHairT,Tdelt,Balance,LeafBLC,LeafH,LeafLE,LeafIR,
     &        GH1,SH1,LH1,E1,ConvertHumidityPa2kgm3,ResSC,IRout1,GH
!----------------------------------------------------

      IF (Ws <= 0) THEN
        Ws1 = 0.001
      ELSE
        Ws1 = Ws
      ENDIF

      ! Air vapor density kg m-3
      HumidAirKgm3 = ConvertHumidityPa2kgm3(HumidairPa, TairK)

      ! Heat convection coefficient (W m-2 K-1) for forced convection. 
      ! Nobel page 366
      GHforced = 0.0259 / (0.004 * ((Llength / Ws)**0.5))

      ! Stomatal resistence s m-1
      StomRes  = ResSC(PPFD)

! REVISE - Replace LeafIRout with LeafIR
!      IRoutairT = LeafIROut(tairK, eps)
!XJ      IRoutairT  = LeafIR(TairK + Tdelt, Eps)
       IRoutairT = LeafIR(TairK, Eps)

      ! Latent heat of vaporization (J Kg-1)
      LatHv = LHV(TairK)

      ! Latent heat flux
      LHairT = LeafLE(TairK,HumidAirKgm3,LatHv,GHforced,StomRes,
     &                TranspireType)

      E1 = (Q + IRin - IRoutairT - LHairT)
      IF (E1 .EQ. 0.) THEN
        E1 = -1.
      ENDIF

      Tdelt = 1.
      Balance = 10.
      DO i = 1, 10
        IF (ABS(Balance) > 2) THEN
          ! Boundary layer conductance
          GH1 = LeafBLC(GHforced, Tdelt, Llength)
          ! Convective heat flux
          SH1 = LeafH(Tdelt, GH1)                
          ! Latent heat of vaporization (J Kg-1)
          LatHv = LHV(TairK + Tdelt)             
          LH = LeafLE(TairK + Tdelt, HumidAirKgm3,      
     &                 LatHv, GH1, StomRes, TranspireType)
          LH1 = LH - LHairT

! REVISE - Replace LeafIROut with LeafIR
!          IRout  = LeafIROut(TairK + Tdelt, Eps)
          IRout  = LeafIR(TairK + Tdelt, Eps)
          IRout1 = IRout - IRoutairT
          Tdelt  = E1 / ((SH1 + LH1 + IRout1) / Tdelt)
          Balance = Q + IRin - IRout - SH1 - LH
        ENDIF
      ENDDO

      If (Tdelt > 10)  Tdelt = 10
      If (Tdelt < -10) Tdelt = -10

      Tleaf = TairK + Tdelt
      GH    = LeafBLC(GHforced, Tleaf - TairK, Llength)
      SH    = LeafH(Tleaf - TairK, GH)
      LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv,         
     &                GH, StomRes, TranspireType)

! REVISE - Replace LeafIROut with LeafIR
!      IRout = LeafIROut(Tleaf, Eps)
      IRout = LeafIR(Tleaf, Eps)

      RETURN
      END SUBROUTINE LeafEB

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Calcbeta
!   Calculates the solar zenith angle
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION Calcbeta(Day, Lat, Hour)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      INTEGER :: Day

      REAL :: Rpi, Hour, Lat, SinDelta,                 
     &          CosDelta, A, B, Sinbeta, Calcbeta
      REAL,PARAMETER :: Pi = 3.14159, Rpi180 = 57.29578
!--------------------------------------------------------------------
      SinDelta = -SIN(0.40907) * COS(6.28 * (Day + 10) / (365))
      CosDelta = (1 - SinDelta**2.)**0.5

      A = SIN(Lat / Rpi180) * SinDelta
      B = COS(Lat / Rpi180) * Cosdelta
      Sinbeta = A + B * COS(2 * Pi * (Hour - 12) / 24)
      Calcbeta = ASIN(Sinbeta) * 57.29578

      END FUNCTION Calcbeta

! REVISE - Delete DIstomata
!!

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION CalcEccentricity
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION CalcEccentricity(Day)

      IMPLICIT NONE
      INTEGER :: Day
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: CalcEccentricity
!----------------------------------------------------------------

      CalcEccentricity = 1 + 0.033 * COS(2*3.14*(Day-10)/365)
    
      END FUNCTION CalcEccentricity

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION WaterVapPres
!
!   Convert water mixing ratio (kg/kg) to water vapor pressure 
!   (Pa or Kpa depending on units of input )
!   Mixing ratio (kg/kg), temp (C), pressure (KPa)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION WaterVapPres(Dens, Pres, WaterAirRatio)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Dens, Pres, WaterVapPres, WaterAirRatio
!----------------------------------------------------------------

      WaterVapPres = (Dens / (Dens + WaterAirRatio)) * Pres

      END FUNCTION WaterVapPres

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Stability
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION Stability(Canopychar, Cantype, Solar, NrCha, NrTyp)

      IMPLICIT NONE
      INTEGER :: Cantype, NrCha, NrTyp
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Solar, Trateboundary, Stability
      REAL,DIMENSION(NrCha, NrTyp)  :: Canopychar
!----------------------------------------------------------------

      Trateboundary = 500

      IF (Solar > Trateboundary) THEN
            ! Daytime temperature lapse rate
        Stability = Canopychar(12, Cantype)
      ELSEIF (Solar > 0) THEN
        Stability = Canopychar(12, Cantype) -                      
     &             ((Trateboundary - Solar) / Trateboundary) *     
     &    (Canopychar(12, Cantype) - Canopychar(13, Cantype))
       ELSE
            ! Nightime temperature lapse rate
         Stability = Canopychar(13, Cantype)
       ENDIF

       END FUNCTION Stability

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION ConvertHumidityPa2kgm3
!
!   Saturation vapor density  (kg/m3)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION ConvertHumidityPa2kgm3(Pa, Tk)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: ConvertHumidityPa2kgm3, Pa, Tk
!--------------------------------------------------------------------

      ConvertHumidityPa2kgm3 = 0.002165 * Pa / Tk

      END FUNCTION ConvertHumidityPa2kgm3

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION ResSC
!
!   Leaf stomatal cond. resistance s m-1
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION ResSC(Par)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Par, SCadj, ResSC
!----------------------------------------------------------------

      SCadj = ((0.0027 * 1.066 * Par) /    
     &        ((1 + 0.0027 * 0.0027 * Par**2.)**0.5))

      IF (SCadj < 0.1) THEN
        ResSC = 2000
      ELSE
        ResSC = 200 / SCadj
      ENDIF

      END FUNCTION ResSC


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafIR
!
!   Calculate IR transfer between leaf and air
!   Added by Alex Guenther and Ling Huang to replace previous
!   MEGAN2.1 IR balance functions
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

       FUNCTION LeafIR(Tk, Eps)

       IMPLICIT NONE
       INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
       REAL :: Eps, Tk, LeafIR
! Stefan-boltzman constant  W m-2 K-4
       REAL,PARAMETER :: Sb = 0.0000000567 
!----------------------------------------------------------------

       LeafIR = Eps * Sb * (2 * (Tk**4.)) 

       END FUNCTION LeafIR


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LHV
!
!   Latent Heat of vaporization(J Kg-1) from Stull p641
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LHV(Tk)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tk, LHV
!----------------------------------------------------------------

! REVISE - Replace 273 with 273.15
!      LHV = 2501000 - (2370 * (Tk - 273))
      LHV = 2501000 - (2370 * (Tk - 273.15))

      END FUNCTION LHV

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafLE
!
!   Latent energy term in Energy balance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tleaf, Ambvap, LatHv, GH, StomRes,                      
     &         TranspireType, SvdTk,LeafRes, Vapdeficit, LeafLE, LE
!----------------------------------------------------------------

      LeafRes    = (1 / (1.075 * (GH / 1231))) + StomRes
      Vapdeficit = (SvdTk(Tleaf) - Ambvap)

! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / 
!                 leaf resistence (s m-1)
      LE = TranspireType * (1 / LeafRes) * LatHv * Vapdeficit
      IF (LE < 0) THEN
        LeafLE = 0
      ELSE
        LeafLE = LE
      ENDIF

      END FUNCTION  LeafLE

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafBLC
!
!   Boundary layer conductance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LeafBLC(GHforced, Tdelta, Llength)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: GHforced, Tdelta, Llength, Ghfree, LeafBLC
!----------------------------------------------------------------

! This is based on Leuning 1995 p.1198 except using molecular 
! conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
! diffusivity so that you end up with a heat convection coefficient 
! (W m-2 K-1) instead of a conductance for free convection

      IF (Tdelta >= 0) THEN
         GhFree = 0.5 * 0.00253 * ((160000000 * Tdelta /            
     &             (Llength**3.))**0.25) / Llength
      ELSE
        GhFree = 0
      ENDIF
      LeafBLC = GHforced + GhFree

      END FUNCTION LeafBLC

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafH
!
!   Convective energy term in Energy balance (W m-2 heat flux 
!      from both sides of leaf)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LeafH(Tdelta, GH)
  
      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tdelta, GH, LeafH
!----------------------------------------------------------------

! 2 sides X conductance X Temperature gradient
      LeafH = 2 * GH * Tdelta
  
      END FUNCTION LeafH

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION SvdTk
!
!   Saturation vapor density  (kg/m3)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      Function SvdTk(Tk)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tk, Svp, SvdTk
!----------------------------------------------------------------

! Saturation vapor pressure (millibars)
      Svp = 10**((-2937.4 / Tk) - (4.9283 * LOG10(Tk)) + 23.5518)  
      SvdTk = 0.2165 * Svp / Tk

      END FUNCTION  SvdTk
