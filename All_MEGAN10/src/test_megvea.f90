      PROGRAM TEST_MEGVEA

! PURPOSE: Calculate Vegetation emission activity (EA) for each emission
!		class as the product of EA for individual drivers
!		calculated by the following functions
!
! Vegetation Emission Activity (EA) algorithm FUNCTIONS
!
!   GAMTLD: EA Temperature response (light dependent emission) 
!   GAMTLI: EA Temperature response (light independent emission)
!   GAMP: EA Light response
!   GAMTP: combines GAMLD, GAMLI, GAMP to get canopy average
!
!   GAMLA: EA leaf age response 
!   GAMBD: EA bidirectional exchange LAI response
!
!   CDEA: Canopy depth emission response
!
!   GAMHW: EA response to high wind storms
!   GAMHT: EA resposne to high temperature
!   GAMLT: EA response to low temperature
!
!   GAMAQ: EA response to air pollution
!
!   GAMCO2: EA CO2 response (only applied to isoprene)
!   GAMSM: EA response to soil moisture (multiplied with LDF)
!
! INCLUDE FILES
!     'PARMS3.EXT'   ! I/O API parameters
!     'IODECL3.EXT'  ! I/O API function declarations
!     'FDESC3.EXT'   ! I/O API file desc. data structures
!     'MEGVEA.EXT'    ! coefficients
!
!  INPUT Files
!	Single value for each location
!		LDF: Light dependent fraction (for categories other than
!		monoterpene, use constant values from MEGVEA.EXT)
!		AQ:  Air Quality indicator
!	Time series for each location
!		LAI: Leaf Area Index 
!		GAMSM: emission activity response to soil moisture
!		MaxT: Daily maximum temperature (K)
!		MinT: Daily minimum temperature (K)
!		MaxWS: Daily mximum wind speed (m/s)
!               D_TEMP: Daily average temperature (K)
!               D_PPFD: Daily averaged PPFD (umol/m2.s)
!	Hourly time series for each canopy layer in each location
!		sunfrac: fraction of leaves that are sunlit
!		SUNT: leaf Temperature of sunlit leaves (K)
!		SUNP: sun leaf PPFD visible light (micromol/m2/s)
!		SHAT: leaf Temperature of shade leaves (K)
!		SHAP: shade leaf PPFD visible light (micromol/m2/s)
!	
!  OUTPUT Files
!	Hourly time series for each location 
!		Emission activity for each of the 20 emission types
!
!
! HISTORY:
!   Based on code initiated by Alex Guenther in 1990s
!   Coded in FORTRAN as
!   MEGAN: Jack Chen 11/04
!   MEGANv2.04: Tan 11/21/06
!   MEGANv2.1: X. Wang 11/04/2007
!       Modified by Julia Lee-Taylor 03/18/2008
!       Modified by Xuemei Wang 09/30/2008
!       Modified by Tan 07/28/2011
!   MEGAN3.0:
!   Alex Guenther and Ling Huang Feb 2017
!   Frank Binkowski and Dongmei Yang 2019     
!       Modulize to make MEGAN3 used in CMAQ inline
!       Remove AQFILE  
!       Remove stress options 

      USE All_MEGAN
      
      IMPLICIT NONE

! INCLUDE FILES
      INCLUDE 'PARMS3.EXT'   ! I/O API parameters
      INCLUDE 'IODECL3.EXT'  ! I/O API function declarations
      INCLUDE 'FDESC3.EXT'   ! I/O API file desc. data structures

!...  EXTERNAL FUNCTIONS and their descriptions:
      INTEGER, EXTERNAL       ::   ENVINT
      INTEGER, EXTERNAL       ::   ENVYN
      LOGICAL      DSCGRID
      EXTERNAL     DSCGRID

!...  Program I/O files: From run script
! Program name
      CHARACTER*16  :: PROGNAME = 'MEGVEA'
! Netcdf file
      CHARACTER*16  :: LAIS46   = 'LAIS46'     ! LAI file logical name
!     CHARACTER*16  :: AQFILE   = 'AQFILE'     ! Air quality index file
!     CHARACTER*16  :: SMFILE   = 'SMFILE'     ! Soil moisture acitvity file
      CHARACTER*16  :: CANMET   = 'CANMET'     ! Canopy meteorology file
      CHARACTER*16  :: DailyMET = 'DailyMET'   ! Daily meteorology file
      CHARACTER*16  :: LDFILE   = 'LDFILE'     ! Light dependent fraction file
                                               ! for certain MT/SQT

! output file
      CHARACTER*16  :: MGNERS = 'MGNERS'       ! Emission activity 

!...  Parameters for file units
      INTEGER  LOGDEV                          ! Logfile unit number

!...  External parameters
! From run script
      INTEGER       SDATE          ! Start date YYYYDDD
      INTEGER       STIME          ! Start time HHMMSS
      INTEGER       RLENG          ! Run length HHMMSS
!     INTEGER       Layers         ! Vertical layers
!     INTEGER       N_MaxT         ! number of past days for maximum temperature
!     INTEGER       N_MinT         ! number of past days for minimum temperature
!     INTEGER       N_MaxWS        ! number of past days for maximum wind speed

!     LOGICAL       GAMBD_YN       ! flag for GAMBD
!     LOGICAL       GAMAQ_YN       ! flag for GAMAQ
!     LOGICAL       GAMCO2_YN      ! flag for GAMCO2
!     LOGICAL       GAMHW_YN       ! flag for GAMHW
!     LOGICAL       GAMHT_YN       ! flag for GAMHT
!     LOGICAL       GAMLT_YN       ! flag for GAMLT
!     LOGICAL       GAMSM_YN       ! flag for GAMSM

! I/O API file parameters
      INTEGER :: JDATE        ! Date YYYYDDD from inpname
      INTEGER :: JTIME        ! Time HHMMSS from inpname
      INTEGER :: NCOLS        ! Number of columns
      INTEGER :: NROWS        ! Number of rows
      INTEGER :: NLAYS        ! Number of vertical layers
      INTEGER :: MXREC        ! Total number of timesteps
      INTEGER :: TSTEP        ! Time step

!...  Internal parameters
! Internal parameters (status and buffer)
      INTEGER       IOS                    ! i/o status
      CHARACTER*256 MESG                   ! message buffer

! Parameter for output species
!      INTEGER, PARAMETER :: NEMIS = NCLASS
                          ! number of emission classes

! Local variables and their descriptions:
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

      REAL, ALLOCATABLE :: ER( :,: )      ! Output emission buffer
      REAL, ALLOCATABLE :: NON_DIMGARMA (:,:,:)

      REAL, ALLOCATABLE :: GAMLA( :,: )      ! EA leaf age response
      REAL, ALLOCATABLE :: GAMSM( :,: )      ! emission activity response to soil moisture
      REAL, ALLOCATABLE :: GAMAQ( :,: )      ! EA response to air pollution
      REAL, ALLOCATABLE :: AQI( :,: )        ! air quality index (i.e.W126)
      REAL, ALLOCATABLE :: GAMBD( :,: )      ! EA bidirectional exchange LAI response
      REAL, ALLOCATABLE :: GAMHT( :,: )      ! EA response to high temperature
      REAL, ALLOCATABLE :: GAMLT( :,: )      ! EA response to low temperature
      REAL, ALLOCATABLE :: GAMHW( :,: )      ! EA response to high wind speed
      REAL, ALLOCATABLE :: GAMCO2( :,: )     ! EA response to CO2
      REAL, ALLOCATABLE :: GAMTP( :,: )      ! combines GAMLD, GAMLI, GAMP to get canopy average
      REAL, ALLOCATABLE :: LDFMAP( :,: )     ! light depenedent fraction map

      REAL, ALLOCATABLE :: LAIp( :,: )    ! Previous time step LAI
      REAL, ALLOCATABLE :: LAIc( :,: )    ! Current time step LAI

!     variable from CANMET
      REAL, ALLOCATABLE :: SunT (:,:,:)   ! Sun leaf temperature (K)
      REAL, ALLOCATABLE :: ShaT (:,:,:)   ! Shade leaf temperature (K)
      REAL, ALLOCATABLE :: SunP (:,:,:)   ! Sun leaf PPFD (umol/m2.s)
      REAL, ALLOCATABLE :: ShaP (:,:,:)   ! Shade leaf PPFD (umol/m2.s)
      REAL, ALLOCATABLE :: SunF (:,:,:)   ! Sun leaf fraction
      REAL, ALLOCATABLE :: CDEA (:,:,:)   ! Emission response to canopy depth

!     variable from DailyMET
      REAL, ALLOCATABLE :: Max_temp( :,:,: )
      REAL, ALLOCATABLE :: Max_wind( :,:,: )
      REAL, ALLOCATABLE :: Min_temp( :,:,: )
      REAL, ALLOCATABLE :: D_TEMP(:,:)    ! Daily averagye temperature (K)
      REAL, ALLOCATABLE :: D_PPFD(:,:)    ! Daily averagye PPFD (umol/m2.s)

      REAL, ALLOCATABLE :: MaxT( :,: )    ! Maximum temperature of
                                          ! previous n days (K)
      REAL, ALLOCATABLE :: MinT( :,: )    ! Minimum temperature of
                                          ! previous n days (K)
      REAL, ALLOCATABLE :: MaxWS( :,: )   ! Maximum wind speed of
                                          ! previous n days (m/s)

      REAL, ALLOCATABLE :: VPGWT(:), Ea1L(:), Ea2L(:)
!      REAL :: GAMP, GAMTLD, GAMTLI
      INTEGER :: LAIp_I,LAIc_I
      INTEGER :: MXLAI

! loop indices
      INTEGER :: IDATE, ITIME
      INTEGER :: S, T, I, J, K

!***********************************************************************

!--=====================================================================
!...  Begin program
!--=====================================================================

!----------------------------------------------------------------
!.....1) File set up and assign I/O parameters
!----------------------------------------------------------------
!...  Initialize log file unit
      LOGDEV = INIT3()
           !  Now I/O API is set up, and LOGUNIT is the unit number
           !  for the log file (or it 6 for st'd output).

!...  Get input parameters from run script
      MESG = 'Model start date (YYYYDDD)'
      SDATE = ENVINT( 'SDATE', MESG, JDATE, IOS )

      MESG = 'Model start time (HHMMSS)'
      STIME = ENVINT( 'STIME', MESG, JTIME, IOS )

      MESG = 'Model run length (HHMMSS)'
      RLENG = ENVINT( 'RLENG', MESG, MXREC*10000, IOS )

!     MESG = 'Canopy vertical layers'
!     Layers = ENVINT( 'Layers', MESG, 5, IOS )

!     MESG = 'Apply air quality stress'
!     GAMAQ_YN = ENVYN ( 'GAMAQ_YN', MESG, .FALSE., IOS )

!     MESG = 'Apply CO2 stress'
!     GAMCO2_YN = ENVYN ( 'GAMCO2_YN', MESG, .FALSE., IOS )

!     MESG = 'Apply emission response to high wind storms '
!     GAMHW_YN = ENVYN ( 'GAMHW_YN', MESG, .FALSE., IOS )
 
!     MESG = 'Apply emission response to high temperature '
!     GAMHT_YN = ENVYN ( 'GAMHT_YN', MESG, .FALSE., IOS )

!     MESG = 'Apply emission response to low temperature '
!     GAMLT_YN = ENVYN ( 'GAMLT_YN', MESG, .FALSE., IOS )

!     MESG = 'Apply emission response to soil moisture '
!     GAMSM_YN = ENVYN ( 'GAMSM_YN', MESG, .FALSE., IOS )

!     MESG = 'Apply emission response to LAI'
!     GAMBD_YN = ENVYN ( 'GAMBD_YN', MESG, .FALSE., IOS )

      !MESG = 'Number of past days for maximum temperature to be used'
      !N_MaxT = ENVINT( 'N_MaxT', MESG, 5, IOS )

      !MESG = 'Number of past days for minimum temperature to be used'
      !N_MinT = ENVINT( 'N_MinT', MESG, 5, IOS )

      !MESG = 'Number of past days for maximum wind speed to be used'
      !N_MaxWS = ENVINT( 'N_MaxWS', MESG, 5, IOS )

      CALL ENVSTR( 'GDNAM3D', MESG, 'ASACA36km', GDNAM, IOS )
      IF( .NOT. DSCGRID( GDNAM, CNAME, GDTYP3D,                  &
                   P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, &
                   XORIG3D, YORIG3D, XCELL3D, YCELL3D,          &
                   NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Open files

! LAI file
      WRITE(MESG,1030) 'Checking up LAI file',0,0,0
      CALL M3MESG( MESG )
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

! AQI file
!     IF ( GAMAQ_YN ) THEN
!       WRITE(MESG,1030) 'Checking up AQ file',0,0,0
!       CALL M3MESG( MESG )
!       IF ( .NOT. OPEN3( AQFILE, FSREAD3, PROGNAME ) ) THEN
!          CALL NAMEVAL (AQFILE, MESG)  ! get input file name and path
!          MESG = 'Could not open file '//TRIM(MESG)
!          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!       ENDIF
!       ! Check grid
!       IF ( .NOT. FILCHK3 ( AQFILE,                               &
!                  GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
!          MESG = 'AQFILE has differenet grid definition'
!          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!       ENDIF
!       IF ( .NOT. DESC3( AQFILE ) ) THEN
!          CALL NAMEVAL (AQFILE, MESG)  ! get input file name and path
!          MESG = 'Could not get description of '//TRIM(MESG)
!          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!       ENDIF
!     ENDIF

! LDF file
      WRITE(MESG,1030) 'Checking up LDF file',0,0,0
      CALL M3MESG( MESG )
      IF ( .NOT. OPEN3( LDFILE, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (LDFILE, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( LDFILE,                              &
                  GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'LDFILE has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( LDFILE ) ) THEN
         CALL NAMEVAL (LDFILE, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

! GAMSM
!     IF ( GAMSM_YN ) THEN
!        WRITE(MESG,1030) 'Checking up GAMSMfile',0,0,0
!        CALL M3MESG( MESG )
!        IF ( .NOT. OPEN3( SMFILE, FSREAD3, PROGNAME ) ) THEN
!           CALL NAMEVAL (SMFILE, MESG)  ! get input file name and path
!           MESG = 'Could not open file '//TRIM(MESG)
!           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!        ENDIF
!        ! Check grid
!        IF ( .NOT. FILCHK3 ( SMFILE,                                &
!                     GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
!           MESG = 'SMFILE has differenet grid definition'
!           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!        ENDIF
!        IF ( .NOT. DESC3( SMFILE ) ) THEN
!           CALL NAMEVAL (SMFILE, MESG)  ! get input file name and path
!           MESG = 'Could not get description of '//TRIM(MESG)
!           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!        ENDIF
!     ENDIF

! DAYMET
      WRITE(MESG,1030) 'Checking up DAYMET',0,0,0
      CALL M3MESG( MESG )
      IF ( .NOT. OPEN3( DailyMET, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (DailyMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( DailyMET,                               &
                    GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'DailyMET has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( DailyMET ) ) THEN
         CALL NAMEVAL (DailyMET, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

! CANMET
      WRITE(MESG,1030) 'Checking up CANMET',0,0,0
      CALL M3MESG( MESG )
      IF ( .NOT. OPEN3( CANMET, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( CANMET,                                 &
                   GRDDED3, NCOLS3D, NROWS3D, Layers, NTHIK3D))  THEN
         MESG = 'CANMET has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( CANMET ) ) THEN
         CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
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

!...  Check start date, start time, end date, end time in CANMET
      IDATE = SDATE; ITIME = STIME
      IF ( .NOT. CHECK3( CANMET, 'SunleafTK', IDATE, ITIME ) ) THEN
         MESG = 'Starting time not on met file'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      CALL NEXTIME ( IDATE, ITIME, RLENG-10000 )
      IF ( .NOT. CHECK3( CANMET, 'SunleafTK', IDATE, ITIME ) ) THEN
         MESG = 'Ending time not on met file'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IDATE = SDATE; ITIME = STIME

!...  Set output parameters that are different from met file and open file
      SDATE3D = SDATE                ! From run-script
      STIME3D = STIME                ! From run-script
      MXREC3D = RLENG / 10000
      MXREC = MXREC3D
      NLAYS3D = 1
      NVARS3D = 20

      DO S = 1, NEMIS
         VNAME3D(S) = TRIM( MGN_SPC( S ) )
         VDESC3D(s) = 'Environmental activity factor for '//  TRIM( MGN_SPC(S) )
         UNITS3D(s) = 'Non-Dimension '
         VTYPE3D(s) = M3REAL
!         print*,'VNAME=', vname3d(s),VDESC3d(s),UNITS3d(s)
      ENDDO

      IF ( .NOT. OPEN3( MGNERS, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!----------------------------------------------------------------
!.....2) Calculate emission activity 
!----------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( NON_DIMGARMA (NEMIS, NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'NON_DIMGARMA', PROGNAME )
      ALLOCATE ( ER ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'ER', PROGNAME )
      ALLOCATE ( LAIc    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAIc',        PROGNAME )
      ALLOCATE ( LAIp    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAIp',        PROGNAME )
      ALLOCATE ( GAMLA   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMLA',       PROGNAME )
      ALLOCATE ( GAMSM   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMSM',       PROGNAME )
      ALLOCATE ( LDFMAP   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LDFMAP',       PROGNAME )
      ALLOCATE ( AQI   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'AQI',       PROGNAME )
      ALLOCATE ( GAMAQ   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMAQ',       PROGNAME )
      ALLOCATE ( GAMTP   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMTP',       PROGNAME )
      ALLOCATE ( GAMBD   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMBD',       PROGNAME )
      ALLOCATE ( GAMHT   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMHT',       PROGNAME )
      ALLOCATE ( GAMLT   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMLT',       PROGNAME )
      ALLOCATE ( GAMHW   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMHW',       PROGNAME )
      ALLOCATE ( GAMCO2  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMCO2',      PROGNAME )
      ALLOCATE ( SunT  ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunT',   PROGNAME )
      ALLOCATE ( ShaT  ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'ShaT', PROGNAME )
      ALLOCATE ( SunP  ( NCOLS, NROWS, Layers ),     STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunP',     PROGNAME )
      ALLOCATE ( ShaP  ( NCOLS, NROWS, Layers ),    STAT = IOS )
      CALL CHECKMEM    ( IOS, 'ShaP',   PROGNAME )
      ALLOCATE ( SunF  ( NCOLS, NROWS, Layers ),      STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunF',  PROGNAME )
      ALLOCATE ( CDEA  ( NCOLS, NROWS, Layers ),      STAT = IOS )
      CALL CHECKMEM    ( IOS, 'CDEA',  PROGNAME )
      ALLOCATE ( MaxT  ( NCOLS, NROWS ),   STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MaxT',     PROGNAME )
      ALLOCATE ( MinT  ( NCOLS, NROWS ),   STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MinT',     PROGNAME )
      ALLOCATE ( MaxWS  ( NCOLS, NROWS ),  STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MaxWS',    PROGNAME )
      !ALLOCATE ( Max_temp ( N_MaxT, NCOLS, NROWS ),   STAT = IOS )
      !CALL CHECKMEM    ( IOS, 'Max_temp',     PROGNAME )
      !ALLOCATE ( Min_temp ( N_MinT, NCOLS, NROWS ),   STAT = IOS )
      !CALL CHECKMEM    ( IOS, 'Min_temp',     PROGNAME )
      !ALLOCATE ( Max_wind ( N_MaxWS, NCOLS, NROWS ),   STAT = IOS )
      CALL CHECKMEM    ( IOS, 'Max_wind',     PROGNAME )
      ALLOCATE ( D_TEMP  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'D_TEMP',   PROGNAME )
      ALLOCATE ( D_PPFD  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'D_PPFD',   PROGNAME )
      ALLOCATE ( VPGWT  ( Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'VPGWT',   PROGNAME )
      ALLOCATE ( Ea1L  ( Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'Ea1L',   PROGNAME )
      ALLOCATE ( Ea2L  ( Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'Ea2L',   PROGNAME )

!...  Read DailyMET

      IDATE = SDATE; ITIME = STIME
      IF ( .NOT. READ3(DailyMET,'D_TEMP', 1,IDATE,ITIME,D_TEMP)) THEN
          MESG = 'Error reading D_TEMP'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. READ3(DailyMET,'D_PPFD', 1,IDATE,ITIME,D_PPFD)) THEN
          MESG = 'Error reading D_PPFD'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
      ENDIF

!...    Read MaxT
!     IF ( GAMHT_YN ) THEN
!        DO T = 1, N_MaxT
!        IDATE = SDATE; ITIME = STIME
!        IF ( .NOT. READ3(DailyMET,'MaxT', 1,IDATE,ITIME,         & 
!                                       Max_temp(T,:,:))) THEN
!            MESG = 'Error reading MaxT'
!            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
!        ENDIF
!        CALL NEXTIME(IDATE,ITIME,-240000*T)
!        ENDDO
!        MaxT = MAXVAL(Max_temp,1)
!      ENDIF
  
!...    Read MinT
!     IF ( GAMLT_YN ) THEN
!        DO T = 1, N_MinT
!        IDATE = SDATE; ITIME = STIME
!        IF ( .NOT. READ3(DailyMET,'MinT', 1,IDATE,ITIME,              &
!                                               Min_temp(T,:,:))) THEN
!            MESG = 'Error reading MinT'
!            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
!        ENDIF
!        CALL NEXTIME(IDATE,ITIME,-240000*T)
!        ENDDO
!        MinT = MINVAL(Min_temp,1)
!     ENDIF

!...    Read MaxWS
!     IF ( GAMHW_YN ) THEN
!        DO T = 1, N_MaxWS
!        IDATE = SDATE; ITIME = STIME
!        IF ( .NOT. READ3(DailyMET,'MaxWS', 1,IDATE,ITIME,           &
!                                          Max_wind(T,:,:))) THEN
!            MESG = 'Error reading MaxWS'
!            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
!        ENDIF
!        CALL NEXTIME(IDATE,ITIME,-240000*T)
!        ENDDO
!        MaxWS = MAXVAL(Max_wind,1)
!      ENDIF
  

! read air quality index
!     IF ( GAMAQ_YN ) THEN
!        IF ( .NOT. READ3(AQFILE,'W126', ALLAYS3,0,0,AQI)) THEN
!            MESG = 'Error reading AQFILE'
!            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
!        ENDIF
!     ENDIF

!...  Start the loop over the time period
      IDATE = SDATE
      ITIME = STIME
      DO T = 1, MXREC

        WRITE(MESG,*) 'Processing at :', IDATE, ITIME
        CALL M3MESG( MESG )

        ! Find LAIc from date
        CALL FINDLAI(IDATE,MXLAI,LAIp_I,LAIc_I)
        WRITE(MESG,*) 'Found LAI current period for YYYYJJJ : ', IDATE,LAIc_I
        CALL M3MESG( MESG )

        WRITE(MESG, '(I0.2)') LAIc_I
        IF ( .NOT. READ3(LAIS46,'LAI'//TRIM(MESG),ALLAYS3,0,0,LAIc)) THEN
          MESG = 'Error reading LAI at current time step'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        WRITE(MESG,'(I0.2)')  LAIp_I
        IF ( .NOT. READ3(LAIS46,'LAI'//TRIM(MESG),ALLAYS3,0,0,LAIp)) THEN
          MESG = 'Error reading LAI at previous time step'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! soil moisture activity factor

!       IF ( GAMSM_YN ) THEN
!          IF ( .NOT. READ3(SMFILE,'GAMSM', 1,IDATE,ITIME,GAMSM)) THEN
!              MESG = 'Error reading GAMSM'
!              CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
!          ENDIF
!       ELSE
           GAMSM = 1.0
!       ENDIF

! read CANMET

        IF ( .NOT. READ3(CANMET,'SunleafTK', ALLAYS3,IDATE,ITIME, SunT)) THEN
            MESG = 'Error reading SunleafTK'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(CANMET,'ShadeleafTK', ALLAYS3,IDATE,ITIME, ShaT)) THEN
            MESG = 'Error reading ShadeleafTK'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(CANMET,'SunPPFD', ALLAYS3,IDATE,ITIME, SunP)) THEN
            MESG = 'Error reading SunPPFD'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(CANMET,'ShadePPFD', ALLAYS3,IDATE,ITIME, ShaP)) THEN
            MESG = 'Error reading ShadePPFD'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(CANMET,'SunFrac', ALLAYS3,IDATE,ITIME,  SunF)) THEN
            MESG = 'Error reading SunFrac'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

!       ! Go over all the emission classes
        DO S = 1, NEMIS


! Read LDF 
        IF ( S .EQ. 3 .OR. S .EQ. 4 .OR. S .EQ. 5 .OR. S .EQ. 6 .OR.      &
                        S .EQ. 9 .OR. S .EQ. 10 ) THEN
! Read LDF from LDFILE for MT_PINE, MT_ACYC, MT_CAMP, MT_SABI, SQT_HR,SQT_LR
          WRITE(MESG,*) 'A', 'Reading LDF from LDFILE for '//VNAME3D(S)
          CALL M3MESG(MESG)
          WRITE(MESG,'(A)') VNAME3D(S)
          WRITE(MESG,'(I0.2)') S
          IF ( .NOT. READ3(LDFILE,'LDF'//TRIM(MESG),                   &
                                          ALLAYS3,0,0,LDFMAP)) THEN
            MESG = 'Error reading LDF for'//VNAME3D(S)
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF
        ELSE
! Read LDF from MEGVEA.EXT
            LDFMAP = LDF(S)
        ENDIF

! leaf age activity factor
          CALL GAMMA_A( NCOLS, NROWS, S,                              &
                         LAIp, LAIc, D_TEMP, GAMLA )

! Emission response to canopy depth
          CALL GAMMA_CD( NCOLS, NROWS, Layers, LAIc, CDEA )

! EA bidirectional exchange LAI response
!         IF ( GAMBD_YN ) THEN
!            CALL GAMMA_LAIbidir(NCOLS, NROWS, LAIc, GAMBD)
!         ELSE
             GAMBD = 1.0
!         ENDIF

! emission activity response to air quality
!         IF ( GAMAQ_YN ) THEN
!            CALL GAMMA_AQ(NCOLS, NROWS, S, AQI, GAMAQ)
!         ELSE
             GAMAQ = 1.0
!         ENDIF

! EA response to high temperature
!         IF ( GAMHT_YN ) THEN
!            CALL GAMMA_HT(NCOLS, NROWS, S, MaxT, GAMHT)
!         ELSE
             GAMHT = 1.0
!         ENDIF

! EA response to low temperature
!         IF ( GAMLT_YN ) THEN
!            CALL GAMMA_LT(NCOLS, NROWS, S, MinT, GAMLT)
!         ELSE
             GAMLT = 1.0
!         ENDIF

! EA response to high wind speed
!         IF ( GAMHW_YN ) THEN 
!            CALL GAMMA_HW(NCOLS, NROWS, S, MaxWS, GAMHW)
!         ELSE
             GAMHW = 1.0
!         ENDIF

! EA response to CO2 only for isoprene
          IF ( S .EQ. 1 ) THEN
!           IF ( GAMCO2_YN ) THEN
!             CALL GAMMA_CO2(NCOLS, NROWS, GAMCO2)
!           ELSE
                 GAMCO2 = 1.0
!           ENDIF
          ENDIF

! EA response to canopy temperature/light

        IF ( Layers .EQ. 5 ) THEN
          VPGWT(1) = 0.1184635
          VPGWT(2) = 0.2393144
          VPGWT(3) = 0.284444444
          VPGWT(4) = 0.2393144
          VPGWT(5) = 0.1184635
        ELSE
          DO K = 1,Layers
          VPGWT(K) =1/Layers
          ENDDO
        ENDIF

          DO I = 1, NCOLS
            DO J = 1, NROWS

              DO K = 1, Layers
                Ea1L(K)  = CDEA(I,J,K) *                             &
                     GAMTLD(SunT(I,J,K),D_TEMP(I,J),S) *             &
                     GAMP(SunP(I,J,K),D_PPFD(I,J)) *  SunF(I,J,K) +  &
                     GAMTLD(ShaT(I,J,K),D_TEMP(I,J),S) *             &
                     GAMP(ShaP(I,J,K),D_PPFD(I,J))                   &
                    * (1-SunF(I,J,K))

                Ea2L(K) = GAMTLI(SunT(I,J,K),S)* SunF(I,J,K)+        &
                   GAMTLI(ShaT(I,J,K),S)*(1-SunF(I,J,K))

              ENDDO   ! ENDDO canopy layers
            GAMTP(I,J)=SUM((Ea1L(:)*LDFMAP(I,J) +                    &
                           Ea2L(:)*(1-LDFMAP(I,J)))* VPGWT(:))
            ENDDO   ! NROWS
          ENDDO ! NCOLS

!          write(MESG,*) 'GAMSM(100,100)',GAMSM(100,100)
!          CALL M3MESG( MESG )
 
! ... Calculate emission activity factor
!       GAMCO2 only applied to isoprene
          IF ( S .EQ. 1 ) THEN
          ER(:,:) = LAIc * GAMTP * GAMCO2 * GAMLA *                 &
              GAMHW * GAMAQ * GAMHT * GAMLT * GAMSM * LDFMAP
!       GAMBD only applied to ethanol and acetaldehyde
          ELSE IF ( S .EQ. 13 ) THEN
          ER(:,:) = LAIc * GAMTP * GAMBD * GAMLA *                  &
              GAMHW * GAMAQ * GAMHT * GAMLT * GAMSM * LDFMAP
!       For CO (S=19), LDFMAP not applied
          ELSE IF ( S .EQ. 19 ) THEN
          ER(:,:) = LAIc * GAMTP * GAMLA *                          &
              GAMHW * GAMAQ * GAMHT * GAMLT * GAMSM
          ELSE
          ER(:,:) = LAIc * GAMTP * GAMLA *                          &
              GAMHW * GAMAQ * GAMHT * GAMLT * GAMSM * LDFMAP
          ENDIF
          WHERE(ER(:,:).GT.0.0)
           NON_DIMGARMA (S,:,:) = ER(:,:)
          ELSEWHERE
           NON_DIMGARMA (S,:,:) = 0.0
          ENDWHERE

        ENDDO  ! End loop of species (S)

!----------------------------------------------------------------
!.....3) Write out the calculated EA 
!----------------------------------------------------------------
!... Write emission to file
        WRITE(MESG,1030) 'Writing emission at ',T,IDATE,ITIME
        CALL M3MESG( MESG )
        DO S = 1,NEMIS
          IF (.NOT. WRITE3(MGNERS, VNAME3D(S),IDATE,ITIME,         &
                   NON_DIMGARMA (S,:,:))) THEN
            CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
            MESG = 'Error writing to file: '//TRIM(MESG)
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

        ENDDO ! End loop for emission species (S)

       CALL NEXTIME( IDATE, ITIME, TSTEP )
      ENDDO ! End loop for time step (T)

!... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

      DEALLOCATE ( ER )
      DEALLOCATE ( NON_DIMGARMA )
      DEALLOCATE ( GAMLA )
      DEALLOCATE ( GAMSM )
      DEALLOCATE ( GAMTP )
      DEALLOCATE ( AQI )
      DEALLOCATE ( GAMAQ )
      DEALLOCATE ( GAMBD )
      DEALLOCATE ( GAMHT )
      DEALLOCATE ( GAMLT )
      DEALLOCATE ( LDFMAP )
      DEALLOCATE ( GAMHW )
      DEALLOCATE ( GAMCO2 )
      DEALLOCATE ( LAIc )
      DEALLOCATE ( LAIp )
      DEALLOCATE ( SunT )
      DEALLOCATE ( ShaT )
      DEALLOCATE ( SunP )
      DEALLOCATE ( ShaP )
      DEALLOCATE ( SunF )
      DEALLOCATE ( CDEA )
      DEALLOCATE ( MaxT )
      DEALLOCATE ( MinT )
      DEALLOCATE ( MaxWS )
      DEALLOCATE ( Max_temp )
      DEALLOCATE ( Max_wind )
      DEALLOCATE ( Min_temp )
      DEALLOCATE ( D_TEMP )
      DEALLOCATE ( D_PPFD )
!
! ... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

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

      END PROGRAM TEST_MEGVEA
