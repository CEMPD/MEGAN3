!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in
!  continuous development by various groups and is based on information
!  from these groups: Federal Government employees, contractors working
!  within a United States Government contract, and non-Federal sources
!  including research institutions.  These groups give the Government
!  permission to use, prepare derivative works of, and distribute copies
!  of their work in the CMAQ system to the public and to permit others
!  to do so.  The United States Environmental Protection Agency
!  therefore grants similar permission to use the CMAQ system software,
!  but users are requested to provide copies of derivative works or
!  products designed to operate in the CMAQ system to the United States
!  Government without restrictions as to use by others.  Software
!  that is used with the CMAQ system but distributed under the GNU
!  General Public License or the GNU Lesser General Public License is
!  subject to their copyright restrictions.
!------------------------------------------------------------------------!
C Adopted from CMAQ files and modified for MEGAN3.1 by Ling Huang
C 2019/07/15

      MODULE BDSNP_MOD
      
      IMPLICIT NONE
      REAL,    ALLOCATABLE, SAVE :: NDEPRATE  ( :,: )  ! "ng N m-2 s-1" deposition rate used to update NDEPRES
      LOGICAL, SAVE     :: BDSNP
      LOGICAL, SAVE     :: EPIC
      
      CONTAINS
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HRNOBDSNP( IDATE, ITIME, TSTEP, NC, NR, COSZEN,
     &            TASFC, SSOLAR, PRES, SOILM, SOILT, RSTYP,LAI,
     &                     SNOCOV, CFRAC, WSPD10, GAREA,BDSNP_NO )

  
      USE CANOPY_NOX_MOD
     
      IMPLICIT NONE

C Includes:
      INCLUDE 'PARMS3.EXT'   ! I/O API parameters
      INCLUDE 'IODECL3.EXT'  ! I/O API function declarations
      INCLUDE 'FDESC3.EXT'   ! I/O API file desc. data structures

C Arguments:
      INTEGER, INTENT( IN )  :: IDATE             ! current simulation date (YYYYDDD)
      INTEGER, INTENT( IN )  :: ITIME             ! current simulation time (HHMMSS)
      INTEGER, INTENT( IN )  :: TSTEP( 3 )        ! time step vector (HHMMSS)
      INTEGER, INTENT( IN ) :: NC      ! no. columns
      INTEGER, INTENT( IN ) :: NR      ! no. rows
      REAL, INTENT( IN ) :: GAREA   ! area of grid cell (m2)
C     These are arrays 
      REAL,    INTENT( IN ) :: COSZEN( NC,NR )        ! cosine of zenith angle
      REAL,    INTENT( IN ) :: TASFC    ( NC,NR )        ! surface air temperature [K]
      REAL,    INTENT( IN ) :: SSOLAR( NC,NR )        ! surface radiation [w/m**2]
      REAL,    INTENT( IN ) :: PRES  ( NC,NR )        ! surface pressure [Pa]
      REAL,    INTENT( IN ) :: SOILM ( NC,NR )         ! soil moisture [m3/m3] (PX)
      REAL,    INTENT( IN ) :: SOILT  ( NC,NR )        ! soil temperature [K] (PX)
      REAL,    INTENT( IN ) :: RSTYP  ( NC,NR )        ! soil type
      REAL,    INTENT( IN ) :: LAI  ( NC,NR )        ! leaf area index (m2/m2)
      REAL,    INTENT( IN ) :: SNOCOV  ( NC,NR )        ! snow cover
      REAL,    INTENT( IN ) :: CFRAC  ( NC,NR )        ! cloud fraction
      REAL,    INTENT( IN ) :: WSPD10  ( NC,NR )        ! wind speed at 10m
      REAL,    INTENT( OUT ) :: BDSNP_NO ( NC,NR ) ! output NO emissions in nanomol/m^2/s

C External Functions
      INTEGER, EXTERNAL       ::   ENVINT
      INTEGER, EXTERNAL       ::   ENVYN
      INTEGER, EXTERNAL       ::   TIME2SEC
      INTEGER, EXTERNAL       ::   SECSDIFF

C Parameters:
C      INTEGER, PARAMETER :: MXRHRS = 24     ! no. of rainfall hours for YL95 algorithm
      INTEGER, PARAMETER :: MAXSTYPES = 11
!                 IF ( SOILCAT .LE. 14 ) THEN  <<< NEW PX METCRO2D ??? 14 = water
C      REAL,    PARAMETER :: CFNODRYFC = ( 1.0 / 3.0 ) * ( 1.0 / 30.0 )
      
      ! Scale factor so that fertilizer emission = 1.8 Tg N/yr (Stehfest and Bouwman, 2006)
      ! before canopy reduction
      REAL*8, PARAMETER :: FERT_SCALE = 0.0068 
      ! Value calculated by running the 2x2.5 GEOS-Chem model
      
      REAL*8,  PARAMETER :: TAU_MONTHS   = 6. ! this is the decay time for dep. N reservoir, fert is 4 months
      REAL*8,  PARAMETER :: SECPERDAY    = 86400.d0
      REAL*8,  PARAMETER :: DAYSPERMONTH = 30.
      REAL*8,  PARAMETER :: TAU_SEC      = TAU_MONTHS * DAYSPERMONTH * SECPERDAY
      
      ! Conversion factor from [ng N/m2/s] to [molec/cm2/s]
      REAL*8,  PARAMETER :: UNITCONV = 4.3d9
      
      ! New soil biomes based on Steinkamp et al., 2011
      INTEGER, PARAMETER :: NSOIL    = 24

      ! Canopy wind extinction coefficients
      ! (cf. Yienger & Levy [1995], Sec 5), now a function of the MODIS/KOPPEN biometype (J.D. Maasakkers)
       REAL*8,  PARAMETER :: SOILEXC(NSOIL)    = (/ 
     &  0.10, 0.50, 0.10, 0.10, 0.10,
     &  0.10, 0.10, 0.10, 0.10, 1.00,
     &  1.00, 1.00, 1.00, 2.00, 4.00,
     &  4.00, 4.00, 4.00, 4.00, 4.00,
     &  4.00, 2.00, 0.10, 2.00                  /)

      ! Steinkamp and Lawrence, 2011 A values, wet biome coefficients
      ! for each of the 24 soil biomes [ng N/m2/s] 
      REAL*8,  PARAMETER :: A_BIOME(NSOIL)          =                (/ 
     &   0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.09, 0.09, 0.01, 
     &   0.84, 0.84, 0.24, 0.42, 0.62, 0.03, 0.36, 0.36, 0.35, 
     &   1.66, 0.08, 0.44, 0.57, 0.57, 0.57     /)
        
C Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
C Pleim-Xiu Land-Surface and PBL Model (PX-LSM)
C See Jacquemin B. and Noilhan J. (1990), Bound.-Layer Meteorol., 52, 93-134.

      REAL :: SATURATION( MAXSTYPES ) =
     &                  (/ 0.395, 0.410, 0.435, 0.485,
     &                     0.451, 0.420, 0.477, 0.476,
     &                     0.426, 0.482, 0.482 /)

C Local Variables:

      CHARACTER( 16 ), SAVE :: SOILINSTATE = 'SOILINSTATE' ! logical name for input NO soil data, restart file
      CHARACTER( 16 ), SAVE :: SOILOUT = 'SOILOUT' ! logical name for output NO soil data - same format as soilinstate
C Land use files for BDSNP: both time independant in CMAQ sense and absolutely - e.g. fertilizer does not vary with year
      CHARACTER( 16 ), SAVE :: CLIMAFILE = 'CLIMAFILE' ! climate, arid 
      CHARACTER( 16 ), SAVE :: CLIMNAFILE = 'CLIMNAFILE' ! climate, nonarid 
      CHARACTER( 16 ), SAVE :: LANDTYPEFILE = 'LANDTYPEFILE' ! biome type 
      CHARACTER( 16 ), SAVE :: FERTRESFILE = 'FERTRESFILE' ! fertilizer reservoir
      CHARACTER( 16 ), SAVE :: EPICRESFILE = 'EPICRESFILE' ! Epic file
      CHARACTER( 16 ), SAVE :: NDEPFILE = 'NDEPFILE' ! Nitrogen deposition file (monthly, generated from CESM2/CAM-chem)

      CHARACTER( 16 ) :: VAR        ! variable name

      REAL,    ALLOCATABLE, SAVE :: SOILMPREV( :,: )  ! soil moisture previous tstep [m3/m3] (PX)
      REAL,    ALLOCATABLE, SAVE :: FERT     ( :,: )  ! "ng N m-2" already - reservoir
      REAL,    ALLOCATABLE, SAVE :: T1_EPIC  ( :,: )  ! "Epic input ng N m-2"
      REAL,    ALLOCATABLE, SAVE :: EPICN    ( :,: )  ! "epic N reservoir ng N m-2"
      REAL,    ALLOCATABLE, SAVE :: NDEPRES  ( :,: )  ! "ng N m-2" deposition reservoir
      REAL,    ALLOCATABLE, SAVE :: DRYPERIOD( :,: )  ! Will determine pulse, hours
      REAL,    ALLOCATABLE, SAVE :: PFACTOR  ( :,: )  ! current pulse factor
C Gridded Canopy NOx reduction factor for BDSNP Soil NO calculations
      REAL,    ALLOCATABLE, SAVE :: CRF   ( :,: )     ! 0-1

C --- diagnostic variables, can be removed in final version      
      REAL,    ALLOCATABLE, SAVE :: CRFAVG   ( :,: )  ! 0-1
      REAL,    ALLOCATABLE, SAVE :: PULSEAVG   ( :,: )  ! 1+
      REAL,    ALLOCATABLE, SAVE :: BASESUM   ( :,: )  ! used in calculating the above two averages
      REAL,    ALLOCATABLE, SAVE :: THETA_DIAG( :,: )  ! diagnositc theta
      REAL,    ALLOCATABLE, SAVE :: WET_DIAG ( :,: )  ! diagnositc wet term
      REAL,    ALLOCATABLE, SAVE :: TEMP_DIAG ( :,: )  ! diagnositc temp term
      REAL,    ALLOCATABLE, SAVE :: A_DIAG ( :,: )  ! diagnositc biome base emissions term
      REAL,    ALLOCATABLE, SAVE :: AFERT_DIAG ( :,: )  ! diagnositc fert emissions term
      REAL,    ALLOCATABLE, SAVE :: NRES_FERT_DIAG ( :,: )  ! diagnositc nres fert
      REAL,    ALLOCATABLE, SAVE :: NRES_DEP_DIAG ( :,: )  ! diagnositc nres  dep    
C ---------------------------------------------------------------------------      
      INTEGER, ALLOCATABLE, Save :: ARID     ( :,: )
      INTEGER, ALLOCATABLE, Save :: NONARID  ( :,: )
      INTEGER, ALLOCATABLE, Save :: LANDTYPE ( :,: )
      REAL,                 SAVE :: EMPOL,EMPOLSUM, EMPOLAVG       ! use to check reasonableness of results, g/hr
      REAL,                 SAVE :: TIMECHECK ! use to output CPU_TIME(TIMECHECK) to see if this section of code is running unreasonably long
      
      INTEGER, SAVE :: COUNTER ! save number of subroutine calls for calculating average CRF or others

      
      INTEGER, SAVE :: EDATE     ! end scenario date
      INTEGER, SAVE :: ETIME     ! end scenario time
      INTEGER, SAVE :: NDATE     ! test date to update rainfall
      INTEGER, SAVE :: NTIME     ! test time to update rainfall
      INTEGER, SAVE :: SDATE     ! scenario start date
      INTEGER, SAVE :: STIME     ! scenario start time
        
      LOGICAL, SAVE :: INITIAL_RUN            ! false: use SOILINSTATE restart file
                                              ! initial_run is set in run.cctm script
C      LOGICAL, SAVE :: INITIAL_DAY = .FALSE.  ! true: 1st 24 hours; no previous data
                                              ! false: previous 24 hours of rainfall
                                              ! are available for HRNO

      LOGICAL, SAVE :: PX_VERSION         ! true: use PX version of MCIP; should always be true

      INTEGER          SOILCAT            ! soil category
      INTEGER          NSTEPS             ! run duration (HHMMSS)
      INTEGER          MONTH              ! month of the year
      INTEGER, SAVE :: MSTEPS             ! run no. of steps
      INTEGER          I, J, K, R, C, L   ! counters
      INTEGER          JDAY, YEAR
      CHARACTER*3      CHARDAY
      CHARACTER*2      CHARMON
      LOGICAL          OK
      INTEGER          IOS                ! IO or memory allocation status
      
C      REAL,    SAVE :: EFAC
       REAL            TEMP_TERM
       REAL            WET_TERM
       REAL            PULSE
       REAL            A_FERT
       REAL            CRF_TERM
       REAL            SOILNOX, FERTDIAG
C      REAL             CFNO               ! NO correction factor
C      REAL             CFNOGRASS          ! NO correction factor for grasslands
C      REAL             TAIR               ! surface temperature
C      REAL             TSOI               ! soil temperature
       REAL             THETA              ! water filled pore space
       REAL             THETAPREV

C      REAL             CFNOWET, CFNODRY, THETA
C      REAL             FAC1, FAC2, FAC3, FAC4

      LOGICAL, SAVE :: USE_SOILT = .TRUE. ! use soil temperature in PX version
                                          ! rather than estimate as in BEIS2

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      CHARACTER( 256 ) :: MESG            ! message buffer
      CHARACTER( 16 )  :: PNAME = 'BDSNPHRNO'  ! procedure name

C-----------------------------------------------------------------------
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         BDSNP = .TRUE.
         ! we need to initialize and allocate:
         ! pulse
         ! length of dry period
         ! soil moisture of previous time step
         ! N reservoir, deposition only
         ! These values can be provided from a restart file. The restart file is 'timeless'.
         ! This means CMAQ isn't checking to see if the restart file is actually from
         ! the immediately prior timstep.

         WRITE( MESG,*) 'FIRSTIME BDSNP',NC,NR
         CALL M3MESG( MESG )
C Determine last timestamp
         SDATE = ENVINT( 'SDATE', 'Scenario Start (YYYYJJJ)', 0,IOS )
         STIME = ENVINT( 'STIME', 'Scenario Start (HHMMSS)', 0, IOS)
         NSTEPS = ENVINT( 'RLENG', 'Scenario Duration (HHMMSS)', 0,IOS )
         MONTH = ENVINT( 'MONTH', 'month', 0,IOS )
         EDATE = SDATE; ETIME = STIME
         CALL NEXTIME( EDATE, ETIME, NSTEPS-10000 )   ! end date & time
         MSTEPS = TIME2SEC( NSTEPS ) / TIME2SEC( TSTEP( 1 ) )

C Check if using PX version of MCIP
         PX_VERSION = ENVYN( 'PX_VERSION', 'MCIP is PX version?',
     &                       .TRUE., IOS )

C     ! make sure it only runs with PX version
         IF( .NOT. PX_VERSION ) THEN
            MESG = "BDSNP Soil NO is only compatible with PX version"
            CALL M3EXIT( PNAME, IDATE, ITIME, MESG, 2)
         END IF

C Allocate memory for data and read
         Write(MESG,*) 'about to allocate memory, NC,NR', NC, NR
         CALL M3MESG( MESG )

         ALLOCATE( SOILMPREV( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'SOILMPREV', PNAME )

         ALLOCATE( DRYPERIOD( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'DRYPERIOD', PNAME )

         ALLOCATE( NDEPRES( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'NDEPRES', PNAME )

         ALLOCATE( NDEPRATE( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'NDEPRATE', PNAME )

         ALLOCATE( PFACTOR( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'PFACTOR', PNAME )

         ALLOCATE( ARID( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'ARID', PNAME )

         ALLOCATE( NONARID( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'NONARID', PNAME )

         ALLOCATE( LANDTYPE( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'LANDTYPE', PNAME )

         ALLOCATE( FERT( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'FERT', PNAME )

         ALLOCATE( T1_EPIC( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'T1_EPIC', PNAME )

         ALLOCATE( EPICN( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'EPICN', PNAME )

         ALLOCATE( CRF( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CRF', PNAME )

         ALLOCATE( CRFAVG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CRFAVG', PNAME )

         ALLOCATE( PULSEAVG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'PULSEAVG', PNAME )

         ALLOCATE( BASESUM( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'BASESUM', PNAME )

C ------ Diagnostics -----------------------------------
         ALLOCATE( THETA_DIAG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'THETA_DIAG', PNAME )

         ALLOCATE( WET_DIAG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'WET_DIAG', PNAME )

         ALLOCATE( TEMP_DIAG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'TEMP_DIAG', PNAME )

         ALLOCATE( A_DIAG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'A_DIAG(', PNAME )

         ALLOCATE( AFERT_DIAG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'AFERT_DIAG', PNAME )

         ALLOCATE( NRES_FERT_DIAG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'NRES_FERT_DIAG', PNAME )

         ALLOCATE( NRES_DEP_DIAG( NC,NR ), STAT=IOS )
         CALL CHECKMEM( IOS, 'NRES_DEP_DIAG', PNAME )

C-----------------------------------------------------------------------------
C Initial run if the model hasn't been run before, otherwise use a restart file
C to determine DRYPERIOD, pulse state, prev. timestep soil moisture, and N reservoir.
         INITIAL_RUN = ENVYN( 'INITIAL_RUN', MESG, .FALSE., IOS )
         WRITE( MESG,* ) 'BDSNP Initial run?', INITIAL_RUN
         CALL M3MESG( MESG )
        
C If initial run, initialize some variables, otherwise get them from file
         IF ( INITIAL_RUN ) THEN

            PFACTOR   = 1d0   ! array
            DRYPERIOD = 0d0   ! array

            SOILMPREV = 0d0   ! array
            FERT      = 0d0   ! array
            NDEPRES   = 0d0   ! array
            EMPOL  = 0d0
            EMPOLSUM  = 0d0
            EMPOLAVG  = 0d0
            BASESUM = 0.0
            CRFAVG = 0.0
            PULSEAVG = 0.0
            !NDEPRES   = 0d0   ! array
            !attempt to use steady state condition to reduce spin up time by setting dN/dt = 0
            !or NDEPRES = Dep rate * tau, the decay time
!
C open nitrogen deposition file
        IF ( .NOT. OPEN3( NDEPFILE, FSREAD3, PNAME ) ) THEN
            CALL NAMEVAL (NDEPFILE, MESG)  ! get input file name and path
            MESG = 'Could not open file '//TRIM(MESG)
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        ENDIF

        WRITE(CHARMON,'(i0.2)') MONTH
        VAR =  'NITROGEN'// CHARMON
        MESG = 'Reading nitrogen deposition for '//VAR
        CALL M3MESG(MESG)

        IF ( .NOT. READ3(NDEPFILE,VAR, ALLAYS3,0,0,
     &                            NDEPRATE )) THEN
          MESG = 'Error reading NDEPRATE from NDEPFILE'
          CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
        ENDIF

            DO R = 1, NR
                DO C = 1, NC
                NDEPRES( C,R ) = NDEPRATE( C,R )*TAU_SEC
C           check for negatives
          IF( NDEPRES(C,R) .lt. 0.0 ) THEN
            Write(MESG,*) 'NDEPRES negative', NDEPRES, ' ',NDEPRATE(C,R)
            CALL M3EXIT( PNAME, IDATE, ITIME, MESG, 2 )
          ELSE IF( NDEPRATE(C,R) .lt. 0.0 ) THEN
            Write(MESG,*) 'NDEPRATE negative', NDEPRES,' ',NDEPRATE(C,R)
            CALL M3EXIT( PNAME, IDATE, ITIME, MESG, 2 )
          END IF
              END DO
           END DO

        ELSE ! SOILINSTATE file available

C Open saved NO soil state data file
         IF ( .NOT. OPEN3( SOILINSTATE, FSREAD3, PNAME ) ) THEN
            CALL NAMEVAL (SOILINSTATE, MESG)  ! get input file name and path
            MESG = 'Could not open file '//TRIM(MESG)
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
         ENDIF

C Get description of NO soil state data file
         IF ( .NOT. DESC3( SOILINSTATE ) ) THEN
            CALL NAMEVAL (SOILINSTATE, MESG)  ! get input file name and path
            MESG = 'Could not get description of SOILINSTATE'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
         ENDIF

C Check that the file start date and time are consistent
            IF ( SDATE3D .NE. SDATE ) THEN
               MESG = 'Starting date is wrong'
               CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF

            IF ( STIME3D .NE. STIME ) THEN
               MESG = 'Starting time is wrong'
               CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF

C Check grid
         IF ( .NOT. FILCHK3 ( SOILINSTATE,
     &                 GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
            MESG = 'SOILINSTATE has differenet grid definition'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
         ENDIF

C Read data from file
         IF ( .NOT. READ3(SOILINSTATE,'DRYPERIOD', ALLAYS3,0,0,
     &                            DRYPERIOD )) THEN
           MESG = 'Error reading DRYPERIOD from SOILINSTATE'
           CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
         ENDIF

         IF ( .NOT. READ3(SOILINSTATE,'PFACTOR', ALLAYS3,0,0,
     &                            PFACTOR )) THEN
           MESG = 'Error reading PFACTOR from SOILINSTATE'
           CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
         ENDIF

         IF ( .NOT. READ3(SOILINSTATE,'SOILMPREV', ALLAYS3,0,0,
     &                            SOILMPREV )) THEN
           MESG = 'Error reading SOILMPREV from SOILINSTATE'
           CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
         ENDIF

         IF ( .NOT. READ3(SOILINSTATE,'NDEPRES', ALLAYS3,0,0,
     &                            NDEPRES )) THEN
           MESG = 'Error reading NDEPRES from SOILINSTATE'
           CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
         ENDIF

         IF ( .NOT. READ3(SOILINSTATE,'NDEPRATE_DIAG', ALLAYS3,
     &                  0,0,NDEPRATE )) THEN
           MESG = 'Error reading NDEPRATE_DIAG from SOILINSTATE'
           CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
         ENDIF

C Close input file
         IF ( .NOT. CLOSE3( SOILINSTATE ) ) THEN
            MESG = 'Could not close file "'//TRIM( SOILINSTATE )//'"'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2)
         END IF

        END IF  ! initial run check

        Write(MESG,*) 'about to read climate file'
        CALL M3MESG(MESG)

C open CLIMAFILE
        IF ( .NOT. OPEN3( CLIMAFILE, FSREAD3, PNAME ) ) THEN
            CALL NAMEVAL (CLIMAFILE, MESG)  ! get input file name and path
            MESG = 'Could not open file '//TRIM(MESG)
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        ENDIF

       IF ( .NOT. READ3(CLIMAFILE,'ARID', ALLAYS3,0,0,
     &                            ARID )) THEN
          MESG = 'Error reading ARID from CLIMAFILE'
          CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
       ENDIF

       IF ( .NOT. OPEN3( CLIMNAFILE, FSREAD3, PNAME ) ) THEN
            CALL NAMEVAL (CLIMNAFILE, MESG)  ! get input file name and path
            MESG = 'Could not open file '//TRIM(MESG)
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        ENDIF

        IF ( .NOT. READ3(CLIMNAFILE,'NONARID', ALLAYS3,0,0,
     &                            NONARID )) THEN
           MESG = 'Error reading NONARID from CLIMNAFILE'
           CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
        ENDIF

C open LANDTYPEFILE
        Write(MESG,*) 'about to read LANDTYPE file'
        CALL M3MESG(MESG)
        IF ( .NOT. OPEN3( LANDTYPEFILE, FSREAD3, PNAME ) ) THEN
            CALL NAMEVAL (LANDTYPE, MESG)  ! get input file name and path
            MESG = 'Could not open file '//TRIM(MESG)
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        ENDIF
         IF ( .NOT. READ3(LANDTYPEFILE,'LANDTYPE', ALLAYS3,0,0,
     &                            LANDTYPE )) THEN
           MESG = 'Error reading LANDTYPE from LANDTYPEFILE'
           CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
         ENDIF
      WRITE(MESG,*) 'finished IF FIRSTIME BDSNP'
      CALL M3MESG(MESG)
      END IF ! FIRSTIME

      IF ( SECSDIFF( IDATE,ITIME,EDATE,ETIME ) .LE. TIME2SEC(TSTEP(2))) 
     &            GO TO 9999

C Choose EPIC or Potter for the fertilizer data
      EPIC = ENVYN( 'EPIC', MESG, .FALSE., IOS )
      WRITE(MESG,*) 'Use EPIC outputs for fertilizer (+ deposition) 
     &                                                data?',EPIC
      CALL M3MESG(MESG)
   
      IF ( EPIC ) THEN
C Read Time EPIC specific data
       Write(MESG,*) 'about to read EPIC data:'
       CALL M3MESG(MESG)

C     read EPIC fertilizer reservoir
       IF ( .NOT. OPEN3( EPICRESFILE, FSREAD3, PNAME ) ) THEN
         CALL NAMEVAL (EPICRESFILE, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PNAME, IDATE, ITIME, MESG, 2 )
       ENDIF
C get the day's fertilizer from EPIC
       IF ( .NOT. READ3(EPICRESFILE,'T1_EPIC',ALLAYS3,IDATE,0,
     &                            T1_EPIC )) THEN
           MESG = 'Error reading '// 'T1_EPIC' //' from EPICRESFILE'
           CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
        ENDIF

        DO R = 1, NR
         DO C = 1, NC
             EPICN( C,R ) =  T1_EPIC( C,R )

          END DO
        END DO

      ELSE ! EPIC logic variable is FALSE
C    read day dependant fertilizer reservoir e.g. Potter et al 2010
        IF ( .NOT. OPEN3( FERTRESFILE, FSREAD3, PNAME ) ) THEN
            CALL NAMEVAL (FERTRESFILE, MESG)  ! get input file name and path
            MESG = 'Could not open file '//TRIM(MESG)
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        ENDIF
C get the day for the file
        YEAR = INT( FLOAT( IDATE ) / 1000.0 )
        JDAY = IDATE - YEAR * 1000
        WRITE(CHARDAY,'(i3.3)') JDAY
        VAR =  'FERT'// CHARDAY
        MESG = 'Reading fertilizer for '//VAR
        CALL M3MESG(MESG)
C get the days' fertilizer
        IF ( .NOT. READ3(FERTRESFILE,VAR, ALLAYS3,0,0,
     &                            FERT )) THEN
           MESG = 'Error reading '// VAR //' from FERTRESFILE'
           CALL M3EXIT(PNAME,0,0,MESG,2)
        ENDIF

      ENDIF ! EPIC or not

C Fertilizer N reservoir already calculated and read from file, update deposition reservoir from dep rate
      DO R = 1, NR
        DO C = 1, NC
            CALL GET_NDEPRES( TSTEP, NDEPRES( C,R ), NDEPRATE(C,R),
     &                                            TAU_SEC, C, R)
        END DO
      END DO
      WRITE(MESG,*) 'READ ALL DATA BDSNP'
      CALL M3MESG(MESG)

C Calculate temporal non-speciated soil NO emissions to EMPOL

C attempt to remove artefact by setting inputs to constant value
C      LANDTYPE = 1
C     If False Don't do any calculations to test - replicate 0 output
      IF( .TRUE. ) THEN
      CALL GET_CANOPY_NOX(IDATE, ITIME, NC, NR, COSZEN,
     & TASFC, SSOLAR, PRES, LANDTYPE, LAI, SNOCOV, CFRAC, WSPD10, CRF)
      DO R = 1, NR
         DO C = 1, NC

            SOILNOX  = 0d0
            FERTDIAG = 0d0

      ! ------Loop Over MODIS/Koppen  LANDTYPE
C            DO K = 1, 24

C            IF ( LANDTYPE( C,R ) .EQ. K ) THEN !Skip LANDTYPEs not pre
            K = LANDTYPE( C,R ) !Skip LANDTYPE not present
            ! Temperature-dependent term of soil NOx emissions
            ! [unitless]
            ! Uses PX soil temperature instead of inferring from air
            ! temperature
            TEMP_TERM = SOILTEMP( SOILT(C,R) )

            ! Use THETA instead of boolean wet/dry climate
            SOILCAT = INT( RSTYP( C,R ) )
            IF ( SOILCAT .LE. MAXSTYPES ) THEN !not water
               THETA = SOILM( C,R ) / SATURATION( SOILCAT )
               THETAPREV = SOILMPREV( C,R ) / SATURATION( SOILCAT )
               ! Soil moisture scaling of soil NOx emissions
               WET_TERM = SOILWET( THETA , ARID( C,R ), NONARID( C,R ))
            ELSE
               WET_TERM = 0d0
               THETA = 0d0
            END IF
            ! Cumulative multiplication factor (over baseline emissions)
            ! that accounts for soil pulsing
C            PFACTOR( C,R ) = PULSING( THETA, TSTEP, THETAPREV,
C     &                       PFACTOR( C,R ), DRYPERIOD( C,R ) )
            PULSE = PULSING( THETA, TSTEP, THETAPREV,
     &                       PFACTOR( C,R ), DRYPERIOD( C,R ) )

      IF ( EPIC ) THEN
C        WRITE(LOGDEV, *) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
C        WRITE(LOGDEV, *) 'Using fertilizer files from EPIC outputs...'
           A_FERT = EPICRESERVOIRN( EPICN( C,R ) ) !adds epic N reservoir (deposition considered)
C            WRITE( LOGDEV,*) 'FERTADD complete ',A_FERT
      ELSE
            A_FERT = FERTADD( FERT( C,R ) , NDEPRES( C,R ) ) !adds reservoirs returns emission rates
      END IF

            ! Canopy reduction factor
            CRF_TERM  = CRF( C,R )
C                  CRF_TERM  = SOILCRF( K, LAI,
C     &                        R_CANOPY(K),
C     &                               WINDSQR, SUNCOS )

         !  SOILNOX includes fertilizer
            SOILNOX   = ( A_BIOME(K) + A_FERT )  !don't forget to check parenthesis when uncommenting
     &             * ( TEMP_TERM * WET_TERM * PULSE )
     &             * ( 1.d0 - CRF_TERM  )

         !  FERTDIAG, only used for the fertilizer diagnostic, note
         !  includes DEP
         ! not actually used for anything at the moment, only
         ! diagnostics
            FERTDIAG  = ( A_FERT )
     &             * ( TEMP_TERM * WET_TERM * PULSE )
     &             * ( 1.d0 - CRF_TERM  )

C            END IF !LANDTYPE check

C            ENDDO

           !scale emissions

           !diagnostics
           EMPOL = SOILNOX * GAREA * 3600.0 * 10.0**-9![ng N/m2/s] * m2 * s/hr * g/ng
           BDSNP_NO(C,R) = SOILNOX / 14 ![nmol/m2/s]
           ! sum various quantities for daily averaging
           EMPOLSUM = EMPOLSUM + EMPOL
           BASESUM(C,R) = BASESUM(C,R) + ( A_BIOME(K) + A_FERT ) !don'tforget check paren when uncommenting
     &             * ( TEMP_TERM * WET_TERM)
     &             * GAREA * 3600.0 * 10.0**-9 ![ng N/m2/s] * m2 * s/hr * g/ng
           PULSEAVG(C,R) = PULSEAVG(C,R) + ( A_BIOME(K) + A_FERT ) !don'tforget check paren when uncommenting
     &             * ( TEMP_TERM * WET_TERM * PULSE )
     &             * GAREA * 3600.0 * 10.0**-9 ![ng N/m2/s] * m2 * s/hr * g/ng
           CRFAVG(C,R) = CRFAVG(C,R) + ( A_BIOME(K) + A_FERT ) !don'tforget check paren when uncommenting
     &             * ( TEMP_TERM * WET_TERM )
     &             * ( 1.d0 - CRF_TERM  )
     &             * GAREA * 3600.0 * 10.0**-9 ![ng N/m2/s] * m2 * s/hr * g/ng

C--------- MORE DIAGNOSTICS  ---------------------------------
           A_DIAG( C,R ) = A_BIOME(K)
           AFERT_DIAG( C,R ) = A_FERT
           NRES_FERT_DIAG( C,R ) = FERT( C,R )
           NRES_DEP_DIAG( C,R )  = NDEPRES( C,R )
           WET_DIAG( C,R ) = WET_TERM
           THETA_DIAG( C,R ) = THETA
           TEMP_DIAG( C,R ) = TEMP_TERM

C -----------------------------------------------------
        END DO ! columns
      END DO ! rows

      ELSE ! add things until it dies
      WRITE( MESG,*) 'BDSNP testing terms'
      CALL M3MESG(MESG)
         DO R = 1, NR
         DO C = 1, NC
            K = LANDTYPE( C,R ) !Skip LANDTYPE not present
            ! Temperature-dependent term of soil NOx emissions
            ! [unitless]
            ! Uses PX soil temperature instead of inferring from air
            ! temperature
            TEMP_TERM = SOILTEMP( SOILT(C,R) )
            ! Use THETA instead of boolean wet/dry climate
            SOILCAT = INT( RSTYP( C,R ) )
            IF ( SOILCAT .LE. MAXSTYPES ) THEN !not water
               THETA = SOILM( C,R ) / SATURATION( SOILCAT )
               THETAPREV = SOILMPREV( C,R ) / SATURATION( SOILCAT )
               ! Soil moisture scaling of soil NOx emissions
               WET_TERM = SOILWET( THETA , ARID( C,R ), NONARID( C,R ))
            ELSE
               WET_TERM = 0d0
               THETA = 0d0
            END IF
            ! Cumulative multiplication factor (over baseline emissions)
            ! that accounts for soil pulsing
C            PFACTOR( C,R ) = PULSING( THETA, TSTEP, THETAPREV,
C     &                       PFACTOR( C,R ), DRYPERIOD( C,R ) )
            PULSE = PULSING( THETA, TSTEP, THETAPREV,
     &                       PFACTOR( C,R ), DRYPERIOD( C,R ) )

         END DO
         END DO
      END IF ! end do nothing test if

      SOILMPREV = SOILM !save soilM array to soilMprev for next time step
      MESG = 'BDSNP calculated emissions'
      CALL M3MESG(MESG)
      CALL CPU_TIME(TIMECHECK)
      WRITE(MESG,*)  'PROCESS TOOK:', TIMECHECK, 'SECONDS'
      CALL M3MESG(MESG)
      EMPOLAVG = EMPOLSUM/FLOAT(NC*NR)
      WRITE( MESG,*) 'average value:', EMPOLAVG
      CALL M3MESG(MESG)
      EMPOLSUM = 0d0 !array
      COUNTER = COUNTER + 1
      RETURN

9999  CONTINUE

C Create soil NO state save file at the end of the run for restart purposes

C Final timestamp
      NDATE = EDATE; NTIME = ETIME

C Build description for, and create/open soil NO emissions output file
      FTYPE3D = GRDDED3
      SDATE3D = NDATE
      STIME3D = NTIME
      TSTEP3D = 0   ! make it a time-independent file
      NCOLS3D = NC
      NROWS3D = NR
      NLAYS3D = 1
      NVARS3D = 14
      MXREC3D = 1

      VNAME3D = ' '
      VNAME3D( 1 ) = 'PFACTOR'
      VNAME3D( 2 ) = 'DRYPERIOD'
      VNAME3D( 3 ) = 'NDEPRES'
      VNAME3D( 4 ) = 'SOILMPREV'
C --- DIAGNOSTICS ---------------------------
      VNAME3D( 5 ) = 'THETA_DIAG'
      VNAME3D( 6 ) = 'WET_TERM_DIAG'
      VNAME3D( 7 ) = 'TEMP_DIAG'
      VNAME3D( 8 ) = 'TEMP_TERM_DIAG'
      VNAME3D( 9 ) = 'A_DIAG'
      VNAME3D( 10 ) = 'NRES_FERT_DIAG'
      VNAME3D( 11 ) = 'AFERT_DIAG'
      VNAME3D( 12 ) = 'NDEPRATE_DIAG'
      VNAME3D( 13 ) = 'CRFAVG'
      VNAME3D( 14 ) = 'PULSEAVG'
c -------------------------------------------
      UNITS3D = ' '
      UNITS3D( 1 ) = 'REAL'
      UNITS3D( 2 ) = 'REAL'
      UNITS3D( 3 ) = 'REAL'
      UNITS3D( 4 ) = 'REAL'
      UNITS3D( 5 ) = 'REAL'
      UNITS3D( 6 ) = 'REAL'
      UNITS3D( 7 ) = 'REAL'
      UNITS3D( 8 ) = 'REAL'
      UNITS3D( 9 ) = 'REAL'
      UNITS3D( 10 ) = 'REAL'
      UNITS3D( 11 ) = 'REAL'
      UNITS3D( 12 ) = 'REAL'
      UNITS3D( 13 ) = 'REAL'
      UNITS3D( 14 ) = 'REAL'


      VDESC3D( 1 ) = 'NO emission current pulse factor'
      VDESC3D( 2 ) = 'length of the dry period in hours'
      VDESC3D( 3 ) = 'soil N reservoir from deposition'
      VDESC3D( 4 ) = 'Soil moisture prev. timestep m3/m3'
      VDESC3D( 5 ) = 'moisture WFPS 0-1'
      VDESC3D( 6 ) = 'moisture scale factore diagnostic'
      VDESC3D( 7 ) = 'temperature diagnostic'
      VDESC3D( 8 ) = 'temperature scale factor diagnostic'
      VDESC3D( 9 ) = 'biome base emission diagnostic'
      VDESC3D( 10 ) = 'NRES fert only diagnostic'
      VDESC3D( 11 ) = 'fertilizer emission factor diagnostic'
      VDESC3D( 12 ) = 'N deposition rate diagnostic'
      VDESC3D( 13 ) = 'canopy reduction factor diagnostic'
      VDESC3D( 14 ) = 'pulse factor diagnostic'

      VTYPE3D = 0
      VTYPE3D( 1:14 ) = M3REAL

      FDESC3D = ' '
      FDESC3D( 1 ) = 'Gridded soil state data for soil NO emissions'
      FDESC3D( 2 ) = '/From/ ' // PNAME
      FDESC3D( 3 ) = '/Version/ MEGAN3.1'

C Open NO soil state save file
      IF ( .NOT. OPEN3( SOILOUT, FSCREA3, PNAME ) ) THEN
         CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PNAME, IDATE, ITIME, MESG, 2 )
      ENDIF

      IF ( .NOT. WRITE3(SOILOUT,'PFACTOR',NDATE,NTIME,PFACTOR )) THEN
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. WRITE3(SOILOUT,'DRYPERIOD',NDATE,NTIME,DRYPERIOD )) THEN
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. WRITE3(SOILOUT,'NDEPRES',NDATE,NTIME,NDEPRES )) THEN
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. WRITE3(SOILOUT,'SOILMPREV',NDATE,NTIME,SOILM )) THEN
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

C ---- DIAGNOSTICS -------------------------------------------------------------------
      IF ( .NOT. WRITE3(SOILOUT,'THETA_DIAG',NDATE,NTIME,THETA_DIAG )) THEN ! diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing THETA_DIAG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. WRITE3(SOILOUT,'WET_TERM_DIAG',NDATE,NTIME,WET_DIAG )) THEN ! diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing WET_TERM_DIAG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. WRITE3(SOILOUT,'TEMP_DIAG',NDATE,NTIME,SOILT)) THEN ! diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing TEMP_DIAG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT.WRITE3(SOILOUT,'TEMP_TERM_DIAG',NDATE,NTIME,TEMP_DIAG)) THEN ! diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing TEMP_TERM_DIAG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT.WRITE3(SOILOUT,'A_DIAG',NDATE,NTIME,A_DIAG)) THEN ! diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing A_DIAG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT.WRITE3(SOILOUT,'AFERT_DIAG',NDATE,NTIME,AFERT_DIAG)) THEN ! diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing AFERT_DIAG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT.WRITE3(SOILOUT,'NRES_FERT_DIAG',NDATE,NTIME,
     &                                  NRES_FERT_DIAG)) THEN ! diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing NRES_FERT_DIAG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT.WRITE3(SOILOUT,'NDEPRATE_DIAG',NDATE,NTIME,
     &                                  NDEPRATE)) THEN !diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing NDEPRATE_DIAG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT.WRITE3(SOILOUT,'CRFAVG',NDATE,NTIME,
     &                                  CRFAVG/BASESUM)) THEN !diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing CRFAVG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT.WRITE3(SOILOUT,'PULSEAVG',NDATE,NTIME,
     &                                  PULSEAVG/BASESUM)) THEN !diagnostic theta
        CALL NAMEVAL (SOILOUT, MESG)  ! get input file name and path
        MESG = 'Error writing PULSEAVG to file: '//TRIM(MESG)
        CALL M3EXIT(PNAME,IDATE,ITIME,MESG,2)
      ENDIF

C close SOILOUT
         IF ( .NOT. CLOSE3( SOILOUT ) ) THEN
            MESG = 'Could not close file "'//TRIM( SOILOUT )//'"'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2)
         END IF

      COUNTER = 0
C ------------------------------------------------------------------------------------
      WRITE( MESG,* )
     &      'Timestep written to', SOILOUT,
     &      'for date and time', NDATE, NTIME

      RETURN

94040 FORMAT( 5X, 3( A, :, 1X ), I8, ":", I6.6)

      END SUBROUTINE HRNOBDSNP
C----------------------------------------------------------------------------------

         REAL FUNCTION PULSING( THETA, TSTEP, THETAPREV,
     &                       PFACTOR, DRYPERIOD )!_____
! !DESCRIPTION: Function PULSING calculates the increase (or "pulse") of
!  soil NOx emission that happens after preciptiation falls on dry soil.
!\\
!\\
!  According to  Yan et al., [2005] , this pulsing process is thought to
!  be due to a release of inorganic nitrogen trapped on top of the dry soil
!  and a subsequent reactivation of water-stressed bacteria, which then
!  metabolize the excess nitrogen. This can happen in seasonally dry
!  grasslands and savannahs or over freshly fertilized fields.
!  Soil NOx emissions consist of baseline emissions plus discrete "pulsing"
!  episodes.  We follow the Yan et al., [2005] algorithm, where the pulse
!  (relative to the flux pre wetting) is determined by the antecedent dry
!  period, with a simple logarithmic relationship,
!
!  PFACTOR = 13.01 ln ( DRYPERIOD ) -  53.6
!
!  ,where PFACTOR is the magnitude of peak flux relative to prewetting flux,
!  and DRYPERIOD  is the length of the antecedent dry period in hours.
!
!  The pulse decays with
!
!  PFACTOR = PFACTOR * EXP( -0.068d0 * TSTEP(HOURS) )

         IMPLICIT NONE

         INTEGER, EXTERNAL       ::   TIME2SEC

C Function arguments:
         INTEGER, INTENT( IN )    :: TSTEP( 3 )        ! time step vector (HHMMSS)
         REAL,    INTENT( IN )    :: THETA, THETAPREV  ! only avilable if PX version
         REAL,    INTENT( INOUT ) :: DRYPERIOD
         REAL,    INTENT( INOUT ) :: PFACTOR
C Local Variables
         REAL MOISTDIFF
         REAL DTHOURS
         DTHOURS = TIME2SEC(TSTEP(2))/3600.0
         ! If soil moisture less than 0.3 and no pulse is taking place
         IF ( THETA < 0.3D0 .and. PFACTOR == 1.D0) THEN

            ! Get change in soil moisture since previous timestep
            MOISTDIFF = ( THETA - THETAPREV )

            ! If change in soil moisture is > 0.01 (rains)
            IF ( MOISTDIFF > 0.01 ) THEN

               !Initialize new pulse factor (dry period hours)
               PFACTOR = 13.01 * LOG( DRYPERIOD ) - 53.6

               ! If dry period < ~3 days then no pulse
               IF ( PFACTOR < 1.0 ) PFACTOR = 1.0

                  ! Reinitialize dry period
                  DRYPERIOD = 0

                ! If no rain (i.e.,  change in soil moisture is < 0.01)
               ELSE
                ! Add one timestep to dry period
                DRYPERIOD = DRYPERIOD + DTHOURS

            ENDIF

         ! If box is already pulsing , then decay pulse one timestep
         ELSEIF ( PFACTOR /= 1.d0) THEN

            ! Decay pulse
            PFACTOR   = PFACTOR * EXP( -0.068d0 * DTHOURS )

            ! Update dry period
            IF ( THETA < 0.3D0 ) DRYPERIOD = DRYPERIOD + DTHOURS

            ! If end of pulse
            IF ( PFACTOR < 1.d0 ) PFACTOR = 1.d0

         ENDIF
         PULSING = PFACTOR
         RETURN

         END FUNCTION PULSING!_____

C---------------------------------------------------------------------------------------------
      
         SUBROUTINE GET_NDEPRES( TSTEP, NDEPRES, NDEPRATE, TAU_SEC,C,R )
C Get the deposition rate of the appropriate species for the appropriate timestep, add to reservoir and decay.
C Return reservoir amount.

         IMPLICIT NONE

         INTEGER, EXTERNAL       ::   TIME2SEC

C Function arguments:

         INTEGER, INTENT( IN )  :: TSTEP( 3 )        ! time step vector (HHMMSS)
         INTEGER, INTENT( IN )  :: C
         INTEGER, INTENT( IN )  :: R
         REAL*8,  INTENT( IN )  :: TAU_SEC
         REAL,    INTENT( IN )  :: NDEPRATE
         REAL,    INTENT( INOUT ) :: NDEPRES

C Local Variables
         CHARACTER( 256 ) :: MESG            ! message buffer
         CHARACTER( 16 )  :: PNAME = 'GET_NDEPRES'  ! procedure name
         REAL*8  :: C1 ! a factor
         REAL*8  :: C2  ! another one
         REAL*8  :: TS_SEC ! time step in seconds
         real NDEPTEMP
C           check for negatives
          IF( NDEPRES < 0.0 ) THEN
          WRITE(MESG,*) 'NDEPRES negative'
          Write(*,*) 'In GET_NDEPRES:'
          Write(*,*) 'NDEPRES negative', NDEPRES,' '
          write(*,*) 'TS, TAU, C1, C2:', TS_SEC, TAU_SEC,C1,C2
          CALL M3EXIT( PNAME, 0, 0, MESG, 2)
          ELSE IF( NDEPRATE < 0.0 ) THEN
          MESG = 'NDEPRATE negative'
          Write(*,*) 'In GET_NDEPRES:'
          Write(*,*) 'NDEPRATE negative', NDEPRATE
          CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
          END IF

         ! takes the NDEPRATE and uses it to update NDEPRES before
         ! clearing it.

         !Do mass balance (see Intro to Atm Chem Chap. 3)
         !m(t) = m(0) * exp(-t/tau) + Source * tau * (1 - exp(-t/tau))
         TS_SEC = TIME2SEC(TSTEP(2))
         C1 = EXP( - TS_SEC / TAU_SEC)
         C2 = 1.d0 - C1
               NDEPTEMP = NDEPRES
               NDEPRES = NDEPRES*C1+NDEPRATE*TAU_SEC*C2
C           check for negatives
          IF( NDEPRES < 0.0 ) THEN
          MESG = 'negative'
          Write(*,*) 'In GET_NDEPRES:'
          Write(*,*) 'NDEPRES negative', NDEPRES
          write(*,*) 'TS, TAU, C1, C2:', TS_SEC, TAU_SEC,C1,C2
          CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
          END IF

         RETURN

         END SUBROUTINE GET_NDEPRES

C -----------------------------------------------------------------------------
         REAL FUNCTION SOILTEMP( SOILT )
C Calculate the soil temperature factor

         IMPLICIT NONE
C Function arguments:
         REAL, INTENT( IN )       :: SOILT !kelvin, soil temperature
C Local Variables
         REAl SOILTC !temperature in degrees celsius
         CHARACTER( 256 ) :: MESG            ! message buffer
         CHARACTER( 16 )  :: PNAME = 'SOILTEMP'  ! procedure name
         SOILTC = SOILT - 273.16

         IF ( SOILTC <= 0d0 ) THEN
         ! No soil emissions if temp below freezing
         SOILTEMP = 0d0
C         BENCHMARKING:
C         MESG = 'temperature less than 0 in august florida?'
C         CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT1 )

         ELSE

         ! Caps temperature response at 30C
         IF ( SOILTC >= 30.d0 ) SOILTC = 30.d0

         SOILTEMP =  EXP( 0.103 * SOILTC )

         ENDIF
         RETURN

         END FUNCTION SOILTEMP

C ---------------------------------------------------------------------------------------------------------
         REAL FUNCTION FERTADD( FERT , DEPN )
C Add fertilizer reservoir to deposition reservoir and create N driven
C emission factor
         IMPLICIT NONE
C Function arguments:
         REAL, INTENT( IN )       :: FERT !fertilizer reservoir [ngN/m2]
         REAL, INTENT( IN )       :: DEPN !deposition reservoir [ngN/m2]

C Local Variables
         REAL*8,  PARAMETER :: SECPERYEAR    = 86400.d0 * 365.
         ! Scale factor so that fertilizer emission = 1.8 Tg N/yr
         ! (Stehfest and Bouwman, 2006)
         ! before canopy reduction
         REAL*8, PARAMETER :: FERT_SCALE = 0.0068
         ! Value calculated by running the 2x2.5 GEOS-Chem model
         ! (J.D. Maasakkers)
         FERTADD = FERT + DEPN
         FERTADD = FERTADD / SECPERYEAR * FERT_SCALE

         RETURN

         END FUNCTION FERTADD

C -------------------------------------------------------------------------------------
         REAL FUNCTION EPICRESERVOIRN( EPICN )
C Add EPIC N reservoir (FERT+DEP) in lieu to FERTADD taking from Potter data

         IMPLICIT NONE
C Functions arguments:
         REAL, INTENT( IN )       :: EPICN !ConvertedTo'ngN/m2'from'kg/ha
C         REAL, INTENT( IN )       :: DEPN !deposition reservoir[ngN/m2]

C Local Variables
         REAL*8,  PARAMETER :: SECPERYEAR_EPIC    = 86400.d0 * 365.
         ! Scale factor so that fertilizer emission = 1.8 Tg N/yr
         ! (Stehfest and Bouwman, 2006)
         ! before canopy reduction
         REAL*8, PARAMETER :: FERT_SCALE_EPIC = 0.0068
         ! Value calculated by running the 2x2.5 GEOS-Chem model
         ! (J.D. Maasakkers)
         EPICRESERVOIRN = EPICN
         EPICRESERVOIRN = EPICRESERVOIRN / SECPERYEAR_EPIC * FERT_SCALE_EPIC


         RETURN
         END FUNCTION EPICRESERVOIRN
C -------------------------------------------------------------------------------------

C Local Variables

         REAL FUNCTION SOILWET( THETA , ARID, NONARID)
C Calculate the soil moisture factor

         IMPLICIT NONE
C Function arguments:
         REAL, INTENT( IN )       :: THETA !0-1 soil moisture
         INTEGER, INTENT( IN )    :: ARID !1 indicates arid cell
         INTEGER, INTENT( IN )    :: NONARID !1 indicates nonarid cell, if both 0 then
C Local Variables

         IF ( ARID .EQ. 1 ) THEN !ARID, Max poison at theta = .2
         SOILWET = 8.24*THETA*EXP(-12.5*THETA*THETA)
         ELSE IF (NONARID .EQ. 1 ) THEN !NONARID Max Poisson at theta =.3
         SOILWET = 5.5*THETA*EXP(-5.55*THETA*THETA)
         ELSE !neither arid nor nonarid, water or non-emitting cell
         SOILWET = 0.0
         END IF

         RETURN

         END FUNCTION SOILWET

C -------------------------------------------------------------------
      END MODULE BDSNP_MOD
      
