      PROGRAM MEGSEA

!***********************************************************************
!   This program computes soil NO emission activity factor and isoprene
!   soil moisture activity using MCIP output variables.
!
!  DESCRIPTION:
!
!     Uses new NO algorithm NO = Normalized*Tadj*Padj*Fadj*Cadj
!     to estimate NO emissions
!     Information needed to estimate NO emissions
!     Julian Day          (integer)    JDATE
!     Surface Temperature (MCIP field) TA    (K)
!     Soil Moisture       (MCIP field) SOILM (M**3/M**3) (LSOIL)
!          (ratio of volume of water per volume of soil)
!     Soil Temperature    (MCIP field) SOILT (K)         (LSOIL)
!     Soil Type           (MCIP field) ISLTYP            (LSOIL)
!
!     saturation values for soil types (constants)       (LSOIL)
!     FOR PX Version, the Temperature adjustment factor accounts for wet
!     and dry soils
!                and  the precipitation adjustment factor accounts for
!                saturated soils
!     FOR the non-PX version, the basic algorithm remains with a
!     temperature adjustment factor (dry soil)
!                     and no adjustment for saturated soils
!
!     The following arrays are updated after each call to SOILNOX
!     PULTYPE   type of NO emission pulse
!     PULSEDATE julian date for the beginning of an NO pulse
!     PULSETIME        time for the beginning of an NO pulse
!
!     The calculation are based on the following paper by J.J. Yienger
!     and H. Levy II
!     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol
!     100,11447-11464,1995
!
!     The Temperature Adjustment Factor is based on section 4.2 for wet
!     and dry soils with the following modification (PX version):
!       Instead of classifying soils as either 'wet' or 'dry', the wet
!       and dry adjustment is calculated at each grid cell.  A linear 
!       interpolation between the wet and dry adjustment factor is made 
!       using the relative amount of soil moisture in the top layer (1cm)
!       as the interpolating factor.  The relative amount of soil moisture 
!       is determined by taking the MCIP soil moisture field and dividing by the
!       saturation value defined for each soil type in the PX version of MCIP
!       the soil temperature is used in PX version
!
!     The Precipation Adjustment factor is based on section 4.1 with the
!     following modifications.
!       The rainrate is computed from the MCIP directly using a 24 hr daily total.
!       THe types of Pulses as described in YL95 were used to estimate
!       the NO emission rate.
!
!    Also see the following paper for more information:
!    Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!    Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!    by Tom Pierce and Lucille Bender
!
!    REFERENCES
!
!    JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
!    J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and 
!      Nitric Oxide Emissions from Agricultural Proc Proceedings of the Air and Waste 
!      Management Association/U.S. Environmental Protection Agency EMission Inventory 
!      Conference, Raleigh October 26-28, 1999 Raleigh NC
!  PRECONDITIONS REQUIRED:
!     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
!     NO emission pulse type, soil moisture from previous time step, julian date
!     of NO emission pulse start, time of NO emission pulse start, soil type, 
!     SOIL TYPES, Land use data
!
!  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
!     FERTILIZER_ADJ computes fertlizer adjustment factor
!     VEG_ADJ        computes vegatation adjustment factor
!     GROWSEASON     computes day of growing season
!
! HISTORY:
!   07/21/11: Imported from SMOKE-BEIS v3.14 for MEGEAN v2.10 (Tan)
!   03/19/17: Make as an indpendent program (MEGSEA) (Ling Huang)
!   03/31/17: Add calculation for soil moisture activity (Ling Huang)
!   06/10/19: Add an option to use BDSNP model to calculate soil NO
!             emissions (Ling Huang)
!*********************************************************************

      USE SOILNOX_FX
     
      USE BDSNP_MOD             ! BDSNP model to calcualte biogenic soil
                                ! microbe NO emissions

      IMPLICIT NONE

! INCLUDE FILES
      INCLUDE 'PARMS3.EXT'   ! I/O API parameters
      INCLUDE 'IODECL3.EXT'  ! I/O API function declarations
      INCLUDE 'FDESC3.EXT'   ! I/O API file desc. data structures
      INCLUDE 'MEGSEA.EXT'   ! wilting point info.

!...  EXTERNAL FUNCTIONS and their descriptions:
      INTEGER, EXTERNAL       ::   ENVINT
      INTEGER, EXTERNAL       ::   ENVYN
      LOGICAL      DSCGRID
      EXTERNAL     DSCGRID

!...  Program I/O files: From run script

! Program name
      CHARACTER*16  :: PROGNAME = 'MEGSEA'

! Netcdf file
      CHARACTER*16  :: MGNMET = 'MGNMET'       ! Meteorology file  
      CHARACTER*16  :: LAIS46 = 'LAIS46'       ! LAI file  
      CHARACTER*16  :: CANTYP = 'CANTYP'     ! canopy type file logical name

! output file
      CHARACTER*16  :: MGNSEA = 'MGNSEA'       ! Emission activity 

!...  Parameters for file units
      INTEGER  LOGDEV                          ! Logfile unit number

!...  External parameters
! From run script
      INTEGER       SDATE          ! Start date YYYYDDD
      INTEGER       STIME          ! Start time HHMMSS
      INTEGER       RLENG          ! Run length HHMMSS

! I/O API file parameters
      INTEGER :: JDATE        ! Date YYYYDDD from inpname
      INTEGER :: JTIME        ! Time HHMMSS from inpname
      INTEGER :: NCOLS        ! Number of columns
      INTEGER :: NROWS        ! Number of rows
      INTEGER :: NLAYS        ! Number of vertical layers
      INTEGER :: MXREC        ! Total number of timesteps
      INTEGER :: TSTEP        ! Time step
      INTEGER :: TSTEP3(3)    ! 3D time step used by CMAQ

!...  Internal parameters
! Internal parameters (status and buffer)
      INTEGER       IOS                    ! i/o status
      CHARACTER*256 MESG                   ! message buffer

! Local variables and their descriptions:
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

!     variable from LAI file
      REAL, ALLOCATABLE :: LAIc( :,: )    ! Current time step LAI
      REAL, ALLOCATABLE :: LAT ( :,: )    ! Latitude
      REAL, ALLOCATABLE :: LON ( :,: )    ! Longitude

!     variable from MGNMET
      REAL, ALLOCATABLE :: TEMP (:,:)   ! Temperautre (K)
      REAL, ALLOCATABLE :: PRES (:,:)   ! Surface pressure [Pa]
      REAL, ALLOCATABLE :: SNOCOV (:,:)   ! Snow cover; needed by BDSNP
      REAL, ALLOCATABLE :: CFRAC (:,:)   ! Cloud fraction; needed by BDSNP
      REAL, ALLOCATABLE :: WSPD10 (:,:)   ! wind speed; needed by BDSNP
      REAL, ALLOCATABLE :: SSOLAR (:,:)   ! Surface radiation [w/m**2]
      REAL, ALLOCATABLE :: PAR (:,:)   ! Photosynthetically Active Radiation [w/m**2]
      REAL, ALLOCATABLE :: COSZEN (:,:)   ! zenith angle
      REAL, ALLOCATABLE :: PRECADJ (:,:)   
      INTEGER, ALLOCATABLE :: SLTYP ( :,: ) ! soil type
      REAL, ALLOCATABLE :: SOILM ( :,: ) ! soil moisture
      REAL, ALLOCATABLE :: SOILT ( :,: ) ! soil temperature
 
      REAL, ALLOCATABLE    :: RSTYP( :,: )

      REAL, ALLOCATABLE :: CTF( :, :, : ) ! Canopy type factor arra

!     output variable
      REAL, ALLOCATABLE :: CFNO ( :,: )      ! Emission activity for crop
      REAL, ALLOCATABLE :: CFNOG ( :,: )     ! Emission activity for grass
      REAL, ALLOCATABLE :: GAMSM ( :,: )     ! Soil moisture activity for isoprene
      REAL, ALLOCATABLE :: GAMNO ( :,: )     ! Final NO emission activity
      REAL, ALLOCATABLE :: BDSNP_NO ( :,: )     ! BDSNP NO emissions(nmol/s/m2)

      INTEGER :: LAIp_I,LAIc_I
      INTEGER :: GDAY, GLEN
      INTEGER :: MXLAI,MXCT
      REAL :: t1,wilt,TMO1,TMO2
      REAL :: GAREA

      LOGICAL :: LSOIL = .TRUE.
      LOGICAL :: LBDSNP

! loop indices
      INTEGER :: IDATE, ITIME
      INTEGER :: T,I,J,I_CT

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

      CALL ENVSTR( 'GDNAM3D', MESG, 'ASACA36km', GDNAM, IOS )
      IF( .NOT. DSCGRID( GDNAM, CNAME, GDTYP3D,
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Open files


! Canopy type file
      WRITE(MESG,1030) 'Checking up canopy type file',0,0,0
      CALL M3MESG( MESG )
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
      WRITE(MESG,1030) 'Checking up LAI file',0,0,0
      CALL M3MESG( MESG )
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

      GAREA = XCELL3D * YCELL3D
      NCOLS = NCOLS3D
      NROWS = NROWS3D
      TSTEP = TSTEP3D
      TSTEP3(1) = TSTEP3D
      TSTEP3(2) = TSTEP3D
      TSTEP3(3) = TSTEP3D

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
      MXREC3D = RLENG / 10000
      MXREC = MXREC3D
      NLAYS3D = 1
      NVARS3D = 2

      LBDSNP = ENVYN( 'BDSNP_YN',MESG, .FALSE., IOS)
      IF ( LBDSNP ) THEN
      VNAME3D(1) = 'BDSNP_NO'
      UNITS3D(1) = 'nmol/m2/s'
      VTYPE3D(1) = M3REAL
      VDESC3D(1) = ' '
      ELSE
      VNAME3D(1) = 'GAMNO'
      UNITS3D(1) = ' '
      VTYPE3D(1) = M3REAL
      VDESC3D(1) = ' '
      ENDIF

      VNAME3D(2) = 'GAMSM'
      UNITS3D(2) = ' '
      VTYPE3D(2) = M3REAL
      VDESC3D(2) = ' '


      IF ( .NOT. OPEN3( MGNSEA, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNSEA, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!----------------------------------------------------------------
!.....2) Calculate emission activity 
!----------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( TEMP    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'TEMP',        PROGNAME )
      ALLOCATE ( PRES    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PRES',        PROGNAME )
      ALLOCATE ( PAR    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PAR',        PROGNAME )
      ALLOCATE ( SSOLAR    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SSOLAR',        PROGNAME )
      ALLOCATE ( COSZEN    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'COSZEN',        PROGNAME )
      ALLOCATE ( PRECADJ    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PRECADJ',        PROGNAME )
      ALLOCATE ( LAIc    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAIc',        PROGNAME )
      ALLOCATE ( LAT    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAT',        PROGNAME )
      ALLOCATE ( LON    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LON',        PROGNAME )
      ALLOCATE ( SOILT    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SOILT',  PROGNAME )
      ALLOCATE ( SOILM ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SOILM',  PROGNAME )
      ALLOCATE ( SLTYP ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SLTYP',  PROGNAME )
      ALLOCATE ( RSTYP ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'RSTYP',  PROGNAME )
      ALLOCATE ( SNOCOV ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SNOCOV',  PROGNAME )
      ALLOCATE ( CFRAC ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'CFRAC',  PROGNAME )
      ALLOCATE ( WSPD10 ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'WSPD10',  PROGNAME )
      ALLOCATE ( CFNO  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'CFNO',   PROGNAME )
      ALLOCATE ( CFNOG ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'CFNOG',  PROGNAME )
      ALLOCATE ( GAMNO ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMNO',  PROGNAME )
      ALLOCATE ( BDSNP_NO ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'BDSNP_NO',  PROGNAME )
      ALLOCATE ( GAMSM ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMSM',  PROGNAME )
      ALLOCATE ( CTF( NRTYP, NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'CTF', PROGNAME )

!...  Start the loop over the time period
      IDATE = SDATE
      ITIME = STIME
      DO T = 1, MXREC

        WRITE(MESG,*) 'Processing at :', IDATE, ITIME
        CALL M3MESG( MESG )
! ...  Initialize hourly variables
        TEMP = 0.
        PRES = 0.
        PAR = 0.
        SSOLAR = 0.
        COSZEN = 0.
        LAIc = 0.
        CFNO = 0.
        CFNOG = 0.
        GAMNO = 0.
        BDSNP_NO = 0.

        IF ( .NOT. READ3(MGNMET,'TEMP2',  ALLAYS3,IDATE,ITIME,TEMP)) THEN
          MESG = 'Error reading temperature'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PRES',  ALLAYS3,IDATE,ITIME,PRES)) THEN
          MESG = 'Error reading surface pressure'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PAR',  ALLAYS3,IDATE,ITIME,PAR)) THEN
          MESG = 'Error reading PAR'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PREC_ADJ',ALLAYS3,IDATE,ITIME,PRECADJ))
     &     THEN
          MESG = 'Error reading precipitation adjustment'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! if using met2mgn_rad45
        SSOLAR =  PAR/0.45 ! convert PAR to solar radiation
! if using default met2mgn
!        SSOLAR =  PAR/0.5 ! convert PAR to solar radiation

!...  Read CANTYP
      DO I_CT = 1, MXCT
         IF ( .NOT. READ3(CANTYP,'CTS',1,0,(I_CT-1)*10000,
     &                                               CTF(I_CT,:,:))) THEN
            MESG = 'Error reading CTS'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
         ENDIF
      ENDDO

! Read LAIS46 file 
! Find LAIc from date
        CALL FINDLAI(IDATE,MXLAI,LAIp_I,LAIc_I)
        WRITE(MESG,*) 'Found LAI current period for YYYYJJJ : ',
     &  IDATE,LAIc_I
        CALL M3MESG( MESG )

        WRITE(MESG,'(I0.2)') LAIc_I
        IF ( .NOT. READ3(LAIS46,'LAI'//TRIM(MESG),ALLAYS3,0,0,LAIc)) THEN
          MESG = 'Error reading LAI at current time step'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(LAIS46,'LAT',1,0,0,LAT)) THEN
          MESG = 'Error reading LAT'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(LAIS46,'LONG',1,0,0,LON)) THEN
          MESG = 'Error reading LON'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF


! ... Calculate emission activity factor for NO
! ... Two options in MEGAN3.1: YL (same model as MEGAN3) or BDSNP (new
!     in MEGAN3.1)
        IF ( LBDSNP ) THEN
           WRITE(MESG,1030) 'Use BDSNP to estimate soil NO emissions'
           CALL M3MESG( MESG )
!      Compute zenith angle
           CALL CZANGLE ( IDATE,ITIME,NCOLS,NROWS,LAT,LON,COSZEN )
! Read soil moisture data
        IF ( .NOT. READ3(MGNMET,'SOIM1',ALLAYS3,IDATE,ITIME,SOILM)) THEN
            MESG = 'Error reading SOIM1'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! Read soil temperature data
        IF ( .NOT. READ3(MGNMET,'SOIT1',ALLAYS3,IDATE,ITIME,SOILT)) THEN
            MESG = 'Error reading SOIT1'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! Read soil type 
        IF ( .NOT. READ3(MGNMET,'SLTYP',ALLAYS3,IDATE,ITIME,RSTYP)) THEN
            MESG = 'Error reading SLTYP'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! Read snow cover
        IF ( .NOT. READ3(MGNMET,'SNOCOV',ALLAYS3,IDATE,ITIME,SNOCOV)) THEN
            MESG = 'Error reading SNOCOV'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! Read cloud fraction
        IF ( .NOT. READ3(MGNMET,'CFRAC',ALLAYS3,IDATE,ITIME,CFRAC)) THEN
            MESG = 'Error reading CFRAC'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! Read wind speed at 10 m
        IF ( .NOT. READ3(MGNMET,'WSPD10',ALLAYS3,IDATE,ITIME,WSPD10)) THEN
            MESG = 'Error reading WSPD10'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

!      Call BDSNP
           CALL HRNOBDSNP( IDATE,ITIME,TSTEP3,NCOLS,NROWS,COSZEN,TEMP,
     &             SSOLAR,PRES,SOILM,SOILT,RSTYP,LAIc,SNOCOV,CFRAC,
     &                                          WSPD10,GAREA,BDSNP_NO)
        ELSE
           WRITE(MESG,1030) 'Use YL algorithm to estimate soil 
     &              NO emissions'
           CALL M3MESG( MESG )
           WRITE(MESG,1030) 'Estimating soil NOx adj: ',T,IDATE,ITIME
           CALL M3MESG( MESG )
           IF ( READ3(MGNMET,'SOIM1', ALLAYS3,IDATE,ITIME,SOILM) .AND.
     &          READ3(MGNMET,'SOIT1', ALLAYS3,IDATE,ITIME,SOILT) .AND.
     &          READ3(MGNMET,'SLTYP', ALLAYS3,IDATE,ITIME,RSTYP) ) THEN

             MESG = 'Using SOIL parameters in NOx adjustment'
             CALL M3MESG( MESG )
             LSOIL = .TRUE.
             SLTYP = INT(RSTYP)
           ELSE
             MESG = 'SOIL parameters are not available'
             CALL M3MESG( MESG )
             LSOIL = .FALSE.
           ENDIF
           CALL SOILNOX(IDATE,ITIME,NCOLS,NROWS,
     &                  TEMP,LSOIL,SLTYP, SOILM, SOILT,
     &                  LAIc, LAT, PRECADJ,
     &                  CFNO, CFNOG )

           DO I = 1,NCOLS
             DO J = 1,NROWS
               CALL GROWSEASON(IDATE,LAT(I,J),GDAY,GLEN)
               IF (GDAY .EQ. 0) THEN
                ! non growing season
                ! CFNOG for everywhere
                  GAMNO(I,J) = CFNOG(I,J)

                ELSE IF (GDAY .GT. 0 .AND. GDAY .LE. 366) THEN
                ! growing season
                ! CFNOG for everywhere except crops
                TMO1 = 0.
                TMO2 = 0.
                DO I_CT = 1,5
                  TMO1 = TMO1 + CTF(I_CT,I,J)
                  TMO2 = TMO2 + CTF(I_CT,I,J) * CFNOG(I,J)
                ENDDO
                ! CFNO for crops
                TMO1 = TMO1 + CTF(6,I,J)
                TMO2 = TMO2 + CTF(6,I,J) * CFNO(I,J)
                IF (TMO1 .EQ. 0.0) THEN
                   GAMNO(I,J) = 0.0
                ELSE
                   GAMNO(I,J) = TMO2 / TMO1
                ENDIF

              ELSE
                ! bad GDAY
                WRITE(MESG,*) 'Bad GDAY ',GDAY
                CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
              ENDIF

              ENDDO  !NROWS
           ENDDO  !NCOLS

        ENDIF ! YL or BDSNP

        WRITE(MESG,1030) 'Finished soil NOx adj: ',T,IDATE,ITIME
        CALL M3MESG( MESG )

! ... Calculate soil moisture activity factor for isoprene emissions
        IF ( READ3(MGNMET,'SOIM1', ALLAYS3,IDATE,ITIME,SOILM) .AND.
     &       READ3(MGNMET,'SOIT1', ALLAYS3,IDATE,ITIME,SOILT) .AND.
     &       READ3(MGNMET,'SLTYP', ALLAYS3,IDATE,ITIME,RSTYP) ) THEN
        WRITE(MESG,1030) 'Calculating soil moisture activity
     &                                    factor: ',T,IDATE,ITIME

        CALL M3MESG( MESG )
        SLTYP = INT(RSTYP)
           DO I = 1, NCOLS
             DO J = 1, NROWS
               wilt = WWLT(SLTYP(I,J))
               t1 = wilt + d1
               IF ( SOILM(I,J) < wilt ) THEN
                   GAMSM(I,J) = 0
               ELSE IF ( SOILM(I,J) >= wilt .AND. SOILM(I,J) < t1 ) THEN
                   GAMSM(I,J) = (SOILM(I,J) - wilt)/d1
               ELSE
                   GAMSM(I,J) = 1
               END IF
             END DO ! NROWS
           END DO ! NCOLS
         ENDIF

!----------------------------------------------------------------
!.....3) Write out the calculated EA 
!----------------------------------------------------------------
!... Write emission to file
        WRITE(MESG,1030) 'Writing emission at ',T,IDATE,ITIME
        CALL M3MESG( MESG )

! #1
      IF ( LBDSNP ) THEN
        vname3d(1) = 'BDSNP NO emission'
        vtype3d(1) =  m3real
        units3d(1) = 'nmol/m2/s'
        vdesc3d(1) = 'NO emission estimated from BDSNP algorithm'
        IF (.NOT. WRITE3(MGNSEA,'BDSNP_NO',IDATE,ITIME,
     &               BDSNP_NO(:,:))) THEN
          CALL NAMEVAL (MGNSEA, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
      ELSE
        vname3d(1) = 'Soil NO activity factor '
        vtype3d(1) =  m3real
        units3d(1) =  'No dimension            '
        vdesc3d(1) = 'Soil NO activity factor'
        IF ( .NOT. WRITE3(MGNSEA,'GAMNO',IDATE,ITIME,
     &               GAMNO(:,:))) THEN
          CALL NAMEVAL (MGNSEA, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
      ENDIF

! #2
        vname3d(2) = 'GAMSM '
        vtype3d(2) =  m3real
        units3d(2) =  'No dimension            '
        vdesc3d(2) = 'soil moisture activity factor'

        IF ( .NOT. WRITE3(MGNSEA,'GAMSM',IDATE,ITIME,
     &               GAMSM(:,:))) THEN
          CALL NAMEVAL (MGNSEA, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

       CALL NEXTIME( IDATE, ITIME, TSTEP )
      ENDDO ! End loop for time step (T)

!... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

      DEALLOCATE ( TEMP    )
      DEALLOCATE ( PRECADJ )
      DEALLOCATE ( LAIc    )
      DEALLOCATE ( LAT     )
      DEALLOCATE ( CFNO )
      DEALLOCATE ( CFNOG )
      DEALLOCATE ( SOILM )
      DEALLOCATE ( SOILT )
      DEALLOCATE ( SLTYP )
      DEALLOCATE ( RSTYP )
      DEALLOCATE ( GAMSM )
      DEALLOCATE ( GAMNO )
      DEALLOCATE ( BDSNP_NO )
      DEALLOCATE ( PAR )
      DEALLOCATE ( SSOLAR )
      DEALLOCATE ( PRES )
      DEALLOCATE ( SNOCOV )
      DEALLOCATE ( WSPD10 )
      DEALLOCATE ( CFRAC )
      DEALLOCATE ( LON )
      DEALLOCATE ( CTF )
!
! ... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

!--=====================================================================
!...  FORMAT
!--=====================================================================
1030  FORMAT (A30,I8,X,I8,X,I8)

!--=====================================================================
!...  End program
!--=====================================================================

      END PROGRAM MEGSEA

