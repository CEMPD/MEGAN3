
      SUBROUTINE SOILNOX( JDATE, JTIME, NX, NY,
     &                 TA, LSOIL, ISLTYP, SOILM, SOILT,
     &                 LAIc, LAT,
     &                 PRECADJ,
     &                 CFNO, CFNOG )

!***********************************************************************
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
!     FOR PX Version, the Temperature adjustment factor accounts for wet and dry soils
!                and  the precipitation adjustment factor accounts for saturated soils
!     FOR the non-PX version, the basic algorithm remains with a temperature adjustment factor (dry soil)
!                     and no adjustment for saturated soils
!
!
!     The following arrays are updated after each call to SOILNOX
!     PULTYPE   type of NO emission pulse 
!     PULSEDATE julian date for the beginning of an NO pulse 
!     PULSETIME        time for the beginning of an NO pulse
!  
!     The calculation are based on the following paper by J.J. Yienger and H. Levy II
!     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!
!     The Temperature Adjustment Factor is based on section 4.2 for wet and dry soils with
!       the following modification (PX version):
!       Instead of classifying soils as either 'wet' or 'dry', the wet and dry adjustment is 
!       calculated at each grid cell.  A linear interpolation between the wet and dry adjustment
!       factor is made using the relative amount of soil moisture in the top layer (1cm)
!       as the interpolating factor.  The relative amount of soil moisture is determined by
!       taking the MCIP soil moisture field and dividing by the saturation value defined for each
!       soil type in the PX version of MCIP
!       the soil temperature is used in PX version
!
!     The Precipation Adjustment factor is based on section 4.1 with the following modifications.
!       The rainrate is computed from the MCIP directly using a 24 hr daily total. 
!       THe types of Pulses as described in YL95 were used to estimate the NO emission
!       rate.  
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
!    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and Nitric Oxide Emissions from Agricultural Processes
!       Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!        Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!
C  PRECONDITIONS REQUIRED:
C     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
C     NO emission pulse type, soil moisture from previous time step, julian date
C     of NO emission pulse start, time of NO emission pulse start,
C     soil type, SOIL TYPES, Land use data
C
C  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
C     FERTILIZER_ADJ computes fertlizer adjustment factor
C     VEG_ADJ        computes vegatation adjustment factor
C     GROWSEASON     computes day of growing season
C     
C  REVISION  HISTORY:
C    10/01 : Prototype by GAP
C    10/03 : modified transition to non growing season for jul-oct of the year
C    08/04 : Converted to SMOKE code style by C. Seppanen
C    07/21/11 : Imported form SMOKE-BEIS v3.14 for MEGAN v2.10
C 
C***********************************************************************

        USE SOILNOX_FX

        IMPLICIT NONE
        
C.........  INCLUDES
        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C.........  ARGUMENTS and their descriptions
        INTEGER, INTENT (IN)  :: JDATE   !  current simulation date (YYYYDDD)
        INTEGER, INTENT (IN)  :: JTIME   !  current simulation time (HHMMSS)
        INTEGER, INTENT (IN)  :: NX      !  no. columns
        INTEGER, INTENT (IN)  :: NY      !  no. rows

        REAL, INTENT (IN)  ::  TA      ( NX, NY )    !  air temperature (K)
        REAL, INTENT (IN)  ::  SOILM   ( NX, NY )    !  soil moisture (m3/m3)
        REAL, INTENT (IN)  ::  SOILT   ( NX, NY )    !  soil temperature (K)
        REAL, INTENT (IN)  ::  PRECADJ ( NX, NY )    !  precip adjustment
        REAL, INTENT (IN)  ::  LAIc    ( NX, NY )    !  soil temperature (K)
        REAL, INTENT (IN)  ::  LAT     ( NX, NY )    !  Latitude
        REAL, INTENT (IN OUT)  ::  CFNO    ( NX, NY )    !  NO correction factor
        REAL, INTENT (IN OUT)  ::  CFNOG   ( NX, NY )    !  NO correction factor for grass
        
        INTEGER, INTENT (IN)  ::  ISLTYP  ( NX, NY )    !  soil type

        LOGICAL, INTENT (IN) :: LSOIL              ! true: using PX version of MCIP
        
C.........  Local ARRAYS
C Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
C       PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
C See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
        INTEGER, PARAMETER :: MAXSTYPES = 11
        REAL SATURATION( MAXSTYPES )
        DATA SATURATION / 0.395, 0.410, 0.435, 0.485,
     &                    0.451, 0.420, 0.477, 0.476,
     &                    0.426, 0.482, 0.482        /

C.........  SCRATCH LOCAL VARIABLES and their descriptions:
        INTEGER         R, C, L      ! counters
        INTEGER         SOILCAT      ! soil category
        
        REAL            CF           ! NO correction factor
        REAL            CFG          ! NO correction factor for grasslands
        REAL            TAIR         ! surface temperature
        REAL            TSOI         ! soil temperature
        REAL            CFNOWET, CFNODRY, RATIO


        CHARACTER(256)  MESG         ! message buffer
        
        CHARACTER(16) :: PROGNAME = 'SOILNOX'   !  program name

C***********************************************************************

 
C.....  Loop through cells
        DO R = 1, NY
        DO C = 1, NX

          TAIR = TA( C, R )         ! unit in degree K

C.......  Check max and min bounds for temperature
          IF (TAIR < 200.0) THEN
             WRITE( MESG, 94010 ) 'TAIR=', TAIR,
     &              'out of range at (C,R)=', C, R
             CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
          END IF

          IF (TAIR > 315.0 ) THEN
              WRITE( MESG, 94020 ) 'TAIR=', TAIR,
     &               'out of range at (C,R)=', C, R,
     &               ' resetting to 315K'
              CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
              TAIR = 315.0
          END IF

C.......  CFNOG
          IF( TAIR > 303.00 ) TAIR = 303.00

          IF ( TAIR > 268.8690 ) THEN  
              CFG = EXP( 0.04686 * TAIR - 14.30579 ) ! grass (from BEIS2)
          ELSE
              CFG = 0.0
          END IF

          CFNOG(C,R) = CFG

C.......  CFNO
          IF( .NOT. LSOIL ) THEN
          ! no soil

             TSOI = 0.72 * TAIR + 82.28
             IF (TSOI <= 273.16) TSOI = 273.16
             IF (TSOI >= 303.16) TSOI = 303.16

             CFNODRY = (1./3.) * (1./30.) * (TSOI-273.16)  ! see YL 1995 Equa 9a p. 11452
             IF (TSOI <= 283.16) THEN         ! linear cold case
                 CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28 ! see YL 1995 Equ 7b
             ELSE                             ! exponential case
                 CFNOWET = EXP(0.103 * (TSOI-273.16)) *
     &                     EXP(-0.103 * 30.0)
             END IF
             CF = 0.5 * CFNOWET + 0.5 * CFNODRY

          ELSE
          ! soil

             TSOI = SOILT( C,R )
             IF (TSOI <= 273.16) TSOI = 273.16
             IF (TSOI >= 303.16) TSOI = 303.16

             CFNODRY = (1./3.)*(1./30.)*(TSOI-273.16)  ! see YL 1995 Equa 9a p. 11452
             IF (TSOI <= 283.16) THEN         ! linear cold case
                CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28 ! see YL 1995 Equ 7b
             ELSE                             ! exponential case
                CFNOWET = EXP(0.103 * (TSOI-273.16)) *
     &                    EXP(-0.103 * 30.0)
             END IF

             SOILCAT = INT( ISLTYP( C,R ) )
             IF( SOILCAT > 0 .AND. SOILCAT <= MAXSTYPES ) THEN
                 RATIO = SOILM( C,R ) / SATURATION( SOILCAT )
                 CF = RATIO*CFNOWET + (1.-RATIO)*CFNODRY
             ELSE
                 CF = 0.
             END IF

          END IF  ! Endif LSOIL

          CFNO(C,R) = CF *
     &                FERTLZ_ADJ( JDATE, LAT(C,R) ) *
     &                VEG_ADJ( LAIc(C,R) ) * PRECADJ(C,R)

        END DO  ! loop over columns
        END DO  ! loop over rows

C******************  FORMAT  STATEMENTS   ******************************
94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020   FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )

        RETURN

        END SUBROUTINE SOILNOX

