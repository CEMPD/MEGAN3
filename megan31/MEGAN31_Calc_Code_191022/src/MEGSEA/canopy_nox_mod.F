      MODULE CANOPY_NOX_MOD
      !Used by BDSNP_MOD to calculate canopy reduction of above canopy Soil NOx Emissions
      !Modified for MEGAN3.1 by Ling Huang 2019/07/14
      
      IMPLICIT NONE

      PUBLIC :: GET_CANOPY_NOX
      
C      PRIVATE :: DIFFG 
      
      CONTAINS
      
C ----------------------------------------------------------------------

      SUBROUTINE GET_CANOPY_NOX(JDATE, JTIME, NC, NR, COSZEN,
     & TASFC, SSOLAR, PRES, LANDTYPE, LAI, SNOCOV, CFRAC, WSPD, CRF) ! called tmpbeis, change called BDSNP, add K argument

      IMPLICIT NONE

C Includes:
      INCLUDE 'PARMS3.EXT'   ! I/O API parameters
      INCLUDE 'IODECL3.EXT'  ! I/O API function declarations
      INCLUDE 'FDESC3.EXT'   ! I/O API file desc. data structures

C Arguments
      INTEGER, INTENT( IN )  :: JDATE             ! current simulation date (YYYYDDD)
      INTEGER, INTENT( IN )  :: JTIME             ! current simulation time (HHMMSS)
      INTEGER, INTENT( IN ) :: NC      ! no. columns
      INTEGER, INTENT( IN ) :: NR      ! no. rows
C     These are arrays 
      REAL,    INTENT( IN ) :: COSZEN( NC,NR )        ! cosine of zenith angle
      REAL,    INTENT( IN ) :: TASFC ( NC,NR )        ! surface air temperature [K]
      REAL,    INTENT( IN ) :: SSOLAR( NC,NR )        ! surface radiation [w/m**2]
      REAL,    INTENT( IN ) :: PRES  ( NC,NR )        ! surface pressure [Pa]
      INTEGER, INTENT( IN ) :: LANDTYPE( NC,NR )     ! the biome type in each cell
      REAL,    INTENT( IN ):: LAI  ( NC,NR )        ! leaf area index (m2/m2)
      REAL,    INTENT( IN ):: SNOCOV  ( NC,NR )        ! snow cover
      REAL,    INTENT( IN ):: CFRAC  ( NC,NR )        ! cloud fraction
      REAL,    INTENT( IN ):: WSPD  ( NC,NR )        ! cloud fraction
      REAL,    INTENT( OUT ):: CRF  ( NC,NR )        ! outputs the canopy reduction factor
      ! !LOCAL VARIABLES:
      !
      CHARACTER( 16 )  :: PNAME = 'CANOPY_NOX'  ! procedure name
      INTEGER          IOS                ! IO or memory allocation status
      CHARACTER( 256 ) :: MESG            ! message buffer
      ! Scalars
      INTEGER :: C, R, K, KK, MY_NCOLS, MY_NROWS
      REAL*8  :: F0,     HSTAR, XMW              
      REAL*8  :: DTMP1,  DTMP2, DTMP3,  DTMP4, GFACT, GFACI
      REAL*8  :: RT,     RAD0,  RIX,    RIXX,  RDC,   RLUXX
      REAL*8  :: RGSX,   RCLX,  TEMPK,  TEMPC, WINDSQR
      REAL*8 :: VFNEW
      
      LOGICAL, SAVE          :: FIRSTCANOPY = .TRUE. 

      ! Arrays
      REAL*8  :: RI  (24)       
      REAL*8  :: RLU (24)      
      REAL*8  :: RAC (24)      
      REAL*8  :: RGSS(24)     
      REAL*8  :: RGSO(24)     
      REAL*8  :: RCLS(24)     
      REAL*8  :: RCLO(24)
! !DEFINED PARAMETERS:
!      
      INTEGER, PARAMETER :: SNIRI(24) = (/9999, 200, 9999, 9999, 9999, 9999, 
     & 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 400, 400, 
     & 200, 200, 200, 9999, 200/)

      INTEGER, PARAMETER :: SNIRLU(24) = (/9999, 9000, 9999, 9999, 9999, 
     & 9999, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 9000, 1000, 
     & 9000, 9000, 9000, 9000, 1000, 9000, 9999, 9000/)

      INTEGER, PARAMETER :: SNIRAC(24) = (/0, 300, 0, 0, 0, 0, 100, 100, 
     & 100, 100, 100, 100, 100, 100, 2000, 2000, 2000, 2000, 2000, 2000, 
     & 2000, 200, 100, 200/)

      INTEGER, PARAMETER :: SNIRGSS(24) = (/0, 0, 100, 1000, 100, 1000, 350, 
     & 350, 350, 350, 350, 350, 350, 350, 500, 200, 500, 500, 500, 500, 
     & 200, 150, 400, 150/)

      INTEGER, PARAMETER :: SNIRGSO(24) = (/2000, 1000, 3500, 400, 3500, 
     & 400, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 
     & 200, 200, 200, 150, 300, 150/)

      INTEGER, PARAMETER :: SNIRCLS(24) = (/9999, 2500, 9999, 9999, 9999, 
     & 9999, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 9999, 
     & 2000, 2000, 2000, 2000, 9999, 2000, 9999, 2000/)
    
      INTEGER, PARAMETER :: SNIRCLO(24) = (/9999, 1000, 1000, 9999, 1000, 
     & 9999, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 9999, 
     & 1000, 1000, 1000, 1000, 9999, 1000, 9999, 1000/)

      INTEGER, PARAMETER :: SNIVSMAX(24) = (/10, 100, 100, 10, 100, 10, 100, 
     & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 
     & 100, 100, 100, 100/)   

      REAL, PARAMETER :: DRYCOEFF(20) = (/-0.358, 3.02, 3.85, -0.0978, -3.66, 
     & 12, 0.252, -7.8, 0.226, 0.274, 1.14, -2.19, 0.261, -4.62, 0.685, 
     & -0.254, 4.37, -0.266, -0.159, -0.206 /)   

      ! Canopy wind extinction coefficients
      ! (cf. Yienger & Levy [1995], Sec 5), now a function of the MODIS/KOPPEN biometype (J.D. Maasakkers)
       REAL*8,  PARAMETER :: SOILEXC(24)    = (/ 
     &  0.10, 0.50, 0.10, 0.10, 0.10,
     &  0.10, 0.10, 0.10, 0.10, 1.00,
     &  1.00, 1.00, 1.00, 2.00, 4.00,
     &  4.00, 4.00, 4.00, 4.00, 4.00,
     &  4.00, 2.00, 0.10, 2.00                  /)     
      

      ! Molecular weight of water [kg]
      REAL*8, PARAMETER :: XMWH2O = 18d-3
      !
      ! Ventilation velocity for NOx, day & night values [m/s]
      REAL*8,  PARAMETER :: VFDAY   = 1.0d-2
      REAL*8,  PARAMETER :: VFNIGHT = 0.2d-2 
      REAL*8, PARAMETER :: PRESS  = 1.5d5

      ! Set physical parameters
      HSTAR = 0.01d0              ! Henry's law constant
      F0    = 0.1d0               ! Reactivity factor for biological oxidation 
      XMW   = 46d-3               ! Molecular wt of NO2 (kg)

      IF( FIRSTCANOPY ) THEN
        FIRSTCANOPY = .FALSE.
      END IF
   
      CRF = 0d0 ! array
      
      ! begin calculating canopy reduction factor
      DO C=1, NC
        DO R=1, NR
          IF(LAI(C,R) > 0.0) THEN
            TEMPC = TASFC(C,R) - 273.15d0 ! convert kelvin to Celsius
      ! Compute bulk surface resistance for gases.    
         !                                  
         !  Adjust external surface resistances for temperature; 
         !  from Wesely [1989], expression given in text on p. 1296.        
            RT = 1000.0D0 * EXP( -TEMPC - 4.0d0 )
         
                
         !--------------------------------------------------------------
         ! Get surface resistances - loop over biome types K
         !
         ! The land types within each grid square are defined using the 
         ! Olson land-type database.  Each of the Olson land types is 
         ! assigned a corresponding "deposition land type" with 
         ! characteristic values of surface resistance components.  
         ! There are 74 Olson land-types but only 11 deposition 
         ! land-types (i.e., many of the Olson land types share the 
         ! same deposition characteristics).  Surface resistance 
         ! components for the "deposition land types" are from Wesely 
         ! [1989] except for tropical forests [Jacob and Wofsy, 1990] 
         ! and for tundra [Jacob et al., 1992].  All surface resistance 
         ! components are normalized to a leaf area index of unity.
         !--------------------------------------------------------------
	!Set biometype
         
            K = LANDTYPE( C,R )
            
            ! Set second loop variable to K to allow snow/ice correction
	     KK = K

            ! If the surface is snow or ice, then set K=3
            IF ( SNOCOV(C,R) .EQ. 1 ) KK = 3

		!USE new MODIS/KOPPEN Biometypes to read data

            ! Read the internal resistance RI (minimum stomatal resistance 
            ! for water vapor, per unit area of leaf) from the IRI array; 
            ! a '9999' value means no deposition to stomata so we impose a 
            ! very large value for RI.
            RI(K) = DBLE( SNIRI(KK) )
            IF ( RI(K) >= 9999.D0 ) RI(K)= 1.D12
            
            ! Cuticular resistances IRLU read in from 'drydep.table'
            ! are per unit area of leaf; divide them by the leaf area index 
            ! to get a cuticular resistance for the bulk canopy.  If IRLU is 
            !'9999' it means there are no cuticular surfaces on which to 
            ! deposit so we impose a very large value for RLU.
            IF ( SNIRLU(KK) >= 9999 .OR. LAI(C,R) <= 0d0 ) THEN
               RLU(K)  = 1.D6
            ELSE
               RLU(K)= DBLE( SNIRLU(KK) ) / LAI(C,R) + RT
            ENDIF

            ! The following are the remaining resistances for the Wesely
            ! resistance-in-series model for a surface canopy
            ! (see Atmos. Environ. paper, Fig.1).  
            RAC(K)  = MAX( DBLE( SNIRAC(KK)  ),      1d0 )
            RGSS(K) = MAX( DBLE( SNIRGSS(KK) ) + RT, 1d0 )
            RGSO(K) = MAX( DBLE( SNIRGSO(KK) ) + RT, 1d0 ) 
            RCLS(K) =      DBLE( SNIRCLS(KK) ) + RT           
            RCLO(K) =      DBLE( SNIRCLO(KK) ) + RT 

            IF (  RAC(K) >= 9999.D0 ) RAC(K)  = 1d12
            IF ( RGSS(K) >= 9999.D0 ) RGSS(K) = 1d12
            IF ( RGSO(K) >= 9999.D0 ) RGSO(K) = 1d12
            IF ( RCLS(K) >= 9999.D0 ) RCLS(K) = 1d12         
            IF ( RCLO(K) >= 9999.D0 ) RCLO(K) = 1d12

            !-------------------------------------------------------------
            ! Adjust stomatal resistances for insolation and temperature:  
            ! 
            ! Temperature adjustment is from Wesely [1989], equation (3).
            ! 
            ! Light adjustment by the function BIOFIT is described by Wang 
            ! [1996].  It combines:
            !
            ! - Local dependence of stomal resistance on the intensity I 
            !   of light impinging the leaf; this is expressed as a 
            !   multiplicative factor I/(I+b) to the stomatal resistance 
            !   where b = 50 W m-2
            !   (equation (7) of Baldocchi et al. [1987])
            ! - Radiative transfer of direct and diffuse radiation in the 
            !   canopy using equations (12)-(16) from Guenther et al. 
            !   [1995]
            ! - Separate accounting of sunlit and shaded leaves using
            !   equation (12) of Guenther et al. [1995]
            ! - Partitioning of the radiation at the top of the canopy 
            !   into direct and diffuse components using a 
            !   parameterization to results from an atmospheric radiative 
            !   transfer model [Wang, 1996]
            !
            ! The dependent variables of the function BIOFIT are the leaf 
            ! area index (XYLAI), the cosine of zenith angle (SUNCOS) and 
            ! the fractional cloud cover (CFRAC).  The factor GFACI 
            ! integrates the light dependence over the canopy depth; so
            ! be scaled by LAI to yield a bulk canopy value because that's 
            ! already done in the GFACI formulation.
            !-------------------------------------------------------------

            ! Radiation @ sfc [W/m2]
            RAD0 = SSOLAR(C,R)
            
            ! Internal resistance
            RIX  = RI(K)

            ! Skip the following block if the resistance RIX is high
            IF ( RIX < 9999d0 ) THEN
               GFACT = 100.0D0

               IF ( TEMPC > 0.D0 .AND. TEMPC < 40.D0) THEN
                  GFACT = 400.D0 / TEMPC / ( 40.0D0 - TEMPC )
               ENDIF

               GFACI = 100.D0

               IF ( RAD0 > 0d0 .AND. LAI(C,R) > 0d0 ) THEN
                  GFACI= 1d0 / 
     &                   BIOFIT( DRYCOEFF,       LAI(C,R),
     &                           COSZEN(C,R), CFRAC(C,R)    )
               ENDIF
            
               RIX = RIX * GFACT * GFACI
            ENDIF
            
            ! Compute aerodynamic resistance to lower elements in lower 
            ! part of the canopy or structure, assuming level terrain - 
            ! equation (5) of Wesely [1989].                     
            RDC = 100.D0*(1.0D0+1000.0D0/(RAD0 + 10.D0))

            ! Loop over species; species-dependent corrections to resistances
            ! are from equations (6)-(9) of Wesely [1989].
            !
            ! NOTE: here we only consider NO2 (bmy, 6/22/09)
            RIXX   = RIX * DIFFG( TASFC(C,R), PRESS, XMWH2O ) /
     &                     DIFFG( TASFC(C,R), PRESS, XMW    )
     &             + 1.D0 / ( HSTAR/3000.D0 + 100.D0*F0  )

            RLUXX  = 1.D12

            IF ( RLU(K) < 9999.D0 ) THEN
               RLUXX = RLU(K) / ( HSTAR / 1.0D+05 + F0 )
            ENDIF
            
            ! To prevent virtually zero resistance to species with huge HSTAR, 
            ! such as HNO3, a minimum value of RLUXX needs to be set. 
            ! The rationality of the existence of such a minimum is 
            ! demonstrated by the observed relationship between Vd(NOy-NOx) 
            ! and Ustar in Munger et al.[1996]; Vd(HNO3) never exceeds 2 cm/s 
            ! in observations. The corresponding minimum resistance is 50 s/m.
            ! was introduced by J.Y. Liang on 7/9/95.
            RGSX = 1d0 / ( HSTAR/1d5/RGSS(K) + F0/RGSO(K) )
            RCLX = 1d0 / ( HSTAR/1d5/RCLS(K) + F0/RCLO(K) )

            ! Get the bulk surface resistance of the canopy
            ! from the network of resistances in parallel and in series 
            ! (Fig. 1 of Wesely [1989])
            DTMP1 = 1.D0 / RIXX
            DTMP2 = 1.D0 / RLUXX
            DTMP3 = 1.D0 / ( RAC(K) + RGSX )
            DTMP4 = 1.D0 / ( RDC      + RCLX )

            ! Save the within canopy depvel of NOx, used in calculating 
            ! the canopy reduction factor for soil emissions [1/s]
            CRF(C,R) = DTMP1 + DTMP2 + DTMP3 + DTMP4
            
            ! Pick proper ventilation velocity for day or night
            IF ( COSZEN( C,R ) > 0d0 ) THEN
               VFNEW = VFDAY              
            ELSE 
               VFNEW = VFNIGHT            
            ENDIF

      ! If the leaf area index and the bulk surface resistance
      ! of the canopy to NOx deposition are both nonzero ...
            IF (CRF(C,R) > 0d0 ) THEN

         ! Adjust the ventilation velocity.  
         ! NOTE: SOILEXC(21) is the canopy wind extinction 
         ! coefficient for the tropical rainforest biome.
              WINDSQR=WSPD(C,R)*WSPD(C,R)
              VFNEW    = (VFNEW * SQRT( WINDSQR/9d0 * 7d0/LAI(C,R)) *
     &                          ( SOILEXC(21)  / SOILEXC(K) ))

         ! Soil canopy reduction factor
              CRF(C,R) = CRF(C,R) / ( CRF(C,R) + VFNEW )
         
C         IF( CRF(C,R) > 1.0 ) THEN
C         write(*,*) 'CANOPY NOX REDUCTION FACTOR TOO HIGH'
C         write(*,*) 'C,R,k,crf,vfnew,dt1,dt2,dt3,dt4'
C         write(*,*) C,R,k, CRF(C,R), VFNEW, DTMP1, DTMP2, DTMP3, DTMP4 
C         write(*,*) 'windsqr, cos, lai'
C         write(*,*) WINDSQR, COSZEN(C,R),LAI(C,R)
C         write(*,*) 'soil21/soilk'
C         write(*,*) SOILEXC(21)/SOILEXC(K) 
C         write(*,*) 'rdc,rclx,rac,rgsx,rluxx,rixx,rclo'
C         write(*,*) RDC, RCLX, RAC(K), RGSX, RLUXX, RIXX, RCLO(K)
C         write(*,*) 'RCLS(K), RGSO(K), RGSS(K), RLU(K)'
C         write(*,*) RCLS(K), RGSO(K), RGSS(K), RLU(K)
C         write(*,*) 'DIFFG(h20)/DIFFG(no2)'        
c         write(*,*) DIFFG( TASFC(C,R), PRESS, XMWH2O ) /
C     &                     DIFFG( TASFC(C,R), PRESS, XMW    )
C        write(*,*) DIFFG( TASFC(C,R), PRESS, XMWH2O )
C         write(*,*) DIFFG( TASFC(C,R), PRESS, XMW    )
C         write(*,*) TASFC(C,R), PRESS, XMW, XMWH2O
C         write(*,*) '1.D0 / ( HSTAR/3000.D0 + 100.D0*F0  )'
C         write(*,*) 1.D0 / ( HSTAR/3000.D0 + 100.D0*F0  )
C         write(*,*) 'RAD0, GFACT, GFACI, BIOFIT'
C         write(*,*) RAD0, GFACT, GFACI, BIOFIT( DRYCOEFF, LAI(C,R),
C     &                           COSZEN(C,R), CFRAC(C,R))
C         write(*,*) 'TEMPC, CFRAC(C,R), DRYCOEFF(K), RT'
C         write(*,*) TEMPC, CFRAC(C,R), DRYCOEFF(K), RT
C         write(*,*) 'RI(K),RLU(K),RAC(K),RGSS(K),RGSO(K),RCLS(K),RCLO(K)'
C         write(*,*) RI(K), RLU(K),RAC(K),RGSS(K),RGSO(K),RCLS(K),RCLO(K)
C         MESG = 'CRF too high'
C         CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT1 )
C         END IF   

            ELSE ! CRF < 0.0
     
         ! Otherwise set the soil canopy reduction factor to zero
              CRF(C,R) = 0d0

            END IF

            IF( CRF(C,R) .LT. 0.0) THEN
            
            MESG = 'CRF Less than 0'
            CALL M3EXIT( PNAME, JDATE, JTIME, MESG, 2 )
            
            ELSE IF( CRF(C,R) .GT. 1.0) THEN
            
            MESG = 'CRF Greater than one'
            CALL M3EXIT( PNAME, JDATE, JTIME, MESG, 2 ) 
            
            END IF
            
          ELSE
            CRF(C,R) = 0.0
          END IF !lai check

        END DO !row loop
      END DO !col loop
            
      END SUBROUTINE GET_CANOPY_NOX
      
      FUNCTION DIFFG( TK, PRESS, XM ) RESULT( DIFF_G )
! !DESCRIPTION: Function DIFFG calculates the molecular diffusivity [m2/s] in 
!  air for a gas X of molecular weight XM [kg] at temperature TK [K] and 
!  pressure PRESS [Pa].
!\\
!\\
!  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
!  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
!  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
!  also gives radii for some other molecules.  Rather than requesting the user
!  to supply a molecular radius we specify here a generic value of 2.E-10 m for
!  all molecules, which is good enough in terms of calculating the diffusivity
!  as long as molecule is not too big.
!
! !INPUT PARAMETERS:
!
      REAL, INTENT(IN) :: TK      ! Temperature [K]
      REAL*8, INTENT(IN) :: PRESS   ! Pressure [Pa]
      REAL*8, INTENT(IN) :: XM      ! Molecular weight of gas [kg]
!
! !RETURN VALUE:
!
      REAL*8             :: DIFF_G  ! Molecular diffusivity [m2/s]
!
! !REVISION HISTORY:
!     22 Jun 2009 - R. Yantosca - Copied from "drydep_mod.f"
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: AIRDEN, Z, DIAM, FRPATH, SPEED            
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER  :: XMAIR  = 28.8d-3 
      REAL*8, PARAMETER  :: RADAIR = 1.2d-10
      REAL*8, PARAMETER  :: PI     = 3.1415926535897932d0
      REAL*8, PARAMETER  :: RADX   = 1.5d-10
      REAL*8, PARAMETER  :: RGAS   = 8.32d0
      REAL*8, PARAMETER  :: AVOGAD = 6.023d23

      !=================================================================
      ! DIFFG begins here!
      !=================================================================

      ! Air density
      AIRDEN = ( PRESS * AVOGAD ) / ( RGAS * TK )

      ! DIAM is the collision diameter for gas X with air.
      DIAM   = RADX + RADAIR

      ! Calculate the mean free path for gas X in air: 
      ! eq. 8.5 of Seinfeld [1986];
      Z      = XM  / XMAIR
      FRPATH = 1d0 /( PI * SQRT( 1d0 + Z ) * AIRDEN*( DIAM**2 ) )

      ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED  = SQRT( 8d0 * RGAS * TK / ( PI * XM ) )

      ! Calculate diffusion coefficient of gas X in air; 
      ! eq. 8.9 of Seinfeld [1986]
      DIFF_G = ( 3d0 * PI / 32d0 ) * ( 1d0 + Z ) * FRPATH * SPEED

      ! Return to calling program
      END FUNCTION DIFFG
      
      SUBROUTINE SUNPARAM(X)

      IMPLICIT NONE


!===============================================
! the sequence is lai,suncos,cloud fraction
!===============================================
!  NN = number of variables (lai,suncos,cloud fraction)
      INTEGER NN
      PARAMETER(NN=3)
!  ND = scaling factor for each variable
      INTEGER ND(NN),I
      DATA ND /55,20,11/
!  X0 = maximum for each variable
      REAL X(NN),X0(NN),XLOW
      DATA X0 /11.,1.,1./

      DO I=1,NN
        X(I)=MIN(X(I),X0(I))
! XLOW = minimum for each variable
        IF (I.NE.3) THEN
          XLOW=X0(I)/REAL(ND(I))
        ELSE
          XLOW= 0.
        END IF
        X(I)=MAX(X(I),XLOW)
        X(I)=X(I)/X0(I)
      END DO

      RETURN
      END SUBROUTINE SUNPARAM
      
      REAL*8 FUNCTION BIOFIT(COEFF1,XLAI1,SUNCOS1,CFRAC1)

      IMPLICIT NONE

!===============================================
! Calculate the light correction
!===============================================
!* BIOFIT and SUNPARAM were written by Y.H. Wang.   
!*             !-------------------------------------------------------------
            ! Adjust stomatal resistances for insolation and temperature:  
            ! 
            ! Temperature adjustment is from Wesely [1989], equation (3).
            ! 
            ! Light adjustment by the function BIOFIT is described by Wang 
            ! [1996].  It combines:
            !
            ! - Local dependence of stomal resistance on the intensity I 
            !   of light impinging the leaf; this is expressed as a 
            !   multiplicative factor I/(I+b) to the stomatal resistance 
            !   where b = 50 W m-2
            !   (equation (7) of Baldocchi et al. [1987])
            ! - Radiative transfer of direct and diffuse radiation in the 
            !   canopy using equations (12)-(16) from Guenther et al. 
            !   [1995]
            ! - Separate accounting of sunlit and shaded leaves using
            !   equation (12) of Guenther et al. [1995]
            ! - Partitioning of the radiation at the top of the canopy 
            !   into direct and diffuse components using a 
            !   parameterization to results from an atmospheric radiative 
            !   transfer model [Wang, 1996]
            !
            ! The dependent variables of the function BIOFIT are the leaf 
            ! area index (XYLAI), the cosine of zenith angle (SUNCOS) and 
            ! the fractional cloud cover (CFRAC).  The factor GFACI 
            ! integrates the light dependence over the canopy depth; so
            ! be scaled by LAI to yield a bulk canopy value because that's 
            ! already done in the GFACI formulation.
!*************************************************************
      INTEGER KK
      PARAMETER (KK=4)
      REAL COEFF1(20),TERM(KK),REALTERM(20)
      REAL XLAI1,SUNCOS1,CFRAC1
      INTEGER K,K1,K2,K3

      TERM(1)=1.
      TERM(2)=XLAI1
      TERM(3)=SUNCOS1
      TERM(4)=CFRAC1
      CALL SUNPARAM(TERM(2))
      K=0
      DO K3=1,KK
        DO K2=K3,KK
          DO K1=K2,KK
            K=K+1
            REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3)
          END DO
        END DO
      END DO
      BIOFIT=0
      DO K=1,20
        BIOFIT=BIOFIT+COEFF1(K)*REALTERM(K)
      END DO
      IF (BIOFIT.LT.0.1) BIOFIT=0.1

      RETURN
      END FUNCTION BIOFIT
      
      END MODULE
      
      
      
!  References:
!  ============================================================================
!  (1 ) Baldocchi, D.D., B.B. Hicks, and P. Camara, "A canopy stomatal
!        resistance model for gaseous deposition to vegetated surfaces",
!        Atmos. Environ. 21, 91-101, 1987.
!  (2 ) Brutsaert, W., "Evaporation into the Atmosphere", Reidel, 1982.
!  (3 ) Businger, J.A., et al., "Flux-profile relationships in the atmospheric 
!        surface layer", J. Atmos. Sci., 28, 181-189, 1971.
!  (4 ) Dwight, H.B., "Tables of integrals and other mathematical data",
!        MacMillan, 1957.
!  (5 ) Guenther, A., and 15 others, A global model of natural volatile
!         organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!  (6 ) Hicks, B.B., and P.S. Liss, "Transfer of SO2 and other reactive
!        gases across the air-sea interface", Tellus, 28, 348-354, 1976.
!  (7 ) Jacob, D.J., and S.C. Wofsy, "Budgets of reactive nitrogen,
!        hydrocarbons, and ozone over the Amazon forest during the wet season",
!        J.  Geophys. Res., 95, 16737-16754, 1990.
!  (8 ) Jacob, D.J., et al, "Deposition of ozone to tundra", J. Geophys. Res., 
!        97, 16473-16479, 1992.
!  (9 ) Levine, I.N., "Physical Chemistry, 3rd ed.", McGraw-Hill, 
!        New York, 1988.
!  (10) Munger, J.W., et al, "Atmospheric deposition of reactive nitrogen 
!        oxides and ozone in a temperate deciduous forest and a sub-arctic 
!        woodland", J. Geophys. Res., in press, 1996.
!  (11) Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, "SO2, sulfate, 
!        and HNO3 deposition velocities computed using regional landuse and
!        meteorological data", Atmos. Environ., 20, 949-964, 1986.
!  (12) Wang, Y.H., paper in preparation, 1996.
!  (13) Wesely, M.L, "Improved parameterizations for surface resistance to
!        gaseous dry deposition in regional-scale numerical models", 
!        Environmental Protection Agency Report EPA/600/3-88/025,
!        Research Triangle Park (NC), 1988.
!  (14) Wesely, M. L., Parameterization of surface resistance to gaseous dry 
!        deposition in regional-scale numerical models.  Atmos. Environ., 23
!        1293-1304, 1989. 
!  (15) Price, H., L. Jaeglé, A. Rice, P. Quay, P.C. Novelli, R. Gammon, 
!        Global Budget of Molecular Hydrogen and its Deuterium Content: 
!        Constraints from Ground Station, Cruise, and Aircraft Observations,
!        submitted to J. Geophys. Res., 2007.      
      
      
