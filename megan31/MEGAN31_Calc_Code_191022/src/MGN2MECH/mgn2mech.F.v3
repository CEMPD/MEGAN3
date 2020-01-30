      PROGRAM MGN2MECH

!***********************************************************************
!   This program does chemical speciation and MECHANISM conversion.
!   The output from megan.f is converted from 20 to 201 species which
!   are then lumped according to the MECHANISM assigned in the run script.  
!   The program loops through all timesteps of the input file.
!
!   Procedure
!   1) File set up and assign I/O parameters
!   2) Conversion from MGN 20 to speciated 201
!   3) Conversion from speciated species to MECHANISM species
!   4) Convert to tonne/hour if needed
!
!   The input file gives variables in units of g-species/sec.
!   All outputs are in mole/sec or tonne/hr depending on assignment.
!
!
!   INPUT:
!           1) MEGVEA output (netCDF-ioapi)
!           2) MEGSEA output (netCDF-ioapi)
!
!   OUTPUT:
!           1) MEGAN speciation or MECHANISM species (netCDF-ioapi)
!
!   Requirement:
!      Requires libnetcdf.a and libioapi.a to compile
!
!      setenv MGERFILE    <DEFANGED_input MEGAN output for emission activity factors>  
!      setenv OUTPFILE    <output speciated emission>
!
!   CALLS:  CHECKMEM
!
!   Originally created by Jack Chen 11/04 for MEGAN v.0
!   For MEGAN v2.0 created by Tan 12/01/06
!   For MEGAN v2.1 created by Xuemei Wang 11/04/07
!   For MEGAN v2.1 to use 150 species created by Xuemei Wang 09/30/09
!   For MEGAN v3 to use 201 species created by Alex Guenther 03/19/17
!
!  History:
!  08/14/07 Tan    - Move to MEGANv2.02 with no update
!  08/29/07 modified by A. Guenther to correct error in assigning   
!           emission factor. This version is called MEGANv2.03
!  10/29/07 modified by A. Guenther to correct omission of diurnal variation    
!           factor. This version is called MEGANv2.04
!  11/04/07 modified by Xuemei Wang to give two options for MAP or lookup table for
!           the emission factors. Also gives options for different chemical MECHANISMs
!           in the code: user modifies the external script to assign MECHANISM.
!           This version is called MEGANv2.1.0
!  06/04/08 modified by J. Lee-Taylor to accept vegetation-dependent speciation factors
!           in table format (RESHAPE tables) rather than from DATA statements.
!  09/30/08  modified by Xuemei Wang to give options for input file and test different mechanisms
!  09/27/11  Tan&Xuemei MEGANv2.10 includes soil NOx adjustment and a lot of updates
!  03/19/17  Ling Huang MEGANv3 (a) updates vegetation emission factors and
!            add additional compounds; (b) emission factor unit from
!            ug/m2/hr to nanomoles/m2/s; (c) fix NO emission bug
!  06/08/17  Ling Huang updated/added new mechanism for MEGAN3
!***********************************************************************

      IMPLICIT NONE

!...  INCLUDES:
      INCLUDE 'PARMS3.EXT'        !  I/O API parameters
      INCLUDE 'IODECL3.EXT'       !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'        !  I/O API file description data structures.

      INCLUDE 'SPC_NOCONVER.EXT'
      INCLUDE 'SPC_CB05.EXT'
      INCLUDE 'SPC_CB6.EXT'
      INCLUDE 'SPC_CBMZ.EXT'
      INCLUDE 'SPC_SAPRC99.EXT'
      INCLUDE 'SPC_SAPRC99Q.EXT'
      INCLUDE 'SPC_SAPRC99X.EXT'
      INCLUDE 'SPC_SAPRCII.EXT'
      INCLUDE 'SPC_SOAX.EXT'
      INCLUDE 'SPC_RADM2.EXT'
      INCLUDE 'SPC_RACM.EXT'
      INCLUDE 'SPC_CB6X.EXT'         ! new in MEGAN3
      INCLUDE 'SPC_CRIv2.EXT'        ! new in MEGAN3
      INCLUDE 'SPC_RACM2.EXT'        ! new in MEGAN3
      INCLUDE 'SPC_SAPRC07.EXT'      ! new in MEGAN3
      INCLUDE 'MAP_CV2CB05.EXT'
      INCLUDE 'MAP_CV2CB6.EXT'      
      INCLUDE 'MAP_CV2CBMZ.EXT'
      INCLUDE 'MAP_CV2SAPRC99.EXT'
      INCLUDE 'MAP_CV2SAPRC99Q.EXT'
      INCLUDE 'MAP_CV2SAPRC99X.EXT'
      INCLUDE 'MAP_CV2SAPRCII.EXT'
      INCLUDE 'MAP_CV2SOAX.EXT'
      INCLUDE 'MAP_CV2RADM2.EXT'
      INCLUDE 'MAP_CV2RACM.EXT'
      INCLUDE 'MAP_CV2CB6X.EXT'
      INCLUDE 'MAP_CV2CRIv2.EXT'
      INCLUDE 'MAP_CV2RACM2.EXT'
      INCLUDE 'MAP_CV2SAPRC07.EXT'

C...........   EXTERNAL FUNCTIONS and their descriptions:
      INTEGER      PROMPTFFILE, STR2INT
      REAL         STR2REAL
      INTEGER      INDEX1
      EXTERNAL     PROMPTFFILE,STR2INT
      EXTERNAL     STR2REAL
      EXTERNAL     INDEX1
      LOGICAL      ENVYN
      EXTERNAL     ENVYN
      LOGICAL      DSCGRID
      EXTERNAL     DSCGRID


!...  Program I/O files
! Program name
      CHARACTER*16 :: PROGNAME = 'MGN2MECH'
! Input MEGAN emission activity factor files
      CHARACTER*16 :: MGNERS   = 'MGNERS'    ! Input MEGAN emission activity factor 
      CHARACTER*16 :: MGNSEA   = 'MGNSEA'    ! Input MEGAN soil NO emission activity factor
! Netcdf file
      CHARACTER*16 :: EFMAPS   = 'EFMAPS'    ! EFMAP input file  name
! Output file
      CHARACTER*16 :: MGNOUT   = 'MGNOUT'    ! Output file logical name
! Parameters for file units
      INTEGER :: logdev                      ! Logfile unit number

!...  Program I/O parameters
!...  External parameters

! from run script
      LOGICAL :: TONPHR      = .FALSE.  ! Output in tons/hr flag
      LOGICAL :: CONVERSION  = .FALSE.  ! Mechanism conversion flag

      CHARACTER*16 :: MECHANISM              ! Mechanism name

! I/O API file parameters
      INTEGER       jdate          ! Looping date YYYYDDD
      INTEGER       jtime          ! Looping time HHMMSS
      INTEGER       ncols          ! Number of columns
      INTEGER       nrows          ! Number of rows
      INTEGER       mxrec          ! Total number of time steps
      INTEGER       sdate          ! Start date YYYYDDD
      INTEGER       stime          ! Start time HHMMSS
      INTEGER       tstep          ! Time step


!...  Internal parameters
! internal paramters (status and buffer)
      INTEGER   ios                          ! i/o status
      CHARACTER*256 :: MESG                  ! Message buffer

! local variables and their descriptions:
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

      INTEGER :: t, s, I, N                   ! Counters
      INTEGER :: nmpmg, nmpsp, nmpmc          ! Counters
      REAL, ALLOCATABLE :: inper(:,:,:)       ! Input emission buffer
      REAL, ALLOCATABLE :: tmper(:,:,:)       ! Temp emission buffer
      REAL, ALLOCATABLE :: outer(:,:,:)       ! Output emission buffer
      REAL, ALLOCATABLE :: EF(:,:,:)          ! Emission factor  
      REAL, ALLOCATABLE :: GAMNO(:,:)         ! NO emission factor
      INTEGER :: n_scon_spc 
      INTEGER :: NVAR,INO
      REAL :: garea 

      INTEGER, ALLOCATABLE ::  spmh_map(:),mech_map(:)   ! speciated species name
      CHARACTER*16, ALLOCATABLE :: mech_spc(:)
      REAL, ALLOCATABLE :: conv_fac(:)        
      REAL,ALLOCATABLE :: mech_mwt(:)
      
      REAL, PARAMETER :: ug2g = 1E-6          ! convert microgram to metric gram
      REAL, PARAMETER :: g2tonne = 1E-6       ! convert microgram to metric ton
      REAL, PARAMETER :: hr2sec = 3600        ! convert hr to second
      REAL, PARAMETER :: n2no   = 2.142857    
      REAL, PARAMETER :: nmol2mol   = 1E-9    ! convert nanomoles to moles
      

C***********************************************************************

!=======================================================================
!...  Begin program
!=======================================================================


!-----------------------------------------------------------------------
!.....1) File set up and assign I/O parameters
!-----------------------------------------------------------------------
!...  Initialize log file unit
      logdev = init3()

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
      IF ( .NOT. OPEN3( EFMAPS, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (EFMAPS, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( EFMAPS,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'EFMAPS has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( EFMAPS ) ) THEN
         CALL NAMEVAL (EFMAPS, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      IF ( .NOT. OPEN3( MGNERS, FSREAD3, PROGNAME ) ) THEN
          MESG = 'Could not open input file '//TRIM(MGNERS)
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( MGNERS,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'MGNERS has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( MGNERS ) ) THEN
          MESG = 'Could not get description of '//TRIM(MGNERS)
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      IF ( .NOT. OPEN3( MGNSEA, FSREAD3, PROGNAME ) ) THEN
          MESG = 'Could not open input file '//TRIM(MGNSEA)
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( MGNSEA,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'MGNSEA has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( MGNSEA ) ) THEN
          MESG = 'Could not get description of '//TRIM(MGNSEA)
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      sdate = sdate3d
      stime = stime3d
      tstep = tstep3d
      mxrec = mxrec3d
      ncols = ncols3d
      nrows = nrows3d
      garea = xcell3d * ycell3d
      jdate = sdate
      jtime = stime
      print *,' Area per grid cell (m2): ',garea

!...  Set output MECHANISM variable variables
! Check to run MECHANISM conversion
      MESG = 'Convert to any ATM MECHANISM'
      CONVERSION = ENVYN ( 'RUN_CONVERSION', MESG, .FALSE., ios )

! Check tonne per hour output
      MESG = 'Output in Tonne per hour'
      TONPHR = ENVYN ( 'SPCTONHR', MESG, .FALSE., ios )

! Set attribute and variables for output
      IF ( CONVERSION ) THEN
        MESG = 'Mechanism Name'
        CALL ENVSTR( 'MECHANISM', MESG, 'SAPRCII', MECHANISM, ios )

        SELECT CASE ( TRIM(MECHANISM) )
          CASE ('SAPRCII')
            n_scon_spc = n_saprcii
            NVARS3D = n_saprcii_spc
          CASE ('SAPRC99')
            n_scon_spc = n_saprc99
            NVARS3D = n_saprc99_spc
          CASE ('RADM2')
            n_scon_spc = n_radm2
            NVARS3D = n_radm2_spc
          CASE ('RACM')
            n_scon_spc = n_racm
            NVARS3D = n_racm_spc
          CASE ('CBMZ')
            n_scon_spc = n_cbmz
            NVARS3D = n_cbmz_spc
          CASE ('SAPRC99X')
            n_scon_spc = n_saprc99_x
            NVARS3D = n_saprc99_x_spc
          CASE ('SAPRC99Q')
            n_scon_spc = n_saprc99_q
            NVARS3D = n_saprc99_q_spc
          CASE ('CB05')
            n_scon_spc = n_cb05
            NVARS3D = n_cb05_spc
          CASE ('CB6')
            n_scon_spc = n_cb6
            NVARS3D = n_cb6_spc
          CASE ('SOAX')
            n_scon_spc = n_soax
            NVARS3D = n_soax_spc
          CASE ('CB6X')
            n_scon_spc = n_cb6x
            NVARS3D = n_cb6x_spc
          CASE ('RACM2')
            n_scon_spc = n_racm2
            NVARS3D = n_racm2_spc
          CASE ('CRIv2')
            n_scon_spc = n_criv2
            NVARS3D = n_criv2_spc
          CASE ('SAPRC07')
            n_scon_spc = n_saprc07
            NVARS3D = n_saprc07_spc
          CASE DEFAULT
            MESG = 'Error: Mechanism conversion, invalid MECHANISM: '
     &            //TRIM(MECHANISM)
            CALL M3EXIT(PROGNAME, 0, 0,MESG,2)
        ENDSELECT
 
        ALLOCATE ( spmh_map(n_scon_spc), STAT = ios )
          CALL CHECKMEM ( ios, 'spmh_map', PROGNAME )
        ALLOCATE ( mech_map(n_scon_spc), STAT = ios )
          CALL CHECKMEM ( ios, 'mech_map', PROGNAME )
        ALLOCATE ( conv_fac(n_scon_spc), STAT = ios )
          CALL CHECKMEM ( ios, 'conv_fac', PROGNAME )
        ALLOCATE ( mech_spc(NVARS3D ), STAT = ios )
          CALL CHECKMEM ( ios, 'mech_spc', PROGNAME )
        ALLOCATE ( mech_mwt(NVARS3D ), STAT = ios )
          CALL CHECKMEM ( ios, 'mech_mwt', PROGNAME )

        SELECT CASE ( TRIM(MECHANISM) )

          CASE ('SAPRCII')
            spmh_map(1:n_scon_spc) = spmh_map_saprcii(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprcii(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprcii(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprcii(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprcii(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SAPRC99')
            spmh_map(1:n_scon_spc) = spmh_map_saprc99(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprc99(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprc99(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprc99(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprc99(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('RADM2')
            spmh_map(1:n_scon_spc) = spmh_map_radm2(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_radm2(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_radm2(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_radm2(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_radm2(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('RACM')
            spmh_map(1:n_scon_spc) = spmh_map_racm(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_racm(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_racm(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_racm(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_racm(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('CBMZ')
            spmh_map(1:n_scon_spc) = spmh_map_cbmz(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_cbmz(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_cbmz(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_cbmz(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_cbmz(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SAPRC99X')
            spmh_map(1:n_scon_spc) = spmh_map_saprc99_X(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprc99_X(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprc99_X(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprc99_X(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprc99_X(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SAPRC99Q')
            spmh_map(1:n_scon_spc) = spmh_map_saprc99_Q(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprc99_Q(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprc99_Q(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprc99_Q(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprc99_Q(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('CB05')
            spmh_map(1:n_scon_spc) = spmh_map_cb05(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_cb05(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_cb05(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_cb05(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_cb05(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('CB6')
            spmh_map(1:n_scon_spc) = spmh_map_cb6(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_cb6(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_cb6(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_cb6(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_cb6(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SOAX')
            spmh_map(1:n_scon_spc) = spmh_map_soax(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_soax(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_soax(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_soax(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_soax(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('CB6X')
            spmh_map(1:n_scon_spc) = spmh_map_cb6x(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_cb6x(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_cb6x(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_cb6x(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_cb6x(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('RACM2')
            spmh_map(1:n_scon_spc) = spmh_map_racm2(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_racm2(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_racm2(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_racm2(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_racm2(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('CRIv2')
            spmh_map(1:n_scon_spc) = spmh_map_criv2(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_criv2(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_criv2(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_criv2(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_criv2(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SAPRC07')
            spmh_map(1:n_scon_spc) = spmh_map_saprc07(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprc07(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprc07(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprc07(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprc07(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
        ENDSELECT
      ELSEIF ( .NOT. CONVERSION ) THEN
        print*,'No conversion'
        NVARS3D = n_spca_spc
        VNAME3D(1:NVARS3D) = spca_spc(1:NVARS3D)
      ELSE
         MESG = 'Error: Conversion flag is not assigned.'
         CALL M3EXIT(PROGNAME, 0, 0,MESG,2)
      ENDIF

! Change the unit according to TONPHR flag
      IF ( TONPHR ) THEN
         UNITS3D(1:NVARS3D) = 'tons/hr'
      ELSE
         UNITS3D(1:NVARS3D) = 'mol/s'
      ENDIF
      VDESC3D(1:NVARS3D) = ' '
      VTYPE3D(1:NVARS3D) = m3real
    
      NVAR = NVARS3D 

      DO s = 1, NVARS3D
         PRINT*,'Output variable:',VNAME3D(s),UNITS3D(s)
      ENDDO
      ftype3d = grdded3
      sdate3d = sdate
      stime3d = stime
      tstep3d = tstep
      mxrec3d = mxrec

      IF ( CONVERSION ) THEN
      MESG = 'Mechanism Name'
      CALL ENVSTR( 'MECHANISM', MESG, 'SAPRCII', MECHANISM, ios )
         SELECT CASE ( TRIM(MECHANISM) )
         CASE ('SAPRCII')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRCII Species:'
         CASE ('SAPRC99')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRC99 Species:'
         CASE ('RADM2')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'RADM2 Species:'
         CASE ('RACM')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'RACM Species:'
         CASE ('CBMZ')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'CBMZ Species:'
         CASE ('SAPRC99X')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRC99X Species:'
         CASE ('SAPRC99Q')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRC99Q Species:'
         CASE ('CB05')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'CB05 Species:'
         CASE ('CB6')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'CB6 Species:'
         CASE ('SOAX')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SOAX Species:'
         CASE ('CB6X')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'CB6X Species:'
         CASE ('RACM2')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'RACM2 Species:'
         CASE ('SAPRC07')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRC07 Species:'
         CASE ('CRIv2')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'CRIv2 Species:'
         CASE DEFAULT
            MESG = 'Error: Mechanism conversion, invalid MECHANISM: '
     &             //TRIM(MECHANISM)
            CALL M3EXIT(PROGNAME,JDATE,JTIME,MESG,2)
         ENDSELECT

      ELSEIF ( .NOT. CONVERSION ) THEN
         FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                  //'MEGAN 201 species'
      ELSE
         MESG = 'Error: Conversion flag is not assigned.'
         CALL M3EXIT(PROGNAME, 0, 0,MESG,2)

      ENDIF

      CALL NAMEVAL ( MGNERS , MESG )  ! get input file name and path
      FDESC3D( 2 ) = 'Input MEGAN file: '//TRIM(MESG)

      IF ( .NOT. OPEN3( MGNOUT, FSCREA3, PROGNAME ) ) THEN
        MESG = 'Could not open file '//TRIM(MGNOUT)
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Allocate memory
      ALLOCATE ( EF   ( N_MGN_SPC, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'EF',    PROGNAME )

      ALLOCATE ( GAMNO ( ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'GAMNO', PROGNAME )


!...  Read EFMAPS

      DO N = 1, N_MGN_SPC
        IF ( .NOT. READ3(EFMAPS,'EF_'//MGN_SPC(N),ALLAYS3,0,0,
     &                                        EF(N,:,:))) THEN
           MESG = 'Error reading '//TRIM('EF_'//MGN_SPC(N))
           CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
        ENDIF
      ENDDO

!.....2) Conversion from MGN 20 to speciated 201
!-----------------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( inper(N_MGN_SPC, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'inper', PROGNAME )
      ALLOCATE ( tmper( n_spca_spc, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'tmper', PROGNAME )
      ALLOCATE ( outer( NVAR, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'outer', PROGNAME )

      INO = INDEX1('NO',N_MGN_SPC,MGN_SPC)
!...  Loop through time
      jdate = sdate
      jtime = stime
      DO t = 1, mxrec
        inper = 0.
        tmper = 0.
        outer = 0.

        DO N = 1, N_MGN_SPC
          IF ( .NOT. READ3(MGNERS,TRIM(MGN_SPC(N)),1,
     &                      jdate,jtime,inper(N,:,:))) THEN
            MESG = 'Error reading '//TRIM(MGN_SPC(nmpmg))
            CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
          ENDIF
        ENDDO

        IF ( .NOT. READ3(MGNSEA,'GAMNO',1,jdate,jtime,GAMNO(:,:))) THEN
          MESG = 'Error reading GAMNO'
          CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
        ENDIF

        DO s = 1, N_SMAP_SPC
            
          nmpmg = mg20_map(s)
          nmpsp = spca_map(s)
          IF (t .EQ. 1) THEN
            MESG='Convert '//MGN_SPC(NMPMG)//' to '//SPCA_SPC(NMPSP)
            CALL M3MESG( MESG )
          ENDIF

         IF ( nmpmg .NE. INO ) then
           !...  Not NO
           MESG = 'Use EFMAPS for '//MGN_SPC(NMPMG)
           CALL M3MESG( MESG )
           tmper(nmpsp,:,:) = inper(nmpmg,:,:) * EF(nmpmg,:,:)
     &                          * effs_all(nmpsp)

         ELSEIF ( nmpmg .EQ. INO ) then
!!-----------------NO Stuff-----------------------
!     GAMNO is emission activity factor
           tmper(nmpsp,:,:) = GAMNO(:,:) * EF(INO,:,:)
     &                           * effs_all(nmpsp)

!-----------------end of NO----------------------
         ENDIF     !IF ( nmpmg .NE. INO ) then

        ENDDO ! End species loop

!-----------------------------------------------------------------------
!.....3) Conversion from speciated species to MECHANISM species
!-----------------------------------------------------------------------
        ! convert from nanomol/m^2/s to mol/m^2/hr
        DO s = 1, n_spca_spc
           tmper(s,:,:) = tmper(s,:,:) * nmol2mol * hr2sec
        ENDDO

        IF ( CONVERSION ) THEN
          ! lumping to MECHANISM species

          DO s = 1, n_scon_spc
            nmpsp = spmh_map(s)         ! Mapping value for SPCA
            nmpmc = mech_map(s)         ! Mapping value for MECHANISM
           IF ( nmpmc .NE. 999 ) THEN
            IF (t .EQ. 1) THEN
              MESG='Convert '//SPCA_SPC(NMPSP)//' to '//MECH_SPC(NMPMC)
              CALL M3MESG( MESG )
            ENDIF
            outer(nmpmc,:,:) = outer(nmpmc,:,:) +
     &                        (tmper(nmpsp,:,:) * conv_fac(s)) *
     &                         garea/hr2sec
            ! units of these species are in mole/s
           ENDIF
          ENDDO ! End species loop
          WRITE(MESG,*) 'conv_fac,outer',conv_fac(1),outer(1,101,101)
          CALL M3MESG(MESG)

        ELSE
          ! get all 201 species into the output array
          outer(:,:,:) = tmper(:,:,:)*garea/hr2sec
          ! units of these species are in mole/s/grid

        ENDIF

!-----------------------------------------------------------------------
!.....4) Convert to tonne/hour if needed
!-----------------------------------------------------------------------
        IF ( TONPHR ) THEN
          IF ( CONVERSION ) THEN
            ! convert from mol/s to ton/hr
            DO s = 1, NVAR
              outer(s,:,:) = outer(s,:,:)*mech_mwt(s)*
     &                       g2tonne*hr2sec
            ENDDO
          ELSEIF ( .NOT. CONVERSION ) THEN
            ! convert from mol/s to ton/hr
            DO s = 1, NVAR
              outer(s,:,:) = outer(s,:,:)*spca_mwt(s)*
     &                       g2tonne*hr2sec
            ENDDO
          ENDIF
        ENDIF

!-----------------------------------------------------------------------
!.....5) Write out the calculated emissions
!-----------------------------------------------------------------------
        print*,'VARS=',NVAR,n_scon_spc
        DO s = 1, NVAR
          IF ( .NOT. WRITE3(MGNOUT,VNAME3D(s),jdate,jtime,
     &                                        outer(s,:,:)    )  ) THEN
             MESG = 'Error writing to file: '//TRIM(MGNOUT)
             CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
          ENDIF
        ENDDO 

!-----------------------------------------------------------------------
        CALL NEXTIME(jdate,jtime,tstep)
      ENDDO ! End time loop

!... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

      DEALLOCATE ( spmh_map )
      DEALLOCATE ( mech_map )
      DEALLOCATE ( conv_fac )
      DEALLOCATE ( mech_spc )
      DEALLOCATE ( mech_mwt )
      DEALLOCATE (EF)
      DEALLOCATE (GAMNO)
      DEALLOCATE (outer)
      DEALLOCATE (inper)
      DEALLOCATE (tmper)

!=======================================================================
!...  FORMAT
!=======================================================================        


!=======================================================================
!...  End program
!=======================================================================
      END PROGRAM MGN2MECH
