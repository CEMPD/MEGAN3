
      PROGRAM TXT2IOAPI
!***********************************************************************
!  This program converts comma deliminated line input to gridded
!    netCDF-IOAPI format.
!  The grid parameter is input from GRIDDESC
!  Program always takes second and third input column as grid-x 
!    and grid-y index.
!  If you export the file from Excel to .csv file, make sure you convert 
!    the DOS text file to Unix text file with "dos2unix" command.
!
!  Requirement:
!     Requires libnetcdf.a and libioapi.a to compile
!     Also, parsline.o and chkmetem.o from SMOKE
!
!     setenv EFSTXTF <ef_map text file>
!     setenv CTTXTF <canopy type    text file>
!     setenv LAITXTF <lai    text file>
!     setenv EFMAPS  <output EFS netCDF-ioapi file>
!     setenv CANTYP  <output CANTYP netCDF-ioapi file>
!     setenv LAIS46  <output LAI netCDF-ioapi file>
!     setenv RUN_EFS <T|F>
!     setenv RUN_LAI <T|F>
!     setenv RUN_CANTYP <T|F>
!     setenv GRIDDESC <grid desc file>
!     setenv W126TXTF <W126 text file> New in MEGAN3
!     setenv LDFTXTF <LDF text file> New in MEGAN3
!Jiang
!     setenv FERTTXTF <FERT text file> New in MEGAN3 
!     setenv ARIDTXTF <ARID text file> New in MEGAN3
!     setenv NONARIDTXTF <NONARID text file> New in MEGAN3
!     setenv NITROGEN <NITROGEN text file> New in MEGAN3
!     setenv LANDTYPE <LANDTYPE text file> New in MEGAN3
!
!  Call:
!     PARSLINE
!     CHECKMEM
!
!  Created by Jack Chen 11/04
!  Modified by Tan 11/13/06 for MEGAN v2.0
!
!  History:
!  08/14/07 Tan    - Move to MEGANv2.02 with no update
!  07/06/11 Tan    - Update for MEGANv2.10
!  03/31/17 Ling Huang - Update for MEGANv3: add three subroutines that
!                  convert LDF, W126 and Canopy type data
!  07/24/19 Xiaoyan Jiang - Add subroutines to convert FERT, ARID, NONARID, NITROGEN, and LANDTYPE data
!***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
      INCLUDE 'PARMS3.EXT'    !  I/O API parameters
      INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

      INCLUDE 'SPC_MGN.EXT'    !  MEGAN species
      INCLUDE 'LAI_MGN.EXT'    !  LAI
      INCLUDE 'CANTYP_MGN.EXT'    !  Canopy type
      INCLUDE 'FERT_MGN.EXT'
C...........   EXTERNAL FUNCTIONS and their descriptions:
      INTEGER      PROMPTFFILE, STR2INT
      REAL         STR2REAL
      LOGICAL      DSCGRID
      LOGICAL      ENVYN
      INTEGER      INDEX1
      EXTERNAL     PROMPTFFILE, DSCGRID
     &             STR2INT, STR2REAL, INDEX1

C.........  LOCAL VARIABLES and their descriptions:
      INTEGER       IOS            ! i/o status
      INTEGER    :: LDEV           ! unit no. for log file
      INTEGER       ITYPE
      CHARACTER*256    MESG        ! message buffer

      LOGICAL       RUN_EFS
      LOGICAL       RUN_CANTYP
      LOGICAL       RUN_LAI
      LOGICAL       RUN_W126
      LOGICAL       RUN_LDF
!Jiang
      LOGICAL       RUN_FERT
      LOGICAL       RUN_ARID
      LOGICAL       RUN_NONARID
      LOGICAL       RUN_NITROGEN
      LOGICAL       RUN_LANDTYPE

      CHARACTER*16  :: EFSTXTF = 'EFSTXTF'   ! Input EF file logical name
      CHARACTER*16  :: CTTXTF = 'CTTXTF'     ! Input canopy type file logical name
      CHARACTER*16  :: LAITXTF = 'LAITXTF'   ! Input LAI file logical name
      CHARACTER*16  :: W126TXTF = 'W126TXTF'   ! Input W126 file logical name
      CHARACTER*16  :: LDFTXTF = 'LDFTXTF'   ! Input W126 file logical name
      CHARACTER*16  :: EFMAPS = 'EFMAPS'     ! Output EF file logical name
      CHARACTER*16  :: CANTYP = 'CANTYP'     ! Output canopy type file logical name
      CHARACTER*16  :: LAIS46 = 'LAIS46'     ! Output LAI file logical name
      CHARACTER*16  :: W126FILE = 'W126FILE'     ! Output W126 file logical name
      CHARACTER*16  :: LDFILE = 'LDFILE'     ! Output LDF file logical name
!Jiang
      CHARACTER*16  :: FERTTXTF = 'FERTTXTF'
      CHARACTER*16  :: FERTRESFILE = 'FERTRESFILE'     ! Output FERT file logical name
      CHARACTER*16  :: ARIDTXTF = 'ARIDTXTF'
      CHARACTER*16  :: ARIDFILE = 'ARIDFILE'     ! Output ARID file logical name
      CHARACTER*16  :: NONARIDTXTF = 'NONARIDTXTF'
      CHARACTER*16  :: NONARIDFILE = 'NONARIDFILE'     ! Output NONARID file logical name
      CHARACTER*16  :: NDEPTXTF = 'NDEPTXTF'
      CHARACTER*16  :: NDEPFILE = 'NDEPFILE'     ! Output NITROGEN file logical name
      CHARACTER*16  :: LANDTYPETXTF = 'LANDTYPETXTF' 
      CHARACTER*16  :: LANDTYPEFILE = 'LANDTYPEFILE'    ! Output LANDTYPE file logical name

!     CHARACTER (LEN=14)   :: LANDTYPETXTF(24)
!     CHARACTER (LEN=14)   :: LANDTYPEFILE(24)

      CHARACTER*16  :: GDNAM

      CHARACTER*16  :: PROGNAME = 'TXT2IOAPI'

!***********************************************************************
!... Begin program
!... Initialize log file unit
      LDEV = INIT3() 

!... Get output file grid parameters
      MESG = 'Coordinate name: '
      CALL ENVSTR( 'GDNAM3D', MESG, 'ASACA36km', GDNAM, IOS )

!... Convert EFS
      RUN_EFS = ENVYN( 'RUN_EFS','Run TXT to IOAPI for EFS',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_EFS'
      ENDIF
      IF (RUN_EFS) THEN
        CALL CNVT_EFS(GDNAM,EFSTXTF,EFMAPS)
      ENDIF

!... Convert Canopy type
      RUN_CANTYP = ENVYN( 'RUN_CANTYP','Run TXT to IOAPI for canopy
     &                              types',.TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_CANTYP'
      ENDIF
      IF (RUN_CANTYP) THEN
        CALL CNVT_CANTYP(GDNAM,CTTXTF,CANTYP)
      ENDIF

!... Convert LAI
      RUN_LAI = ENVYN( 'RUN_LAI','Run TXT to IOAPI for LAI',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_LAI'
      ENDIF
      IF (RUN_LAI) THEN
        CALL CNVT_LAI(GDNAM,LAITXTF,LAIS46)
      ENDIF

!... Convert W126
      RUN_W126 = ENVYN( 'RUN_W126','Run TXT to IOAPI for W126',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_W126'
      ENDIF
      IF (RUN_W126) THEN
        CALL CNVT_W126(GDNAM,W126TXTF,W126FILE)
      ENDIF

!... Convert LDF
      RUN_LDF = ENVYN( 'RUN_LDF','Run TXT to IOAPI for LDF',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_LDF'
      ENDIF
      IF (RUN_LDF) THEN
        CALL CNVT_LDF(GDNAM,LDFTXTF,LDFILE)
      ENDIF

!Jiang
!... Convert FERT
      RUN_FERT = ENVYN( 'RUN_FERT','Run TXT to IOAPI for FERT',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_FERT'
      ENDIF
      IF (RUN_FERT) THEN
        CALL CNVT_FERT(GDNAM,FERTTXTF,FERTRESFILE)
      ENDIF

!... Convert ARID
      RUN_ARID = ENVYN( 'RUN_ARID','Run TXT to IOAPI for ARID',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_ARID'
      ENDIF
      IF (RUN_ARID) THEN
        CALL CNVT_ARID(GDNAM,ARIDTXTF,ARIDFILE)
      ENDIF

!... Convert NONARID
      RUN_NONARID = ENVYN( 'RUN_NONARID','Run TXT to IOAPI for NONARID',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_NONARID'
      ENDIF
      IF (RUN_NONARID) THEN
        CALL CNVT_NONARID(GDNAM,NONARIDTXTF,NONARIDFILE)
      ENDIF

!... Convert NITROGEN
      RUN_NITROGEN = ENVYN( 'RUN_NITROGEN','Run TXT to IOAPI for NITROGEN',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_NITROGEN'
      ENDIF
      IF (RUN_NITROGEN) THEN
        CALL CNVT_NITROGEN(GDNAM,NDEPTXTF,NDEPFILE)
      ENDIF

!... Convert LANDTYPE
      RUN_LANDTYPE = ENVYN( 'RUN_LANDTYPE','Run TXT to IOAPI for LANDTYPE',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_LANDTYPE'
      ENDIF
      IF (RUN_LANDTYPE) THEN

!Jiang
         CALL CNVT_LANDTYPE(GDNAM,LANDTYPETXTF,LANDTYPEFILE)

      ENDIF

!... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

      ENDPROGRAM
