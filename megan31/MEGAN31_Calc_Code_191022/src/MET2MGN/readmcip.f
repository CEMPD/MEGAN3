c     readmcip.f
c
c     Read met variables from MCIP files and write 
c     as variables to the netCDF file 
c
c     Jeremiah Johnson 
c     ENVIRON
c
cccccccccccccccccccccccccccccccccccccc

      subroutine readmcip(sdate,stime,edate,etime,r_sdate,r_stime)

      include 'IODECL3.EXT'
      include 'PARMS3.EXT'
      include 'FDESC3.EXT' 

      integer l,n,nx,ny,IOS,isoilm,isoilt,isoiltyp
      integer jdate,jtime,sdate,stime,edate,etime,ejdate,ejtime
      integer epsdate,r_sdate,r_stime,r_jdate,r_jtime
      integer ENVINT,istat

      character*256 MESG
      character*16, parameter :: ONAME = 'OUTFILE'
      character*16, parameter :: PNAME = 'PFILE'
      character*16, parameter :: METCRO2Dnam1 = 'METCRO2Dfile1'
      character*16, parameter :: METCRO2Dnam2 = 'METCRO2Dfile2'
      character*16, parameter :: METCRO3Dnam = 'METCRO3Dfile'
      character*16, parameter :: METDOT3Dnam = 'METDOT3Dfile'
      character*16, parameter :: SOICROnam1 = 'SOICROfile1'
      character*16, parameter :: SOICROnam2 = 'SOICROfile2'
      character*16, parameter :: pgname = 'TPAR2IOAPI'
      character*16 TMCIP 

      integer, parameter :: MAXSTYPES = 11 

      logical MCIP_rad,SAT_par,ENVYN,INITIAL_HOUR,ifsoil
      logical SOICRO
      real, allocatable::  srad(:,:)
      real, allocatable::  par(:,:)
      real, allocatable::  reftemp(:,:)
      real, allocatable::  pres(:,:)
      real, allocatable::  snocov(:,:)
      real, allocatable::  wspd10(:,:)
      real, allocatable::  cfrac(:,:)
      real, allocatable::  soilm(:,:)
      real, allocatable::  soilt(:,:)
      real, allocatable::  soiltyp(:,:)
      real, allocatable::  qv(:,:)
      real, allocatable::  u_wind(:,:)
      real, allocatable::  v_wind(:,:)
      real, allocatable::  avewind(:,:)
      real, allocatable::  windspd(:,:)
      real, allocatable::  rn(:,:)
      real, allocatable::  rc(:,:)
      real, allocatable::  rain(:,:,:)
      real, allocatable::  rain_acc24(:,:,:)
      real, allocatable::  PTYPE(:,:)
      real, allocatable::  PULSEDATE(:,:)
      real, allocatable::  PULSETIME(:,:)
      real, allocatable::  PREC_ADJ(:,:)
c      real, allocatable::  precip_adjfac(:,:)
      integer :: PULTYPE,PULDATE,PULTIME,soilcat
      real :: precip_adjfac 

      real,    external :: PRECIPFACT
      integer, external :: PULSETYPE


c-----get variables from output file
      nx = ncols3d 
      ny = nrows3d

      istime = stime/100
      ietime = etime/100
      if (ietime.eq.0) then
        ietime = ietime+24
      endif
      jdate = sdate+2000000
      jtime = stime*100
      ejdate = 2000000+edate
      ejtime = etime*100

      r_jdate = r_sdate+2000000
      r_jtime = r_stime*100

      call GETENV ('TMCIP',TMCIP)

      MESG = 'Coordinate name: '
      call ENVSTR( 'GDNAM3D', MESG, 'TX_36km', GDNAM3D, IOS )
      MCIP_rad = ENVYN( 'MCIPRAD', MESG, .FALSE., IOS )
      SAT_par  = ENVYN ( 'SATPAR', MESG, .FALSE., IOS )
      SOICRO  = ENVYN ( 'SOICRO_YN', MESG, .FALSE., IOS )

      epsdate = ENVINT('EPISODE_SDATE','Episode Start Date',0,istat)
      if (istat.gt.0) then
        MESG = 'Bad value for EPISODE_SDATE'
        call M3EXIT('met2mgn',0,0,MESG,2)
      endif

      ifsoil = .true.
                                                       
c-----open MCIP files 

      if(.not. open3(METCRO2Dnam1,fsread3,pgname)) then
        call m3err( 'readmcip', sdate3d, stime3d,
     &    'Could not open or create '//METCRO2Dnam1//' file',.TRUE.)
      else if(.not. desc3(METCRO2Dnam1)) then
        call m3err( 'readmcip', sdate3d, stime3d,
     &   'Could not get description for '//METCRO2Dnam1//' file',.TRUE.)
      endif

      if ( SOICRO ) then
        if(.not. open3(SOICROnam1,fsread3,pgname)) then
          call m3err( 'readmcip', sdate3d, stime3d,
     &    'Could not open or create '//SOICROnam1//' file',.TRUE.)
        else if(.not. desc3(SOICROnam1)) then
          call m3err( 'readmcip', sdate3d, stime3d,
     &   'Could not get description for '//SOICROnam1//' file',.TRUE.)
        endif
      else
        isoilm = INDEX1('SOIM1',NVARS3D,VNAME3D)
        isoilt = INDEX1('SOIT1',NVARS3D,VNAME3D)
        isoiltyp = INDEX1('SLTYP',NVARS3D,VNAME3D)
        if (isoilm.le.0 .or. isoilt.le.0 .or. isoiltyp.le.0) then
          ifsoil = .false.
        endif  
      endif
      
      if(r_sdate.ne.sdate) then
        if(.not. open3(METCRO2Dnam2,fsread3,pgname)) then
          call m3err( 'readmcip', sdate3d, stime3d,
     &      'Could not open or create '//METCRO2Dnam2//' file',.TRUE.)
        endif
        if ( SOICRO ) then
          if(.not. open3(SOICROnam2,fsread3,pgname)) then
          call m3err( 'readmcip', sdate3d, stime3d,
     &      'Could not open or create '//SOICROnam2//' file',.TRUE.)
          endif
        endif

        if(.not. open3(PNAME,fsread3,pgname)) then
          call m3err( 'readmcip', sdate3d, stime3d,
     &      'Could not open or create '//PNAME//' file',.TRUE.)
        endif
      endif

      if(.not. open3(METCRO3Dnam,fsread3,pgname)) then
        call m3err( 'readmcip', sdate3d, stime3d,
     &    'Could not open or create '//METCRO3Dnam//' file',.TRUE.)
      endif

      if(.not. open3(METDOT3Dnam,fsread3,pgname)) then
        call m3err( 'readmcip', sdate3d, stime3d,
     &    'Could not open or create '//METDOT3Dnam//' file',.TRUE.)
      endif

      allocate (reftemp(nx,ny))
      allocate (pres(nx,ny))
      allocate (cfrac(nx,ny))
      allocate (snocov(nx,ny))
      allocate (wspd10(nx,ny))
      allocate (soilm(nx,ny))
      allocate (soilt(nx,ny))
      allocate (soiltyp(nx,ny))
      allocate (qv(nx,ny))
      allocate (u_wind(nx+1,ny+1))
      allocate (v_wind(nx+1,ny+1))
      allocate (avewind(nx+1,ny+1))
      allocate (windspd(nx,ny))
      if(MCIP_rad) then
        allocate (srad(nx,ny))
        allocate (par(nx,ny))
      endif
      allocate (rn(nx,ny))
      allocate (rc(nx,ny))
      allocate (rain(nx,ny,100))
      allocate (rain_acc24(nx,ny,100))
      allocate (PTYPE(nx,ny))
      allocate (PULSEDATE(nx,ny))
      allocate (PULSETIME(nx,ny))
      allocate (PREC_ADJ(nx,ny))

c---Store total precip (RN+RC) in array starting from 23 hours prior to
c   target start time and ending at the target end time   

      do n=1,100
       
c-------Target start time already reached
        if(r_jdate.ge.jdate.and.r_jtime.ge.jtime) then
          if (r_sdate.eq.sdate) then
c
c---First day (only one MCIP file to read in)
c
            if (.not.read3(METCRO2Dnam1,'RC',1,r_jdate,r_jtime,rc)) then
              call m3err( 'readmcip', r_jdate, r_jtime,
     &             'Could not read RC from '//METCRO2Dnam1//' file', .TRUE.)
            endif
            if (.not.read3(METCRO2Dnam1,'RN',1,r_jdate,r_jtime,rn)) then
              call m3err( 'readmcip', r_jdate, r_jtime,
     &             'Could not read RC from '//METCRO2Dnam1//' file', .TRUE.)
            endif
          else 
c
c---Past first day (read in second MCIP file)
c
            if (.not.read3(METCRO2Dnam2,'RC',1,r_jdate,r_jtime,rc)) then
              call m3err( 'readmcip', r_jdate, r_jtime,
     &             'Could not read RC from '//METCRO2Dnam1//' file', .TRUE.)
            endif
            if (.not.read3(METCRO2Dnam2,'RN',1,r_jdate,r_jtime,rn)) then
              call m3err( 'readmcip', r_jdate, r_jtime,
     &             'Could not read RC from '//METCRO2Dnam1//' file', .TRUE.)
            endif
          endif

c-------Before target start time
        else
          if (.not.read3(METCRO2Dnam1,'RC',1,r_jdate,r_jtime,rc)) then
c            call m3err( 'readmcip', r_jdate, r_jtime,
c     &           'Could not read RC at '//METCRO2Dnam1//' file', .TRUE.)
             MESG = 'Could not read RC from '//METCRO2Dnam1//' file'
             call m3warn('readmcip', r_jdate, r_jtime, MESG)
          endif
          if (.not.read3(METCRO2Dnam1,'RN',1,r_jdate,r_jtime,rn)) then
c           call m3err( 'readmcip', r_jdate, r_jtime,
c     &           'Could not read RC at '//METCRO2Dnam1//' file', .TRUE.)
             MESG = 'Could not read RN from '//METCRO2Dnam1//' file'
             call m3warn('readmcip', r_jdate, r_jtime, MESG)
          endif
        endif

        do j=1,ny
        do i=1,nx
          rain(i,j,n)=rn(i,j)+rc(i,j)
        enddo
        enddo

        if(r_jtime.eq.ejtime.and.r_jdate.eq.ejdate) then
          deallocate(rc)
          deallocate(rn)
          goto 700
        endif
        call nextime(r_jdate, r_jtime, 10000)
      enddo

700   do n=1,100

      if(r_sdate.eq.sdate) then  !First day
        if (.not.read3(METCRO2Dnam1,TMCIP,1,jdate,jtime,
     &         reftemp)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO2Dnam1//' file', .TRUE.)
        endif

        if (.not.read3(METCRO2Dnam1,'SNOCOV',1,jdate,jtime,
     &         snocov)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO2Dnam1//' file', .TRUE.)
        endif

        if (.not.read3(METCRO2Dnam1,'CFRAC',1,jdate,jtime,
     &         cfrac)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO2Dnam1//' file', .TRUE.)
        endif

        if (.not.read3(METCRO2Dnam1,'WSPD10',1,jdate,jtime,
     &         wspd10)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO2Dnam1//' file', .TRUE.)
        endif

        if (ifsoil) then 
          if ( SOICRO ) then
            if (.not.read3(SOICROnam1,'SOIM3D',1,jdate,jtime,soilm)) then
              MESG = 'Could not read SOIM3D from '//SOICROnam1//' file'
              call m3mesg(MESG)
              ifsoil = .false.
            endif

            if (.not.read3(SOICROnam1,'SOIT3D',1,jdate,jtime,soilt)) then
              MESG = 'Could not read SOIT3D from '//SOICROnam1//' file'
              call m3mesg(MESG)
              ifsoil = .false.
            endif
          else
            if (.not.read3(METCRO2Dnam1,'SOIM1',1,jdate,jtime,soilm)) then
              MESG = 'Could not read SOIM1 from '//METCRO2Dnam1//' file'
              call m3mesg(MESG)
              ifsoil = .false.
            endif

            if (.not.read3(METCRO2Dnam1,'SOIT1',1,jdate,jtime,soilt)) then
              MESG = 'Could not read SOIT1 from '//METCRO2Dnam1//' file'
              call m3mesg(MESG)
              ifsoil = .false.
            endif
          endif

          if (.not.read3(METCRO2Dnam1,'SLTYP',1,jdate,jtime,soiltyp)) then
              MESG = 'Could not read SLTYP from '//METCRO2Dnam1//' file'
              call m3mesg(MESG)
              ifsoil = .false.
          endif
        endif 

        if(MCIP_rad) then
          if (.not.read3(METCRO2Dnam1,'RGRND',1,jdate,jtime,srad)) then
             call m3err( 'readmcip', jdate, jtime,
     &            'Could not read '//METCRO2Dnam1//' file', .TRUE.)
          endif
          par(:,:)=srad(:,:)/2
        endif

      else  !Not first day
c
c-----If first hour of the day (and not the first day of the episode), read in 
c     pulse variables from current hour (last hour of previous day should 
c     be the same as first hour of current day)  
c
        if((sdate+2000000).eq.jdate.and.(stime*100).eq.jtime) then
          lfirst = .true.
        else
          lfirst = .false.
        endif 
        if (lfirst) then
          if (.not.read3(PNAME,'PTYPE',1,jdate,jtime,
     &           PTYPE)) then
             call m3err( 'readmcip', jdate, jtime,
     &            'Could not read '//PNAME//' file', .TRUE.)
          endif
          if (.not.read3(PNAME,'PULSEDATE',1,jdate,jtime,
     &           PULSEDATE)) then
             call m3err( 'readmcip', jdate, jtime,
     &            'Could not read '//PNAME//' file', .TRUE.)
          endif
          if (.not.read3(PNAME,'PULSETIME',1,jdate,jtime,
     &           PULSETIME)) then
             call m3err( 'readmcip', jdate, jtime,
     &            'Could not read '//PNAME//' file', .TRUE.)
          endif
          if (.not.close3(PNAME)) then
             call m3err( 'readmcip', jdate, jtime,
     &            'Could not clsoe '//PNAME//' file', .TRUE.)
          endif
        endif
        if (.not.read3(METCRO2Dnam2,TMCIP,1,jdate,jtime,
     &         reftemp)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO2Dnam2//' file', .TRUE.)
        endif

        if (.not.read3(METCRO2Dnam2,'SNOCOV',1,jdate,jtime,
     &         snocov)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO2Dnam2//' file', .TRUE.)
        endif


        if (.not.read3(METCRO2Dnam2,'CFRAC',1,jdate,jtime,
     &         cfrac)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO2Dnam2//' file', .TRUE.)
        endif

        if (.not.read3(METCRO2Dnam2,'WSPD10',1,jdate,jtime,
     &         wspd10)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO2Dnam2//' file', .TRUE.)
        endif

        if (ifsoil) then
          if ( SOICRO ) then 
            if (.not.read3(SOICROnam2,'SOIM3D',1,jdate,jtime,soilm)) then
              MESG = 'Could not read SOIM3D from '//SOICROnam2//' file'
              call m3mesg(MESG)
              ifsoil = .false.
            endif

            if (.not.read3(SOICROnam2,'SOIT3D',1,jdate,jtime,soilt)) then
              MESG = 'Could not read SOIMTD from '//SOICROnam2//' file'
              call m3mesg(MESG)
              ifsoil = .false.
            endif
          else ! read soil moisture and soil temperature from METCRO2D
            if (.not.read3(METCRO2Dnam2,'SOIM1',1,jdate,jtime,soilm)) then
              MESG = 'Could not read SOIM1 from '//METCRO2Dnam1//' file'
              call m3mesg(MESG)
              ifsoil = .false.
            endif

            if (.not.read3(METCRO2Dnam2,'SOIT1',1,jdate,jtime,soilt)) then
              MESG = 'Could not read SOIT1 from '//METCRO2Dnam1//' file'
              call m3mesg(MESG)
              ifsoil = .false.
            endif
          endif

          if (.not.read3(METCRO2Dnam2,'SLTYP',1,jdate,jtime,soiltyp)) then
              MESG = 'Could not read SLTYP from '//METCRO2Dnam1//' file'
              call m3mesg(MESG)
              ifsoil = .false.
          endif
        endif

        if(MCIP_rad) then
          if (.not.read3(METCRO2Dnam2,'RGRND',1,jdate,jtime,srad)) then
             call m3err( 'readmcip', jdate, jtime,
     &            'Could not read '//METCRO2Dnam2//' file', .TRUE.)
          endif
          par(:,:)=srad(:,:)/2
        endif
      endif

        if (.not.read3(METCRO3Dnam,'PRES',1,jdate,jtime,pres)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO3Dnam//' file', .TRUE.)
        endif

        if (.not.read3(METCRO3Dnam,'QV',1,jdate,jtime,qv)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METCRO3Dnam//' file', .TRUE.)
        endif

        if (.not.read3(METDOT3Dnam,'UWIND',1,jdate,jtime,u_wind)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METDOT3Dnam//' file', .TRUE.)
        endif

        if (.not.read3(METDOT3Dnam,'VWIND',1,jdate,jtime,v_wind)) then
           call m3err( 'readmcip', jdate, jtime,
     &          'Could not read '//METDOT3Dnam//' file', .TRUE.)
        endif


c----Calculate 24-hr rain accumulation
c
      if(r_sdate.eq.sdate) then
        do l=1,n
          do j=1,ny
          do i=1,nx
            rain_acc24(i,j,n)=rain_acc24(i,j,n)+rain(i,j,l)
          enddo
          enddo
        enddo
      else
        do l=n,n+23,1
          do j=1,ny
          do i=1,nx
            rain_acc24(i,j,n)=rain_acc24(i,j,n)+rain(i,j,l)
          enddo
          enddo
        enddo
      endif
c
c---Calculate adjustment to precip from "pulse" 
c
      if((r_sdate+2000000).eq.jdate.and.(r_stime*100).eq.jtime) then
        INITIAL_HOUR = .TRUE.
      else
        INITIAL_HOUR = .FALSE. 
      endif
c
c-----If initial hour (first hour of first day of episode), 
c     then set pulse variables to zero
c
c      if(r_sdate.eq.sdate) then  !first day
        if(INITIAL_HOUR) then !first day *and* first hour
          do j=1,ny
          do i=1,nx
            PTYPE(i,j)=0.
            PULSEDATE(i,j)=0.
            PULSETIME(i,j)=0.
            soilcat=soiltyp(i,j)
            if(ifsoil) then
              if( soilcat > 0 .and. soilcat <= MAXSTYPES ) then 
                PREC_ADJ(i,j) = 2.
              else
                PREC_ADJ(i,j) = 1.
              endif 
            else
              PREC_ADJ(i,j)=1.
            endif
          enddo
          enddo
        else              !first day but not first hour
c
c-----Calculate adjustment factor for all other hours 
c
          do j=1,ny
          do i=1,nx
            PULTYPE=PTYPE(i,j)
            PULDATE=PULSEDATE(i,j)
            PULTIME=PULSETIME(i,j)
            if(PULSETYPE(rain_acc24(i,j,n)).le.PULTYPE) then
              ! no new rain 
              precip_adjfac=PRECIPFACT(PULTYPE,jdate,jtime,
     &                               PULDATE,PULTIME)
            else 
              ! rain 
              PULDATE = jdate
              PULTIME = jtime
              PULTYPE = PULSETYPE(rain_acc24(i,j,n))
              precip_adjfac=PRECIPFACT(PULTYPE,jdate,jtime,
     &                               PULDATE,PULTIME)
            endif 
            PTYPE(i,j)=PULTYPE 
            PULSEDATE(i,j)=PULDATE
            PULSETIME(i,j)=PULTIME
            PREC_ADJ(i,j)=precip_adjfac
          enddo
          enddo
        endif

c----Interpolate winds at cell corners to cell centers
c
        avewind(:,:)=sqrt(u_wind(:,:)**2+v_wind(:,:)**2)

        do j=1,ny 
        do i=1,nx 
           windspd(i,j)=(avewind(i,j)+avewind(i,j+1)+avewind(i+1,j)
     &                   +avewind(i+1,j+1))/4.
        enddo
        enddo

        sdate3d=sdate+2000000
        stime3d=stime*100
        nlays3d=1
        nthik3d=1
        mxrec3d=0

!-----Initialize and define output variables and attributes                                                                    
      if (ifsoil) then
        nvars3d = 13
      else
        nvars3d = 10
      endif

      if (ifsoil) then
      vname3d(1) =  'SOIM1           '
      vtype3d(1) =  m3real           
      units3d(1) =  'M**3/M**3       '
      vdesc3d(1) =  'volumetric soil moisture in top cm'
                                               
      vname3d(2) =  'SOIT1           '
      vtype3d(2) =  m3real           
      units3d(2) =  'K               '
      vdesc3d(2) =  'soil temperature in top cm'
                                              
      vname3d(3) =  'SLTYP           '
      vtype3d(3) =  m3real           
      units3d(3) =  'CATEGORY        '
      vdesc3d(3) =  'soil texture type by USDA category'
      endif

      vname3d(nvars3d-9) =  'TEMP2           '
      vtype3d(nvars3d-9) =  m3real
      units3d(nvars3d-9) =  'K               '

      if (TMCIP.eq.'TEMP1P5') then
         vdesc3d(nvars3d-9) =  'temperature at 1.5 m'
      else                 
         vdesc3d(nvars3d-9) =  'temperature at 2 m'
      endif
                                    
      vname3d(nvars3d-8) =  'PRES            '
      vtype3d(nvars3d-8) =  m3real           
      units3d(nvars3d-8) =  'Pa              '
      vdesc3d(nvars3d-8) =  'pressure        '
                                     
      vname3d(nvars3d-7) =  'SNOCOV          '
      vtype3d(nvars3d-7) =  m3real
      units3d(nvars3d-7) =  'DECIMAL         '
      vdesc3d(nvars3d-7) =  'snow cover     '

      vname3d(nvars3d-6) =  'CFRAC           '
      vtype3d(nvars3d-6) =  m3real
      units3d(nvars3d-6) =  'FRACTION        '
      vdesc3d(nvars3d-6) =  'total cloud fraction'

      vname3d(nvars3d-5) =  'WSPD10          '
      vtype3d(nvars3d-5) =  m3real
      units3d(nvars3d-5) =  'm/s        '
      vdesc3d(nvars3d-5) =  'wind speed at 10m'

      vname3d(nvars3d-4) =  'QV              '
      vtype3d(nvars3d-4) =  m3real            
      units3d(nvars3d-4) =  'KG/KG           '
      vdesc3d(nvars3d-4) =  'water vapor mixing ratio'
                                     
      vname3d(nvars3d-3) =  'WINDSPD         '
      vtype3d(nvars3d-3) =  m3real           
      units3d(nvars3d-3) =  'm/s             '
      vdesc3d(nvars3d-3) =  'Cell centered Windspeed'

      vname3d(nvars3d-2) =  'RAIN_ACC24       '
      vtype3d(nvars3d-2) =  m3real            
      units3d(nvars3d-2) =  'cm               '
      vdesc3d(nvars3d-2) =  '24-hour accumulated rain'

      vname3d(nvars3d-1) =  'PREC_ADJ         '
      vtype3d(nvars3d-1) =  m3real            
      units3d(nvars3d-1) =  'No dimension     '
      vdesc3d(nvars3d-1) =  'Precip adjustment factor'

      vname3d(nvars3d) =  'PAR             '
      vtype3d(nvars3d) =  m3real           
      units3d(nvars3d) =  'WATTS/M**2      '
      vdesc3d(nvars3d) =  'Photosynthetically Active Radiation'

        if(.not. open3(ONAME,fsunkn3, pgname)) then
            call m3err( 'readmcip', sdate3d, stime3d,
     &           'Could not open or create '//ONAME//' file',.TRUE.)
        endif

        if (ifsoil) then
          if(.not.write3(ONAME,vname3d(1),jdate,jtime,soilm)) then
                call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
          endif
  
          if(.not.write3(ONAME,vname3d(2),jdate,jtime,soilt)) then
                call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
          endif
  
          if(.not.write3(ONAME,vname3d(3),jdate,jtime,soiltyp)) then
                call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
          endif
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-9),jdate,jtime,
     &        reftemp)) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-8),jdate,jtime,pres)) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-7),jdate,jtime,snocov)) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-6),jdate,jtime,cfrac)) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-5),jdate,jtime,wspd10)) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-4),jdate,jtime,qv)) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-3),jdate,jtime,windspd)) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-2),jdate,jtime,
     &                 rain_acc24(:,:,n))) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif

        if(.not.write3(ONAME,vname3d(nvars3d-1),jdate,jtime,
     &                 PREC_ADJ)) then
              call m3err( 'readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
        endif
       
        if(MCIP_rad) then
          if(.not.write3(ONAME,vname3d(nvars3d),jdate,jtime,par)) then
              call m3err('readmcip', jdate, jtime,
     &             'Could not write '//ONAME//' file', .TRUE.)
          
          endif
        endif

c---Write pulse variables to temporary file (PNAME) 
c
        if(jtime.eq.ejtime.and.jdate.eq.ejdate) then 

        if (SAT_par) then 
          call readpar(sdate,stime,edate,etime)
        endif

        nvars3d = 3 
        vname3d(1) =  'PTYPE           '
        vtype3d(1) =  m3real            
        units3d(1) =  'No dimension     '
        vdesc3d(1) =  'type of NO pulse        '
  
        vname3d(2) = 'PULSEDATE        '
        vtype3d(2) =  m3real            
        units3d(2) =  'No dimension     '
        vdesc3d(2) =  'Date of last NO pulse'

        vname3d(3) = 'PULSETIME        '
        vtype3d(3) =  m3real            
        units3d(3) =  'No dimension            '
        vdesc3d(3) =  'Time of last NO pulse'
        sdate3d=jdate
        stime3d=jtime
        mxrec3d=1
c
c--------Write out last hour's pulse vars to intermediate file
c
          if(.not. open3(PNAME,fscrea3, pgname)) then
                call m3err( 'readmcip', sdate3d, stime3d,
     &               'Could not open or create '//PNAME//' file',.TRUE.)
          endif
          if(.not.write3(PNAME,vname3d(1),jdate,jtime,
     &                   PTYPE)) then
                call m3err( 'readmcip', jdate, jtime,
     &               'Could not write '//PNAME//' file', .TRUE.)
          endif
          if(.not.write3(PNAME,vname3d(2),jdate,jtime,
     &                   PULSEDATE)) then
                call m3err( 'readmcip', jdate, jtime,
     &               'Could not write '//PNAME//' file', .TRUE.)
          endif
          if(.not.write3(PNAME,vname3d(3),jdate,jtime,
     &                   PULSETIME)) then
                call m3err( 'readmcip', jdate, jtime,
     &               'Could not write '//PNAME//' file', .TRUE.)
          endif

          deallocate(reftemp)
          deallocate(pres)
          deallocate(wspd10)
          deallocate(snocov)
          deallocate(cfrac)
          deallocate(soilm)
          deallocate(soilt)
          deallocate(soiltyp)
          deallocate(qv)
          deallocate(u_wind)
          deallocate(v_wind)
          deallocate(windspd)
          deallocate(rain_acc24)
          deallocate(PTYPE)
          deallocate(PULSEDATE)
          deallocate(PULSETIME)
          deallocate(PREC_ADJ)
          if(MCIP_rad) then
            deallocate(srad)
            deallocate(par)
          endif
          goto 999 
        endif 
        call nextime(jdate, jtime, 10000)
      enddo
 999  end
      

      REAL FUNCTION PRECIPFACT( PULTYPE, JDATE, JTIME, 
     &                          ADATE, ATIME )

C***********************************************************************
C  DESCRIPTION
C    This internal function computes a precipitation adjustment
C    factor from YL 1995 based on a rain rate. The pulse type is
C    and integer ranging from 0 to 3 indicating the type of 
C    rainfall rate.
C
C  CALL
C    SECSDIFF  -  IOAPI
C
C  HISTORY:
C    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
C***********************************************************************


      IMPLICIT NONE
            
C...  Function arguments
      INTEGER, INTENT(IN OUT) :: PULTYPE
      INTEGER, INTENT(IN) :: JDATE
      INTEGER, INTENT(IN) :: JTIME
      INTEGER, INTENT(IN) :: ADATE
      INTEGER, INTENT(IN) :: ATIME
            
C...  External functions
      INTEGER, EXTERNAL :: SECSDIFF

C...  Local variables
      INTEGER  HRDIFF
            
C-----------------------------------------------------------------------------

      HRDIFF = SECSDIFF( ADATE, ATIME, JDATE, JTIME ) / 3600.
       
      SELECT CASE( PULTYPE )
      CASE( 0 )
          PRECIPFACT = 1.
      CASE( 1 )
          IF( ( HRDIFF / 24. ) < 2. ) THEN
              PRECIPFACT = 11.19 * EXP(-0.805*(HRDIFF+24)/24.)
          ELSE
              PULTYPE = 0
              PRECIPFACT = 1.
          ENDIF
      CASE( 2 )
          IF( ( HRDIFF / 24. ) < 6. ) THEN
              PRECIPFACT = 14.68 * EXP(-0.384*(HRDIFF+24)/24.)
          ELSE
              PULTYPE = 0
              PRECIPFACT = 1.
          ENDIF
      CASE DEFAULT
          IF( ( HRDIFF / 24. ) < 13. ) THEN
              PRECIPFACT = 18.46 * EXP(-0.208*(HRDIFF+24)/24.)
          ELSE
              PULTYPE = 0
              PRECIPFACT = 1.
          ENDIF
      END SELECT

      RETURN

      END FUNCTION PRECIPFACT
!=======================================================================
!=======================================================================
    

!=======================================================================
!=======================================================================
      INTEGER FUNCTION PULSETYPE( RRATE )

C.....This internal function computes the pulse type from a rainfall rate.
C     (See YL 1995).

      IMPLICIT NONE
      
C.....Function arguments
      REAL, INTENT(IN) :: RRATE
            
!-----------------------------------------------------------------------
      IF( RRATE < 0.1 ) THEN
          PULSETYPE = 0
      ELSE IF( RRATE < 0.5 ) THEN
          PULSETYPE = 1
      ELSE IF( RRATE < 1.5 ) THEN
          PULSETYPE = 2
      ELSE
          PULSETYPE = 3
      ENDIF
      RETURN
      END FUNCTION PULSETYPE
