      program met2mgn 
!-----
!     met2mgn produces a Models3-IO/API file with various meteorological 
!     variables to be read into the MEGAN biogenics emissions model. 
!-----
      implicit none

      include 'IODECL3.EXT'
      include 'PARMS3.EXT'
      include 'FDESC3.EXT'

      include 'param.inc'
      include 'fields.inc'
  
      integer idate,sdate,edate,stime,etime,r_sdate,r_stime,epsdate
      integer istat,m1,p1
      integer addday,subday,logunit
      character*16  CNAME
      character*80  runmsg
      character*256 MESG
      logical DSCGRID,MM5_met,WRF_met,MCIP_met,ENVYN
      logical SAT_par,MCIP_rad,MM5_rad,WRF_rad,SOICRO

      character*16, parameter  :: ONAME = 'OUTFILE'
      character*16, parameter  :: pgname = 'met2mgn'
      integer ENVINT,IOS

      data runmsg /'MEGAN meteorology preprocessor'/

      logunit=init3()
      MESG = 'Coordinate name: '
      CALL ENVSTR( 'GDNAM3D', MESG, 'RPO_36km', GDNAM3D, IOS )
      IF( .NOT. DSCGRID( GDNAM3D, CNAME, GDTYP3D,
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( 'met2mgn', 0, 0, MESG, 2 )
      ENDIF
      gdtyp3d=GDTYP3D    

      if (gdtyp3d.ne.2.and.gdtyp3d.ne.5.and.gdtyp3d.ne.1) then
        MESG = 'Only LCP, UTM, and LATLON supported at this time.'
        CALL M3EXIT ( 'met2mgn', 0, 0, MESG, 2 )
      endif

      nlays3d=1
      ftype3d = grdded3
      tstep3d = 10000

      fdesc3d(1) = runmsg
      gdnam3d = GDNAM3D
      vgtyp3d = 2
      vgtop3d = 100.0        !millibars

!-----Get env variables and initialize
!     
      MCIP_met = ENVYN ( 'MCIPMET', MESG, .FALSE., IOS ) 
      MM5_met  = ENVYN ( 'MM5MET', MESG, .FALSE., IOS ) 
!      WRF_met  = ENVYN ( 'WRFMET', MESG, .FALSE., IOS ) 
      WRF_met  = .false.

      SAT_par  = ENVYN ( 'SATPAR', MESG, .FALSE., IOS ) 
      MCIP_rad = ENVYN ( 'MCIPRAD', MESG, .FALSE., IOS ) 
      MM5_rad  = ENVYN ( 'MM5RAD', MESG, .FALSE., IOS ) 
      SOICRO = ENVYN ( 'SOICRO_YN', MESG, .FALSE., IOS ) 
!      WRF_rad  = ENVYN ( 'WRFRAD', MESG, .FALSE., IOS ) 
      WRF_rad  = .false.
!
!-----Error handling for environment variable selection
!
      if(.not.MCIP_met .and. .not.MM5_met
     &   .and. .not.WRF_met) then
        MESG = 'One meteorology input option must be selected' 
        call M3EXIT('met2mgn',0,0,MESG,2)
      endif 

      m1=0
      if(MCIP_met) m1=m1+1
      if(MM5_met)  m1=m1+1
      if(WRF_met)  m1=m1+1

      if(m1.ge.2) then 
        MESG = 'Only one meteorology input option can be selected'
        call M3EXIT('met2mgn',0,0,MESG,2)
      endif

      if(.not.SAT_par .and. .not.MCIP_rad .and. .not.MM5_rad
     &   .and. .not.WRF_rad) then
        MESG = 'One PAR/RAD option must be selected'
        call M3EXIT('met2mgn',0,0,MESG,2)
      endif 

      p1=0
      if(SAT_par)  p1=p1+1
      if(MCIP_rad) p1=p1+1
      if(MM5_rad)  p1=p1+1
      if(WRF_rad)  p1=p1+1

      if(p1.ge.2) then 
        MESG = 'Only one PAR/RAD option can be selected'
        call M3EXIT('met2mgn',0,0,MESG,2)
      endif

      sdate = ENVINT('STDATE','Output start date',0,istat)
      if (istat.gt.0) then
        MESG = 'Bad value for STDATE'
        call M3EXIT('met2mgn',0,0,MESG,2)
      endif

      edate = ENVINT('ENDATE','Output end date',0,istat)
      if (istat.gt.0) then
        MESG = 'Bad value for ENDATE'
        call M3EXIT('met2mgn',0,0,MESG,2)
      endif

      epsdate = ENVINT('EPISODE_SDATE','Episode Start Date',0,istat)
      if (istat.gt.0) then
        MESG = 'Bad value for EPISODE_SDATE'
        call M3EXIT('met2mgn',0,0,MESG,2)
      endif

      stime=mod(sdate,100)
      etime=mod(edate,100)
      sdate=sdate/100
      edate=edate/100
      epsdate=epsdate-2000000

!---Calculate rain start date/time (23 hours before target start time)
      if (epsdate.eq.sdate) then
        r_stime=stime
        r_sdate=sdate
      else
        r_stime=stime+1
        r_sdate=subday(sdate)
      endif

      if(etime.gt.23) then
         etime=etime-24
         edate=addday(edate)
      endif
      if(stime.gt.23) then
         stime=stime-24
         sdate=addday(sdate)
      endif
      if(r_stime.gt.23) then
         r_stime=r_stime-24
         r_sdate=addday(r_sdate)
      endif

      stime=100*stime
      etime=100*etime
      r_stime=100*r_stime

      if (MCIP_met) then
        call readmcip(sdate,stime,edate,etime,r_sdate,r_stime)
      elseif (MM5_met) then 
        call procmm5(sdate,stime,edate,etime) 
      endif

      write(*,*) 'Normal completion by met2mgn'

      stop 
      end
!
!-----Date functions
!
      integer function addday(idate)
      implicit none
      integer idate,iyr,idy
      iyr = idate/1000
      idy = idate - iyr*1000
      if ((mod(iyr,4).eq.0 .and. idy.eq.366) .or.
     &    (mod(iyr,4).ne.0 .and. idy.eq.365)) then
        iyr = iyr + 1
        if (iyr.gt.99) iyr = 0
        addday = iyr*1000 + 1
      else
        addday = idate + 1
      endif
      end
!
      integer function subday(idate)
      implicit none
      integer idate,iyr,idy
      iyr = idate/1000
      idy = idate - iyr*1000
      if (idy.eq.1) then
        iyr = iyr - 1
        if (iyr.lt.0) iyr = 99
        if (mod(iyr,4).eq.0) then
          idy = 366
        else
          idy = 365
        endif
        subday = iyr*1000 + idy
      else
        subday = idate - 1
      endif
      end
