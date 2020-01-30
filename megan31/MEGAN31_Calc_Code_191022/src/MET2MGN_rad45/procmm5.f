      subroutine procmm5 (sdate,stime,edate,etime) 
c-----
c     procmm5 processes multiple MM5 files and extracts temperature and
c     PAR (SRAD/2)
c-----
      include 'IODECL3.EXT'
      include 'PARMS3.EXT'
      include 'FDESC3.EXT'

      include 'param.inc'
      include 'fields.inc'
   
      parameter(mxx=300,mxy=300,mxz=30)
      integer kz1(mxz),kz2(mxz),kzin(mxz)
      real tsurf(mxx,mxy),uwind(mxx,mxy,0:mxz),vwind(mxx,mxy,mxz)
      real height(mxx,mxy,mxz),press(mxx,mxy,mxz),wind
      real temp(mxx,mxy,mxz),dummy(mxx,mxy),rad(mxx,mxy)
      real temp1,press1,uwind1,vwind1,fsurf(mxx,mxy,mxz)
      real deltaz,temp0,press0,z0i 
      integer i,j,k,l,m
      integer nx,ny,nz,jdate,jtime,idate,month,itzon
      integer numdays,iunit,ierr,sdate,edate,stime,etime,hourn
      integer odate,ohour,istat,ioffset,joffset,dtout
      integer jhr,mm5hr,mm5date,nxc,nyc,nzc
      integer addday,subday,itime,zhr,zdate,izon,logunit
      real xorg,yorg,dx,dy,tr,rgrndz,clonin,clatin,tlat1in,tlat2in
      real dxcamx,dycamx,x0camx,y0camx,sumf,sumff 
      real phic,xlonc,tlat1,tlat2,alat,alon,xloc,yloc,deltax 
      real cellon(mxx,mxy),cellat(mxx,mxy)
      real, allocatable :: temp_final(:,:,:)
      real, allocatable :: par_final(:,:,:)
      character*80  runmsg,pgname
      character*16  TMCIP,MM5file(10),gnamex
      character*256 MESG
      character*10 project,kvmeth
      character*16  :: ONAME = 'OUTFILE'
      integer :: ENVINT,IOS,numMM5
c
      logical lstagw,lfirst,ltopo
      logical DSCGRID,MM5_met,MCIP_met,MM5c_met,ENVYN
      logical SAT_par,MCIP_rad,MM5_rad
      data lfirst /.true./
      data ltopo /.true./
      data project /'LCP       '/
      data kvmeth /'OB70       '/
      
      data pgname /'met2mgn'/
      data runmsg /'Met PAR and Temp bin2ioapi converter'/
      character*256 CAMXTP1,CAMXZP1,CAMXUV1,fname,mm5in
      character*256 CAMXTP2,CAMXZP2,CAMXUV2,CAMXLU
c
c      character*100 ifile,ifile_par,ifile_mcip
c
c-----Surface roughness (m) as a function of 11 landuse categories
c     and 5 seasons; based on AERMET model (ref EPA SCRAM website)
c
      real z0lu(11,5)
      data z0lu
     & /1.0,0.20,0.100,1.3,1.3,1.30,0.0001,0.002,0.20,0.150,0.30,
     &  1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,
     &  1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,
     &  1.0,0.01,0.001,0.5,1.3,0.90,0.0001,0.002,0.05,0.006,0.15,
     &  1.0,0.03,0.050,1.0,1.3,1.15,0.0001,0.002,0.20,0.040,0.30/
c
c-----Season indices by month and latitude band
c     Season Indices            Latitude Bands
c     1 = summer                1 = <20    Tropical
c     2 = autumn                2 = 20-35  Sub-tropical
c     3 = winter w/o snow       3 = 35-50  Temperate
c     4 = winter w/ snow        4 = 50-75  Cool
c     5 = spring                5 = >75    Polar
c                    Latitude Band
      integer iseason(5,12)
      data iseason / 1, 3, 3, 3, 3, ! Jan
     &               1, 5, 3, 3, 3, ! Feb
     &               1, 5, 5, 3, 3, ! Mar
     &               1, 5, 5, 5, 3, ! Apr
     &               1, 1, 5, 5, 3, ! May
     &               1, 1, 1, 1, 5, ! Jun
     &               1, 1, 1, 1, 1, ! Jul
     &               1, 1, 1, 1, 2, ! Aug
     &               1, 1, 2, 2, 3, ! Sep
     &               1, 2, 2, 2, 3, ! Oct
     &               1, 2, 2, 3, 3, ! Nov
     &               1, 2, 3, 3, 3/ ! Dec
        MM5file(1)  = 'MM5file1'
        MM5file(2)  = 'MM5file2'
        MM5file(3)  = 'MM5file3'
        MM5file(4)  = 'MM5file4'
        MM5file(5)  = 'MM5file5'
        MM5file(6)  = 'MM5file6'
        MM5file(7)  = 'MM5file7'
        MM5file(8)  = 'MM5file8'
        MM5file(9)  = 'MM5file9'
        MM5file(10) = 'MM5file10'
c
      MM5_rad  = ENVYN ( 'MM5RAD', MESG, .FALSE., IOS )
c
c-----
c-----Initialize variables
c
      jdate = 0
      jhr = 0
      ioffset = -1
      joffset = -1
      dtout = 60

      phic=YCENT3D
      xlonc=XCENT3D
      tlat1=P_ALP3D
      tlat2=P_BET3D
      nxc=NCOLS3D
      nyc=NROWS3D
      nzc=2
      xorg=XORIG3D/1000.
      yorg=YORIG3D/1000.
      dx=XCELL3D/1000.
      dy=YCELL3D/1000.
      nlays3d=1

      mxrec3d = 25
      ftype3d = grdded3
      tstep3d = 10000

      fdesc3d(1) = runmsg
      gdtyp3d = lamgrd3
      gdnam3d = gnamex
      vgtyp3d = 2
      vgtop3d = 100.0        !millibars

c-----adjust starting date and time

      zhr=stime         !Initialize GMT hr
      zdate=sdate       !Initialize GMT date

      allocate (temp_final(nxc,nyc,25))
      allocate (par_final(nxc,nyc,25))

c-----Get lat/lon coords for each grid cell

        do j = 1,nyc
          yloc = yorg + dy*(float(j) - 0.5)
          do  i = 1,nxc
            xloc = xorg + dx*(float(i) - 0.5)
            call lcpgeo(1,phic,xlonc,tlat1,tlat2,xloc,yloc,alon,alat)
           
            cellon(i,j) = alon
            cellat(i,j) = alat
          enddo
        enddo

 70   kzin(1) = 1
      kzin(2) = 2
      do k = 1,2
        kz2(k) = kzin(k)
        if (k.eq.1) then
          kz1(k) = 1
        else
          kz1(k) = kz2(k-1) + 1
        endif
      enddo

        dxcamx=dx
        dycamx=dy
        x0camx=xorg
        y0camx=yorg
        clonin=xlonc
        clatin=phic
        tlat1in=tlat1
        tlat2in=tlat2
        itzon=0 
        numMM5 = ENVINT('numMM5','number of MM5 files',1,istat)
        do 200 nf=1,numMM5
           mm5in = MM5file(nf)
           call GETENV(mm5in,fname)
           write(*,*) fname 
        iunit = 20 + nf
        ierr = 0
        open(unit=iunit,file=fname,status='old',form='unformatted')
        write(*,*)
        write(*,*)'Opened input MM5 file: ',fname
c
c-----Read raw MM5 data and convert to CAMx index/units convention
c
 100    call readmm5(iunit,lfirst,project,nxc,nyc,nzc,
     &               ioffset,joffset,dxcamx,dycamx,x0camx,y0camx,izone,
     &               clonin,clatin,tlat1in,tlat2in,mm5date,mm5hr,itzon,
     &               dtout,nx,ny,nz,deltax,ierr)
c
        if (ierr.eq.1) goto 200

        if (mm5date.lt.jdate .or.
     &     (mm5date.eq.jdate .and. mm5hr.le.jhr)) goto 100
        jdate = mm5date
        jhr = mm5hr
        hr = float(jhr)

        if (jdate.lt.sdate .or.
     &     (jdate.eq.sdate .and. jhr.lt.stime)) goto 100
          lstagw = .true.
          call interp_lcp(nx,ny,nz,nxc,nyc,nzc,kz1,kz2,
     &                    ioffset,joffset,deltax,dxcamx,kvmeth)

c---------
 500  itime=hr+(itzon*100)
      idate=jdate
      if (itime.gt.2300) then
        itime=itime-2400
        idate=addday(idate)
      endif

      odate=jdate+2000000
      ohour=hr*100
      write(*,800) odate, ohour

      if(itime.ne.zhr.or.idate.ne.zdate) goto 100  
c-------
c   Determine roughness length, based on season and lat
c-------
      call caldate(idate)
      month = (idate - 10000*int(idate/10000.))/100  
      hourn=(zhr/100)+1

      do j= 1,nyc
        do i = 1,nxc
          temp_final(i,j,hourn)=tac(i,j,1)
          par_final(i,j,hourn)=rground(i,j)/2.
        enddo
      enddo

c-------  
c     Write out temp (and optionally rad) to netCDF file
c
c-------
      sdate3d = 2000000+zdate  !current UTC date
      stime3d = zhr*100        !current UTC time
      xcell3d = 1000*dx
      ycell3d = 1000*dy
      xcent3d = xlonc
      ycent3d = phic

c----write temp to outfile
       if(.not. open3(ONAME,fsunkn3, pgname)) then
            call m3err( 'met2mgn', sdate3d, stime3d,
     &           'Could not open or create '//ONAME//' file',.TRUE.)
       endif

       if(.not.write3(ONAME,vname3d(1),sdate3d,stime3d,
     &         temp_final(:,:,hourn))) then
            call m3err( 'met2mgn', sdate3d, stime3d,
     &         'Could not write '//ONAME//' file', .TRUE.)
       end if

       if (MM5_rad) then
          if(.not.write3(ONAME,vname3d(2),sdate3d,stime3d,
     &           par_final(:,:,hourn))) then
            call m3err( 'met2mgn', sdate3d, stime3d,
     &           'Could not write '//ONAME//' file', .TRUE.)
          end if
       endif

c-----increment hour/check end date/time
      if(zhr.eq.etime.and.zdate.eq.edate) goto 700  

      zhr=zhr+100
      if(zhr.gt.2300) then
         zhr=zhr-2400
         zdate=addday(zdate)
      endif

      goto 100
 200  continue    

 700  deallocate (temp_final)
      deallocate (par_final)
      write(*,*) 'Reach end date/hour; run complete'

      stop 
 800  format(5x,'Reading met for ',12x,i9.7,':',i6.6)
c 900  format(' nx ny nz dx dy ',3i5,2f6.2)
      
      end
