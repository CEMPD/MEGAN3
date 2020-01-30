c     Read PAR satellite data (W/m**2) from GCIP/SRB little-endian 
c     binary formatted files
c
c     Jeremiah Johnson 
c     ENVIRON
c
cccccccccccccccccccccccccccccccccccccc

      subroutine readpar(sdate,stime,edate,etime)

      include 'IODECL3.EXT'
      include 'PARMS3.EXT'
      include 'FDESC3.EXT' 

      parameter   nlat=61,nlon=121
      integer     year,jday,day,month,idate,hour,irec
      integer     lrec,outhr,lat,lon,nx,ny
      integer     num_MCIP,nf,IOS,daynum,istat,idx
      integer     jdate,jtime,sdate,stime,edate,etime,ejdate,ejtime
      integer :: ENVINT,numSATPAR

      real dfdt1,dfdt2,dfdt22,dfdt23,avg3,xorg,yorg
      real lat1,lon1,reslat,reslon,alon,alat
      real dx,dy,phic,xlonc,tlat1,tlat2
      real fluxh(nlon,nlat,24)
      real glat(nlat),glon(nlon)
      data reslat,reslon,lat1,lon1 /0.5,0.5,24.0,-126.0/

      character*256 ifile_par1,ifile_par2,MESG
      character*16, parameter :: ONAME = 'OUTFILE'
 
      real, allocatable::  frad(:,:,:)
      real, allocatable::  frad1(:,:,:)
      real, allocatable::  frad2(:,:,:)
      real, allocatable::  radt(:,:,:)
      real, allocatable::  frad_fin(:,:,:)
      real, allocatable::  flux(:,:)
      real, allocatable::  xg(:,:)
      real, allocatable::  yg(:,:)

      lrec = 4*nlon

c-----get variables from output file
      nx = ncols3d 
      ny = nrows3d

      allocate (frad(nx,ny,25))
      allocate (frad1(nx,ny,25))
      allocate (frad2(nx,ny,25))
      allocate (radt(nx,ny,24))
      allocate (flux(nx,ny))
      allocate (xg(nx,ny))
      allocate (yg(nx,ny))

      if (gdtyp3d.eq.2) then
        phic=  ycent3d
        xlonc= p_gam3d
        tlat1= p_alp3d
        tlat2= p_bet3d
        xorg=  xorig3d/1000
        yorg=  yorig3d/1000
        dx=    xcell3d/1000
        dy=    ycell3d/1000
      elseif (gdtyp3d.eq.5) then
        iutmzon= p_alp3d
        utm_fe= p_bet3d
        utm_fn= p_gam3d
        xorg=  xorig3d/1000
        yorg=  yorig3d/1000
        dx=    xcell3d/1000
        dy=    ycell3d/1000
      elseif (gdtyp3d.eq.1) then
        xorg=  xorig3d
        yorg=  yorig3d
        dx=    xcell3d
        dy=    ycell3d
      endif
c----
       istime = stime/100
       ietime = etime/100
       if (ietime.eq.0) then
         ietime = ietime+24
       endif
       jdate = 2000000+sdate
       jtime = stime*100
       ejdate = 2000000+edate
       ejtime = etime*100
     
       numSATPAR = ENVINT('numSATPAR','number of PAR files',1,istat)
 
c----open input PAR files 

       call GETENV('SATPARFILE1',ifile_par1)               
       write(*,*) 'PAR input file 1: ',ifile_par1

      open(8,file=ifile_par1,status='old',form='unformatted',
     &     access='direct',recl=lrec,err=997)

       if (numSATPAR.eq.2) then
         call GETENV('SATPARFILE2',ifile_par2)               
         write(*,*) 'PAR input file 2: ',ifile_par2

      open(9,file=ifile_par2,status='old',form='unformatted',
     &     access='direct',recl=lrec,err=997)
       endif

      allocate (frad_fin(nx, ny, 25))
c
c --- Calculate latitude and longitude coordinates
c     of GCIP/SRB dat grid cell centers ---
c
      do lat = 1,nlat
        glat(lat) = (lat - 1)*reslat + lat1
      enddo
      do lon = 1,nlon
        glon(lon) = (lon - 1)*reslon + lon1
      enddo

        do j = 1,ny
          yloc = yorg + dy*(float(j) - 0.5)
          do  i = 1,nx
           xloc = xorg + dx*(float(i) - 0.5)
           if (gdtyp3d.eq.2) then  
             call lcpgeo(1,phic,xlonc,tlat1,tlat2,xloc,yloc,alon,alat)
           elseif (gdtyp3d.eq.5) then
             call utmgeo(1,iutmzon,xloc,yloc,alon,alat)
           elseif (gdtyp3d.eq.1) then
             alon=xloc
             alat=yloc
           endif 
           xg(i,j) = alon
           yg(i,j) = alat
c
c --- Check coordinates of modeling cell centers
c     Trap out of bounds coordinates to edge of GCIP/SRB grid ---
c
          if (xg(i,j).le.lon1) then
            xg(i,j) = lon1 + 0.01*reslon
          elseif (xg(i,j).ge.(nlon-1)*reslon+lon1) then
            xg(i,j) = (nlon-1)*reslon+lon1 - 0.01*reslon
          endif
          if (yg(i,j).le.lat1) then
            yg(i,j) = lat1 + 0.01*reslat
          elseif (yg(i,j).ge.(nlat-1)*reslat+lat1) then
            yg(i,j) = (nlat-1)*reslat+lat1 - 0.01*reslat
          endif
        enddo
      enddo

      idx=0
      daynum = 1
 700  if (daynum.eq.1) sdate=sdate-1

      if (daynum.eq.2.and.numSATPAR.eq.2) idx=1
c
c --- Read the GCIP/SRB data file ---
c
      idate = sdate
      call caldate(idate)
      month = (idate - 10000*int(idate/10000.))/100
      day = (idate - 10000*int(idate/10000.))-(month*100)

      write(*,*)
      write(*,*) '    Reading GCIP/SRB file'
      write(*,*)

      irec = 24*nlat*(day-1)
      do hour = 1,24
        do lat = 1,nlat
          irec = irec + 1
          read(8+idx,rec=irec,err=995) (fluxh(lon,lat,hour),lon=1,nlon)
        enddo
      enddo
c
c --- Loop over hours ---
c
 400  do 300 n = 1,24
        write(*,800) 2000000+sdate,n*10000
c
c --- Find the GCIP/SRB data to interpolate ---
c
        do 200 j = 1,ny
          do 200 i = 1,nx
            do l = 1,nlon
              lon = l
              if (glon(lon).gt.xg(i,j)) goto 101
            enddo
            write(*,*) 'ERROR: accessing data west of GCIP/SRB grid'
c
 101        do l = 1,nlat
              lat = l
              if (glat(lat).gt.yg(i,j)) goto 102
            enddo
            write(*,*) 'ERROR: accessing data north of GCIP/SRB grid'
            stop
c
 102        if (lat.lt.2) then
             write(*,*) 'ERROR: accessing data south of GCIP/SRB grid'
              stop
            endif
            if (lon.lt.2) then
              write(*,*) 'ERROR: accessing data west of GCIP/SRB grid'
              stop
            endif
c
c --- Adjust for time zones used in GCIP/SRB data ---
c
            iadj1 = nint(glon(lon-1)/15.0) + 8
            iadj2 = nint(glon(lon)/15.0)   + 8
            n1 = n + iadj1
            n2 = n + iadj2
c
c --- Skip early/late hours where a time zone
c     adjustment spills the GCIP/SRB data arrays ---
c
            if (n1.lt.1.or.n2.lt.1.or.n1.gt.24.or.n2.gt.24) then
              flux(i,j) = -999.
              goto 200
            endif
c
            if (fluxh(lon-1,lat-1,n1).eq.-999. .or.
     &          fluxh(lon,lat-1,n2)  .eq.-999. .or.
     &          fluxh(lon-1,lat,n1)  .eq.-999. .or.
     &          fluxh(lon,lat,n2)    .eq.-999.) then
              flux(i,j) = -999.
            else
              dfdx1 = (fluxh(lon,lat-1,n2) - fluxh(lon-1,lat-1,n1))/
     &                reslon
              dfdx2 = (fluxh(lon,lat,n2) - fluxh(lon-1,lat,n1))/reslon
              f1 = fluxh(lon-1,lat-1,n1) + dfdx1*(xg(i,j) - glon(lon-1))
              f2 = fluxh(lon-1,lat,n1) + dfdx2*(xg(i,j) - glon(lon-1))
              dfdy = (f2 - f1)/reslat
              flux(i,j) = f1 + dfdy*(yg(i,j) - glat(lat-1))
            endif
 200    continue
c
c --- Write gridded PAR to netcdf format file
c     replacing -999 codes with zero ---
c
        do i = 1,nx
          do j = 1,ny
               radt(i,j,n) = amax1(0.0,flux(i,j))
          enddo
        enddo
 300  enddo  ! end hr loop

c-----Interpolate hourly averages to instantaneous rad values

      do i=1,nx
        do j=1,ny

c-----Interpolate rad at 0h and 24h using trends

          avg3=(radt(i,j,3)-radt(i,j,1))/2.
          if(avg3.gt.0) then
             frad(i,j,1)=radt(i,j,1)-(avg3/2.)
          else
             frad(i,j,1)=radt(i,j,1)+(avg3/2.)
          endif
          frad(i,j,1)=amax1(0.0,frad(i,j,1))

          avg22=(radt(i,j,22)-radt(i,j,24))/2.
         if(avg22.gt.0) then
             frad(i,j,25)=radt(i,j,24)-(avg22/2.)
          else
             frad(i,j,25)=radt(i,j,24)+(avg22/2.)
          endif
          frad(i,j,25)=amax1(0.0,frad(i,j,25))
c------
          do outhr=2,24
             frad(i,j,outhr)=(radt(i,j,outhr-1)+
     &                       radt(i,j,outhr))/2.
          enddo
        enddo
      enddo

      if(daynum.eq.1) then

         do i=1,nx
           do j=1,ny
             do n=1,25
               frad1(i,j,n)=frad(i,j,n)
             enddo
           enddo
         enddo
         daynum = 2
         sdate=sdate+1
         goto 700
      else
         do i=1,nx
           do j=1,ny
             do n=1,25
               frad2(i,j,n)=frad(i,j,n)
             enddo
           enddo
         enddo
      endif
c-----
c  check for missing data (valid data exists for hr before and after)
c  and fill if necessary
c----
      do i=1,nx
        do j=1,ny
          do n=17,23
            if (frad1(i,j,n)   .eq. 0 .and.   
     &          frad1(i,j,n-1) .ne. 0 .and.
     &          frad1(i,j,n+1) .ne. 0) then
                frad1(i,j,n) = 0.5*(frad1(i,j,n-1) +
     &                                 frad1(i,j,n+1))
                write(*,*) 'Missing PAR data filled at hr: ', n-16
            endif
          enddo
          do n=2,17
            if (frad2(i,j,n)   .eq. 0 .and.   
     &          frad2(i,j,n-1) .ne. 0 .and.
     &          frad2(i,j,n+1) .ne. 0) then
                frad2(i,j,n) = 0.5*(frad2(i,j,n-1) +
     &                                 frad2(i,j,n+1))
                write(*,*) 'Missing PAR data filled at hr: ', n+8
            endif
          enddo
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          do n=1,8
            n2=16+n
            frad_fin(i,j,n)=frad1(i,j,n2)
          enddo
          do n=9,25
            n2=n-8
            frad_fin(i,j,n)=frad2(i,j,n2)
          enddo
        enddo
      enddo

      sdate3d=jdate
      stime3d=jtime 
      nlays3d=1
      nthik3d=1

      if(.not. open3(ONAME,fsunkn3, pgname)) then
            call m3err( 'readpar', sdate3d, stime3d,
     &        'Could not open or create '//ONAME//' file',.TRUE.)
      else if(.not. desc3(ONAME)) then 
            call m3err( 'readpar', sdate3d, stime3d,
     &        'Could not get description for '//ONAME//' file',.TRUE.)
      endif
c      vname3d(nvars3d) =  'PAR             '
c      vtype3d(nvars3d) =  m3real
c      units3d(nvars3d) =  'WATTS/M**2      '
c      vdesc3d(nvars3d) =  'Photosynthetically Active Radiation'

c-----write sat PAR to outfile
       do n=istime+1,25
          if(.not.write3(ONAME,vname3d(nvars3d),jdate,jtime,
     &         frad_fin(:,:,n))) then
               call m3err( 'readpar.f', jdate, jtime,
     &           'Could not write '//ONAME//' file', .TRUE.)
           endif

          if(jtime.eq.ejtime.and.jdate.eq.ejdate) goto 999        

          call nextime(jdate, jtime, 10000)
       enddo

800   format(5x,'Reading PAR data (LST) for',1x,i9.4,':',i6.6)
995   write(*,*) 'Error reading PAR input file1'
      stop
996   write(*,*) 'Error reading PAR input file2'
      stop
997   write(*,*) 'Error opening PAR input file'
      stop
998   write(*,*) 'Error reading PAR input file'
      stop
999   deallocate (frad_fin) 
      deallocate (frad) 
      deallocate (frad1) 
      deallocate (frad2) 
      deallocate (radt) 
      deallocate (flux) 
      deallocate (xg) 
      deallocate (yg) 
      
      end
