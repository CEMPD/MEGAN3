      program ioapi2uam
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Copyright (C) 2007-2007  ENVIRON
c
c
c   This program is free software; you can redistribute it and/or
c   modify it under the terms of the GNU General Public License
c   as published by the Free Software Foundation; either version 2
c   of the License, or (at your option) any later version.
c
c   This program is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c   GNU General Public License for more details.
c
c   To obtain a copy of the GNU General Public License
c   write to the Free Software Foundation, Inc.,
c   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
c
c
c   For comments and questions, send to bkoo@environcorp.com
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     IOAPI2UAM converts CMAQ 1-D emissions files (I/O API) to CAMx
c     low-level emissions files (UAM-IV). It only converts file format
c     with no species mapping. If a CMAQ variable name is longer than
c     10 characters, it truncates the name. Emission rate is converted
c     from mol/s (or g/s) to mol/hr (or g/hr). It also shifts time-zone
c     from GMT to user-selected local time.
c
c     INPUT ENVIRONMENTAL VARIABLES:
c       INFILE1       - Logical name for input file 1 (current day)
c       INFILE2       - Logical name for input file 2 (next day)
c                       Required only if additional data is needed
c                       due to time zone shifting
c                       Map projection consistency won't be checked
c       OUTFILE       - Logical name for output file
c       TZONE         - Output time zone (8 for PST, etc.)
c       SDATE         - Output start date (YYJJJ)
c       STIME         - Output start time in TZONE (HHMMSS)
c       RLENG         - Output run length (HHMMSS)
c                       240000 for a CAMx daily emissions input
c
c     HISTORY:
c       created by bkoo (06/20/2007)
c
      include 'PARMS3.EXT'
      include 'IODECL3.EXT'
      include 'FDESC3.EXT'

      integer :: LOGUNIT
      integer :: ENVINT,TIME2SEC,JSTEP3
      integer :: JDATE,JTIME,CDATE,CTIME,RLENG
      integer, parameter :: MUNIT = 99

      character(16), parameter :: INFILE1 = 'INFILE1'
      character(16), parameter :: INFILE2 = 'INFILE2'
      character(16), parameter :: OUTFILE = 'OUTFILE'
      character(16), parameter :: PGNAME  = 'IOAPI2UAM'
      character(16) :: INFILE
      character(256) :: MESG

      character(4), dimension(10) :: name
      character(4), dimension(60) :: note
      integer :: iutm,nspec,ibdate,iedate,nx,ny,nz
      real :: btime,etime,xorg,yorg,delx,dely
      integer, parameter :: nseg = 1, iseg = 1, ixseg = 0, iyseg = 0
      integer, parameter :: nzlowr = 0, nzuppr = 0
      real, parameter :: refx = 0.0, refy = 0.0
      real, parameter :: htsur = 0.0, htlow = 0.0, htupp = 0.0
      real, parameter :: sec2hr = 3600.
      character(10), parameter :: namecamx = 'EMISSIONS '
      character(10) :: tmpnam

      character(4), allocatable :: mspec(:,:)
      real, allocatable :: BUFF(:,:,:,:), EM(:,:,:,:)
      integer :: istat

      integer :: itzone,nstep,jchk,jstep
      integer :: i,j,l,n
c
c     Initialize I/O-API
c
      LOGUNIT = INIT3()
c
c     Open input file 1
c
      INFILE = INFILE1
      if (.not.OPEN3(INFILE,FSREAD3,PGNAME)) then
        MESG = 'Cannot open ' // TRIM(INFILE)
        call M3EXIT(PGNAME,0,0,MESG,1)
      endif
c
c     Get description of input file 1
c
      if (.not.DESC3(INFILE)) then
        MESG = 'Cannot get description of ' // TRIM(INFILE)
        call M3EXIT(PGNAME,0,0,MESG,1)
      endif
c
c     Check header info
c
      if (GDTYP3D.ne.LAMGRD3 .and. GDTYP3D.ne.UTMGRD3 .and.
     &    GDTYP3D.ne.LATGRD3 .and. GDTYP3D.ne.POLGRD3) then
        MESG = 'Grid type of ' // TRIM(INFILE) // ' is not supported'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif
      iutm = 0
      if (GDTYP3D.eq.UTMGRD3) iutm = INT(P_ALP3D)
      if (FTYPE3D.ne.GRDDED3) then
        MESG = 'Data type of ' // TRIM(INFILE) // ' must be GRDDED3'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif
c
c     Open output file
c
      CALL NAMEVAL( OUTFILE, MESG )
      open(MUNIT,file=MESG,status='NEW',form='UNFORMATTED')
c
c     Write output header
c
      do n = 1, 10
        name(n) = namecamx(n:n)
      enddo
      do n = 1, 60
        note(n) = FDESC3D(1)(n:n) ! First line of CMAQ file description
      enddo

      itzone = ENVINT('TZONE','Output Time Zone',0,istat)
      if (istat.gt.0) then
        MESG = 'Bad value for TZONE'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif
      JDATE = ENVINT('SDATE','Output Start Date',SDATE3D,istat)
      if (istat.gt.0) then
        MESG = 'Bad value for SDATE'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif
      JTIME = ENVINT('STIME','Output Start Time',STIME3D,istat)
      if (istat.gt.0) then
        MESG = 'Bad value for STIME'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif
      RLENG = ENVINT('RLENG','Output Run Length',240000,istat)
      if (istat.gt.0) then
        MESG = 'Bad value for RLENG'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif

      CDATE = JDATE
      CTIME = JTIME

      ibdate = MOD(CDATE,100000)
      btime = REAL(CTIME/10000) + REAL(TIME2SEC(MOD(CTIME,10000)))/3600.

      call NEXTIME( CDATE, CTIME, RLENG )

      iedate = MOD(CDATE,100000)
      etime = REAL(CTIME/10000) + REAL(TIME2SEC(MOD(CTIME,10000)))/3600.

      write(*,'(/,A)')'Output period (start date/time & end date/time):'
      write(*,'(2(i,f))') ibdate,btime,iedate,etime

      nstep = RLENG / TSTEP3D
      xorg  = XORIG3D
      yorg  = YORIG3D
      delx  = XCELL3D
      dely  = YCELL3D
      nx    = NCOLS3D
      ny    = NROWS3D
      nz    = 1       ! 1-D emissions
      nspec = NVARS3D

      allocate (mspec(10,nspec), stat = istat)
      if (istat.ne.0) then
        MESG = 'Memory allocation failed: MSPEC'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif

      do l = 1, nspec
        read(VNAME3D(l),'(10a1)') (mspec(n,l),n=1,10)
      enddo

      write(MUNIT)name,note,nseg,nspec,ibdate,btime,iedate,etime
      write(MUNIT)refx,refy,iutm,xorg,yorg,delx,dely,nx,ny,nz,
     &            nzlowr,nzuppr,htsur,htlow,htupp
      write(MUNIT)ixseg,iyseg,nx,ny
      write(MUNIT)((mspec(n,l),n=1,10),l=1,nspec)
c
c     Allocate buffer memnory
c
      allocate (BUFF(nx,ny,nz,nspec), stat = istat)
      if (istat.ne.0) then
        MESG = 'Memory allocation failed: BUFF'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif
      allocate (EM(nx,ny,nz,nspec), stat = istat)
      if (istat.ne.0) then
        MESG = 'Memory allocation failed: EM'
        call M3EXIT(PGNAME,0,0,MESG,2)
      endif
c
c     Read/write time-dependent data
c
      CDATE = JDATE
      CTIME = JTIME
      call NEXTIME( CDATE, CTIME, itzone*10000 )

      jchk = JSTEP3(CDATE,CTIME,SDATE3D,STIME3D,TSTEP3D)
      if ( jchk.lt.1 .or. jchk.gt.MXREC3D ) then
        MESG = 'Cannot find start date/time in ' // TRIM(INFILE)
        call M3EXIT(PGNAME,CDATE,CTIME,MESG,2)
      endif
ccc      write(*,'(a,2i)')'Reading ',CDATE,CTIME
      if (.not.READ3(INFILE,ALLVAR3,1,CDATE,CTIME,EM)) then
        MESG = 'Cannot read data from ' // TRIM(INFILE)
        call M3EXIT(PGNAME,CDATE,CTIME,MESG,1)
      endif

      do jstep = 1, nstep

        call NEXTIME( CDATE, CTIME, TSTEP3D )
        jchk = JSTEP3(CDATE,CTIME,SDATE3D,STIME3D,TSTEP3D)
        if ( jchk.lt.1 .or. jchk.gt.MXREC3D ) then
          if (INFILE.eq.INFILE2) then
            MESG = 'Cannot find the following date/time in ' //
     &              TRIM(INFILE)
            call M3EXIT(PGNAME,CDATE,CTIME,MESG,2)
          endif
c
c     Open input file 2
c
          INFILE = INFILE2
          if (.not.OPEN3(INFILE,FSREAD3,PGNAME)) then
            MESG = 'Cannot open ' // TRIM(INFILE)
            call M3EXIT(PGNAME,0,0,MESG,1)
          endif
c
c     Get description of input file 2
c
          if (.not.DESC3(INFILE)) then
            MESG = 'Cannot get description of ' // TRIM(INFILE)
            call M3EXIT(PGNAME,0,0,MESG,1)
          endif
c
c     Check file type and dimensions
c
          if (.not.FILCHK3(INFILE,GRDDED3,nx,ny,ALLAYS3,1)) then
            MESG = 'Inconsistent file type and dimension between ' //
     &              TRIM(INFILE) // ' and ' // TRIM(INFILE1)
            call M3EXIT(PGNAME,0,0,MESG,2)
          endif
c
c     Check species order
c
          if (NVARS3D.ne.nspec) then
            MESG = 'Different number of species between ' //
     &              TRIM(INFILE) // ' and ' // TRIM(INFILE1)
            call M3EXIT(PGNAME,0,0,MESG,2)
          endif
          do l = 1, nspec
            write(tmpnam,'(10a1)') (mspec(n,l),n=1,10)
            if (VNAME3D(l)(1:10).ne.tmpnam) then
              write(*,'(/,A)')
     &             ' No. OUTFILE     INFILE2 (first 10 characters only)'
              do i = 1, nspec
                write(*,'(i3,2x,10a1,2x,a10)') i,(mspec(n,i),n=1,10),
     &                                         VNAME3D(i)
              enddo
              MESG = 'Inconsistent species list between ' //
     &              TRIM(INFILE) // ' and ' // TRIM(INFILE1)
              call M3EXIT(PGNAME,0,0,MESG,2)
            endif
          enddo
        endif ! jchk

ccc        write(*,'(a,2i)')'Reading ',CDATE,CTIME
        if (.not.READ3(INFILE,ALLVAR3,1,CDATE,CTIME,BUFF)) then
          MESG = 'Cannot read data from ' // TRIM(INFILE)
          call M3EXIT(PGNAME,CDATE,CTIME,MESG,1)
        endif

        EM = 0.5 * ( EM + BUFF ) * sec2hr ! hourly averages (unit conversion)

        ibdate = MOD(JDATE,100000)
        btime = REAL(JTIME/10000) +
     &          REAL(TIME2SEC(MOD(JTIME,10000)))/3600.

        call NEXTIME( JDATE, JTIME, TSTEP3D )

        iedate = MOD(JDATE,100000)
        etime = REAL(JTIME/10000) +
     &          REAL(TIME2SEC(MOD(JTIME,10000)))/3600.

        write(MUNIT) ibdate,btime,iedate,etime
        write(*,'(a,2(i,f))')'Writing ',ibdate,btime,iedate,etime

        do l = 1, nspec
          write(MUNIT) iseg,(mspec(n,l),n=1,10),
     &                                  ((EM(i,j,1,l),i=1,nx),j=1,ny)
        enddo

        EM = BUFF

      enddo
c
c     Close output files
c
      close(MUNIT)

      MESG = 'Successful completion of ' // PGNAME
      call M3EXIT(PGNAME,0,0,MESG,0)
	
      end

