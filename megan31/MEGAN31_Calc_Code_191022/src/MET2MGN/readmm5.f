      subroutine readmm5(iomm5,lfirst,project,nxc,nyc,nzc,
     &                   ioffset,joffset,dxcamx,dycamx,x0camx,y0camx,
     &                   izone,clonin,clatin,tlat1in,tlat2in,mm5date,
     &                   mm5hr,itzon,dtout,nx,ny,nz,deltax,ierr)
c
c-----READMM5 reads hourly MM5 data from file, and converts gridded data into
c     conventions for CAMx
c
      implicit none
      include 'param.inc'
      include 'fields.inc'
c
      integer addday,subday
c
      integer iomm5,nxc,nyc,nzc,ioffset,joffset,izone,mm5date,nx,ny,nz,
     &        ierr,itzon,dtout
      real dxcamx,dycamx,x0camx,y0camx,deltax,clonin,clatin,tlat1in,
     &     tlat2in
      logical lfirst
      character*10 project
c
      integer lucats,i,j,k,n,nxcrs,nycrs,nest,ic,jc,mshfac,mm5yr,mm5mo,
     &        mm5day,mm5hr,mm5mn,mm5sc,lu,kbot,ktop,inuf,nfrc
      parameter(lucats=24)
      integer maplu(lucats)
      logical lnewhr,lconv
      real rlu(mnx,mny),z0cat(lucats)
      real pp(mnz),cc(mnz),tt(mnz),qq(mnz),pr(mnz),qp(mnz)
      real rd,g,dxcrs,clat,clon,tlat1,tlat2,xcen,ycen,dxrat,riloc,
     &     rjloc,x0mm5,y0mm5,xfmm5,yfmm5,rlon,rlat,xsw,ysw,xse,yse,
     &     xne,yne,xnw,ynw,xmin,xmax,ymin,ymax,hrfrac,ps0,ts0,tlp,
     &     p0,term1,term2,rtmpr,rtmpc,tv,rho,tausum,dz,c1t,qcmean,
     &     precip,depth,drpdia,drpmas,pi,drpvel,cscav,zmin,cliq,time,
     &     xoffset,yoffset
c
c-----MM5 big header
c
      integer bhi(50,20)
      character*80 bhic(50,20)
      real bhr(20,20)
      character*80 bhrc(20,20)
      integer iflag
c
c-----MM5 subheader
c
      integer ndim,istart(4),iend(4)  
      real*4  xtime,sigmid(mnz)
      character*4 cstag, ctype
      character*24 curdate
      character*8 cfield
      character*25 cunit
      character*46 cdescrip

      data rd /287./
      data g /9.8/
c
c-----Data statements to map MM5 24 category USGS landuse to 11 category
c     CAMx landuse
c
      data maplu(1)  / 1/ !Urban
      data maplu(2)  /10/ !Dry Crop/Pasture
      data maplu(3)  /10/ !Irrigated Crop/Pasture
      data maplu(4)  /10/ !Mixed Dry/Irrigated Crop/Pasture
      data maplu(5)  /10/ !Crop/Grassland Mosaic
      data maplu(6)  / 4/ !Crop/Woodland Mosaic
      data maplu(7)  / 3/ !Grassland
      data maplu(8)  / 3/ !Shrubland
      data maplu(9)  / 3/ !Mix Shrub/Grass
      data maplu(10) / 3/ !Savanna
      data maplu(11) / 4/ !Deciduous Broadleaf Forest
      data maplu(12) / 4/ !Deciduous Needleleaf Forest
      data maplu(13) / 5/ !Evergreen Broadleaf Forest
      data maplu(14) / 5/ !Evergreen Needleleaf Forest
      data maplu(15) / 6/ !Mixed Forest
      data maplu(16) / 7/ !Water
      data maplu(17) / 9/ !Herbaceous Wetland
      data maplu(18) / 5/ !Wooded Wetland
      data maplu(19) / 8/ !Barren, Sparse Vegetation
      data maplu(20) / 3/ !Herbaceous Tundra
      data maplu(21) / 4/ !Wooded Tundra
      data maplu(22) / 6/ !Mixed Tundra
      data maplu(23) /11/ !Bare Ground Tundra
      data maplu(24) / 8/ !Snow or Ice
c
      lnewhr = .true.
      pi = acos(-1.0)
c
c-----Read an MM5 header flag
c
  10  continue
      read(iomm5,end=999) iflag
c
c-----This is the "big header" at the top of the file
c
      if (iflag.eq.0) then
        read(iomm5) bhi,bhr,bhic,bhrc
c
c-----First time through, initialize cloud and calculate grid parameters
c
        if (lfirst) then
c
c          do k = 1,nz
cc            do j = 1,ny-1
c              do i = 1,nx-1
c                ql(i,j,k) = 0.
c                qi(i,j,k) = 0.
c                qr(i,j,k) = 0.
c                qs(i,j,k) = 0.
c                qg(i,j,k) = 0.
c                qc(i,j,k) = 0.
c                qcold(i,j,k) = 0.
c                qpr(i,j,k) = 0.
c                qps(i,j,k) = 0.
c                qpg(i,j,k) = 0.
c              enddo
c            enddo
c          enddo
c
c-----MM5 coarse grid
c
          nycrs = bhi(5,1)
          nxcrs = bhi(6,1)
          dxcrs = bhr(1,1)/1000.
          clat  = bhr(2,1)
          clon  = bhr(3,1)
          tlat1 = bhr(5,1)
          tlat2 = bhr(6,1)
          xcen = 1 + float(nxcrs - 1)/2.
          ycen = 1 + float(nycrs - 1)/2.
          write(*,*)
          write(*,*)'Grid parameters for the Coarse MM5 domain'
          write(*,'(a,2i10)')  '            NX,NY:',nxcrs,nycrs
          write(*,'(a,f10.3)') '               DX:',dxcrs
          write(*,'(a,2f10.3)')'  Central Lat/lon:',clat,clon
          write(*,'(a,2f10.3)')'   True Latitudes:',tlat1,tlat2
          write(*,'(a,2f10.3)')'Grid center (I,J):',xcen,ycen
c
c-----MM5 nest to be processed
c
          nest  = bhi(13,1)
          ny    = bhi(16,1)
          nx    = bhi(17,1)
          nz    = bhi(12,11)
          dxrat = float(bhi(20,1))
          rjloc = bhr(10,1)
          riloc = bhr(11,1)
          deltax = dxcrs/dxrat
          x0mm5 = dxcrs*(riloc - xcen)
          y0mm5 = dxcrs*(rjloc - ycen)
          xfmm5 = x0mm5 + float(nx-1)*deltax
          yfmm5 = y0mm5 + float(ny-1)*deltax
          write(*,*)
          write(*,*)'Grid parameters for the input MM5 domain'
          write(*,'(a,i10)')   '          NEST ID:',nest
          write(*,'(a,3i10)')  '         NX,NY,NZ:',nx,ny,nz
          write(*,'(a,f10.3)') '               DX:',deltax
          write(*,'(a,2f10.3)')'Location of (1,1):',riloc,rjloc
          write(*,'(a,2f10.3)')'    SW x/y corner:',x0mm5,y0mm5
          write(*,'(a,2f10.3)')'    NE x/y corner:',xfmm5,yfmm5
          if (nx.gt.mnx .or. ny.gt.mny .or. nz.gt.mnz) then
            write(*,*)'MM5 dimensions too large for arrays'
            write(*,*)'Increase array dimensions in param.inc ',
     &                'and recompile'
            stop
          endif
          if (project.eq.'LCP       ') then
            if (clat.eq.clatin .and. clon.eq.clonin .and.
     &         ((tlat1.eq.tlat1in .and. tlat2.eq.tlat2in) .or. 
     &          (tlat1.eq.tlat2in .and. tlat2.eq.tlat1in))) then
              project = 'LCP1      '
              write(*,'(a)')'Input LCP identical to MM5 LCP'
              xoffset = (x0camx - x0mm5)/deltax
              yoffset = (y0camx - y0mm5)/deltax
              ioffset = nint((x0camx - x0mm5)/deltax)
              joffset = nint((y0camx - y0mm5)/deltax)
              if (abs(xoffset-float(ioffset)).gt.0.001 .or.
     &            abs(yoffset-float(joffset)).gt.0.001) then
                write(*,'(a)') 
     &              'Grid corners do not align with MM5 dot points'
                stop
              endif
              write(*,'(a,2i10)')'      I,J offsets:',ioffset,joffset
            else
              project = 'LCP2      '
              write(*,'(a)')'Interpolating to input grid LCP'
            endif
          endif
          if (project.ne.'LCP1      ') then
            do i = 1,nx
              xdot(i) = x0mm5 + deltax*(i - 1)
              xcrs(i) = x0mm5 + deltax*(i - 0.5)
            enddo
            do j = 1,ny
              ydot(j) = y0mm5 + deltax*(j - 1)
              ycrs(j) = y0mm5 + deltax*(j - 0.5)
            enddo
            do j = 1,ny
              do i = 1,nx
                call lcpgeo(1,clat,clon,tlat1,tlat2,xdot(i),ydot(j),
     &                      rlon,rlat)
                if (project.eq.'LATLON    ') then
                  xdprj(i,j) = rlon
                  ydprj(i,j) = rlat
                elseif (project.eq.'UTM       ') then
                  call utmgeo(0,izone,xdprj(i,j),ydprj(i,j),rlon,rlat)
                elseif (project.eq.'LCP2      ') then
                  call lcpgeo(0,clatin,clonin,tlat1in,tlat2in,
     &                        xdprj(i,j),ydprj(i,j),rlon,rlat)
                endif
              enddo
            enddo
          endif
c
c-----CAMx grid
c
          if (project.ne.'LCP1      ') then
            do ic = 1,nxc
              xc(ic) = x0camx + dxcamx*(ic - 0.5)
            enddo
            do jc = 1,nyc
              yc(jc) = y0camx + dycamx*(jc - 0.5)
            enddo
            do jc = 1,nyc
              do ic = 1,nxc
                if (project.eq.'UTM       ') then
                  call utmgeo(1,izone,xc(ic),yc(jc),rlon,rlat)
                elseif (project.eq.'LCP2      ') then
                  call lcpgeo(1,clatin,clonin,tlat1in,tlat2in,xc(ic),
     &                        yc(jc),rlon,rlat)
                elseif (project.eq.'LATLON    ') then
                  rlon = xc(ic)
                  rlat = yc(jc)
                endif
                call lcpgeo(0,clat,clon,tlat1,tlat2,xclcp(ic,jc),
     &                      yclcp(ic,jc),rlon,rlat)
              enddo
            enddo
            xsw = xclcp(1,1)
            ysw = yclcp(1,1)
            xse = xclcp(nxc,1)
            yse = yclcp(nxc,1)
            xne = xclcp(nxc,nyc)
            yne = yclcp(nxc,nyc)
            xnw = xclcp(1,nyc)
            ynw = yclcp(1,nyc)
          else
            deltax = anint(deltax*1000.)/1000.
            dxcamx = anint(dxcamx*1000.)/1000.
            dycamx = anint(dycamx*1000.)/1000.
            if (deltax.eq.dxcamx) then
              xsw = x0mm5 + deltax*(ioffset + 0.5)
              ysw = y0mm5 + deltax*(joffset + 0.5)
            else
              xsw = x0mm5 + deltax*ioffset + 0.5*dxcamx
              ysw = y0mm5 + deltax*joffset + 0.5*dxcamx
            endif
            xne = xsw + dxcamx*(nxc - 1.)
            yne = ysw + dycamx*(nyc - 1.)
            xse = xne
            yse = ysw
            xnw = xsw
            ynw = yne
          endif
          write(*,*)
          write(*,*)'Grid parameters for the input domain in LCP space'
          write(*,'(a,3i10)')  '         NX,NY,NZ:',nxc,nyc,nzc
          write(*,'(a,2f10.3)')'            DX,DY:',dxcamx,dycamx
          write(*,'(a,2f10.3)')' SW cell midpoint:',xsw,ysw
          write(*,'(a,2f10.3)')' SE cell midpoint:',xse,yse
          write(*,'(a,2f10.3)')' NE cell midpoint:',xne,yne
          write(*,'(a,2f10.3)')' NW cell midpoint:',xnw,ynw
          if (project.eq.'LCP1      ') then
            if (dxcamx.ne.dycamx) then
              write(*,*)
     &           'For matching LCP projections, input DX must equal DY'
              stop
            endif
            if (deltax.eq.dxcamx .and.
     &          (nxc.gt.nx-1 .or. nyc.gt.ny-1 .or. nzc.gt.nz)) then
              write(*,*)'Input grid dimensions exceed MM5 grid',
     &                  ' dimensions'
              write(*,*)'It must range from 1 through nx-1 and ny-1'
              stop
            elseif (deltax.ne.dxcamx) then
              if (amod(deltax,dxcamx).gt. 0.001) then
                write(*,*)'For matching LCP projections,'
                write(*,*)'MM5 DX must be integer multiple of input 
     &                     grid DX'
                stop
              endif
              mshfac = nint(deltax/dxcamx)
              if (ioffset.lt.1 .or. joffset.lt.1 .or. 
     &            ioffset+nxc/mshfac.gt.nx-1 .or. 
     &            joffset+nyc/mshfac.gt.ny-1) then
                write(*,*)'For matching LCP projections,'
                write(*,*)'Finer input grid mesh cannot span entire 
     &                      MM5 mesh'
                write(*,*)'It must range from 2 through nx-1 and ny-1'
                stop
              endif
            endif
          endif
c
          if (project.ne.'LCP1      ') then
            xmin = x0mm5 + deltax/2.
            xmax = xfmm5 - deltax/2.
            ymin = y0mm5 + deltax/2.
            ymax = yfmm5 - deltax/2.
            if (xsw.lt.xmin .or. xnw.lt.xmin .or.
     &          ysw.lt.ymin .or. yse.lt.ymin .or.
     &          xse.gt.xmax .or. xne.gt.xmax .or.
     &          ynw.gt.ymax .or. yne.gt.ymax) then
              write(*,*)'Input grid ranges outside MM5 grid'
              stop
            endif
            do jc = 1,nyc
              do ic = 1,nxc
                do i = 1,nx
                  if (xdot(i).gt.xclcp(ic,jc)) then
                    idot(ic,jc) = i
                    goto 101
                  endif
                enddo
                write(*,*)'Did not find idot',ic,jc,xclcp(ic,jc)    
                stop
 101            do i = 1,nx
                  if (xcrs(i).gt.xclcp(ic,jc)) then
                    icrs(ic,jc) = i
                    goto 102
                  endif
                enddo
                write(*,*)'Did not find icrs',ic,jc,xclcp(ic,jc)    
                stop
 102            do j = 1,ny
                  if (ydot(j).gt.yclcp(ic,jc)) then
                    jdot(ic,jc) = j
                    goto 103
                  endif
                enddo
                write(*,*)'Did not find jdot',ic,jc,yclcp(ic,jc)    
                stop
 103            do j = 1,ny
                  if (ycrs(j).gt.yclcp(ic,jc)) then
                    jcrs(ic,jc) = j
                    goto 104
                  endif
                enddo
                write(*,*)'Did not find jcrs',ic,jc,yclcp(ic,jc)    
                stop
 104            continue
              enddo
            enddo
          endif

        endif
        goto 10
c
c-----This is a data "sub-header"; read data fields.
c
      elseif (iflag.eq.1) then
        read (iomm5) ndim,istart,iend,xtime,cstag,ctype,curdate,
     &                 cfield,cunit,cdescrip
        if (lnewhr) then
          lnewhr = .false.
          read(curdate(3:4),'(i2)') mm5yr
          read(curdate(6:7),'(i2)') mm5mo
          read(curdate(9:10),'(i2)') mm5day
          read(curdate(12:13),'(i2)') mm5hr
          read(curdate(15:16),'(i2)') mm5mn
          read(curdate(18:19),'(i2)') mm5sc
c
c-----Convert MM5 date to Julian day and convert hour to
c     user-specified time zone
c
          mm5mn = nint(float(mm5mn) + anint(float(mm5sc)/60.))
          nfrc = 60/dtout
          do n = 0,nfrc 
            if (abs(mm5mn - n*dtout).lt.2) then
              mm5mn = n*dtout
              goto 20
            endif
          enddo
          write(*,'(a,a)')'MM5 clock is nuts! Cannot sychronize to',
     &                    ' user-specified output inverval'
          stop
            
 20       mm5date = mm5yr*10000 + mm5mo*100 + mm5day

          call juldate(mm5date)
          if (mm5mn.gt.59) then
            mm5mn = mm5mn - 60
            mm5hr = mm5hr + 1
            if (mm5hr.gt.23) then
              mm5hr = mm5hr - 24
              mm5date = addday(mm5date)
            endif
          endif
          mm5hr = mm5hr - itzon
          if (mm5hr.lt.0) then
            mm5hr = mm5hr + 24
            mm5date = subday(mm5date)
          endif
          if (mm5hr.gt.23) then
            mm5hr = mm5hr - 24
            mm5date = addday(mm5date)
          endif
          mm5hr = 100*mm5hr + mm5mn
c          write(*,'(/,a,t30,i6.5,i5.4)') ' MM5 date/time (YYJJJ HHMM):',
c     &                                   mm5date,mm5hr
        endif

        if (cfield.eq.'U       ') then
          read(iomm5) (((ua(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
        elseif (cfield.eq.'V       ') then
          read(iomm5) (((va(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
        elseif (cfield.eq.'T       ') then
          read(iomm5) (((ta(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
        elseif (cfield.eq.'Q       ') then
          read(iomm5) (((qa(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
C        elseif (cfield.eq.'CLW     ') then
C          read(iomm5) (((ql(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
C        elseif (cfield.eq.'ICE     ') then
C          read(iomm5) (((qi(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
C        elseif (cfield.eq.'RNW     ') then
C          read(iomm5) (((qr(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
C        elseif (cfield.eq.'SNOW    ') then
C          read(iomm5) (((qs(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
C        elseif (cfield.eq.'GRAUPEL ') then
C          read(iomm5) (((qg(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
C        elseif (cfield.eq.'RAIN CON') then
C          read(iomm5) ((rainc(i,j),j=1,ny),i=1,nx)
C        elseif (cfield.eq.'RAIN NON') then
C          read(iomm5) ((rainr(i,j),j=1,ny),i=1,nx)
C        elseif (cfield.eq.'TKE     ') then
C          read(iomm5) (((tke(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
        elseif (cfield.eq.'PP      ') then
          read(iomm5) (((pa(i,j,nz-k+1),j=1,ny),i=1,nx),k=1,nz)
        elseif (cfield.eq.'PSTARCRS') then
          read(iomm5) ((psax(i,j),j=1,ny),i=1,nx)
        elseif (cfield.eq.'GROUND T') then
          read(iomm5) ((tsrf(i,j),j=1,ny),i=1,nx)
        elseif (cfield.eq.'SWDOWN  ') then
          read(iomm5) ((rgrnd(i,j),j=1,ny),i=1,nx)
C        elseif (cfield.eq.'PBL HGT ') then
C          read(iomm5) ((pbl(i,j),j=1,ny),i=1,nx)
        elseif (cfield.eq.'TERRAIN ') then
          read(iomm5) ((topo(i,j),j=1,ny),i=1,nx)
        elseif (cfield.eq.'LAND USE') then
          read(iomm5) ((rlu(i,j),j=1,ny),i=1,nx)
        elseif (cfield.eq.'SFZ0    ') then
          read(iomm5,err=30) (z0cat(i),i=1,lucats)
          goto 10
 30       write(*,*)'Error reading MM5 surface roughness for 24 landuse'
          write(*,*)'categories. You need to edit READMM5.F to'
          write(*,*)'adjust the number of categories and to alter the'
          write(*,*)'landuse mapping algorithm.'
        elseif (cfield.eq.'SIGMAH  ') then
          read(iomm5) (sigmid(nz-k+1),k=1,nz)
        else
          read(iomm5)
        endif
        goto 10
c
c-----Finished reading data for the hour
c
      elseif (iflag.eq.2) then
C        write(*,*)
C        write(*,*) 'Finished reading MM5 output for ',curdate
      endif
c
c-----Map surface roughness to grid according to landuse map
c     and convert from cm to m.  Map MM5 landuse to input grid landuse.
c
      do j = 1,ny-1
        do i = 1,nx-1
          lu = int(rlu(i,j))
          z0(i,j) = z0cat(lu)/100.
          do n = 1,11
            clu(i,j,n) = 0.
          enddo
          clu(i,j,maplu(lu)) = 1.
        enddo
      enddo
c
c-----Calculate cartesian height from sigma levels, convert height MSL
c     to height AGL (m), and determine P* at dot points
c
      if (lfirst) then
        ptop = bhr(2,2)
        ps0 = bhr(2,5)
        ts0 = bhr(3,5)
        tlp = bhr(4,5)
        sigma(0) = 1.
        do k = 1,nz
          sigma(k) = 2.*sigmid(k) - sigma(k-1)
        enddo
        do j = 1,ny-1
          do i = 1,nx-1
            p0 = sigma(0)*psax(i,j) + ptop
            term1 = rd*tlp/(2.*g)*(alog(p0/ps0))**2
            term2 = rd*ts0/g * alog(p0/ps0)
            trn(i,j) = -(term1+term2)
          enddo
        enddo
        do k = 1,nz
          do j = 1,ny-1
            do i = 1,nx-1
              p0 = sigma(k)*psax(i,j) + ptop
              term1 = rd*tlp/(2.*g)*(alog(p0/ps0))**2
              term2 = rd*ts0/g * alog(p0/ps0)
              zh(i,j,k) = -(term1+term2)
              zh(i,j,k) = zh(i,j,k) - trn(i,j)
              if (zh(i,j,k).lt.0.) then
                write(*,*)'ZH<0: at (i,j,k):',i,j,k,zh(i,j,k)
                stop
              endif
            enddo
          enddo
        enddo

        call xtod(psax,psad,nx,ny)

      endif
c
c-----Convert from pressure perturbation to pressure (mb) by the relation
c     p = pert + (p*)x(sigma) + ptop 
c
      do k = 1,nz
        do j = 1,ny-1
          do i = 1,nx-1
            pa(i,j,k) = (pa(i,j,k) + psax(i,j)*sigmid(k) + ptop)/100.
          enddo
        enddo
      enddo
 100  continue
c
      if (lfirst) lfirst = .false.
      return
c
999   continue
      ierr = 1
      return
      end
