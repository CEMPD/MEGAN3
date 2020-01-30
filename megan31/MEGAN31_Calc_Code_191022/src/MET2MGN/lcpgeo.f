      subroutine lcpgeo(iway,phic,xlonc,truelat1,truelat2,xloc,yloc,
     &                  xlon,ylat)
c
c----CAMx v4.40 061025
c
c     LCPGEO performs Lambert Conformal to geodetic (lat/lon) translation
c
c     Code based on the TERRAIN preprocessor for MM5 v2.0,
c     developed by Yong-Run Guo and Sue Chen, National Center for
c     Atmospheric Research, and Pennsylvania State University
c     10/21/1993
c
c     Copyright 1996-2006
c     ENVIRON International Corporation
c
c     Modifications:
c        1/3/06              Now handles case with only 1 true lat
c
c     Input arguments:
c        iway                Conversion type
c                            0 = geodetic to Lambert Conformal
c                            1 = Lambert Conformal to geodetic
c        phic                Central latitude (deg, neg for southern hem)
c        xlonc               Central longitude (deg, neg for western hem)
c        truelat1            First true latitute (deg, neg for southern hem)
c        truelat2            Second true latitute (deg, neg for southern hem)
c        xloc/yloc           Projection coordinates (km)
c        xlon/ylat           Longitude/Latitude (deg)
c
c     Output arguments:
c        xloc/yloc           Projection coordinates (km)
c        xlon/ylat           Longitude/Latitude (deg)
c
c     Routines called:
c        none
c
c     Called by:
c        GRDPREP
c
      implicit none
c
      integer iway
      real phic,xlonc,truelat1,truelat2,xloc,yloc,xlon,ylat
c
      real conv,a,sign,pole,xn,psi1,psi0,xc,yc,flp,flpp,r,cell,
     &     rxn,cel1,cel2,psx,ylon
c
      data conv/57.29578/, a/6370./
c
c-----Entry Point
c
      if (phic.lt.0) then
        sign = -1.
      else
        sign = 1.
      endif
      pole = 90.
      if (abs(truelat1).gt.90.) then
        truelat1 = 60.
        truelat2 = 30.
        truelat1 = sign*truelat1
        truelat2 = sign*truelat2
      endif
      if (truelat1.eq.truelat2) then
        xn = sin(truelat2/conv)
      else
        xn = alog10(cos(truelat1/conv)) - alog10(cos(truelat2/conv))
        xn = xn/(alog10(tan((45. - sign*truelat1/2.)/conv)) -
     &           alog10(tan((45. - sign*truelat2/2.)/conv)))
      endif
      psi1 = 90. - sign*truelat1
      psi1 = psi1/conv
      if (phic.lt.0.) then
        psi1 = -psi1
        pole = -pole
      endif
      psi0 = (pole - phic)/conv
      xc = 0.
      yc = -a/xn*sin(psi1)*(tan(psi0/2.)/tan(psi1/2.))**xn
c
c-----Calculate lat/lon of the point (xloc,yloc)
c
      if (iway.eq.1) then
        xloc = xloc + xc
        yloc = yloc + yc
        if (yloc.eq.0.) then
          if (xloc.ge.0.) flp = 90./conv
          if (xloc.lt.0.) flp = -90./conv
        else
          if (phic.lt.0.) then
            flp = atan2(xloc,yloc)
          else
            flp = atan2(xloc,-yloc)
          endif
        endif
        flpp = (flp/xn)*conv + xlonc
        if (flpp.lt.-180.) flpp = flpp + 360.
        if (flpp.gt. 180.) flpp = flpp - 360.
        xlon = flpp
c
        r = sqrt(xloc*xloc + yloc*yloc)
        if (phic.lt.0.) r = -r
        cell = (r*xn)/(a*sin(psi1))
        rxn  = 1.0/xn
        cel1 = tan(psi1/2.)*cell**rxn
        cel2 = atan(cel1)
        psx  = 2.*cel2*conv
        ylat = pole - psx
c
c-----Calculate x/y from lat/lon
c
      else
        ylon = xlon - xlonc
        if (ylon.gt. 180.) ylon = ylon - 360.
        if (ylon.lt.-180.) ylon = ylon + 360.
        flp = xn*ylon/conv
        psx = (pole - ylat)/conv
        r = -a/xn*sin(psi1)*(tan(psx/2.)/tan(psi1/2.))**xn
        if (phic.lt.0.) then
          xloc = r*sin(flp)
          yloc = r*cos(flp)
        else
          xloc = -r*sin(flp)
          yloc =  r*cos(flp)
        endif
      endif
c
      xloc = xloc - xc
      yloc = yloc - yc
c
      return
      end
