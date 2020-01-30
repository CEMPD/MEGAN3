      subroutine utmgeo(iway,iutmzon,rx4,ry4,rlon4,rlat4)
c  
c-----UTMGEO performs UTM to geodetic (lat/lon) translation:
c
c     This is a Fortran version of the BASIC program "Transverse Mercator
c     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
c     Based on algorithm taken from "Map Projections Used by the USGS"
c     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
c
c     Input arguments:  
c        iway                Conversion type
c                            0 = geodetic to UTM 
c                            1 = UTM to geodetic
c        iutmzon             UTM zone
c        rx4                 UTM easting (km) 
c        ry4                 UTM northing (km) 
c        rlon4               Longitude (deg, negative for W)
c        rlat4               Latitude (deg)
c              
c     Output arguments:  
c        rx4                 UTM easting (km) 
c        ry4                 UTM northing (km) 
c        rlon4               Longitude (deg)
c        rlat4               Latitude (deg)
c              
      implicit real*8 (a-h,o-z)
c
      real   rx4,ry4,rlon4,rlat4
      real*8 north
      logical lsouth
c      
      parameter(pi=3.14159265358979)
      parameter(degrad=pi/180., raddeg=1./degrad)
      parameter(semimaj=6378206.4, semimin=6356583.8)
c     parameter(e2=1.0-(semimin/semimaj)**2.)
c     parameter(e4=e2*e2, e6=e2*e4, ep2=e2/(1.-e2))
      parameter(scfa=.9996)
      parameter(north=0., east=500000.)
c
      e2=1.0-(semimin/semimaj)**2.0
      e4=e2*e2
      e6=e2*e4
      ep2=e2/(1.-e2)
c
c-----Entry point
c
c-----Set Zone parameters
c
      lsouth = .false.
      if( iutmzon .lt. 0 ) lsouth = .true.
      zone = abs(iutmzon)
      cm = zone*6.0 - 183.0
      cmr = cm*degrad
c
c-----Convert inputs from single to double precision
c
      if (iway.eq.1) then
        xx = 1000.*rx4
        yy = 1000.*ry4
        if (lsouth) yy = yy - 1.D7
      else
        dlat = rlat4
        dlon = rlon4
      endif
c
c-----Lat/Lon to UTM conversion
c
      if (iway.eq.0) then
	rlat = degrad*dlat
	rlon = degrad*dlon

	delam = dlon - cm
	if (delam.lt.-180.) delam = delam + 360.
	if (delam.gt.180.) delam = delam - 360.
	delam = delam*degrad
	
        f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat 
        f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024. 
        f2 = f2*sin(2.*rlat) 
        f3 = 15.*e4/256.*45.*e6/1024. 
        f3 = f3*sin(4.*rlat) 
        f4 = 35.*e6/3072. 
        f4 = f4*sin(6.*rlat) 
        rm = semimaj*(f1 - f2 + f3 - f4) 
        if (dlat.eq.90. .or. dlat.eq.-90.) then 
          xx = 0. 
          yy = scfa*rm 
        else 
          rn = semimaj/sqrt(1. - e2*sin(rlat)**2) 
          t = tan(rlat)**2 
          c = ep2*cos(rlat)**2 
          a = cos(rlat)*delam 
           
          f1 = (1. - t + c)*a**3/6. 
          f2 = 5. - 18.*t + t**2 + 72.*c - 58.*ep2 
          f2 = f2*a**5/120. 
          xx = scfa*rn*(a + f1 + f2) 
          f1 = a**2/2. 
          f2 = 5. - t + 9.*c + 4.*c**2 
          f2 = f2*a**4/24. 
          f3 = 61. - 58.*t + t**2 + 600.*c - 330.*ep2 
          f3 = f3*a**6/720. 
          yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3)) 
        endif
	xx = xx + east
	yy = yy + north
c
c-----UTM to Lat/Lon conversion
c
      else
        xx = xx - east 
        yy = yy - north 
        e1 = sqrt(1. - e2) 
        e1 = (1. - e1)/(1. + e1) 
        rm = yy/scfa 
        u = 1. - e2/4. - 3.*e4/64. - 5.*e6/256. 
        u = rm/(semimaj*u) 
         
        f1 = 3.*e1/2. - 27.*e1**3./32. 
        f1 = f1*sin(2.*u) 
        f2 = 21.*e1**2/16. - 55.*e1**4/32. 
        f2 = f2*sin(4.*u) 
        f3 = 151.*e1**3./96. 
        f3 = f3*sin(6.*u) 
        rlat1 = u + f1 + f2 + f3 
        dlat1 = rlat1*raddeg 
        if (dlat1.ge.90. .or. dlat1.le.-90.) then 
          dlat1 = dmin1(dlat1,dble(90.) ) 
          dlat1 = dmax1(dlat1,dble(-90.) ) 
          dlon = cm 
        else 
          c1 = ep2*cos(rlat1)**2. 
          t1 = tan(rlat1)**2. 
          f1 = 1. - e2*sin(rlat1)**2. 
          rn1 = semimaj/sqrt(f1) 
          r1 = semimaj*(1. - e2)/sqrt(f1**3) 
          d = xx/(rn1*scfa) 
           
          f1 = rn1*tan(rlat1)/r1 
          f2 = d**2/2. 
          f3 = 5.*3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2 
          f3 = f3*d**2*d**2/24. 
          f4 = 61. + 90.*t1 + 298.*c1 + 45.*t1**2. - 252.*ep2 - 3.*c1**2
          f4 = f4*(d**2)**3./720. 
          rlat = rlat1 - f1*(f2 - f3 + f4) 
          dlat = rlat*raddeg 
           
          f1 = 1. + 2.*t1 + c1 
          f1 = f1*d**2*d/6. 
          f2 = 5. - 2.*c1 + 28.*t1 - 3.*c1**2 + 8.*ep2 + 24.*t1**2. 
          f2 = f2*(d**2)**2*d/120. 
          rlon = cmr + (d - f1 + f2)/cos(rlat1) 
          dlon = rlon*raddeg 
          if (dlon.lt.-180.) dlon = dlon + 360. 
          if (dlon.gt.180.) dlon = dlon - 360. 
        endif 
      endif
c
c-----Convert precision of outputs
c
      if (iway.eq.1) then
        rlat4 = REAL(dlat)
        rlon4 = REAL(dlon)
      else
        rx4 = REAL(xx/1000.)
        if (lsouth) yy = yy + 1.D7
        ry4 = REAL(yy/1000.)
      endif
c
      return
      end
