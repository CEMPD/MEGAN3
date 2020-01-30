      subroutine micromet(temp,temp0,press,press0,uwind,vwind,deltaz,z0,
     &                    tr)
c 
c-----MICROMET calculates surface layer micro-meteorological flux-gradient
c     relationships and variables based on Louis (1979), and diagnoses
c     new wind and temperatures at standard probe heights.
c  
c     Modifications:
c        none
c             
c     Input arguments:  
c        temp                Layer 1 temperature (K) 
c        temp0               Surface temperature (K)
c        press               Layer 1 pressure (mb)
c        press0              Surface pressure (mb)
c        uwind               Layer 1 midpoint U-component speed (m/s)
c        vwind               Layer 1 midpoint V-component speed (m/s)
c        deltaz              Layer 1 midpoint height (m)
c        z0                  Surface roughness length (m)
c             
c     Output arguments:  
c        tr                  temperature at 2 m reference height (K)
c 
      real theta,theta0,thetar,dtheta,thetabar,wind,ufrac,vfrac,ri
      real zscale,cm,ch,fm,fh,ustar,ustar2,thstar,el
      real w2,pr       
      data vk/0.4/, g/9.8/, gamma/0.286/
c
c-----Calculate potential temperature and richardson number
c
      theta = temp*(1000./press)**gamma
      theta0 = temp0*(1000./press0)**gamma
      dtheta = theta - theta0
      thetabar = (theta + theta0)/2.
      wind = sqrt(uwind**2 + vwind**2)
      ufrac = uwind/wind
      vfrac = vwind/wind
      ri = (g/thetabar)*deltaz*dtheta/wind**2
c
c-----Determine stability functions
c
      zscale = vk/alog(deltaz/z0)
      if (ri.lt.0.) then
        cm    = 69.56*sqrt(deltaz/z0)*zscale**2
        ch    = 49.82*sqrt(deltaz/z0)*zscale**2
        fm    = 1. - 9.4*ri/(1. + cm*sqrt(abs(ri)))
        fh    = 1. - 9.4*ri/(1. + ch*sqrt(abs(ri)))
      else
        fm = 1./((1. + 4.7*ri)**2)
        fh = fm
      endif
c
c-----Calculate micromet variables
c
      ustar2 = fm*(wind*zscale)**2
      ustar2 = amax1(1.e-10,ustar2)
      ustar = sqrt(ustar2)
      thstar = 1.35*zscale**2*wind*dtheta*fh/ustar
      el = ustar2*temp/(vk*g*thstar + 1.e-10)
c
c-----Determine T at 2m
c
      w2 = ustar/(vk*sqrt(fm))*alog(2/z0)
      thetar = theta0 + ustar*thstar/(1.35*w2*fh*(vk/alog(2/z0))**2)  
      pr = press0 + (2/deltaz)*(press - press0)                       
      tr = thetar*(1000./pr)**(-gamma)                                  
      if (thstar.gt.0 .and. tr.gt.temp) then
        tr = temp0 + (2/deltaz)*(temp - temp0)
      elseif (thstar.lt.0 .and. tr.lt.temp) then
        tr = temp0 + (2/deltaz)*(temp - temp0)
      endif

      return
      end
