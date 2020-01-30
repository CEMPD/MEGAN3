      subroutine interp_lcp(nx,ny,nz,nxc,nyc,nzc,kz1,kz2,ioffset,
     &                      joffset,deltax,dxcamx,kvmeth)
c
c-----INTERP horizontally interpolates, then vertically aggregates,
c     MM5 data on the met model grid to the CAMx grid
c
c     NOTE: CAMx grid is simply a windowed area of the MM5 Lambert Conformal
c           grid, and the CAMx physical height grid is a coarser set of the
c           MM5 sigma-p coordinate system. This version allows simple
c           interpolation from coarse to fine resolution, assuming that the 
c           fine CAMx grid exactly aligns with the coarser MM5 grid 
c           (i.e., the CAMx grid is effectively a nest of the MM5 grid).
c
      implicit none
      include 'param.inc'
      include 'fields.inc'
c
      integer nx,ny,nz,nxc,nyc,nzc,kz1(mnzc),kz2(mnzc),ioffset,joffset
      real deltax,dxcamx
      character*10 kvmeth
c
      integer i,j,k,n,istd,jstd,istx,jstx,iid,jjd,iix,jjx,ic,jc,ii,jj
      real mshfac,du1,du2,u1,u2,du,dv1,dv2,v1,v2,dv,pstarx,pstary
c
c-----Couple MM5 variables to P*
c
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ua(i,j,k) = ua(i,j,k)*psad(i,j)
            va(i,j,k) = va(i,j,k)*psad(i,j)
          enddo
        enddo
        do j = 1,ny-1
          do i = 1,nx-1
            ta(i,j,k) = ta(i,j,k)*psax(i,j)
            pa(i,j,k) = pa(i,j,k)*psax(i,j)
            qa(i,j,k) = qa(i,j,k)*psax(i,j)
          enddo
        enddo
      enddo
c
c-----Compute the starting indices of MM5 grid points to use for
c     the CAMx grid
c
      istd = ioffset + 2
      jstd = joffset + 2
      istx = istd - 1 
      jstx = jstd - 1 
      mshfac = float(nint(deltax/dxcamx))
c
c-----Interpolate wind components from Arakawa B grid (MM5) to
c     Arakawa C grid (CAMx)
c
      do k = 1,nz
        if (mshfac.eq.1.) then
c
c-----Simple case of windowing the MM5 grid with 1:1 meshing
c
          do j = 1,nyc
            do i = 1,nxc
              iid = istd + i - 1
              jjd = jstd + j - 1
              utmp(i,j,k) = (ua(iid,jjd,k) + ua(iid,jjd-1,k))/2.
              vtmp(i,j,k) = (va(iid,jjd,k) + va(iid-1,jjd,k))/2.
              iix = istx + i - 1
              jjx = jstx + j - 1
              ttmp(i,j,k)  = ta(iix,jjx,k)
              ptmp(i,j,k)  = pa(iix,jjx,k)
              qtmp(i,j,k)  = qa(iix,jjx,k)
              cwtmp(i,j,k) = qc(iix,jjx,k)
              prtmp(i,j,k) = qpr(iix,jjx,k)
              pstmp(i,j,k) = qps(iix,jjx,k)
              pgtmp(i,j,k) = qpg(iix,jjx,k)
              odtmp(i,j,k) = tau(iix,jjx,k)
              if (kvmeth.eq.'TKE') tktmp(i,j,k) = tke(iix,jjx,k)
              ztmp(i,j,k) = zh(iix,jjx,k)
              if (k.eq.1) then
                pstar(i,j) = psax(iix,jjx)
                tsfc(i,j)  = tsrf(iix,jjx)
                topcx(i,j) = topo(iix,jjx)
                rground(i,j) = rgrnd(iix,jjx)

                do n = 1,11
                  lucx(i,j,n)  = clu(iix,jjx,n)
                enddo
                if (kvmeth.eq.'OB70' .or. kvmeth.eq.'CMAQ') then
                  pblc(i,j) = pbl(iix,jjx)
                  z0c(i,j)  = z0(iix,jjx)
                endif
              endif
            enddo
          enddo
        else
c
c-----Complex case of windowing and interpolating to finer resolution
c     Winds first
c
          jc = 0
          do j = jstd,jstd+nyc/int(mshfac)-1
            ic = 0
            do i = istd,istd+nxc/int(mshfac)-1

              du1 = ua(i,j-1,k) - ua(i-1,j-1,k)
              du2 = ua(i,j,k) - ua(i-1,j,k)
              do ii = 1,int(mshfac)
                u1 = ua(i-1,j-1,k) + du1*ii/mshfac
                u2 = ua(i-1,j,k) + du2*ii/mshfac
                du = u2 - u1
                do jj = 1,int(mshfac)
                  utmp(ic+ii,jc+jj,k) = u1 + du*
     &                                  (1. + 2.*(jj - 1))/(2.*mshfac)
                enddo
              enddo

              dv1 = va(i-1,j,k) - va(i-1,j-1,k)
              dv2 = va(i,j,k) - va(i,j-1,k)
              do jj = 1,int(mshfac)
                v1 = va(i-1,j-1,k) + dv1*jj/mshfac
                v2 = va(i,j-1,k) + dv2*jj/mshfac
                dv = v2 - v1
                do ii = 1,int(mshfac)
                  vtmp(ic+ii,jc+jj,k) = v1 + dv*
     &                                  (1. + 2.*(ii - 1))/(2.*mshfac)
                enddo
              enddo

              ic = ic + int(mshfac)
            enddo
            jc = jc + int(mshfac)
          enddo
c
c-----Now cell-centered variables
c
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 zh(1,1,k),ztmp(1,1,k))
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 ta(1,1,k),ttmp(1,1,k))
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 pa(1,1,k),ptmp(1,1,k))
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 qa(1,1,k),qtmp(1,1,k))
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 qc(1,1,k),cwtmp(1,1,k))
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 qpr(1,1,k),prtmp(1,1,k))
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 qps(1,1,k),pstmp(1,1,k))
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 qpg(1,1,k),pgtmp(1,1,k))
          call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                 tau(1,1,k),odtmp(1,1,k))
          if (kvmeth.eq.'TKE') then
            call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                   tke(1,1,k),tktmp(1,1,k))
          elseif ((kvmeth.eq.'OB70' .or. kvmeth.eq.'CMAQ')
     &            .and. k.eq.1) then
            call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                   pbl,pblc)
            call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                   z0,z0c)
          endif
          if (k.eq.1) then
            call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                   psax,pstar)
            call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                   tsrf,tsfc)
            call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                   topo,topcx)
            call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                   rgrnd,rground)
            do n = 1,11
              call finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,
     &                     clu(1,1,n),lucx(1,1,n))
            enddo
          endif

        endif
      enddo
c
c-----Map momentum and thermodynamic variables onto the CAMx vertical
c     grid structure
c
      call vertmap(nxc,nyc,nzc,kz1,kz2,sigma,utmp,uac)
      call vertmap(nxc,nyc,nzc,kz1,kz2,sigma,vtmp,vac)
      call vertmap(nxc,nyc,nzc,kz1,kz2,sigma,ttmp,tac)
      call vertmap(nxc,nyc,nzc,kz1,kz2,sigma,ptmp,pac)
      call vertmap(nxc,nyc,nzc,kz1,kz2,sigma,qtmp,qac)
c
c-----Decouple vertically interpolated variables from P*
c
       do k = 1,nzc
         do j = 1,nyc
           do i = 1,nxc
             if (i.eq.nxc) then
               uac(i,j,k) = uac(i,j,k)/pstar(i,j)
             else
               pstarx = (pstar(i,j) + pstar(i+1,j))/2.
               uac(i,j,k) = uac(i,j,k)/pstarx
             endif
             if (j.eq.nyc) then
               vac(i,j,k) = vac(i,j,k)/pstar(i,j)
             else
               pstary = (pstar(i,j) + pstar(i,j+1))/2.
               vac(i,j,k) = vac(i,j,k)/pstary
             endif
           enddo
         enddo
         do j = 1,nyc
           do i = 1,nxc
             tac(i,j,k) = tac(i,j,k)/pstar(i,j)
             pac(i,j,k) = pac(i,j,k)/pstar(i,j)
             qac(i,j,k) = qac(i,j,k)/pstar(i,j)
           enddo
         enddo
       enddo
c
c-----Map layer interface heights and TKE to CAMx vertical grid
c
      do j = 1,nyc 
        do i = 1,nxc 
          do k = 1,nzc 
            zhc(i,j,k) = ztmp(i,j,kz2(k))
            if (kvmeth.eq.'TKE') tkc(i,j,k) = tktmp(i,j,kz2(k))
          enddo 
        enddo 
      enddo 
c
      return
      end

c-----------------------------------------------------------------------

      subroutine finelcp(istx,jstx,nx,ny,nxc,nyc,mshfac,cc,ff)
c
c-----FINELCP horizontally interpolates an input field on LCP cross points
c     to a finer LCP grid
c
      implicit none
      include 'param.inc'
c
      integer istx,jstx,nx,ny,nxc,nyc
      real cc(mnx,mny),ff(mnxc,mnyc),mshfac
c
      integer i,j,ic,jc,ii,jj
      real dc1,dc2,c1,c2,dc

      jc = -1
      do j = jstx,jstx+nyc/int(mshfac)
        ic = -1
        do i = istx,istx+nxc/int(mshfac)

          dc1 = cc(i,j-1) - cc(i-1,j-1)
          dc2 = cc(i,j) - cc(i-1,j)
          do ii = 1,int(mshfac)
            c1 = cc(i-1,j-1) + dc1*ii/mshfac
            c2 = cc(i-1,j) + dc2*ii/mshfac
            dc = c2 - c1
            do 100 jj = 1,int(mshfac)
              if (ic+ii.lt.1 .or. jc+jj.lt.1 .or.
     &            ic+ii.gt.nxc .or. jc+jj.gt.nyc) goto 100
              ff(ic+ii,jc+jj) = c1 + dc*jj/mshfac
 100        enddo
          enddo

          ic = ic + int(mshfac)
        enddo
        jc = jc + int(mshfac)
      enddo

      return
      end
