      subroutine vertmap(nx,ny,nz,kz1,kz2,sigma,xin,xout)
c
c-----VERTMAP vertically aggregates MM5 data on the met model grid
c     to the CAMx grid 
c 
c     NOTE: the CAMx physical height grid is a coarser set of the 
c           MM5 sigma-p coordinate system  
c
      implicit none
      include 'param.inc'
c
      integer nx,ny,nz,kz1(mnzc),kz2(mnzc)
      real xin(mnxc,mnyc,mnz),xout(mnxc,mnyc,mnzc)
      real sigma(0:mnz)
c
      integer i,j,k,kk
      real sum,dsigma
c
      do j = 1,ny
        do i = 1,nx
          do k = 1,nz
            sum = 0.
            do kk = kz1(k),kz2(k) 
              dsigma = sigma(kk-1) - sigma(kk)
              sum = sum + xin(i,j,kk)*dsigma
            enddo
            dsigma = sigma(kz1(k)-1) - sigma(kz2(k))
            xout(i,j,k) = sum/dsigma
          enddo
        enddo
      enddo
c
      return
      end
