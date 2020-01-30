      subroutine xtod(arrin,arrout,nx,ny)

c-----Interpolates a 2-D array from cross (x) points to dot (d) points

      implicit none
      include 'param.inc'
c
      integer nx,ny,i,j
      real arrin(mnx,mny),arrout(mnx,mny)

      arrout(1,1) = arrin(1,1)
      arrout(nx,ny) = arrin(nx-1,ny-1)
      arrout(1,ny) = arrin(1,ny-1)
      arrout(nx,1) = arrin(nx-1,1)

      do j = 2,ny-1
        do i = 2,nx-1
          arrout(i,j) = 0.25*(arrin(i-1,j-1) + arrin(i,j-1) + 
     &                        arrin(i,j) + arrin(i-1,j))
        enddo
      enddo

      do i = 2,nx-1
        arrout(i,1) = 0.5*(arrin(i-1,1) + arrin(i,1))
        arrout(i,ny) = 0.5*(arrin(i-1,ny-1) + arrin(i,ny-1))
      enddo

      do j = 2,ny-1
        arrout(1,j) = 0.5*(arrin(1,j-1) + arrin(1,j))
        arrout(nx,j) = 0.5*(arrin(nx-1,j-1) + arrin(nx-1,j))
      enddo

      return
      end
