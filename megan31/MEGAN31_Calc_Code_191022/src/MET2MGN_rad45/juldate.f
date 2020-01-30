      subroutine juldate(idate)
c 
c-----JULDATE converts date from calender (YYMMDD) format to Julian
c     (YYJJJ) format
c
      implicit none
c
      integer idate
      integer nday(12),iyear,imonth,iday,mday,n,jday
c
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
c
c-----Entry point
c
      iyear = idate/10000
      imonth = (idate - iyear*10000)/100
      iday = idate - iyear*10000 - imonth*100
c
      nday(2) = 28
      if (mod(iyear,4).eq.0) nday(2) = 29
      mday = 0
      do 10 n = 1,imonth-1
        mday = mday + nday(n)
 10   continue
      jday = mday + iday
      idate = iyear*1000 + jday
c
      return
      end
