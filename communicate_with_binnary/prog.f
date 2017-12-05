      program prog
      implicit none

      real*4 h,dcl,dcr
      integer LIN,LER,LOT
      parameter(LER=0,LIN=5,LOT=6)

      read(LIN,*) h,dcl,dcr
      write(LOT,10) h, dcl, dcr

10    format(F4.1,F4.1,F4.1)
      end program prog
