      program prog
      implicit none

      real*4 h,dcl,dcr
      integer LIN,LER,LOT
      parameter(LER=0,LIN=5,LOT=6)

      read(LIN,*) h,dcl,dcr
      write(LOT,*) "got", h, dcl, dcr

!     10    format(F3.1,F3.1,F3.1)
      end program prog
