!*==LGSTR.spg  processed by SPAG 6.72Dc at 15:08 on  5 Dec 2017
      FUNCTION LGSTR(Str)
!-----
!       function to find the length of a string
!       this will only be used with file system path names
!       thus the first blank
!       indicates the end of the string
!-----
      IMPLICIT NONE
!*--LGSTR10
      CHARACTER*(*) Str
      INTEGER(kind=4) LGSTR
      INTEGER n , i
      n = LEN(Str)
      LGSTR = 1
      DO i = n , 1 , -1
         LGSTR = i
         IF ( Str(i:i)/=' ' ) GOTO 99999
      ENDDO
99999 END

