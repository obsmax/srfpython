!*==SETDSP.spg  processed by SPAG 6.72Dc at 15:08 on  5 Dec 2017
      SUBROUTINE SETDSP(Nmdisp)
      IMPLICIT NONE
!*--SETDSP4
!*** Start of declarations inserted by SPAG
      REAL dval , obserr , per , val
      INTEGER ilvry , imode , iporg , LIN , LOT
!*** End of declarations inserted by SPAG
!-----
!       Changes
!       21 JUL 2002 - ensure proper format for GAMMA
!-----
!       interactively set up dispersion file
!-----
      PARAMETER (LIN=5,LOT=6)
      CHARACTER Nmdisp*(*)
      CHARACTER*1 lorr(2) , porg(3)
      CHARACTER ostr*12
      DATA lorr/'L' , 'R'/ , porg/'C' , 'U' , 'G'/
 
      WRITE (LOT,*) ' Interactively setting up dispersion file:' ,      &
                  & Nmdisp
!-----
!       open file
!-----
      OPEN (2,FILE=Nmdisp,STATUS='unknown',FORM='formatted',            &
           &ACCESS='sequential')
      REWIND 2
!-----
!       NEW UNITS ARE ALWAYS km AND DATA IS ALWAYS PERIOD
!-----
!       prompt for dispersion information
!-----
      WRITE (LOT,*) ' Enter ilvry,iporg,imode,per,val,dval'
      WRITE (LOT,*) ' ilvry=1(Love)'
      WRITE (LOT,*) '      =2(Rayleigh)'
      WRITE (LOT,*) ' iporg=1 (phase velocity km/s)'
      WRITE (LOT,*) '      =2 (group velocity km/s)'
      WRITE (LOT,*) '      =3 (gamma 1/km)'
      WRITE (LOT,*) ' imode (mode number) e.g., 0=fundamental, 1=first'
      WRITE (LOT,*) ' per=the period '
      WRITE (LOT,*) '    val=dispersion value, velocity or gamma'
      WRITE (LOT,*) '   dval=error in dispersion value'
      WRITE (LOT,*) '       (Enter 1.0 if stderr from residuals)'
      WRITE (LOT,*)                                                     &
              &'   NOTE: Enter all zeros or negative to terminate input'
 100  READ (LIN,*,END=200,ERR=200) ilvry , iporg , imode , per , val ,  &
                                 & dval
      IF ( ilvry.GT.0 .AND. imode.GE.0 ) THEN
         ostr(1:7) = 'SURF96 '
         ostr(8:8) = lorr(ilvry)
         ostr(9:9) = ' '
         ostr(10:10) = porg(iporg)
         ostr(11:12) = ' X'
         obserr = 0.0
         IF ( iporg.EQ.1 .OR. iporg.EQ.2 ) THEN
            WRITE (2,'(a12,1x,i3,1x,3f11.4)') ostr , imode , per , val ,&
                 & dval
         ELSEIF ( iporg.EQ.3 ) THEN
            WRITE (2,'(a12,1x,i3,1x,f11.4,2e11.3)') ostr , imode , per ,&
                 & val , dval
         ENDIF
         GOTO 100
      ENDIF
 200  CLOSE (2)
      END

