!*==SETMOD.spg  processed by SPAG 6.72Dc at 15:05 on  5 Dec 2017
      SUBROUTINE SETMOD(Nmmodl,Mmax)
      IMPLICIT NONE
!*--SETMOD4
!*** Start of declarations inserted by SPAG
      REAL A , B , D , DLAm , ETAp , ETAs , FREfp , FREfs , QA , QAQb , &
         & QB , REFdep , RHO , v1 , v2 , v3 , v4 , v5 , v6
      INTEGER icnvel , idimen , iflsph , iiso , INVdep , ITYpe , iunit ,&
            & LGSTR , LIN , LOT , lt , Mmax
!*** End of declarations inserted by SPAG
!-----
!       COMMENTS
!       07 MAY 2002 - eliminated a write(2 to an unopened file
!-----
!       interactively set up model file
!-----
      PARAMETER (LIN=5,LOT=6)
      COMMON /PARAM / QAQb , ITYpe , DLAm , INVdep
      CHARACTER Nmmodl*(*)
 
      COMMON /MODTIT/ TITle
      CHARACTER TITle*80
      INTEGER NL
      PARAMETER (NL=100)
      COMMON /ISOMOD/ D(NL) , A(NL) , B(NL) , RHO(NL) , QA(NL) , QB(NL) &
                    & , ETAp(NL) , ETAs(NL) , FREfp(NL) , FREfs(NL)
      COMMON /DEPREF/ REFdep
 
      WRITE (LOT,*) ' Interactively setting up initial model file:' ,   &
                  & Nmmodl
!-----
!       get model format
!-----
      WRITE (LOT,*) 'Is model flat (0) or spherical (1)'
      READ (LIN,*) iflsph
      WRITE (LOT,*) 'Enter descriptive title for this model'
      READ (LIN,'(a)') TITle
 
      WRITE (LOT,*) ' Enter d,a,b,rho,qa,qb'
      WRITE (LOT,*) '  d=0.0 or EOF  indicates halfspace ' ,            &
                   &'and end of input'
!-----
!       get model data
!-----
      Mmax = 0
 100  READ (LIN,*,ERR=1000,END=200) v1 , v2 , v3 , v4 , v5 , v6
      WRITE (LOT,*) v1 , v2 , v3 , v4 , v5 , v6
!-----
!           safety check
!-----
      IF ( v2.LT.v3 ) THEN
         WRITE (LOT,*) ' Error: P Velocity < S velocity:' , v2 , '<' ,  &
                     & v3
         WRITE (LOT,*) '        Reenter this layer'
         GOTO 100
      ENDIF
      IF ( v4.LE.0.0 ) THEN
         WRITE (LOT,*) 'Error:  Density <= 0:' , v4
         WRITE (LOT,*) '        Reenter this layer'
         GOTO 100
      ENDIF
      IF ( v5.LT.0.0 .OR. v6.LT.0.0 ) THEN
         WRITE (LOT,*) 'Error: qa or qb not >= 0' , v5 , v5
         WRITE (LOT,*) '        Reenter this layer'
         GOTO 100
      ENDIF
      Mmax = Mmax + 1
      D(Mmax) = v1
      A(Mmax) = v2
      B(Mmax) = v3
      RHO(Mmax) = v4
      IF ( v5.GT.1.0 ) v5 = 1.0/v5
      IF ( v6.GT.1.0 ) v6 = 1.0/v6
      QA(Mmax) = v5
      QB(Mmax) = v6
      ETAp(Mmax) = 0.0
      ETAs(Mmax) = 0.0
      FREfp(Mmax) = 1.0
      FREfs(Mmax) = 1.0
 
      IF ( v1.NE.0.0 ) GOTO 100
 200  D(Mmax) = 0.0
 
!-----
!       lun I*4 - logical unit for writing model file. This
!                 unit is released after the use of this routine
!       nmmodl  C*(*)   - model name
!       mmax    I*4 - number of layers in the model, last layer is
!                    halfspace
!       title   C*(*)   - title of the model file
!       iunit   I*4 - 0 Kilometer, Gram, Sec
!       iiso    I*4 - 0 isotropic
!                 1 transversely anisotropic
!                 2 general anisotropic
!       iflsph  I*4 - 0 flat earth model
!                 1 spherical earth model
!       idimen  I*4 - 1 1-D
!               - 2 2-D
!               - 3 3-D
!       icnvel  I*4 - 0 constant velocity
!                 1 variable velocity
!------
      iunit = 0
      iiso = 0
      idimen = 1
      icnvel = 0
      lt = LGSTR(TITle)
      CALL PUTMOD(2,Nmmodl,Mmax,TITle(1:lt),iunit,iiso,iflsph,idimen,   &
                & icnvel,.TRUE.)
1000  END

