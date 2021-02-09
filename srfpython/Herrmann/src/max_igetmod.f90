!*==GETMOD.spg  processed by SPAG 6.72Dc at 15:05 on  5 Dec 2017
      SUBROUTINE GETMOD(Rlun,Mname,Mmax,Title,Iunit,Iiso,Iflsph,Idimen, &
                      & Icnvel,Ierr,Listmd)
!-----
!       HISTORY
!
!       09 08 2000  gave ierr an initial default value for g77
!       01 13 2001  put in close(lun) if file is not model file
!       03 MAY 2002     Modify to permit read from standard input
!       06 JUL 2005 moved inquire to permit use of STDIN
!
!-----
!       General purpose model input
!       This model specification is designed to be as
!           general as possible
!
!       Input lines
!       Line 01: MODEL
!       Line 02: Model Name
!       Line 03: ISOTROPIC or ANISOTROPIC or
!           TRANSVERSELY ANISOTROPIC
!       Line 04: Model Units, First character
!           is length (k for kilometer
!           second is mass (g for gm/cc), third is time (s for time)
!       Line 05: FLAT EARTH or SPHERICAL EARTH
!       Line 06: 1-D, 2-D or 3-D
!       Line 07: CONSTANT VELOCITY
!       Line 08: open for future use
!       Line 09: open for future use
!       Line 10: open for future use
!       Line 11: open for future use
!       Lines 12-end:   These are specific to the model
!           For ISOTROPIC the entries are
!           Layer Thickness, P-velocity, S-velocity, Density, Qp, Qs,
!           Eta-P, Eta S (Eta is frequency dependence),
!           FreqRefP, FreqRefP
!-----
!MODEL
!TEST MODEL.01
!ISOTROPIC
!KGS
!FLAT EARTH
!1-D
!CONSTANT VELOCITY
!LINE08
!LINE09
!LINE10
!LINE11
! H  VP  VS   RHO   QP  QS   ETAP   ETAS REFP  REFS
!1.0    5.0 3.0 2.5 0.0 0.0 0.0 0.0 1.0 1.0
!2.0    5.1 3.1 2.6 0.0 0.0 0.0 0.0 1.0 1.0
!7.0    6.0 3.5 2.8 0.0 0.0 0.0 0.0 1.0 1.0
!10.0   6.5 3.8 2.9 0.0 0.0 0.0 0.0 1.0 1.0
!20.0   7.0 4.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0
!40.0   8.0 4.7 3.3 0.0 0.0 0.0 0.0 1.0 1.0
!-----
!-----
!       rlun    I*4 - logical unit for reading model file. This
!                 unit is released after the use of this routine
!       mname   C*(*)   - model name - if this is stdin or
!           STDIN just read
!                 from standard input
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
!       ierr    I*4 - 0 model file correctly read in
!               - -1 file does not exist
!               - -2 file is not a model file
!                 -3 error in the model file
!       listmd  L   - .true. list the model
!------
 
      IMPLICIT NONE
!*--GETMOD85
      CHARACTER Mname*(*) , Title*(*)
      INTEGER Rlun
      INTEGER(kind=4) Mmax, Iunit, Iiso, Iflsph, Idimen, Icnvel
      INTEGER(kind=4) Ierr
      CHARACTER string*80
      LOGICAL Listmd
!-----
!       STDIN I*4 - logical unit for standard input
!       STDOUT I*4 - logical unit for standard output
!-----
      INTEGER STDIN , STDOUT , LER
      PARAMETER (STDIN=5,STDOUT=6,LER=0)
 
      INTEGER NL
      PARAMETER (NL=100)
      COMMON /ISOMOD/ D(NL), A(NL), B(NL), RHO(NL), QA(NL), QB(NL), &
              &ETAp(NL), ETAs(NL), FREfp(NL), FREfs(NL)
      REAL D, A, B, RHO, QA, QB, ETAp, ETAs, FREfp, FREfs
      COMMON /DEPREF/ REFdep
      REAL REFdep
 
      LOGICAL ext
      CHARACTER ftype*80
      INTEGER lun , j , i , irefdp
 
!-----
!       test to see if the file exists
!-----
      Ierr = 0
!-----
!       test for input
!-----
      IF ( Mname(1:5) == 'stdin' .OR. Mname(1:5) == 'STDIN' ) THEN
!-----
!           do not open anything, use standard output
!-----
         lun = STDIN
      ELSE
         lun = Rlun
         INQUIRE (FILE=Mname,EXIST=ext)
         IF ( .NOT.ext ) THEN
            Ierr = -1
            WRITE (LER,*) 'Model file does not exist'
            RETURN
         ENDIF
!-----
!           open the file
!-----
         OPEN (lun,FILE=Mname,STATUS='old',FORM='formatted',            &
              &ACCESS='sequential')
         REWIND lun
      ENDIF
!-----
!       verify the file type
!-----
!-----
!       LINE 01
!-----
      READ (lun,'(a)') ftype
      IF ( ftype(1:5) /= 'model' .AND. ftype(1:5) /= 'MODEL' ) THEN
         Ierr = -2
         WRITE (LER,*) 'Model file is not in model format'
         CLOSE (lun)
         RETURN
      ENDIF
!-----
!       LINE 02
!-----
      READ (lun,'(a)') Title
!-----
!       LINE 03
!-----
      READ (lun,'(a)') string
      IF ( string(1:3) == 'ISO' .OR. string(1:3) == 'iso' ) THEN
         Iiso = 0
      ELSEIF ( string(1:3) == 'TRA' .OR. string(1:3) == 'tra' ) THEN
         Iiso = 1
      ELSEIF ( string(1:3) == 'ANI' .OR. string(1:3) == 'ani' ) THEN
         Iiso = 2
      ENDIF
!-----
!       LINE 04
!-----
      READ (lun,'(a)') string
      IF ( string(1:3) == 'KGS' .OR. string(1:3) == 'kgs' ) Iunit = 0
!-----
!       LINE 05
!-----
      READ (lun,'(a)') string
      IF ( string(1:3) == 'FLA' .OR. string(1:3) == 'fla' ) THEN
         Iflsph = 0
      ELSEIF ( string(1:3) == 'SPH' .OR. string(1:3) == 'sph' ) THEN
         Iflsph = 1
      ENDIF
!-----
!       LINE 06
!-----
      READ (lun,'(a)') string
      IF ( string(1:3) == '1-d' .OR. string(1:3) == '1-D' ) THEN
         Idimen = 1
      ELSEIF ( string(1:3) == '2-d' .OR. string(1:3) == '2-D' ) THEN
         Idimen = 2
      ELSEIF ( string(1:3) == '3-d' .OR. string(1:3) == '3-D' ) THEN
         Idimen = 3
      ENDIF
!-----
!       LINE 07
!-----
      READ (lun,'(a)') string
      IF ( string(1:3) == 'CON' .OR. string(1:3) == 'con' ) THEN
         Icnvel = 0
      ELSEIF ( string(1:3) == 'VAR' .OR. string(1:3) == 'var' ) THEN
         Icnvel = 1
      ENDIF
!-----
!       get lines 8 through 11
!-----
      DO i = 8 , 11
         READ (lun,'(a)') string
      ENDDO
!-----
!       get model specifically for 1-D flat isotropic
!-----
!-----
!       get comment line
!-----
      READ (lun,'(a)') string
      Mmax = 0
      REFdep = 0.0
      irefdp = 0
      IF ( Iiso == 0 ) THEN
 50      j = Mmax + 1
         READ (lun,*,ERR=100,END=100) D(j) , A(j) , B(j) , RHO(j) ,     &
                                    & QA(j) , QB(j) , ETAp(j) , ETAs(j) &
                                    & , FREfp(j) , FREfs(j)
         IF ( D(j) <  0.0 ) THEN
            D(j) = -D(j)
            REFdep = REFdep + D(j)
            irefdp = j
         ENDIF
         Mmax = j
!add max : break if layer is 0, so the rest of the file is ignored
         IF ( D(j) == 0.0 ) GOTO 100
         GOTO 50
      ENDIF
 100  IF ( Mmax <= 0 ) THEN
         Ierr = -3
         WRITE (LER,*) 'Error in model file'
      ELSEIF ( Listmd ) THEN
         Ierr = 0
         WRITE (STDOUT,99001)
99001    FORMAT (' LAYER             H      P-VEL     S-VEL   DENSITY  '&
               & )
         DO i = 1 , Mmax
            WRITE (STDOUT,99002) i , D(i) , A(i) , B(i) , RHO(i)
99002       FORMAT (' ',i5,5x,4F10.3)
            IF ( i == irefdp ) WRITE (STDOUT,99003)
99003       FORMAT (' ','-SURFACE ','- - - - - ','- - - - - ',          &
                   &'- - - - - ','- - - - - -')
         ENDDO
      ENDIF
      IF ( lun /= STDIN ) CLOSE (lun)
      END

