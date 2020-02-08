!*==PUTMOD.spg  processed by SPAG 6.72Dc at 15:07 on  5 Dec 2017
      SUBROUTINE PUTMOD(Wlun,Mname,Mmax,Title,Iunit,Iiso,Iflsph,Idimen, &
                      & Icnvel,Lverby)
!-----
!       CHANGES
!       03 MAY 2002 permit write to standard output
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
!       Line 04: Model Units, First character is length (k for kilometer
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
!       wlun    I*4 - logical unit for writing model file. This
!                 unit is released after the use of this routine
!       mname   C*(*)   - model name
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
!       lverby  L   - .false. quiet output
!------
      IMPLICIT NONE
!*--PUTMOD72
 
      CHARACTER Mname*(*)
      CHARACTER Title*(*)
      INTEGER Wlun
      INTEGER*4 Mmax , Iunit , Iiso , Iflsph , Idimen , Icnvel
      LOGICAL Lverby
!-----
!       STDIN I*4 - logical unit for standard input
!       STDOUT I*4 - logical unit for standard output
!-----
      INTEGER lun
      INTEGER STDOUT , LER
      PARAMETER (STDOUT=6,LER=0)
 
      INTEGER NL
      PARAMETER (NL=200)
      COMMON /ISOMOD/ D(NL) , A(NL) , B(NL) , RHO(NL) , QA(NL) , QB(NL) &
                    & , ETAp(NL) , ETAs(NL) , FREfp(NL) , FREfs(NL)
      REAL D , A , B , RHO , QA , QB , ETAp , ETAs , FREfp , FREfs
      COMMON /DEPREF/ REFdep
      REAL REFdep
 
      INTEGER LGSTR , lt
      INTEGER j
      REAL curdep , dout
 
      LOGICAL ext
      CHARACTER cmnt*110
      cmnt(1:11) = '      H(KM)    '
      cmnt(12:22) = '   VP(KM/S)'
      cmnt(23:33) = '   VS(KM/S)'
      cmnt(34:44) = ' RHO(GM/CC)'
      cmnt(45:55) = '     QP    '
      cmnt(56:66) = '     QS    '
      cmnt(67:77) = '   ETAP    '
      cmnt(78:88) = '   ETAS    '
      cmnt(89:99) = '  FREFP    '
      cmnt(100:110) = '  FREFS    '
 
      lt = LGSTR(Title)
!-----
!       test to see if the file exists
!-----
      IF ( Mname(1:6) == 'stdout' .OR. Mname(1:6) == 'STDOUT' ) THEN
!-----
!           do not open anything, use standard output
!-----
         lun = STDOUT
      ELSE
         INQUIRE (FILE=Mname,EXIST=ext)
         IF ( ext .AND. Lverby ) WRITE (LER,*)                          &
                                      &'Overwriting Existing model File'
         lun = Wlun
!-----
!           open the file
!-----
         OPEN (lun,FILE=Mname,STATUS='unknown',FORM='formatted',        &
              &ACCESS='sequential')
         REWIND lun
      ENDIF
!-----
!       verify the file type
!-----
!-----
!       LINE 01
!-----
      WRITE (lun,'(a)') 'MODEL.01'
!-----
!       LINE 02
!-----
      WRITE (lun,'(a)') Title(1:lt)
!-----
!       LINE 03
!-----
      IF ( Iiso == 0 ) THEN
         WRITE (lun,'(a)') 'ISOTROPIC'
      ELSEIF ( Iiso == 1 ) THEN
         WRITE (lun,'(a)') 'TRANSVERSE ANISOTROPIC'
      ELSEIF ( Iiso == 2 ) THEN
         WRITE (lun,'(a)') 'ANISOTROPIC'
      ENDIF
!-----
!       LINE 04
!-----
      WRITE (lun,'(a)') 'KGS'
!-----
!       LINE 05
!-----
      IF ( Iflsph == 0 ) THEN
         WRITE (lun,'(a)') 'FLAT EARTH'
      ELSEIF ( Iflsph == 1 ) THEN
         WRITE (lun,'(a)') 'SPHERICAL EARTH'
      ENDIF
!-----
!       LINE 06
!-----
      IF ( Idimen == 1 ) THEN
         WRITE (lun,'(a)') '1-D'
      ELSEIF ( Idimen == 2 ) THEN
         WRITE (lun,'(a)') '2-D'
      ELSEIF ( Idimen == 3 ) THEN
         WRITE (lun,'(a)') '3-D'
      ENDIF
!-----
!       LINE 07
!-----
      IF ( Icnvel == 0 ) THEN
         WRITE (lun,'(a)') 'CONSTANT VELOCITY'
      ELSEIF ( Icnvel == 1 ) THEN
         WRITE (lun,'(a)') 'VARIABLE VELOCITY'
      ENDIF
!-----
!       put lines 8 through 11
!-----
      WRITE (lun,'(a)') 'LINE08'
      WRITE (lun,'(a)') 'LINE09'
      WRITE (lun,'(a)') 'LINE10'
      WRITE (lun,'(a)') 'LINE11'
!-----
!       put model specifically for 1-D flat isotropic
!-----
!-----
!       put comment line
!-----
      WRITE (lun,'(a)') cmnt(1:110)
      curdep = 0.0
 
      DO j = 1 , Mmax
         curdep = curdep + ABS(D(j))
         IF ( curdep.LE.REFdep ) THEN
            dout = -D(j)
         ELSE
            dout = D(j)
         ENDIF
 
         WRITE (lun,'(4f11.4,6g11.3)') dout , A(j) , B(j) , RHO(j) ,    &
                                     & QA(j) , QB(j) , ETAp(j) , ETAs(j)&
                                     & , FREfp(j) , FREfs(j)
      ENDDO
      IF ( lun /= STDOUT ) CLOSE (lun)
      END

