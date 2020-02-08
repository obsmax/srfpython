!*==SRFPRE96.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017
!  programme preliminaire pour srfdis96
!  entrees d'origine
!     s'il existe, lit le fichier sobs.d, les modeles de profondeur et dispersion qui y figurent
!     sinon : creation interactive (i.e. STDIN) du fichier sobs.d, du modele de profondeur et du fichier de dispersion
!     le modele de profondeur sert a encadrer les vitesses de phase
!     le modele de dispersion sert a donner les type d'onde et modes desires
!     les parametres de dispersion donnent le nmbre max de modes a calculer si ils apparaissent dans le fichier de dispersion
!  sorties d'origine (merdier)
!     le fichier sobs.d
!     le fichier modele si cree, le fichier de dispersion si cree
!     tmpsrfi.03   contient les infos necessaires a srfdis96
!     tmpsrfi.17   copie inutile du fichier model
!     tmpsrfi.07, tmpsrfi.04, tmpsrfi.12, tmpsrfi.00, tmpmod96.000 tous inutiles pour srfdis96
!***********
!  mes modifs
!   => le programme original ecrit sobs.d (sil nexiste pas) pour le relire juste apres
!   => il ecrit les fichier modele prof et dispersion pour les relire juste apres
!   => j'ai vire ca, mnt les parametres du fichier sobs.d ainsi que le modele de profondeur
!   => et le fichier de dispersion sont lus dans stdin
!   desormais toutes les entrees passent par stdin, et les sorties par stdout a destination de srfdis96
!   nvx inputs :
!      ligne 1:
!      h dcl dcr
!      avec : h,    float, increment en periode pour la conversion phase>groupe (0.005 est raisonable)
!             dcl,  float, increment en vitesse de phase pour la recherche des racines de love
!             dcr,  float, comme dcl pour rayleigh
!      lignes suivantes :
!      le modele de profondeur au format mod96, ne pas changer les lignes de header,
!      j'ai modifie igetmod pour que la couche ac H=0 soit le signal de fin de lecture
!      exemple
!MODEL.01
!modeltitle
!ISOTROPIC
!KGS
!FLAT EARTH
!1-D
!CONSTANT VELOCITY
!LINE08
!LINE09
!LINE10
!LINE11
!      H(KM)   VP(KM/S)   VS(KM/S) RHO(GM/CC)     QP         QS       ETAP       ETAS      FREFP      FREFS
!     0.5000     1.0000     0.8000     2.5000   0.00       0.00       0.00       0.00       1.00       1.00
!     0.5000     1.5000     1.2000     2.6000   0.00       0.00       0.00       0.00       1.00       1.00
!     0.0000     1.6000     1.3000     2.6500   0.00       0.00       0.00       0.00       1.00       1.00
!      ligne suivante:
!      NLC NLU NRC NRU
!      avec : NLC,  int, nombre max de mode a calculer pour Love-phase, 0 = aucun, 1 = fondamental seulement
!             NLU,  int, comme NLC pour le groupe
!             NRC,  int, comme NLC pour le rayleigh
!             NRU,  int, comme NRC pour le groupe
!      lignes suivantes :
!      le fichier de dispersion de reference au format surf96 qui donne les periodes, ondes, types et modes
!      attention NLC,NLU,NRC,NRU n'ont d effet que si les donnees existent dans le fichier de ref.
!      exemple: pas de header, cle SURF96, WAVE, TYPE, FLAG (X), MODE, PERIOD(s), VITESSE(km/s), INCERT(km/s)
!SURF96 R C X   0      1.0000     1.0000     0.1000
!SURF96 R C X   0      1.5000     1.5000     0.1000
!SURF96 R C X   0      1.6000     1.6000     0.1000
!SURF96 R U X   0      1.0000     1.4000     0.1000
!SURF96 R U X   0      1.5000     1.5000     0.1000
!SURF96 R U X   0      1.6000     1.6000     0.1000
!SURF96 R C X   1      1.0000     1.0000     0.1000
!SURF96 R C X   1      1.5000     1.5000     0.1000
!SURF96 R C X   1      1.6000     1.6000     0.1000
!SURF96 R U X   1      1.0000     1.4000     0.1000
!SURF96 R U X   1      1.5000     1.5000     0.1000
!SURF96 R U X   1      1.6000     1.6000     0.1000
!SURF96 L C X   0      1.0000     1.0000     0.1000
!SURF96 L C X   0      1.5000     1.5000     0.1000
!SURF96 L C X   0      1.6000     1.6000     0.1000
!SURF96 L U X   0      1.0000     1.4000     0.1000
!SURF96 L U X   0      1.5000     1.5000     0.1000
!SURF96 L U X   0      1.6000     1.6000     0.1000
!SURF96 L C X   1      1.0000     1.0000     0.1000
!SURF96 L C X   1      1.5000     1.5000     0.1000
!SURF96 L C X   1      1.6000     1.6000     0.1000
!SURF96 L U X   1      1.0000     1.4000     0.1000
!SURF96 L U X   1      1.5000     1.5000     0.1000
!SURF96 L U X   1      1.6000     1.6000     0.1000
!      nouvelles sorties:
!      le fichier modele qui doit etre relu par srfdis96 (fait chier..)
!      les parametres de dispersion et pre-calculs indispensables a srfdis96
!      a savoir
!      nper, nper, earthflat (nombre de periodes distinctes, nombre de periodes distinctes, 0 pour terre plate)
      PROGRAM SRFPRE96
      IMPLICIT NONE
!----------------------------------------------------------------------c
!                                                                    c
!      COMPUTER PROGRAMS IN SEISMOLOGY                               c
!      VOLUME IV                                                     c
!                                                                    c
!      PROGRAM: SRFPRE96                                             c
!                                                                    c
!      COPYRIGHT 1986, 1991, 2001                                    c
!      D. R. Russell, R. B. Herrmann                                 c
!      Department of Earth and Atmospheric Sciences                  c
!      Saint Louis University                                        c
!      221 North Grand Boulevard                                     c
!      St. Louis, Missouri 63103                                     c
!      U. S. A.                                                      c
!                                                                    c
!----------------------------------------------------------------------c
!     CHANGES
!     19 JAN 2002 - jobsyn = 1 for observed, 2 for synthetic to
!             be compatible with   f96subf.f
!     05 MAR 2002 - observed = 'X' synthetic = 'T'
!     12 MAR 2003 - increased NM= in getdsp
!
!
!     This program checks the input control file 'sobs.d' and
!     converts the input data into unformatted binary files
!     to be used in other programs.  The unformatted files
!     are labeled as 'tmpsrfi.xx' where xx is a number from
!     0 to 14.
!
!     Developed by David R. Russell, St. Louis University, Jan. 1984.
!     Restructured Input R. B. Herrmann 8 August 1986
!
!     Restructured August 1986 to simplify input format RBH
!     Modified  September 1986 to include gamma values  RBH
!     Modified  November 2001 to use SURF96 dispersion format and
!     model96 earth model format
!     Modified  January  2002 to internally use the model96 files
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!

      INTEGER iunit,NL,STDIN,STDOUT,i
      REAL onel,oner
      REAL*4 h,dcl,dcr
      PARAMETER (STDIN=5,STDOUT=6)
      PARAMETER (NL=100)

      INTEGER NLG,NLC,NLU,NRG,NRC,NRU,nlayer,earthflat

!-----
!     STDIN   - unit for FORTRAN read from terminal
!     STDOUT   - unit for FORTRAN write to terminal
!     STDERR   - unit for FORTRAN error output to terminal
!     NL    - number of layers in model
!     NL2   - number of columns in model (first NL2/2 are
!           - velocity parameters, second NL2/2 are Q values)
!     h     = percentage period change for group velocity partial
!     dcl   = Love Wave Phase Velocity Search Increment
!     onel  = Love Wave Increment for Backtracking on Step
!     dcr   = Rayleigh Wave Phase Velocity Search Increment
!     oner  = Rayleigh Wave Increment for Backtracking on Step

!     nf() is the control array
!        not used nf(1) = 1 estimated stdev computed from residuals
!        not used         0 no scaling by residuals
!     NLG   = nf(2) = TOTAL number of Love Wave Gamma Modes
!             0 DO NOT PROCESS Love Wave Gamma Data for Q
!     NLC   = nf(3) = Maximum Number of Love Wave Phase Velocity Modes
!             0 DO NOT PROCESS Love Wave Phase Velocity
!     NLU   = nf(4) = Maximum Number of Love Wave Group Velocity Modes
!           0 DO NOT PROCESS Love Wave Group Velocity
!     NRG   = nf(5) = Maximum Number of Rayleigh Wave Gamma Modes
!             0 DO NOT PROCESS Rayleigh Wave Gamma Data for Q
!     NRC   = NRC = Maximum Number of Rayleigh Wave Phase Velocity Modes
!             0 DO NOT PROCESS Rayleigh Wave Phase Velocity
!     NRU   = nf(7) = Maximum Number of Rayleigh Wave Group Velocity Modes
!             0 DO NOT PROCESS Rayleigh Wave Group Velocity
!         not used nf(8) = Model Weighting
!         not used         0 No Weighting
!         not used 1 Read In Layer Velocity Weights
!     nlayer=nf(9) = Number of Layers in Starting Model (from model file)
!         not used nf(10)= Input Format (from model file)
!         not used         0 - Inversion a, rho fixed
!         not used         1 - Inversion Poisson Ratio Fixed, Rho computed from Vp
!         not used nf(11)= Type of Smoothing
!         not used         0 None
!         not used         1 Differential
!     earthflat = nf(12)= Earth Flattening
!             0 No Earth Flattening - Model is Flat Earth Model (always)
!             1 Earth Flattening - Model is Spherical Earth Model (never)
!         not used ifrper = nf(13)= Dispersion Abscissa in Input Dispersion Data
!         not used         (from disp file)
!         not used         0 period (always)
!         not used         1 frequency (never)
!-----
      REAL thicknesses(100),density_values(100),vp_values(100),vs_values(100)
      read(STDIN,"(I3)") nlayer
      read(STDIN,*) thicknesses(1:nlayer-1)
      read(STDIN,*) vp_values(1:nlayer)
      read(STDIN,*) vs_values(1:nlayer)
      read(STDIN,*) density_values(1:nlayer)

      write(STDOUT,"(I3)") nlayer

      do i = 1, nlayer-1
        write(STDOUT,'(F8.3)', advance="no") thicknesses(i)
      end do
      write(STDOUT,'(A1)') ""

      do i = 1, nlayer
        write(STDOUT,'(F7.3)', advance="no") vp_values(i)
      end do
      write(STDOUT,'(A1)') ""

      do i = 1, nlayer
        write(STDOUT,'(F7.3)', advance="no") vs_values(i)
      end do
      write(STDOUT,'(A1)') ""

      do i = 1, nlayer
        write(STDOUT,'(F6.3)', advance="no") density_values(i)
      end do
      write(STDOUT,'(A1)') ""

      READ (STDIN,*) h,dcl,dcr
      oner   = 0.0
      onel   = 0.0
      NLG    = 0
      NRG    = 0
      earthflat = 0

      READ (STDIN,*) NLC,NLU,NRC,NRU

      WRITE(STDOUT,"(F7.4)") h
      CALL GETDSP(&
              & NLG,NLC,NLU,NRG,NRC,NRU,nlayer, &
              & earthflat,dcl,dcr,onel,oner)
      END


! #################################################################
      SUBROUTINE GETDSP(NLG,NLC,NLU,NRG,NRC,NRU,nlayer, &
              & earthflat,Dcl,Dcr,Onel,Oner)
      IMPLICIT NONE

      REAL c,cper,dc,Dcl,Dcr,f,obs,obserr,one,    &
         & Onel,Oner,pper,sd
      INTEGER i,idat,igr,ilorr,ilvry,imode,iobs, &
            & iobsyn,iporg,Iunitd,j,k,kmax,LGSTR,STDIN ,&
            & lnobl,STDOUT,ls
      INTEGER lsep,mm,n,nlr,nlrr,NM,nmgr,nmod,nmph ,&
            & NP,nper,nx
      PARAMETER (NM=5000,STDOUT=6,NP=512)
      PARAMETER (STDIN=5)
      INTEGER NLG,NLC,NLU,NRG,NRC,NRU,nlayer,earthflat
      REAL*4 tper(NM),vel(NM),dvel(NM)
      INTEGER*4 lorr(NM),porg(NM),mode(NM),modemx(2,3)
      REAL*4 per(NP),tmp(NM)
      INTEGER*4 key(NM),imap(NM),jmap(NM)
      CHARACTER instr*132
      CHARACTER ic*1
      DATA modemx/0,0,0,0,0,0/
!-----
!     MEANING OF VARIABLES
!     STDIN   - unit for FORTRAN read from terminal
!     STDOUT   - unit for FORTRAN write to terminal
!     STDERR   - unit for FORTRAN error output to terminal
!     NL    - number of layers in model
!     NL2   - number of columns in model (first NL2/2 are
!           - velocity parameters, second NL2/2 are Q values)
!     NP    - number of unique periods
!     NM    - maximum number of observations
!     lorr, porg, mode, per, vel, dvel
!           input parameters for each observation
!           (L/R) (C/U) MODE FREQ VEL SDVEL
!     idat - number of observations
!     per - array of unique periods
!     nper - number of unique periods
!     imap - mapping of observation into unique period
!     modemx(1,1) = max modes love, phase vel, gamma
!     modemx(1,2) = max modes love, group vel
!     modemx(2,1) = max modes rayl, phase vel, gamma
!     modemx(2,2) = max modes rayl, group vel
!
!     to incorporate gamma, we need phase velocity partials
!           for the inversion, so the gamma input will
!           be considered phase for the phase velocity
!           determination
!
!     tmp, key arrays for sort algorithm
!     jmap - temporary mapping
!-----
!     read in data, store in arrays
!     get units and data type
!     NEW units are always km and period(sec)
!-----
      idat = 0
      Iunitd = 0

 100  READ (STDIN,'(a)',END=200) instr
!        WRITE(0,*)idat,' ',instr
      ls = LGSTR(instr)
!-----
!         do the parsing
!-----
      IF ( instr(1:6).EQ.'SURF96' .OR. instr(1:6).EQ.'surf96' ) THEN
!-----
!             now get wave type
!-----
         lsep = 6
         CALL GETBLNK(instr,lsep,ls,lnobl)
         ic = instr(lnobl:lnobl)
         IF ( ic(1:1).EQ.'R' .OR. ic(1:1).EQ.'r' ) THEN
            ilorr = 2
         ELSEIF ( ic(1:1).EQ.'L' .OR. ic(1:1).EQ.'l' ) THEN
            ilorr = 1
         ELSEIF ( ic(1:1).EQ.'A' .OR. ic(1:1).EQ.'a' ) THEN
            ilorr = 0
         ELSE
            GOTO 100
         ENDIF
!-----
!             now get observation type
!-----
         lsep = lnobl
         CALL GETBLNK(instr,lsep,ls,lnobl)
         ic = instr(lnobl:lnobl)
         IF ( ic(1:1).EQ.'C' .OR. ic(1:1).EQ.'c' ) THEN
            iobs = 1
         ELSEIF ( ic(1:1).EQ.'U' .OR. ic(1:1).EQ.'u' ) THEN
            iobs = 2
         ELSEIF ( ic(1:1).EQ.'G' .OR. ic(1:1).EQ.'g' ) THEN
            iobs = 3
         ELSE
            GOTO 100
         ENDIF
!-----
!             now get whether observed or synthetic
!-----
         lsep = lnobl
         CALL GETBLNK(instr,lsep,ls,lnobl)
         ic = instr(lnobl:lnobl)
         IF ( ic(1:1).EQ.'T' .OR. ic(1:1).EQ.'t' ) THEN
            iobsyn = 2
         ELSEIF ( ic(1:1).EQ.'F' .OR. ic(1:1).EQ.'f' ) THEN
            iobsyn = 1
         ELSEIF ( ic(1:1).EQ.'X' .OR. ic(1:1).EQ.'x' ) THEN
            iobsyn = 1
         ELSE
            GOTO 100
         ENDIF
!-----
!-----
!             now get the values using list directed IO
!-----
         lsep = lnobl
         CALL GETBLNK(instr,lsep,ls,lnobl)
         READ (instr(lnobl:ls),*) imode,pper,obs,obserr
      ENDIF


      n = imode + 1
      f = pper
      c = obs
      dc = obserr
!-----
!     ilorr     = Love (1) or Rayleigh
!     iobs     = Phase (1) Group (2) Gamma(3)
!     n     = mode Fund = 1, First = 2, etc
!     f     = Frequency or Period
!     c     = Velocity of Gamma depending on iobs
!     dc    = Error in Velocity of Gamma, depending on iobs
!-----
!     before adding to the data set, ensure that the
!     data are to be used
!-----
!     increment n for internal use
!-----
!      IF ( ilorr.NE.1 .OR. iobs.NE.1 .OR. n.LE.NLC ) THEN
!         IF ( ilorr.NE.1 .OR. iobs.NE.2 .OR. n.LE.NLU ) THEN
!            IF ( ilorr.NE.1 .OR. iobs.NE.3 .OR. n.LE.NLG ) THEN
!               IF ( ilorr.NE.2 .OR. iobs.NE.1 .OR. n.LE.NRC ) THEN
!                  IF ( ilorr.NE.2 .OR. iobs.NE.2 .OR. n.LE.NRU ) THEN
!                     IF ( ilorr.NE.2 .OR. iobs.NE.3 .OR. n.LE.NRG ) THEN
!                     ENDIF
!                  ENDIF
!               ENDIF
!            ENDIF
!         ENDIF
!      ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ilorr.NE.1 <=> ilorr.EQ.2  and inversly
      IF (( ilorr.EQ.2 .OR. iobs.NE.1 .OR. n.LE.NLC ) .AND. &
         &( ilorr.EQ.2 .OR. iobs.NE.2 .OR. n.LE.NLU ) .AND. &
         &( ilorr.EQ.2 .OR. iobs.NE.3 .OR. n.LE.NLG ) .AND. &
         &( ilorr.EQ.1 .OR. iobs.NE.1 .OR. n.LE.NRC ) .AND. &
         &( ilorr.EQ.1 .OR. iobs.NE.2 .OR. n.LE.NRU ) .AND. &
         &( ilorr.EQ.1 .OR. iobs.NE.3 .OR. n.LE.NRG )) THEN

            idat = idat + 1
            lorr(idat) = ilorr
            porg(idat) = iobs
            mode(idat) = n
            !-----
            !     SURF96 input is always period!!!
            !     SURF96 DISPERSION UNITS ARE ALWAYS km/sec and 1/km
            !-----
            tper(idat) = f
            vel(idat) = c
            IF ( dc.EQ.0.0 ) dc = 1.0
            dvel(idat) = dc
            key(idat) = idat
            tmp(idat) = tper(idat)
            !!-----
            !!     make gamma seem to be phase data
            !!-----
            !                        mm = iobs
            !                        IF ( mm.EQ.3 ) mm = 1
            !                        IF ( n.GT.modemx(ilorr,mm) ) modemx(ilorr,mm) = n
            IF ( n.GT.modemx(ilorr,iobs) ) modemx(ilorr,iobs) = n
      ENDIF

      GOTO 100

 200  CALL SORT(tmp,key,idat)
      CALL UNIQ(per,tper,key,idat,nper,imap)
!-----
!     now look at Love/Rayleigh
!     followed by Mode
!-----
!     this is perhaps an excessive linear search, but it
!     will work
!-----
!     fix up period count
!-----
!      WRITE (STDOUT,*) nper,nper,earthflat
      WRITE (STDOUT,'(I4)') nper
!-----
!     adjust NLC, NLU, NRC, NRU for internal use
!-----
      NLC = modemx(1,1)
      NLU = modemx(1,2)
      NRC = modemx(2,1)
      NRU = modemx(2,2)
!     here order output for several control files
!
!     For observed data, we order as
!           LOVE-RAYL
!                 PHASE(or GAMMA)  - GROUP
!                       MODE
!                             write ilvry,iporg,nmod,per,val,dval
!     For disperion computation order the output as
!           LOVE-RAYL
!                 MODE
!                       range of periods to compute
!                       write range
!-----
      DO ilvry = 1,2
         IF ( ilvry.EQ.1 ) THEN
            nmph = modemx(1,1) ! NLC
            nmgr = modemx(1,2) ! NLU
            one = Onel
            dc = Dcl
         ELSE
            nmph = modemx(2,1) ! NRC
            nmgr = modemx(2,2) ! NRU
            one = Oner
            dc = Dcr
         ENDIF
!-----
!                 ENFORCE USER MODE LIMITS
!-----
         kmax = nper ! < you son of a bitch
         IF ( nmgr.EQ.0 ) igr = 0
         IF ( nmph.EQ.0 ) igr = 1
!           if(nmgr.gt.0 .and. nmph.gt.0 .and. nmgm.gt.0)igr=2
         IF ( nmgr.GT.0 .AND. nmph.GT.0 ) igr = 2
         nx = MAX(nmph,nmgr)
!         WRITE (STDOUT,*) kmax,nx,dc,one,igr!,H
         WRITE (STDOUT,"(I4,I4,F7.4,F7.4,I4)") kmax,nx,dc,one,igr!,H
!        WRITE (STDOUT,*, advance="no") (per(i),i=1,kmax)
         do i = 1, kmax
             WRITE (STDOUT,"(F8.3)", advance="no") per(i)
         end do
         WRITE (STDOUT,"(A1)") ""

         DO iporg = 1,2
            DO nmod = 1,modemx(ilvry,iporg)
               nlr = 0
               DO i = 1,idat
                  IF ( lorr(i).EQ.ilvry .AND. MOD(porg(i),2)            &
                     & .EQ.MOD(iporg,2) .AND. mode(i).EQ.nmod ) THEN
                     nlr = nlr + 1
                     tmp(nlr) = tper(i)
                     key(nlr) = nlr
                     jmap(nlr) = i
                  ENDIF
               ENDDO
               IF ( nlr.GT.0 ) THEN
                  CALL SORT(tmp,key,nlr)
                  nlrr = nlr
                  DO i = 1,nlr
                     j = jmap(key(i))
                     k = imap(j)

                  ENDDO
               ELSE
                  key(1) = 1
                  nlrr = 1
                  c = 0
                  sd = 1
                  cper = 0.0

               ENDIF
            ENDDO
         ENDDO
!-----
!     for Love or Rayleigh find the period limits
!     for each mode, so that an inclusive comb is constructed
!     for phase velocity search
!-----
         CALL GETLIM(modemx,idat,lorr,mode,imap,ilvry)
      ENDDO
!        close(4,status='keep')
      END
!*==GETLIM.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017

!#################################################################
      SUBROUTINE GETLIM(Modemx,Idat,Lorr,Mode,Imap,Ilvry)
      IMPLICIT NONE
!*--GETLIM572
!*** Start of declarations inserted by SPAG
      INTEGER i,Ilvry,im,j,STDOUT,md,n
!*** End of declarations inserted by SPAG
      INTEGER*4 Modemx(2,2),Idat,Lorr(*),Imap(*),Mode(*)
!-----
!     get limits on dispersion periods for dispersion program
!     to speed determination of higher modes, we develop
!     an inclusive comb of periods to be evaluated for
!     each mode such that the number of periods at
!     the next higher mode is always within the
!     range of the previous mode
!
!     to do this trick, we work backwords and then output the
!     desired results
!
!-----
      PARAMETER (STDOUT=6)
      INTEGER*4 is(100),ie(100)
      DATA is/100*0/,ie/100*0/
      md = 0
      DO i = 1,2
         IF ( Modemx(Ilvry,i).GT.md ) md = Modemx(Ilvry,i)
      ENDDO
!-----
!     perform linear searches for simplicity
!-----
      DO n = md,1,-1
         DO j = 1,Idat
            IF ( Mode(j).EQ.n .AND. Lorr(j).EQ.Ilvry ) THEN
               im = Imap(j)
               IF ( is(n).EQ.0 .OR. is(n).GT.im ) is(n) = im
               IF ( ie(n).EQ.0 .OR. ie(n).LT.im ) ie(n) = im
            ENDIF
         ENDDO
      ENDDO
!-----
!     fill out comb
!-----
      DO n = md,2,-1
         IF ( is(n).LT.is(n-1) ) is(n-1) = is(n)
         IF ( is(n-1).EQ.0 ) is(n-1) = is(n)
         IF ( ie(n).GT.ie(n-1) ) ie(n-1) = ie(n)
      ENDDO
!-----
!     output starting with the first mode
!-----
      DO n = 1,md
         WRITE (STDOUT,"(I4,I4)") is(n),ie(n)
      ENDDO
      END

!#################################################################
      SUBROUTINE UNIQ(Y,X,Key,Nx,Ny,Imap)
      IMPLICIT NONE
!*--UNIQ629
!*** Start of declarations inserted by SPAG
      INTEGER i,j,Nx,Ny
!*** End of declarations inserted by SPAG

!-----
!     this subroutine takes a sorted list, x(key(i))
!     and determines the unique elements y() in it
!     and returns the unique number ny
!     imap(i) = ny maps original into unique
!-----
      REAL*4 Y(*),X(*)
      INTEGER*4 Key(*),Imap(*)
!        WRITE(0,*)'nx,ny,imap:',nx,ny,(imap(j),j=1,nx)
!        WRITE(0,*)'x:',(x(j),j=1,nx)
!        WRITE(0,*)'key:',(key(j),j=1,nx)
!-----
!     the first element is unique
!-----
      Ny = 1
      Y(Ny) = X(Key(1))
      Imap(Key(1)) = Ny
      DO i = 1,Nx
         j = Key(i)
         IF ( Y(Ny).LT.X(j) ) THEN
            Ny = Ny + 1
            Y(Ny) = X(j)
         ENDIF
         Imap(j) = Ny
      ENDDO
      END
!*==SORT.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017

!#################################################################
      SUBROUTINE SORT(X,Key,N)
      IMPLICIT NONE
!*--SORT665
!*** Start of declarations inserted by SPAG
      INTEGER i,j,ktmp
      REAL tmp
!*** End of declarations inserted by SPAG
!-----
!     Starting with x(1) ,,, x(n)
!     return   the xarray sorted in increasing order
!     also return the pointers key to the initial array.
!     For example given x = [ 3, 1, 2 ]
!     the returned values are
!                       x = [ 1, 2, 3 ]
!                     key = [ 2, 3, 1 ]
!-----
!        Reference: http://en.wikipedia.org/wiki/Bubble_sort
!-----
      INTEGER N
      REAL X(N)
      INTEGER Key(N)
      DO i = 1,N
         Key(i) = i
      ENDDO
      DO i = N,1,-1
         DO j = 1,i - 1
            IF ( X(j).GT.X(j+1) ) THEN
               tmp = X(j)
               X(j) = X(j+1)
               X(j+1) = tmp
               ktmp = Key(j)
               Key(j) = Key(j+1)
               Key(j+1) = ktmp
            ENDIF
         ENDDO
      ENDDO


      END
!*==GETBLNK.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017
!#################################################################
      SUBROUTINE GETBLNK(Instr,Lsep,Ls,Lnobl)
      IMPLICIT NONE
!*--GETBLNK706
!*** Start of declarations inserted by SPAG
      INTEGER i,igotit
!*** End of declarations inserted by SPAG
!-----
!     determine first non-blank character
!
!     instr   Ch* Character string to be parsed
!     lsep    I*4 index of last non blank character
!     ls  I*4 length of input string
!     lnobl   I*4 index of first non blank character
!-----
      CHARACTER Instr*(*)
      INTEGER Lsep,Ls,Lnobl
      CHARACTER tab*1
      tab = CHAR(9)
      Lnobl = Lsep + 1
      igotit = 0
      DO i = Lsep + 1,Ls
         IF ( igotit.EQ.0 ) THEN
            IF ( Instr(i:i).NE.' ' .AND. Instr(i:i).NE.tab ) THEN
               Lnobl = i
               igotit = 1
            ENDIF
         ENDIF
      ENDDO
      END

