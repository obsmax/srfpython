!*==SRFPRE96.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017
!  programme preliminaire pour srfdis96
!  entrees d'origine
!     s'il existe, lit le fichier sobs.d, les modeles de profondeur et dispersion qui y figurent
!     sinon : creation interactive (i.e. STDIN) du fichier sobs.d, du modele de profondeur et du fichier de dispersion
!     le modele de profondeur sert a encadrer les vitesses de phase
!     le modele de dispersion sert a donner les type d'onde et modes desires
!     les parametres de dispersion donnent le nmbre max de mode a calculer ssi ils apparaissent dans le fichier de dispersion
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
      PROGRAM SRFPRE96
      IMPLICIT NONE
!*--SRFPRE9686
!*** Start of declarations inserted by SPAG
      INTEGER iflsph , iunit , iunitd , LIN , m , m2 , NL
      REAL onel , oner
!*** End of declarations inserted by SPAG
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
      PARAMETER (LIN=5)
      PARAMETER (NL=200)
!-----
!     LIN   - unit for FORTRAN read from terminal
!     LOT   - unit for FORTRAN write to terminal
!     LER   - unit for FORTRAN error output to terminal
!     NL    - number of layers in model
!     NL2   - number of columns in model (first NL2/2 are
!           - velocity parameters, second NL2/2 are Q values)
!-----
!        integer nf10(NL)
!        data nf10/NL*1/
      INTEGER nf(13)
!        real dd(NL2)
!        logical wc(NL2)
!        dimension a(NL),b(NL),d(NL),r(NL),rat(NL),qbinv(NL)
!        character*80 nmmodl, nmdisp
      CHARACTER*80 nmdisp
!        logical ext
!        common/param/qaqb,itype,dlam,invdep
!        data qbinv/NL *0.0/
      REAL*4 h , dcl , dcr
!-----
!     h     = percentage period change for group velocity partial
!     dcl   = Love Wave Phase Velocity Search Increment
!     onel  = Love Wave Increment for Backtracking on Step
!     dcr   = Rayleigh Wave Phase Velocity Search Increment
!     oner  = Rayleigh Wave Increment for Backtracking on Step
!-----
      READ (LIN,*) h , dcl , dcr
      oner = 0.0
      onel = 0.0
!-----
!     nf() is the control array
!     nf(1) = 1 estimated stdev computed from residuals
!             0 no scaling by residuals
!     nf(2) = TOTAL number of Love Wave Gamma Modes
!             0 DO NOT PROCESS Love Wave Gamma Data for Q
!     nf(3) = Maximum Number of Love Wave Phase Velocity Modes
!             0 DO NOT PROCESS Love Wave Phase Velocity
!     nf(4) = Maximum Number of Love Wave Group Velocity Modes
!           0 DO NOT PROCESS Love Wave Group Velocity
!     nf(5) = Maximum Number of Rayleigh Wave Gamma Modes
!             0 DO NOT PROCESS Rayleigh Wave Gamma Data for Q
!     nf(6) = Maximum Number of Rayleigh Wave Phase Velocity Modes
!             0 DO NOT PROCESS Rayleigh Wave Phase Velocity
!     nf(7) = Maximum Number of Rayleigh Wave Group Velocity Modes
!             0 DO NOT PROCESS Rayleigh Wave Group Velocity
!     nf(8) = Model Weighting
!             0 No Weighting
!             1 Read In Layer Velocity Weights
!     nf(9) = Number of Layers in Starting Model (from model file)
!     nf(10)= Input Format (from model file)
!             0 - Inversion a, rho fixed
!             1 - Inversion Poisson Ratio Fixed, Rho computed from Vp
!     nf(11)= Type of Smoothing
!             0 None
!             1 Differential
!     nf(12)= Earth Flattening
!             0 No Earth Flattening - Model is Flat Earth Model
!             1 Earth Flattening - Model is Spherical Earth Model
!     nf(13)= Dispersion Abscissa in Input Dispersion Data
!             (from disp file)
!             0 period
!             1 frequency
!-----
      nf(1) = 0
      nf(2) = 0
      nf(5) = 0
      nf(8) = 0
      nf(11) = 0
      nf(12) = 0
!-----
!     nmmodl= name of file containing model information
!-----
!        call getmdl(nmmodl,nf(9),nf(10),iunit,iflsph)
      CALL GETMDL('STDIN',nf(9),nf(10),iunit,iflsph)
      nf(12) = iflsph
      m = nf(9)
      m2 = m + m
      READ (LIN,*) nf(3) , nf(4) , nf(6) , nf(7)
!-----
!     nmdisp= name of file containing dispersion data
!-----
!
!     main loop...check input data and write to unformatted files
!
      CALL GETDSP(nmdisp,nf,dcl,dcr,onel,oner,h,iunitd)
      END
!*==GETMDL.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017
!#################################################################
      SUBROUTINE GETMDL(Nmmodl,Mmax,Nfmod,Iunit,Iflsph)
      IMPLICIT NONE
!*--GETMDL223
!*** Start of declarations inserted by SPAG
      REAL DL , ETAp , ETAs , FREfp , FREfs , QA , QB , REFdep , RHO ,  &
         & VA , VB
      INTEGER icnvel , idimen , ierr , Iflsph , iiso , Iunit , LGSTR ,  &
            & lt , Mmax , Nfmod
!*** End of declarations inserted by SPAG
!-----
!     igetmod common
!-----
      INTEGER NL , NLAY
      PARAMETER (NL=200,NLAY=200)

      COMMON /ISOMOD/ DL(NLAY) , VA(NLAY) , VB(NLAY) , RHO(NLAY) ,      &
                    & QA(NLAY) , QB(NLAY) , ETAp(NLAY) , ETAs(NLAY) ,   &
                    & FREfp(NLAY) , FREfs(NLAY)
      COMMON /DEPREF/ REFdep
      CHARACTER Nmmodl*(*)
      COMMON /MODTIT/ TITle
      CHARACTER*80 TITle
!-----
!     open the model file. These values will be saved
!     in internal files so that the original model is not
!     modified in any way
!-----
      CALL GETMOD(2,Nmmodl,Mmax,TITle,Iunit,iiso,Iflsph,idimen,icnvel,  &
                & ierr,.FALSE.)
!----
!     mmax    nf9
!
!     nfmod   nf10    1 invert for S - P, Rho calculated from S
!             0 invert for S but P, Rho fixed
!-----
      Nfmod = 1
      Iunit = 0
      lt = LGSTR(TITle)
      CALL PUTMOD(2,'stdout',Mmax,TITle(1:lt),Iunit,iiso,Iflsph,idimen, &
                & icnvel,.FALSE.)
!        write(LOT,*) mmax,title(1:lt),iunit,iiso,iflsph,
!     1      idimen,icnvel,.false.
      END
!*==GETDSP.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017


!#################################################################
      SUBROUTINE GETDSP(Nmdisp,Nf,Dcl,Dcr,Onel,Oner,H,Iunitd)
      IMPLICIT NONE
!*--GETDSP270
!*** Start of declarations inserted by SPAG
      REAL c , cper , dc , Dcl , Dcr , f , H , obs , obserr , one ,     &
         & Onel , Oner , pper , sd
      INTEGER i , idat , ifrper , igr , ilorr , ilvry , imode , iobs ,  &
            & iobsyn , iporg , Iunitd , j , k , kmax , l , LGSTR , LIN ,&
            & lnobl , LOT , ls
      INTEGER lsep , m , mm , n , nlr , nlrr , NM , nmgr , nmod , nmph ,&
            & NP , nper , nx
!*** End of declarations inserted by SPAG
!-----
!     nmdisp - file containing dispersion data
!     nf - integer array containing control flags
!-----
      PARAMETER (NM=5000,LOT=6,NP=512)
!-----
!     LIN   - unit for FORTRAN read from terminal
!     LOT   - unit for FORTRAN write to terminal
!     LER   - unit for FORTRAN error output to terminal
!     NL    - number of layers in model
!     NL2   - number of columns in model (first NL2/2 are
!           - velocity parameters, second NL2/2 are Q values)
!     NP    - number of unique periods
!     NM    - maximum number of observations
!-----
      PARAMETER (LIN=5)
      CHARACTER Nmdisp*(*)
      INTEGER Nf(13)
      REAL*4 tper(NM) , vel(NM) , dvel(NM)
      INTEGER*4 lorr(NM) , porg(NM) , mode(NM) , modemx(2,3)
      REAL*4 per(NP) , tmp(NM)
      INTEGER*4 key(NM) , imap(NM) , jmap(NM)
      CHARACTER instr*132
      CHARACTER ic*1
      DATA modemx/0 , 0 , 0 , 0 , 0 , 0/
!-----
!     MEANING OF VARIABLES
!
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
!-----
!     read in data, store in arrays
!-----
      idat = 0
!       open(1,file=nmdisp,status='old',form='formatted',
!    1            access='sequential')
!       rewind 1
!        open(4,file='tmpsrfi.03',form='unformatted',access='sequential')
!        rewind 4
!-----
!     get units and data type
!     NEW units are always km and period(sec)
!-----
      Iunitd = 0
      ifrper = 0
      Nf(13) = ifrper
!           read(1,'(a)',end=1001)instr
 100  READ (LIN,'(a)',END=200) instr
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
         READ (instr(lnobl:ls),*) imode , pper , obs , obserr
      ENDIF
      l = ilorr
      m = iobs
      n = imode + 1
      f = pper
      c = obs
      dc = obserr
!-----
!     l     = Love (1) or Rayleigh
!     m     = Phase (1) Group (2) Gamma(3)
!     n     = mode Fund = 1, First = 2, etc
!     f     = Frequency or Period
!     c     = Velocity of Gamma depending on m
!     dc    = Error in Velocity of Gamma, depending on m
!-----
!     before adding to the data set, ensure that the
!     data are to be used
!-----
!     increment n for internal use
!-----
      IF ( l.NE.1 .OR. m.NE.1 .OR. n.LE.Nf(3) ) THEN
         IF ( l.NE.1 .OR. m.NE.2 .OR. n.LE.Nf(4) ) THEN
            IF ( l.NE.1 .OR. m.NE.3 .OR. n.LE.Nf(2) ) THEN
               IF ( l.NE.2 .OR. m.NE.1 .OR. n.LE.Nf(6) ) THEN
                  IF ( l.NE.2 .OR. m.NE.2 .OR. n.LE.Nf(7) ) THEN
                     IF ( l.NE.2 .OR. m.NE.3 .OR. n.LE.Nf(5) ) THEN
                        idat = idat + 1
                        lorr(idat) = l
                        porg(idat) = m
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
!-----
!     make gamma seem to be phase data
!-----
                        mm = m
                        IF ( mm.EQ.3 ) mm = 1
                        IF ( n.GT.modemx(l,mm) ) modemx(l,mm) = n
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      GOTO 100
!       close (1)
!        WRITE(0,*)'idat:',idat
 200  CALL SORT(tmp,key,idat)
!        WRITE(0,*)'idat:',idat
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
      WRITE (LOT,*) nper , nper , Nf(12)
!        write(4) nper,nper,nf(12)
!-----
!     adjust nf(3) nf(4) nf(5) nf(6) for internal use
!-----
      Nf(3) = modemx(1,1)
      Nf(4) = modemx(1,2)
      Nf(6) = modemx(2,1)
      Nf(7) = modemx(2,2)
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
      DO ilvry = 1 , 2
         IF ( ilvry.EQ.1 ) THEN
            nmph = modemx(1,1)
            nmgr = modemx(1,2)
            one = Onel
            dc = Dcl
         ELSE
            nmph = modemx(2,1)
            nmgr = modemx(2,2)
            one = Oner
            dc = Dcr
         ENDIF
!-----
!                 ENFORCE USER MODE LIMITS
!-----
         kmax = nper
         IF ( nmgr.EQ.0 ) igr = 0
         IF ( nmph.EQ.0 ) igr = 1
!           if(nmgr.gt.0 .and. nmph.gt.0 .and. nmgm.gt.0)igr=2
         IF ( nmgr.GT.0 .AND. nmph.GT.0 ) igr = 2
         nx = MAX(nmph,nmgr)
         WRITE (LOT,*) kmax , nx , dc , one , igr , H
         WRITE (LOT,*) (per(i),i=1,kmax)
!              write(4) kmax,nx,dc,one,igr,h
!              write(4) (per(i),i=1,kmax)
         DO iporg = 1 , 2
            DO nmod = 1 , modemx(ilvry,iporg)
               nlr = 0
               DO i = 1 , idat
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
                  DO i = 1 , nlr
                     j = jmap(key(i))
                     k = imap(j)
!                                   write(LOT,*)vel(j),
!     1                                       dvel(j),k,
!     2                                       per(k)
!     3                                       ,tper(j),porg(j)
!     4                                          ,ilvry,iporg,
!     5                                          nmod
!                                      write(8)ilvry,porg(j),
!     1                                       nmod,tper(j),
!     2                                       vel(j),dvel(j)
                  ENDDO
               ELSE
                  key(1) = 1
                  nlrr = 1
                  c = 0
                  sd = 1
                  cper = 0.0
!                                      write(8)ilvry,iporg,
!     1                                     nmod,cper,
!     2                                     c,sd
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
      INTEGER i , Ilvry , im , j , LOT , md , n
!*** End of declarations inserted by SPAG
      INTEGER*4 Modemx(2,2) , Idat , Lorr(*) , Imap(*) , Mode(*)
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
      PARAMETER (LOT=6)
      INTEGER*4 is(100) , ie(100)
      DATA is/100*0/ , ie/100*0/
      md = 0
      DO i = 1 , 2
         IF ( Modemx(Ilvry,i).GT.md ) md = Modemx(Ilvry,i)
      ENDDO
!-----
!     perform linear searches for simplicity
!-----
      DO n = md , 1 , -1
         DO j = 1 , Idat
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
      DO n = md , 2 , -1
         IF ( is(n).LT.is(n-1) ) is(n-1) = is(n)
         IF ( is(n-1).EQ.0 ) is(n-1) = is(n)
         IF ( ie(n).GT.ie(n-1) ) ie(n-1) = ie(n)
      ENDDO
!-----
!     output on unit 4 starting with the first mode
!-----
      DO n = 1 , md
         WRITE (LOT,*) is(n) , ie(n)
!              write(4)is(n),ie(n)
      ENDDO
      END
!*==UNIQ.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017

!#################################################################
      SUBROUTINE UNIQ(Y,X,Key,Nx,Ny,Imap)
      IMPLICIT NONE
!*--UNIQ629
!*** Start of declarations inserted by SPAG
      INTEGER i , j , Nx , Ny
!*** End of declarations inserted by SPAG

!-----
!     this subroutine takes a sorted list, x(key(i))
!     and determines the unique elements y() in it
!     and returns the unique number ny
!     imap(i) = ny maps original into unique
!-----
      REAL*4 Y(*) , X(*)
      INTEGER*4 Key(*) , Imap(*)
!        WRITE(0,*)'nx,ny,imap:',nx,ny,(imap(j),j=1,nx)
!        WRITE(0,*)'x:',(x(j),j=1,nx)
!        WRITE(0,*)'key:',(key(j),j=1,nx)
!-----
!     the first element is unique
!-----
      Ny = 1
      Y(Ny) = X(Key(1))
      Imap(Key(1)) = Ny
      DO i = 1 , Nx
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
      INTEGER i , j , ktmp
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
      DO i = 1 , N
         Key(i) = i
      ENDDO
      DO i = N , 1 , -1
         DO j = 1 , i - 1
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
      INTEGER i , igotit
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
      INTEGER Lsep , Ls , Lnobl
      CHARACTER tab*1
      tab = CHAR(9)
      Lnobl = Lsep + 1
      igotit = 0
      DO i = Lsep + 1 , Ls
         IF ( igotit.EQ.0 ) THEN
            IF ( Instr(i:i).NE.' ' .AND. Instr(i:i).NE.tab ) THEN
               Lnobl = i
               igotit = 1
            ENDIF
         ENDIF
      ENDDO
      END

