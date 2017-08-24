c  programme preliminaire pour srfdis96
c  entrees d'origine
c     s'il existe, lit le fichier sobs.d, les modeles de profondeur et dispersion qui y figurent
c     sinon : creation interactive (i.e. STDIN) du fichier sobs.d, du modele de profondeur et du fichier de dispersion 
c     le modele de profondeur sert a encadrer les vitesses de phase
c     le modele de dispersion sert a donner les type d'onde et modes desires
c     les parametres de dispersion donnent le nmbre max de mode a calculer ssi ils apparaissent dans le fichier de dispersion
c  sorties d'origine (merdier)
c     le fichier sobs.d
c     le fichier modele si cree, le fichier de dispersion si cree
c     tmpsrfi.03   contient les infos necessaires a srfdis96
c     tmpsrfi.17   copie inutile du fichier model
c     tmpsrfi.07, tmpsrfi.04, tmpsrfi.12, tmpsrfi.00, tmpmod96.000 tous inutiles pour srfdis96
c***********
c  mes modifs
c   => le programme original ecrit sobs.d (sil nexiste pas) pour le relire juste apres
c   => il ecrit les fichier modele prof et dispersion pour les relire juste apres
c   => j'ai vire ca, mnt les parametres du fichier sobs.d ainsi que le modele de profondeur
c   => et le fichier de dispersion sont lus dans stdin
c   desormais toutes les entrees passent par stdin, et les sorties par stdout a destination de srfdis96
c   nvx inputs :
c      ligne 1:
c      h dcl dcr
c      avec : h,    float, increment en periode pour la conversion phase>groupe (0.005 est raisonable)
c             dcl,  float, increment en vitesse de phase pour la recherche des racines de love
c             dcr,  float, comme dcl pour rayleigh
c      lignes suivantes :
c      le modele de profondeur au format mod96, ne pas changer les lignes de header, 
c      j'ai modifie igetmod pour que la couche ac H=0 soit le signal de fin de lecture
c      exemple
cMODEL.01
cmodeltitle
cISOTROPIC
cKGS
cFLAT EARTH
c1-D
cCONSTANT VELOCITY
cLINE08
cLINE09
cLINE10
cLINE11
c      H(KM)   VP(KM/S)   VS(KM/S) RHO(GM/CC)     QP         QS       ETAP       ETAS      FREFP      FREFS    
c     0.5000     1.0000     0.8000     2.5000   0.00       0.00       0.00       0.00       1.00       1.00    
c     0.5000     1.5000     1.2000     2.6000   0.00       0.00       0.00       0.00       1.00       1.00    
c     0.0000     1.6000     1.3000     2.6500   0.00       0.00       0.00       0.00       1.00       1.00 
c      ligne suivante:
c      NLC NLU NRC NRU      
c      avec : NLC,  int, nombre max de mode a calculer pour Love-phase, 0 = aucun, 1 = fondamental seulement       
c             NLU,  int, comme NLC pour le groupe 
c             NRC,  int, comme NLC pour le rayleigh 
c             NRU,  int, comme NRC pour le groupe
c      lignes suivantes :
c      le fichier de dispersion de reference au format surf96 qui donne les periodes, ondes, types et modes
c      attention NLC,NLU,NRC,NRU n'ont d effet que si les donnees existent dans le fichier de ref.
c      exemple: pas de header, cle SURF96, WAVE, TYPE, FLAG (X), MODE, PERIOD(s), VITESSE(km/s), INCERT(km/s)
cSURF96 R C X   0      1.0000     1.0000     0.1000
cSURF96 R C X   0      1.5000     1.5000     0.1000
cSURF96 R C X   0      1.6000     1.6000     0.1000
cSURF96 R U X   0      1.0000     1.4000     0.1000
cSURF96 R U X   0      1.5000     1.5000     0.1000
cSURF96 R U X   0      1.6000     1.6000     0.1000
cSURF96 R C X   1      1.0000     1.0000     0.1000
cSURF96 R C X   1      1.5000     1.5000     0.1000
cSURF96 R C X   1      1.6000     1.6000     0.1000
cSURF96 R U X   1      1.0000     1.4000     0.1000
cSURF96 R U X   1      1.5000     1.5000     0.1000
cSURF96 R U X   1      1.6000     1.6000     0.1000
cSURF96 L C X   0      1.0000     1.0000     0.1000
cSURF96 L C X   0      1.5000     1.5000     0.1000
cSURF96 L C X   0      1.6000     1.6000     0.1000
cSURF96 L U X   0      1.0000     1.4000     0.1000
cSURF96 L U X   0      1.5000     1.5000     0.1000
cSURF96 L U X   0      1.6000     1.6000     0.1000
cSURF96 L C X   1      1.0000     1.0000     0.1000
cSURF96 L C X   1      1.5000     1.5000     0.1000
cSURF96 L C X   1      1.6000     1.6000     0.1000
cSURF96 L U X   1      1.0000     1.4000     0.1000
cSURF96 L U X   1      1.5000     1.5000     0.1000
cSURF96 L U X   1      1.6000     1.6000     0.1000
c      nouvelles sorties:
c      le fichier modele qui doit etre relu par srfdis96 (fait chier..)
c      les parametres de dispersion et pre-calculs indispensables a srfdis96
        program srfpre96
c----------------------------------------------------------------------c
c                                                                    c
c      COMPUTER PROGRAMS IN SEISMOLOGY                               c
c      VOLUME IV                                                     c
c                                                                    c
c      PROGRAM: SRFPRE96                                             c
c                                                                    c
c      COPYRIGHT 1986, 1991, 2001                                    c
c      D. R. Russell, R. B. Herrmann                                 c
c      Department of Earth and Atmospheric Sciences                  c
c      Saint Louis University                                        c
c      221 North Grand Boulevard                                     c
c      St. Louis, Missouri 63103                                     c
c      U. S. A.                                                      c
c                                                                    c
c----------------------------------------------------------------------c
c     CHANGES
c     19 JAN 2002 - jobsyn = 1 for observed, 2 for synthetic to 
c             be compatible with   f96subf.f
c     05 MAR 2002 - observed = 'X' synthetic = 'T'
c     12 MAR 2003 - increased NM= in getdsp
c
c
c     This program checks the input control file 'sobs.d' and
c     converts the input data into unformatted binary files
c     to be used in other programs.  The unformatted files
c     are labeled as 'tmpsrfi.xx' where xx is a number from
c     0 to 14.
c
c     Developed by David R. Russell, St. Louis University, Jan. 1984.
c     Restructured Input R. B. Herrmann 8 August 1986
c
c     Restructured August 1986 to simplify input format RBH
c     Modified  September 1986 to include gamma values  RBH
c     Modified  November 2001 to use SURF96 dispersion format and
c     model96 earth model format
c     Modified  January  2002 to internally use the model96 files
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c

        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200,NL2=NL+NL)
c-----
c     LIN   - unit for FORTRAN read from terminal
c     LOT   - unit for FORTRAN write to terminal
c     LER   - unit for FORTRAN error output to terminal
c     NL    - number of layers in model
c     NL2   - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        integer nf10(NL)
        data nf10/NL*1/
        integer nf(13)
        real dd(NL2)
        logical wc(NL2)
        dimension a(NL),b(NL),d(NL),r(NL),rat(NL),qbinv(NL)
        character*80 nmmodl, nmdisp
        logical ext
        common/param/qaqb,itype,dlam,invdep
        data qbinv/NL *0.0/
       
        real*4 h,dcl,dcr
c-----
c     machine dependent initialization
c-----
        call mchdep()
c-----
c     test for the existence of sobs.d
c     if it does not exist interactively set it up
c-----
!        inquire(file='sobs.d',exist=ext)
!        if(.not.(ext))    call setobs(a,b,d,r,rat,qbinv)
c-----
!        open(1,file='sobs.d',status='old',access='sequential')
!        rewind 1
!          read(1,*) h,dcl,onel,dcr,oner
           read(LIN,*) h,dcl,dcr
           oner = 0.0
           onel = 0.0
C         write(LOT,*) h,dcl,onel,dcr,oner
c-----
c     h     = percentage period change for group velocity partial
c     dcl   = Love Wave Phase Velocity Search Increment
c     onel  = Love Wave Increment for Backtracking on Step
c     dcr   = Rayleigh Wave Phase Velocity Search Increment
c     oner  = Rayleigh Wave Increment for Backtracking on Step
c-----
!        read(1,*) (nf(i),i=1,8),nf(11),nf(12)
C       write(LOT,*) (nf(i),i=1,8),nf(11),nf(12)
c        read(LIN,*) (nf(i),i=1,8),nf(11),nf(12)
         nf(1) = 0 
         nf(2) = 0
         nf(5) = 0
         nf(8) = 0
         nf(11) = 0
         nf(12) = 0
c         read(LIN,*) nf(3),nf(4),nf(6),nf(7) 
c-----
c     nf() is the control array
c     nf(1) = 1 estimated stdev computed from residuals
c             0 no scaling by residuals
c     nf(2) = TOTAL number of Love Wave Gamma Modes
c             0 DO NOT PROCESS Love Wave Gamma Data for Q
c     nf(3) = Maximum Number of Love Wave Phase Velocity Modes
c             0 DO NOT PROCESS Love Wave Phase Velocity
c     nf(4) = Maximum Number of Love Wave Group Velocity Modes
c           0 DO NOT PROCESS Love Wave Group Velocity
c     nf(5) = Maximum Number of Rayleigh Wave Gamma Modes
c             0 DO NOT PROCESS Rayleigh Wave Gamma Data for Q
c     nf(6) = Maximum Number of Rayleigh Wave Phase Velocity Modes
c             0 DO NOT PROCESS Rayleigh Wave Phase Velocity
c     nf(7) = Maximum Number of Rayleigh Wave Group Velocity Modes
c             0 DO NOT PROCESS Rayleigh Wave Group Velocity
c     nf(8) = Model Weighting
c             0 No Weighting
c             1 Read In Layer Velocity Weights
c     nf(9) = Number of Layers in Starting Model (from model file)
c     nf(10)= Input Format (from model file)
c             0 - Inversion a, rho fixed
c             1 - Inversion Poisson Ratio Fixed, Rho computed from Vp
c     nf(11)= Type of Smoothing
c             0 None
c             1 Differential
c     nf(12)= Earth Flattening
c             0 No Earth Flattening - Model is Flat Earth Model
c             1 Earth Flattening - Model is Spherical Earth Model
c     nf(13)= Dispersion Abscissa in Input Dispersion Data 
c             (from disp file)
c             0 period
c             1 frequency
c-----
!        read(1,'(a)')nmmodl
!        write(LOT,*)nmmodl
c-----
c     nmmodl= name of file containing model information
c-----
!        call getmdl(nmmodl,nf(9),nf(10),iunit,iflsph)
        call getmdl('STDIN',nf(9),nf(10),iunit,iflsph)
        nf(12) = iflsph
        m=nf(9)
        m2 = m + m
         read(LIN,*) nf(3),nf(4),nf(6),nf(7) 
c-----
c     nmdisp= name of file containing dispersion data
c-----
!        read(1,'(a)')nmdisp
C       write(LOT,*)nmdisp
!        close (1)
c-----
c     initialize model weights
c-----
!        do 2 i=1,m
!            dd(i) = 1.0
!            dd(i+m) = dd(i)
!            wc(i) = .true.
!            wc(i+m) = .true.
!            if(i.eq.m)then
!                wc(i) = .false.
!                wc(i+m) = .false.
!            else
!                wc(i) = .true.
!                wc(i+m) = .true.
!            endif
!    2   continue
!        open(8,file='tmpsrfi.12',form='unformatted',access='sequential')
!        rewind 8
!        write(8) nf(8),m
!        write(8)(dd(i),i=1,m2)
!        write(8)(wc(i),i=1,m2)
!        close(8,status='keep')
!        open(8,file='tmpsrfi.04',form='unformatted',access='sequential')
!        rewind 8
!        write(8) nf(8),m
!        write(8)(dd(i),i=1,m2)
!        write(8)(wc(i),i=1,m2)
c
c     main loop...check input data and write to unformatted files
c
        call getdsp(nmdisp,nf,dcl,dcr,onel,oner,h,iunitd)
              nf34=max(nf(3),nf(4))
              nf67=max(nf(6),nf(7))
c-----
c     close the file with the internal dispersion values
c-----
!       close(8,status='keep')
c-----
c     nfilt   0   no weighting no smoothing
c         1      weighting no smoothing
c         2       no weighting    smoothing
c         3          weighting    smoothing
c----

        if(nf(8).eq.0.and.nf(11).eq.0) nfilt=0
        if(nf(8).eq.1.and.nf(11).eq.0) nfilt=1
        if(nf(8).eq.0.and.nf(11).eq.1) nfilt=2
        if(nf(8).eq.1.and.nf(11).eq.1) nfilt=3
              nf34 = max(nf(3),nf(4))
              nf67 = max(nf(6),nf(7))
c-----
c     this value of nup forces a run of dispersion first
c-----
              nup=2
              nzer=0
              qaqb = 2.25
              invcsl = 0
              wref = 1.0
              lstinv = -1
        invdep = 1
        itype = 0
        dlam = 1.0
        twnmin = -5.0
        twnmax = 20.0
        iter = 0
        nurftn = -1
        pval = 0.5
        sigv = 0.05
        sigr = 0.05
        sigg = 0.00005
        idtwo = 0
        idum2 = 0
        idum3 = 0
        idum4 = 0
        idum5 = 0
        rdum1 = 0.0
        rdum2 = 0.0
        rdum3 = 0.0
        rdum4 = 0.0
        rdum5 = 0.0
c-----
c     iprog is a binary OR: 2= rftn, 1=surf
c-----
!        iprog = 0 + 1
!        call pttmp0(iprog,nzer,nf(1),nf(2),nf34,nf(5),nf67,nf10,nfilt,
!     1            nup,dlam,qaqb,wref,invcsl,lstinv,
!     2      twnmin,twnmax,iter,nurftn,invdep,pval,sigv,sigr,sigg,
!     3      idtwo,idum2,idum3,idum4,idum5,
!     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        end

        subroutine pttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,
     1      nfilt,nup,dlam,qaqb,wref,invcsl,lstinv,
     2      twnmin,twnmax,iter,nurftn,invdep,pval,sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        integer NL
        parameter (NL=200)
        integer nf10(NL)
c-----
c     update control file
c-----
!        open(1,file='tmpsrfi.00',form='unformatted'
!     1            ,access='sequential')
!        rewind 1
!        write(1) iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
!     1      nup,dlam,qaqb,wref,invcsl,lstinv,
!     2      twnmin,twnmax,iter,nurftn,invdep,pval,sigv,sigr,sigg,
!     3      idtwo,idum2,idum3,idum4,idum5,
!     4      rdum1,rdum2,rdum3,rdum4,rdum5
!        close(1,status='keep')
        return
        end

        subroutine psrho(a,b,r,ratio,nswtch,iunit)
c-----
c relations between p and s velocities, poissons ratio and density
c if(nswtch.eq.0)  given a and b find ratio
c              1   given a and ratio find b
c              2   given b and ratio find a
c              3   given a find r
c              4   given b and ratio find a and r
c              5   given a and ratio find b and r
c              6   given r find a
c iunit   0 km,km/sec,gm/cc
c         1 ft,ft/sec,gm/cc
c         2  m, m/sec,gm/cc
c-----
c-----
c     convert units to km/sec if necessary
c-----
        if(iunit.eq.1)then
              a = a * 0.3048006 / 1000.0
              b = b * 0.3048006 / 1000.0
        elseif(iunit.eq.2)then
              a = a / 1000.0
              b = b / 1000.0
        endif
        if(nswtch.eq.0) then
c-----
c determine poissons ratio
c-----
              a2 = a**2
              b2 = b**2
              ratio = (2.*b2-a2)/(2.*b2-2.*a2)
        elseif(nswtch.eq.1 .or. nswtch.eq.5)then
c-----
c determine s-wave velocity
c-----
              b = a*sqrt((1.-2.*ratio)/(2.-2.*ratio))
        elseif(nswtch.eq.2 .or. nswtch.eq.4)then
c-----
c determine p-wave velocity
c-----
              if(ratio.eq..5.or.b.eq.0.0) then
                    a = 1.52
              else
                    a = b/sqrt((1.-2.*ratio)/(2.-2.*ratio))
              endif
        endif
c-----
c     determine densities from P-velocity
c-----
        if(nswtch.eq.3 .or. nswtch.eq.4 .or. nswtch.eq.5)then
              call gtrofa(r,a)
        endif
c-----
c     determin P-velocity from density
c-----
        if(nswtch.eq.6)then
              call gtaofr(r,a)
        endif
c-----
c     convert km/sec units back to original units
c-----
        if(iunit.eq.1)then
              a = a * 1000.0 / 0.3048006
              b = b * 1000.0 / 0.3048006
        elseif(iunit.eq.2)then
              a = a * 1000.0
              b = b * 1000.0
        endif
        return
        end

        subroutine gtrofa(r,a)
c-----
c     obtain density (r,gm/cc) as a function of P-vel (a,km/sec)
c-----
        dimension rp(18)
        data rp/1.65,1.85,2.05,2.15,2.23,2.32,2.39,2.50,2.60,2.70,2.85,
     # 2.98,3.14,3.31,3.49,3.66,3.86,4.05/
c-----
c determine density
c-----
        if(a.gt.1.5 .and. a.lt.10.0)then
              xa = (a-1.5)*2.
              na = xa+1
              dr = rp(na+1)-rp(na)
              da=a-na*.5-1.0
              r = rp(na)+dr*da*2.
        elseif(a.le.1.5)then
              r = 1.0
        elseif(a.ge.10.0)then
              r = 4.10
        endif
        return
c-----
c     determine P-vel (a,km/sec) as a function of density (r,gm/cc)
c-----
        entry gtaofr(r,a)
        if(r.lt.1.65)then
              a = 1.5
        elseif(r.gt.4.05)then
              a = 10.0
        elseif(r.ge.1.65 .and. r.le.4.05)then
              do 10 i=2,18
                    if(rp(i).gt.r) go to 20
10          continue
20          xa = i*.5+.5
              dr = rp(i)-rp(i-1)
              a=((r-rp(i-1))/dr)*.5+xa
        endif
        return
        end

        subroutine getmdl(nmmodl,mmax,nfmod,iunit,iflsph)
c-----
c     igetmod common
c-----
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        character nmmodl*(*)
        common/modtit/title
        character*80 title
c-----
c     open the model file. These values will be saved
c     in internal files so that the original model is not
c     modified in any way
c-----
        call getmod(2,nmmodl,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
c----
c     mmax    nf9
c     
c     nfmod   nf10    1 invert for S - P, Rho calculated from S
c             0 invert for S but P, Rho fixed
c-----
        nfmod = 1
        iunit = 0
        LT = LGSTR(TITLE)
!        call putmod(2,'tmpsrfi.17',mmax,title(1:lt),iunit,iiso,iflsph,
!     1      idimen,icnvel,.false.)
        call putmod(2,'stdout',mmax,title(1:lt),iunit,iiso,iflsph,
     1      idimen,icnvel,.false.)
!        call putmod(2,'tmpmod96.000',mmax,title(1:lt),iunit,iiso,iflsph,
!     1      idimen,icnvel,.false.)
        return
        end

        subroutine getdsp(nmdisp,nf,dcl,dcr,onel,oner,h,iunitd)
c-----
c     nmdisp - file containing dispersion data
c     nf - integer array containing control flags
c-----
        parameter (NM=5000, LOT=6,NP=512)
c-----
c     LIN   - unit for FORTRAN read from terminal
c     LOT   - unit for FORTRAN write to terminal
c     LER   - unit for FORTRAN error output to terminal
c     NL    - number of layers in model
c     NL2   - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c     NP    - number of unique periods
c     NM    - maximum number of observations
c-----
        parameter(LER=0,LIN=5)
        character nmdisp*(*)
        integer nf(13)
        real*4 tper(NM),vel(NM),dvel(NM)
        integer*4 lorr(NM), porg(NM), mode(NM),modemx(2,3)
        real*4 per(NP), tmp(NM)
        integer*4 key(NM), imap(NM), jmap(NM)
        character instr*132
        character ic*1
        data modemx/0,0,0,0,0,0/
c-----
c     MEANING OF VARIABLES
c
c     lorr, porg, mode, per, vel, dvel
c           input parameters for each observation
c           (L/R) (C/U) MODE FREQ VEL SDVEL
c     idat - number of observations
c     per - array of unique periods
c     nper - number of unique periods
c     imap - mapping of observation into unique period
c     modemx(1,1) = max modes love, phase vel, gamma
c     modemx(1,2) = max modes love, group vel
c     modemx(2,1) = max modes rayl, phase vel, gamma
c     modemx(2,2) = max modes rayl, group vel
c
c     to incorporate gamma, we need phase velocity partials
c           for the inversion, so the gamma input will
c           be considered phase for the phase velocity
c           determination
c
c     tmp, key arrays for sort algorithm
c     jmap - temporary mapping
c-----
c-----
c     read in data, store in arrays
c-----
        idat = 0
!       open(1,file=nmdisp,status='old',form='formatted',
!    1            access='sequential')
!       rewind 1
!        open(4,file='tmpsrfi.03',form='unformatted',access='sequential')
!        rewind 4
c-----
c     get units and data type
c     NEW units are always km and period(sec)
c-----
        iunitd = 0
        ifrper = 0
        nf(13)=ifrper
 1000       continue
!           read(1,'(a)',end=1001)instr
           read(LIN,'(a)',end=1001)instr
C        WRITE(0,*)idat,' ',instr
            ls = lgstr(instr)
c-----
c         do the parsing
c-----
            if(instr(1:6).eq.'surf96' .or.
     1          instr(1:6).eq.'SURF96')then
c-----
c             now get wave type
c-----      
                lsep = 6
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'R'. or. ic(1:1).eq.'r')then
                    ilorr = 2
                else if(ic(1:1).eq.'L'. or. ic(1:1).eq.'l')then
                    ilorr = 1
                else if(ic(1:1).eq.'A'. or. ic(1:1).eq.'a')then
                    ilorr = 0
                else
                    go to 1000
                endif
c-----
c             now get observation type
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'C'. or. ic(1:1).eq.'c')then
                    iobs = 1
                else if(ic(1:1).eq.'U'. or. ic(1:1).eq.'u')then
                    iobs = 2
                else if(ic(1:1).eq.'G'. or. ic(1:1).eq.'g')then
                    iobs = 3
                else
                    go to 1000
                endif
c-----
c             now get whether observed or synthetic
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'T'. or. ic(1:1).eq.'t')then
                    iobsyn = 2
                else if(ic(1:1).eq.'F'. or. ic(1:1).eq.'f')then
                    iobsyn = 1
                else if(ic(1:1).eq.'X'. or. ic(1:1).eq.'x')then
                    iobsyn = 1
                else
                    go to 1000
                endif
c-----
c-----
c             now get the values using list directed IO
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                read(instr(lnobl:ls),*) imode,pper,obs,obserr
            endif
            l = ilorr
            m = iobs
            n = imode + 1
            f = pper
            c = obs
            dc = obserr
c-----
c     l     = Love (1) or Rayleigh
c     m     = Phase (1) Group (2) Gamma(3)
c     n     = mode Fund = 1, First = 2, etc
c     f     = Frequency or Period
c     c     = Velocity of Gamma depending on m
c     dc    = Error in Velocity of Gamma, depending on m
c-----
c     before adding to the data set, ensure that the
c     data are to be used
c-----
c     increment n for internal use
c-----
        if(l.eq.1 .and. m.eq.1 .and. n.gt.nf(3))goto 1000
        if(l.eq.1 .and. m.eq.2 .and. n.gt.nf(4))goto 1000
        if(l.eq.1 .and. m.eq.3 .and. n.gt.nf(2))goto 1000
        if(l.eq.2 .and. m.eq.1 .and. n.gt.nf(6))goto 1000
        if(l.eq.2 .and. m.eq.2 .and. n.gt.nf(7))goto 1000
        if(l.eq.2 .and. m.eq.3 .and. n.gt.nf(5))goto 1000
        idat = idat+1
        lorr(idat)=l
        porg(idat)=m
        mode(idat)=n
c-----
c     SURF96 input is always period!!!
c     SURF96 DISPERSION UNITS ARE ALWAYS km/sec and 1/km
c-----
              tper(idat) = f
        vel(idat)=c
        if(dc.eq.0.0)dc=1.0
        dvel(idat)=dc
        key(idat)=idat
        tmp(idat)=tper(idat)
c-----
c     make gamma seem to be phase data
c-----
        mm = m
        if(mm.eq.3)mm=1
        if(n.gt.modemx(l,mm))modemx(l,mm)=n
        goto 1000
 1001       continue
!       close (1)
C        WRITE(0,*)'idat:',idat
        call sort(tmp,key,idat)
C        WRITE(0,*)'idat:',idat
        call uniq(per,tper,key,idat,nper,imap)
c-----
c     now look at Love/Rayleigh
c     followed by Mode
c-----
c     this is perhaps an excessive linear search, but it
c     will work
c-----
c     fix up period count
c-----
        write(LOT,*)nper,nper,nf(12)
!        write(4) nper,nper,nf(12)
c-----
c     adjust nf(3) nf(4) nf(5) nf(6) for internal use
c-----
        nf(3) = modemx(1,1)
        nf(4) = modemx(1,2)
        nf(6) = modemx(2,1)
        nf(7) = modemx(2,2)
c     here order output for several control files
c
c     For observed data, we order as
c           LOVE-RAYL
c                 PHASE(or GAMMA)  - GROUP
c                       MODE
c                             write ilvry,iporg,nmod,per,val,dval
c     For disperion computation order the output as
c           LOVE-RAYL
c                 MODE
c                       range of periods to compute
c                       write range
c-----
        do 2000 ilvry = 1,2
              if(ilvry.eq.1)then
                    nmph = modemx(1,1)
                    nmgr = modemx(1,2)
                    one = onel
                    dc = dcl
              else
                    nmph = modemx(2,1)
                    nmgr = modemx(2,2)
                    one = oner
                    dc = dcr
              endif
c-----
c                 ENFORCE USER MODE LIMITS
c-----
              kmax = nper
              if(nmgr.eq.0)igr = 0
              if(nmph.eq.0 )igr = 1
c           if(nmgr.gt.0 .and. nmph.gt.0 .and. nmgm.gt.0)igr=2
              if(nmgr.gt.0 .and. nmph.gt.0 )igr=2
              nx = max(nmph,nmgr)
              write(LOT,*) kmax,nx,dc,one,igr,h
              write(LOT,*) (per(i),i=1,kmax)
!              write(4) kmax,nx,dc,one,igr,h
!              write(4) (per(i),i=1,kmax)
              do 100  iporg=1,2
                    do 200 nmod=1,modemx(ilvry,iporg)
                          nlr = 0
                          do 300 i=1,idat
                                if(lorr(i).eq.ilvry .and.
     1                                 mod(porg(i),2).eq.mod(iporg,2)
     2                                 .and.mode(i).eq.nmod)then
                                      nlr = nlr +1
                                      tmp(nlr)=tper(i)
                                      key(nlr)=nlr
                                      jmap(nlr)=i
                                endif
  300                         continue
                                if(nlr.gt.0)then
                                      call sort(tmp,key,nlr)
                                      nlrr = nlr
                                do 400 i=1,nlr
                                      j=jmap(key(i))
                                      k = imap(j)
C                                   write(LOT,*)vel(j),
C     1                                       dvel(j),k,
C     2                                       per(k)
C     3                                       ,tper(j),porg(j)
C     4                                          ,ilvry,iporg,
C     5                                          nmod
!                                      write(8)ilvry,porg(j),
!     1                                       nmod,tper(j),
!     2                                       vel(j),dvel(j)
  400                               continue
                                else
                                      key(1) = 1
                                      nlrr = 1
                                      c = 0
                                      sd = 1
                                      cper = 0.0
!                                      write(8)ilvry,iporg,
!     1                                     nmod,cper,
!     2                                     c,sd
                                endif
  200                   continue
  100             continue
c-----
c     for Love or Rayleigh find the period limits
c     for each mode, so that an inclusive comb is constructed
c     for phase velocity search
c-----
        call getlim(modemx,idat,lorr,mode,imap,ilvry)
 2000       continue
!        close(4,status='keep')
        return
        end

        subroutine getlim(modemx,idat,lorr,mode,imap,ilvry)
        integer*4 modemx(2,2), idat,lorr(*),imap(*),mode(*)
c-----
c     get limits on dispersion periods for dispersion program
c     to speed determination of higher modes, we develop
c     an inclusive comb of periods to be evaluated for
c     each mode such that the number of periods at
c     the next higher mode is always within the
c     range of the previous mode
c
c     to do this trick, we work backwords and then output the
c     desired results
c
c-----
        parameter (LIN=5, LOT=6, LER=0)
        integer*4 is(100),ie(100)
        data is/100*0/,ie/100*0/
        md = 0
        do 100 i=1,2
              if(modemx(ilvry,i).gt.md)md=modemx(ilvry,i)
  100 continue
c-----
c     perform linear searches for simplicity
c-----
        do 200 n=md,1,-1
              do 250 j=1,idat
                    if(mode(j).eq.n.and.lorr(j).eq.ilvry)then
                          im = imap(j)
                          if(is(n).eq.0.or.is(n).gt.im)is(n)=im
                          if(ie(n).eq.0.or.ie(n).lt.im)ie(n)=im
                    endif
  250       continue
  200 continue
c-----
c     fill out comb
c-----
        do 300 n=md,2,-1
              if(is(n).lt.is(n-1))is(n-1)=is(n)
              if(is(n-1).eq.0)is(n-1)=is(n)
              if(ie(n).gt.ie(n-1))ie(n-1)=ie(n)
  300 continue
c-----
c     output on unit 4 starting with the first mode
c-----
        do 400 n=1,md
              write(LOT,*)is(n),ie(n)
!              write(4)is(n),ie(n)
  400 continue
        return
        end

        subroutine uniq(y,x,key,nx,ny,imap)
        
c-----
c     this subroutine takes a sorted list, x(key(i))
c     and determines the unique elements y() in it
c     and returns the unique number ny
c     imap(i) = ny maps original into unique
c-----
        real*4 y(*), x(*)
        integer*4 key(*), imap(*)
C        WRITE(0,*)'nx,ny,imap:',nx,ny,(imap(j),j=1,nx)
C        WRITE(0,*)'x:',(x(j),j=1,nx)
C        WRITE(0,*)'key:',(key(j),j=1,nx)
c-----
c     the first element is unique
c-----
        ny = 1
        y(ny) = x(key(1))
        imap(key(1)) = ny
        do 100 i=1,nx
              j = key(i)
              if(y(ny).lt.x(j))then
                    ny = ny + 1
                    y(ny) = x(j)
              endif
              imap(j) = ny
  100       continue
        return
        end

       subroutine sort(x,key,n)
c-----
c     Starting with x(1) ,,, x(n)
c     return   the xarray sorted in increasing order
c     also return the pointers key to the initial array. 
c     For example given x = [ 3, 1, 2 ]
c     the returned values are
c                       x = [ 1, 2, 3 ]        
c                     key = [ 2, 3, 1 ]
c-----
c        Reference: http://en.wikipedia.org/wiki/Bubble_sort
c-----
       integer n
       real x(n)
       integer key(n)

       do i=1,n
           key(i) = i
       enddo
       do i = n, 1, -1
           do j = 1 , i -1
               if(x(j) .gt. x(j+1))then
                   tmp = x(j)
                   x(j) = x(j+1)
                   x(j+1) = tmp
                   ktmp = key(j)
                   key(j) = key(j+1)
                   key(j+1) = ktmp
                endif
           enddo
       enddo
       return
       end

        subroutine setobs(a,b,d,r,rat,qbinv)
c-----
c     set up the control file sobs.d interactively
c-----
        parameter (LIN=5, LOT=6, LER=0)
        parameter (NL=200,NL2=NL+NL)
        character*80 nmmodl, nmdisp
        real a(NL),b(NL),d(NL),r(NL),rat(NL),qbinv(NL)
        integer nf(12)
        logical ext
c-----
c     open the file sobs.d
c-----
        open(1,file='sobs.d',status='unknown',form='formatted',
     1            access='sequential')
        rewind 1
c-----
c     get control for dispersion search
c-----
        write(LOT,*)
     1 ' Enter h,dcl,dcr'
        write(LOT,*)
     1 '      h = fraction change in period to get group vel'
        write(LOT,*)
     1 '            (0.005 is reasonable)'
        write(LOT,*)
     1 '      dcl, dcr are phase velocity increment in root'
        write(LOT,*)
     1 '            search for Love and Rayl respectively'
        read(LIN,*)h,dcl,dcr
        onel = 0.0
        oner = 0.0
        write(1,*)h,dcl,onel,dcr,oner
c-----
c     get second line of sobs.d
c-----
        write(LOT,*)
     1 ' Enter 1 if SW variance based on residuals or'
        write(LOT,*)
     1 '       0 if SW variance based on obs std err'
        read(LIN,*)nf(1)
        write(LOT,*)
     1 ' Enter maximum number of Love gamma modes to process'
        write(LOT,*)
     1 '      0 means DO NO PROCESS LOVE gamma data'
        read(LIN,*)nf(2)
        write(LOT,*)
     1 ' Enter maximum number of Love Phvel modes to process'
        write(LOT,*)
     1 '      0 means DO NO PROCESS LOVE phase vel data'
        read(LIN,*)nf(3)
        write(LOT,*)
     1 ' Enter maximum number of Love Gpvel modes to process'
        write(LOT,*)
     1 '      0 means DO NO PROCESS LOVE group vel data'
        read(LIN,*)nf(4)
        write(LOT,*)
     1 ' Enter maximum number of Rayl gamma modes to process'
        write(LOT,*)
     1 '      0 means DO NO PROCESS RAYL gamma data'
        read(LIN,*)nf(5)
        write(LOT,*)
     1 ' Enter maximum number of Rayl Phvel modes to process'
        write(LOT,*)
     1 '      0 means DO NO PROCESS RAYL phase vel data'
        read(LIN,*)nf(6)
        write(LOT,*)
     1 ' Enter maximum number of Rayl Gpvel modes to process'
        write(LOT,*)
     1 '      0 means DO NO PROCESS RAYL group vel data'
        read(LIN,*)nf(7)
C       write(LOT,*)
C     1 ' Model Layer Weighting'
C       write(LOT,*)
C     1 '      0 No weighting'
C       write(LOT,*)
C     1 '      1 Read in weights for layers'
C       read(LIN,*)nf(8)
        nf(8) = 0
C       write(LOT,*)' Enter inversion technique'
C       write(LOT,*)'    0   invert for Vs :Va,rho fixed'
C       write(LOT,*)'    1 : invert for Vs :Poisson fixed, rho from Vp'
C       read(LIN,*)nf(10)
C       write(LOT,*)
C     1 ' Enter Type of Smoothing'
C       write(LOT,*)
C     1 '      0 None'
C       write(LOT,*)
C     1 '      1 Differential'
C       read(LIN,*)nf(11)
        nf(11) = 1
c-----
c     nf(12) on flat/spherical model is now build into the 
c            earth model file
c-----
        nf(12) = 0
        write(1,'(10i5)')(nf(i),i=1,8),nf(11),nf(12)
c-----
c     get model file information
c-----
        write(LOT,*)
     1 ' Enter name of model file'
        read(LIN,'(a)')nmmodl
        write(1,'(a)')nmmodl
c-----
c     see whether the model file exists
c-----
        inquire(file=nmmodl,exist=ext)
c-----
c     if it does not, interactively create it
c-----
        if(.not.(ext))call setmod(nmmodl,mmax)
c-----
c     get dispersion file
c-----
        write(LOT,*)
     1 ' Enter name of dispersion file'
        read(LIN,'(a)')nmdisp
        write(1,'(a)')nmdisp
c-----
c     if the file does not exist, then interactively build it
c-----
        inquire(file=nmdisp,exist=ext)
        if(.not.(ext))call setdsp(nmdisp)
        return
        end

        subroutine perr(str)
        character str*(*)
        parameter(LER=0,LIN=5,LOT=6)
        dimension nf(10)
        common/param/qaqb,itype,dlam,invdep
        integer NL
        parameter (NL=200)
        integer nf10(NL)
        data nf10/NL*1/
        data nf/10*0/
            nzer=-1
            nf34 = 0
            nfilt = 0
            nup=0
            wref = 1.0
            invcsl = 0
            lstinv = -1
            twnmin = -5.0
            twnmax =  20.0
            iter = 0 
            nurftn = 0 
            invdep = 1
            pval = 0.5
        sigv = 0.05
        sigr = 0.05
        sigg = 0.00005
        idtwo = 0
        idum2 = 0
        idum3 = 0
        idum4 = 0
        idum5 = 0
        rdum1 = 0.0
        rdum2 = 0.0
        rdum3 = 0.0
        rdum4 = 0.0
        rdum5 = 0.0
        write(LER,*)str
c-----
c     iprog is a binary OR: 2= rftn, 1=surf
c-----
        iprog = 0 + 1
        call pttmp0(iprog,nzer,nf(1),nf(2),nf34,nf(5),
     1      nf67,nf10,nfilt,nup,
     2      dlam,qaqb,wref,invcsl,lstinv,
     3      twnmin,twnmax,iter,nurftn,invdep,pval,sigv,sigr,sigg,
     4      idtwo, idum2, idum3, idum4, idum5,
     5      rdum1, rdum2, rdum3, rdum4, rdum5)
        stop
        end

        subroutine getblnk(instr,lsep,ls,lnobl)
c-----
c     determine first non-blank character
c
c     instr   Ch* Character string to be parsed
c     lsep    I*4 index of last non blank character
c     ls  I*4 length of input string
c     lnobl   I*4 index of first non blank character
c-----
        character instr*(*)
        integer lsep,ls,lnobl
        character tab*1
        tab=char(9)
        lnobl = lsep+1
        igotit = 0
        do 1000 i=lsep+1,ls
            if(igotit.eq.0)then
            if(instr(i:i).ne.' ' .and. instr(i:i).ne.tab)then
                lnobl = i
                igotit = 1
            endif
            endif
 1000   continue
        return
        end

