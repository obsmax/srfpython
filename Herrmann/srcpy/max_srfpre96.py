import sys, numpy as np

# -----------------------------------------------------------------
def unpackmod96(string):
    """read files at mod96 format (see Herrmann's doc)
    """
    string = [line.strip() for line in string.split('\n')]
    string.remove('')
    # assert string[0].strip() == "MODEL.01"
    title = string[1].strip()
    isotropic = string[2].strip().upper() == "ISOTROPIC"
    kgs = string[3].strip().upper() == "KGS"
    flatearth = string[4].strip().upper() == "FLAT EARTH"
    # oned = string[5].strip().upper() == "1-D"
    # cstvelo = string[6].strip().upper() == "CONSTANT VELOCITY"
    # assert string[7].strip() == "LINE08"
    # assert string[8].strip() == "LINE09"
    # assert string[9].strip() == "LINE10"
    # assert string[10].strip() == "LINE11"
    # header = string[11]

    nlayer = len(string) - 12
    #H, VP, VS, RHO, QP, QS, ETAP, ETAS, FREFP, FREFS = [np.empty(nlayer, float) for _ in xrange(10)]
    DAT = np.empty((nlayer, 10), float)

    for n in xrange(nlayer):
        DAT[n, :] = np.asarray(string[12 + n].split(), float)

    H, VP, VS, RHO, QP, QS, ETAP, ETAS, FREFP, FREFS = [DAT[:, j] for j in xrange(10)]

    assert not H[0]
    assert H[:-1].all()
    Z = np.concatenate(([0.], H[:-1].cumsum()))
    return title, isotropic, kgs, flatearth, nlayer, Z, H, VP, VS, RHO, QP, QS, ETAP, ETAS, FREFP, FREFS
# -----------------------------------------------------------------
def readmod96(filename):
    with open(filename, 'r') as fid:
        L = fid.readlines()
    return unpackmod96("".join(L))
# -----------------------------------------------------------------
def unpacksurf96(string):
    string = [line.strip() for line in string.split('\n')]
    string.remove('')
    npoints = len(string)


    datatypes = ['|S1', '|S1', '|S1', int, float, float, float]
    WAVE, TYPE, FLAG, MODE, PERIOD, VALUE, DVALUE = [np.empty(npoints, dtype = d) for d in datatypes]
    NLC, NLU, NRC, NRU = 0, 0, 0, 0
    for n in xrange(npoints):
        l = string[n].split()
        WAVE[n], TYPE[n], FLAG[n] = np.asarray(l[1:4], "|S1")
        MODE[n] = int(l[4])
        PERIOD[n], VALUE[n], DVALUE[n] = np.asarray(l[5:], float)
        if   WAVE[n] == "L":
            if   TYPE[n] == "C": NLC += 1
            elif TYPE[n] == "U": NLU += 1
        elif WAVE[n] == "R":
            if   TYPE[n] == "C": NRC += 1
            elif TYPE[n] == "U": NRU += 1
        else: raise
    return WAVE, TYPE, FLAG, MODE, PERIOD, VALUE, DVALUE, NLC, NLU, NRC, NRU
# -----------------------------------------------------------------
def readsurf96(filename):
    with open(filename, 'r') as fid:
        L = fid.readlines()
    return unpacksurf96("".join(L))
# -----------------------------------------------------------------
if __name__ == "__main__":
    title, isotropic, kgs, flatearth, nlayer, Z, H, VP, VS, RHO, QP, QS, ETAP, ETAS, FREFP, FREFS = readmod96("bidon0.mod96")
    WAVE, TYPE, FLAG, MODE, PERIOD, VALUE, DVALUE, NLC, NLU, NRC, NRU = readsurf96("bidon0.surf96")





# #!#################################################################
#       SUBROUTINE GETDSP(Nmdisp,Nf,Dcl,Dcr,Onel,Oner,H,Iunitd)
#       IMPLICIT NONE
# #!*--GETDSP270
# #!*** Start of declarations inserted by SPAG
#       REAL c , cper , dc , Dcl , Dcr , f , H , obs , obserr , one ,     &
#          & Onel , Oner , pper , sd
#       INTEGER i , idat , ifrper , igr , ilorr , ilvry , imode , iobs ,  &
#             & iobsyn , iporg , Iunitd , j , k , kmax , l , LGSTR , LIN ,&
#             & lnobl , LOT , ls
#       INTEGER lsep , m , mm , n , nlr , nlrr , NM , nmgr , nmod , nmph ,&
#             & NP , nper , nx
# #!*** End of declarations inserted by SPAG
# #!-----
# #!     nmdisp - file containing dispersion data
# #!     nf - integer array containing control flags
# #!-----
#       PARAMETER (NM=5000,LOT=6,NP=512)
# #!-----
# #!     LIN   - unit for FORTRAN read from terminal
# #!     LOT   - unit for FORTRAN write to terminal
# #!     LER   - unit for FORTRAN error output to terminal
# #!     NL    - number of layers in model
# #!     NL2   - number of columns in model (first NL2/2 are
# #!           - velocity parameters, second NL2/2 are Q values)
# #!     NP    - number of unique periods
# #!     NM    - maximum number of observations
# #!-----
#       PARAMETER (LIN=5)
#       CHARACTER Nmdisp*(*)
#       INTEGER Nf(13)
#       REAL*4 tper(NM) , vel(NM) , dvel(NM)
#       INTEGER*4 lorr(NM) , porg(NM) , mode(NM) , modemx(2,3)
#       REAL*4 per(NP) , tmp(NM)
#       INTEGER*4 key(NM) , imap(NM) , jmap(NM)
#       CHARACTER instr*132
#       CHARACTER ic*1
#       DATA modemx/0 , 0 , 0 , 0 , 0 , 0/
# #!-----
# #!     MEANING OF VARIABLES
# #!
# #!     lorr, porg, mode, per, vel, dvel
# #!           input parameters for each observation
# #!           (L/R) (C/U) MODE FREQ VEL SDVEL
# #!     idat - number of observations
# #!     per - array of unique periods
# #!     nper - number of unique periods
# #!     imap - mapping of observation into unique period
# #!     modemx(1,1) = max modes love, phase vel, gamma
# #!     modemx(1,2) = max modes love, group vel
# #!     modemx(2,1) = max modes rayl, phase vel, gamma
# #!     modemx(2,2) = max modes rayl, group vel
# #!
# #!     to incorporate gamma, we need phase velocity partials
# #!           for the inversion, so the gamma input will
# #!           be considered phase for the phase velocity
# #!           determination
# #!
# #!     tmp, key arrays for sort algorithm
# #!     jmap - temporary mapping
# #!-----
# #!-----
# #!     read in data, store in arrays
# #!-----
#       idat = 0
# #!       open(1,file=nmdisp,status='old',form='formatted',
# #!    1            access='sequential')
# #!       rewind 1
# #!        open(4,file='tmpsrfi.03',form='unformatted',access='sequential')
# #!        rewind 4
# #!-----
# #!     get units and data type
# #!     NEW units are always km and period(sec)
# #!-----
#       Iunitd = 0
#       ifrper = 0
#       Nf(13) = ifrper
# #!           read(1,'(a)',end=1001)instr
#  100  READ (LIN,'(a)',END=200) instr
# #!        WRITE(0,*)idat,' ',instr
#       ls = LGSTR(instr)
# #!-----
# #!         do the parsing
# #!-----
#       IF ( instr(1:6).EQ.'SURF96' .OR. instr(1:6).EQ.'surf96' ) THEN
# #!-----
# #!             now get wave type
# #!-----
#          lsep = 6
#          CALL GETBLNK(instr,lsep,ls,lnobl)
#          ic = instr(lnobl:lnobl)
#          IF ( ic(1:1).EQ.'R' .OR. ic(1:1).EQ.'r' ) THEN
#             ilorr = 2
#          ELSEIF ( ic(1:1).EQ.'L' .OR. ic(1:1).EQ.'l' ) THEN
#             ilorr = 1
#          ELSEIF ( ic(1:1).EQ.'A' .OR. ic(1:1).EQ.'a' ) THEN
#             ilorr = 0
#          ELSE
#             GOTO 100
#          ENDIF
# #!-----
# #!             now get observation type
# #!-----
#          lsep = lnobl
#          CALL GETBLNK(instr,lsep,ls,lnobl)
#          ic = instr(lnobl:lnobl)
#          IF ( ic(1:1).EQ.'C' .OR. ic(1:1).EQ.'c' ) THEN
#             iobs = 1
#          ELSEIF ( ic(1:1).EQ.'U' .OR. ic(1:1).EQ.'u' ) THEN
#             iobs = 2
#          ELSEIF ( ic(1:1).EQ.'G' .OR. ic(1:1).EQ.'g' ) THEN
#             iobs = 3
#          ELSE
#             GOTO 100
#          ENDIF
# #!-----
# #!             now get whether observed or synthetic
# #!-----
#          lsep = lnobl
#          CALL GETBLNK(instr,lsep,ls,lnobl)
#          ic = instr(lnobl:lnobl)
#          IF ( ic(1:1).EQ.'T' .OR. ic(1:1).EQ.'t' ) THEN
#             iobsyn = 2
#          ELSEIF ( ic(1:1).EQ.'F' .OR. ic(1:1).EQ.'f' ) THEN
#             iobsyn = 1
#          ELSEIF ( ic(1:1).EQ.'X' .OR. ic(1:1).EQ.'x' ) THEN
#             iobsyn = 1
#          ELSE
#             GOTO 100
#          ENDIF
# #!-----
# #!-----
# #!             now get the values using list directed IO
# #!-----
#          lsep = lnobl
#          CALL GETBLNK(instr,lsep,ls,lnobl)
#          READ (instr(lnobl:ls),*) imode , pper , obs , obserr
#       ENDIF
#       l = ilorr
#       m = iobs
#       n = imode + 1
#       f = pper
#       c = obs
#       dc = obserr
# #!-----
# #!     l     = Love (1) or Rayleigh
# #!     m     = Phase (1) Group (2) Gamma(3)
# #!     n     = mode Fund = 1, First = 2, etc
# #!     f     = Frequency or Period
# #!     c     = Velocity of Gamma depending on m
# #!     dc    = Error in Velocity of Gamma, depending on m
# #!-----
# #!     before adding to the data set, ensure that the
# #!     data are to be used
# #!-----
# #!     increment n for internal use
# #!-----
#       IF ( l.NE.1 .OR. m.NE.1 .OR. n.LE.Nf(3) ) THEN
#          IF ( l.NE.1 .OR. m.NE.2 .OR. n.LE.Nf(4) ) THEN
#             IF ( l.NE.1 .OR. m.NE.3 .OR. n.LE.Nf(2) ) THEN
#                IF ( l.NE.2 .OR. m.NE.1 .OR. n.LE.Nf(6) ) THEN
#                   IF ( l.NE.2 .OR. m.NE.2 .OR. n.LE.Nf(7) ) THEN
#                      IF ( l.NE.2 .OR. m.NE.3 .OR. n.LE.Nf(5) ) THEN
#                         idat = idat + 1
#                         lorr(idat) = l
#                         porg(idat) = m
#                         mode(idat) = n
# #!-----
# #!     SURF96 input is always period!!!
# #!     SURF96 DISPERSION UNITS ARE ALWAYS km/sec and 1/km
# #!-----
#                         tper(idat) = f
#                         vel(idat) = c
#                         IF ( dc.EQ.0.0 ) dc = 1.0
#                         dvel(idat) = dc
#                         key(idat) = idat
#                         tmp(idat) = tper(idat)
# #!-----
# #!     make gamma seem to be phase data
# #!-----
#                         mm = m
#                         IF ( mm.EQ.3 ) mm = 1
#                         IF ( n.GT.modemx(l,mm) ) modemx(l,mm) = n
#                      ENDIF
#                   ENDIF
#                ENDIF
#             ENDIF
#          ENDIF
#       ENDIF
#       GOTO 100
# #!       close (1)
# #!        WRITE(0,*)'idat:',idat
#  200  CALL SORT(tmp,key,idat)
# #!        WRITE(0,*)'idat:',idat
#       CALL UNIQ(per,tper,key,idat,nper,imap)
# #!-----
# #!     now look at Love/Rayleigh
# #!     followed by Mode
# #!-----
# #!     this is perhaps an excessive linear search, but it
# #!     will work
# #!-----
# #!     fix up period count
# #!-----
#       WRITE (LOT,*) nper , nper , Nf(12)
# #!        write(4) nper,nper,nf(12)
# #!-----
# #!     adjust nf(3) nf(4) nf(5) nf(6) for internal use
# #!-----
#       Nf(3) = modemx(1,1)
#       Nf(4) = modemx(1,2)
#       Nf(6) = modemx(2,1)
#       Nf(7) = modemx(2,2)
# #!     here order output for several control files
# #!
# #!     For observed data, we order as
# #!           LOVE-RAYL
# #!                 PHASE(or GAMMA)  - GROUP
# #!                       MODE
# #!                             write ilvry,iporg,nmod,per,val,dval
# #!     For disperion computation order the output as
# #!           LOVE-RAYL
# #!                 MODE
# #!                       range of periods to compute
# #!                       write range
# #!-----
#       DO ilvry = 1 , 2
#          IF ( ilvry.EQ.1 ) THEN
#             nmph = modemx(1,1)
#             nmgr = modemx(1,2)
#             one = Onel
#             dc = Dcl
#          ELSE
#             nmph = modemx(2,1)
#             nmgr = modemx(2,2)
#             one = Oner
#             dc = Dcr
#          ENDIF
# #!-----
# #!                 ENFORCE USER MODE LIMITS
# #!-----
#          kmax = nper
#          IF ( nmgr.EQ.0 ) igr = 0
#          IF ( nmph.EQ.0 ) igr = 1
# #!           if(nmgr.gt.0 .and. nmph.gt.0 .and. nmgm.gt.0)igr=2
#          IF ( nmgr.GT.0 .AND. nmph.GT.0 ) igr = 2
#          nx = MAX(nmph,nmgr)
#          WRITE (LOT,*) kmax , nx , dc , one , igr , H
#          WRITE (LOT,*) (per(i),i=1,kmax)
# #!              write(4) kmax,nx,dc,one,igr,h
# #!              write(4) (per(i),i=1,kmax)
#          DO iporg = 1 , 2
#             DO nmod = 1 , modemx(ilvry,iporg)
#                nlr = 0
#                DO i = 1 , idat
#                   IF ( lorr(i).EQ.ilvry .AND. MOD(porg(i),2)            &
#                      & .EQ.MOD(iporg,2) .AND. mode(i).EQ.nmod ) THEN
#                      nlr = nlr + 1
#                      tmp(nlr) = tper(i)
#                      key(nlr) = nlr
#                      jmap(nlr) = i
#                   ENDIF
#                ENDDO
#                IF ( nlr.GT.0 ) THEN
#                   CALL SORT(tmp,key,nlr)
#                   nlrr = nlr
#                   DO i = 1 , nlr
#                      j = jmap(key(i))
#                      k = imap(j)
# #!                                   write(LOT,*)vel(j),
# #!     1                                       dvel(j),k,
# #!     2                                       per(k)
# #!     3                                       ,tper(j),porg(j)
# #!     4                                          ,ilvry,iporg,
# #!     5                                          nmod
# #!                                      write(8)ilvry,porg(j),
# #!     1                                       nmod,tper(j),
# #!     2                                       vel(j),dvel(j)
#                   ENDDO
#                ELSE
#                   key(1) = 1
#                   nlrr = 1
#                   c = 0
#                   sd = 1
#                   cper = 0.0
# #!                                      write(8)ilvry,iporg,
# #!     1                                     nmod,cper,
# #!     2                                     c,sd
#                ENDIF
#             ENDDO
#          ENDDO
# #!-----
# #!     for Love or Rayleigh find the period limits
# #!     for each mode, so that an inclusive comb is constructed
# #!     for phase velocity search
# #!-----
#          CALL GETLIM(modemx,idat,lorr,mode,imap,ilvry)
#       ENDDO
# #!        close(4,status='keep')
#       END
# #!*==GETLIM.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017
#
# #!#################################################################
#       SUBROUTINE GETLIM(Modemx,Idat,Lorr,Mode,Imap,Ilvry)
#       IMPLICIT NONE
# #!*--GETLIM572
# #!*** Start of declarations inserted by SPAG
#       INTEGER i , Ilvry , im , j , LOT , md , n
# #!*** End of declarations inserted by SPAG
#       INTEGER*4 Modemx(2,2) , Idat , Lorr(*) , Imap(*) , Mode(*)
# #!-----
# #!     get limits on dispersion periods for dispersion program
# #!     to speed determination of higher modes, we develop
# #!     an inclusive comb of periods to be evaluated for
# #!     each mode such that the number of periods at
# #!     the next higher mode is always within the
# #!     range of the previous mode
# #!
# #!     to do this trick, we work backwords and then output the
# #!     desired results
# #!
# #!-----
#       PARAMETER (LOT=6)
#       INTEGER*4 is(100) , ie(100)
#       DATA is/100*0/ , ie/100*0/
#       md = 0
#       DO i = 1 , 2
#          IF ( Modemx(Ilvry,i).GT.md ) md = Modemx(Ilvry,i)
#       ENDDO
# #!-----
# #!     perform linear searches for simplicity
# #!-----
#       DO n = md , 1 , -1
#          DO j = 1 , Idat
#             IF ( Mode(j).EQ.n .AND. Lorr(j).EQ.Ilvry ) THEN
#                im = Imap(j)
#                IF ( is(n).EQ.0 .OR. is(n).GT.im ) is(n) = im
#                IF ( ie(n).EQ.0 .OR. ie(n).LT.im ) ie(n) = im
#             ENDIFNLC, NLU, NRC, NRU
#          ENDDO
#       ENDDO
# #!-----
# #!     fill out comb
# #!-----
#       DO n = md , 2 , -1
#          IF ( is(n).LT.is(n-1) ) is(n-1) = is(n)
#          IF ( is(n-1).EQ.0 ) is(n-1) = is(n)
#          IF ( ie(n).GT.ie(n-1) ) ie(n-1) = ie(n)
#       ENDDO
# #!-----
# #!     output on unit 4 starting with the first mode
# #!-----
#       DO n = 1 , md
#          WRITE (LOT,*) is(n) , ie(n)
# #!              write(4)is(n),ie(n)
#       ENDDO
#       END
# #!*==UNIQ.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017
#
# #!#################################################################
#       SUBROUTINE UNIQ(Y,X,Key,Nx,Ny,Imap)
#       IMPLICIT NONE
# #!*--UNIQ629
# #!*** Start of declarations inserted by SPAG
#       INTEGER i , j , Nx , Ny
# #!*** End of declarations inserted by SPAG
#
# #!-----
# #!     this subroutine takes a sorted list, x(key(i))
# #!     and determines the unique elements y() in it
# #!     and returns the unique number ny
# #!     imap(i) = ny maps original into unique
# #!-----
#       REAL*4 Y(*) , X(*)
#       INTEGER*4 Key(*) , Imap(*)
# #!        WRITE(0,*)'nx,ny,imap:',nx,ny,(imap(j),j=1,nx)
# #!        WRITE(0,*)'x:',(x(j),j=1,nx)
# #!        WRITE(0,*)'key:',(key(j),j=1,nx)
# #!-----
# #!     the first element is unique
# #!-----
#       Ny = 1
#       Y(Ny) = X(Key(1))
#       Imap(Key(1)) = Ny
#       DO i = 1 , Nx
#          j = Key(i)
#          IF ( Y(Ny).LT.X(j) ) THEN
#             Ny = Ny + 1
#             Y(Ny) = X(j)
#          ENDIF
#          Imap(j) = Ny
#       ENDDO
#       END
# #!*==SORT.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017
#
# #!#################################################################
#       SUBROUTINE SORT(X,Key,N)
#       IMPLICIT NONE
# #!*--SORT665
# #!*** Start of declarations inserted by SPAG
#       INTEGER i , j , ktmp
#       REAL tmp
# #!*** End of declarations inserted by SPAG
# #!-----
# #!     Starting with x(1) ,,, x(n)
# #!     return   the xarray sorted in increasing order
# #!     also return the pointers key to the initial array.
# #!     For example given x = [ 3, 1, 2 ]
# #!     the returned values are
# #!                       x = [ 1, 2, 3 ]
# #!                     key = [ 2, 3, 1 ]
# #!-----
# #!        Reference: http://en.wikipedia.org/wiki/Bubble_sort
# #!-----
#       INTEGER N
#       REAL X(N)
#       INTEGER Key(N)
#       DO i = 1 , N
#          Key(i) = i
#       ENDDO
#       DO i = N , 1 , -1
#          DO j = 1 , i - 1
#             IF ( X(j).GT.X(j+1) ) THEN
#                tmp = X(j)
#                X(j) = X(j+1)
#                X(j+1) = tmp
#                ktmp = Key(j)
#                Key(j) = Key(j+1)
#                Key(j+1) = ktmp
#             ENDIF
#          ENDDO
#       ENDDO
#
#
#       END
# #!*==GETBLNK.spg  processed by SPAG 6.72Dc at 15:03 on  5 Dec 2017
# #!#################################################################
#       SUBROUTINE GETBLNK(Instr,Lsep,Ls,Lnobl)
#       IMPLICIT NONE
# #!*--GETBLNK706
# #!*** Start of declarations inserted by SPAG
#       INTEGER i , igotit
# #!*** End of declarations inserted by SPAG
# #!-----
# #!     determine first non-blank character
# #!
# #!     instr   Ch* Character string to be parsed
# #!     lsep    I*4 index of last non blank character
# #!     ls  I*4 length of input string
# #!     lnobl   I*4 index of first non blank character
# #!-----
#       CHARACTER Instr*(*)
#       INTEGER Lsep , Ls , Lnobl
#       CHARACTER tab*1
#       tab = CHAR(9)
#       Lnobl = Lsep + 1
#       igotit = 0
#       DO i = Lsep + 1 , Ls
#          IF ( igotit.EQ.0 ) THEN
#             IF ( Instr(i:i).NE.' ' .AND. Instr(i:i).NE.tab ) THEN
#                Lnobl = i
#                igotit = 1
#             ENDIF
#          ENDIF
#       ENDDO
#       END
#
