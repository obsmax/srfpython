! #################################################################
      SUBROUTINE GETDSP(NLC,NLU,NRC,NRU)
      IMPLICIT NONE

      REAL c,cper,f,obs,obserr,one,    &
         & Onel,Oner,pper,sd
      INTEGER i,idat,igr,ilorr,ilvry,imode,iobs, &
            & iobsyn,iporg,Iunitd,j,k,LGSTR,STDIN ,&
            & lnobl,STDOUT,ls
      INTEGER lsep,mm,n,nlr,nlrr,NM,nmgr,nmod,nmph ,&
            & NP,nper,nx
      INTEGER NLC,NLU,NRC,NRU

      PARAMETER (NM=1000,STDOUT=6,NP=512)
      PARAMETER (STDIN=5)

      INTEGER(kind=4) lorr(NM),porg(NM),mode(NM),modemx(2,3)
      INTEGER(kind=4) key(NM),imap(NM),jmap(NM)

      REAL(kind=4) tper(NM),vel(NM)
      REAL(kind=4) per(NP),tmp(NM)

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
      oner   = 0.0
      onel   = 0.0
      idat = 0
      Iunitd = 0

 100  READ (STDIN,'(a)',END=200) instr
!        WRITE(0,*)idat,' ',instr
      ls = LGSTR(instr)
!-----
!         do the parsing
!-----
      IF ( instr(1:6) == 'SURF96' .OR. instr(1:6) == 'surf96' ) THEN
!-----
!             now get wave type
!-----
         lsep = 6
         CALL GETBLNK(instr,lsep,ls,lnobl)
         ic = instr(lnobl:lnobl)
         IF ( ic(1:1) == 'R' .OR. ic(1:1) == 'r' ) THEN
            ilorr = 2
         ELSEIF ( ic(1:1) == 'L' .OR. ic(1:1) == 'l' ) THEN
            ilorr = 1
         ELSEIF ( ic(1:1) == 'A' .OR. ic(1:1) == 'a' ) THEN
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
         IF ( ic(1:1) == 'C' .OR. ic(1:1) == 'c' ) THEN
            iobs = 1
         ELSEIF ( ic(1:1) == 'U' .OR. ic(1:1) == 'u' ) THEN
            iobs = 2
         ELSEIF ( ic(1:1) == 'G' .OR. ic(1:1) == 'g' ) THEN
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
         IF ( ic(1:1)=='T' .OR. ic(1:1)=='t' ) THEN
            iobsyn = 2
         ELSEIF ( ic(1:1) == 'F' .OR. ic(1:1) == 'f' ) THEN
            iobsyn = 1
         ELSEIF ( ic(1:1) == 'X' .OR. ic(1:1) == 'x' ) THEN
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
!-----
!     ilorr     = Love (1) or Rayleigh
!     iobs     = Phase (1) Group (2) Gamma(3)
!     n     = mode Fund = 1, First = 2, etc
!     f     = Frequency or Period
!     c     = Velocity or Gamma depending on iobs
!     dc    = Error in Velocity or Gamma, depending on iobs
!-----
!     before adding to the data set, ensure that the
!     data are to be used
!-----
!     increment n for internal use
!-----
      ! ilorr /= 1 <=> ilorr == 2  and inversly
!      write(STDOUT,*) ilorr, iobs, n
      IF (      ( ilorr == 1 .AND. iobs == 1 .AND. n <= NLC )&
         & .OR. ( ilorr == 1 .AND. iobs == 2 .AND. n <= NLU )&
         & .OR. ( ilorr == 2 .AND. iobs == 1 .AND. n <= NRC )&
         & .OR. ( ilorr == 2 .AND. iobs == 2 .AND. n <= NRU )&
              &) THEN
!            write(STDOUT, *) 'ok'
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
            key(idat) = idat
            tmp(idat) = tper(idat)
            !!-----
            !!     make gamma seem to be phase data
            !!-----
            IF ( n > modemx(ilorr,iobs) ) modemx(ilorr,iobs) = n
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
         IF ( ilvry == 1 ) THEN
            nmph = modemx(1,1) ! NLC
            nmgr = modemx(1,2) ! NLU
            one = Onel
         ELSE
            nmph = modemx(2,1) ! NRC
            nmgr = modemx(2,2) ! NRU
            one = Oner
         ENDIF
!-----
!                 ENFORCE USER MODE LIMITS
!-----

         IF ( nmgr == 0 ) igr = 0
         IF ( nmph == 0 ) igr = 1
!           if(nmgr > 0 .and. nmph > 0 .and. nmgm > 0)igr=2
         IF ( nmgr > 0 .AND. nmph > 0 ) igr = 2
         nx = MAX(nmph,nmgr)

         WRITE (STDOUT,"(I4,I4,I4)") nper,nx,igr!,H

         do i = 1, nper
             WRITE (STDOUT,"(F8.3)", advance="no") per(i)
         end do
         WRITE (STDOUT,"(A1)") ""

         DO iporg = 1,2
            DO nmod = 1,modemx(ilvry,iporg)
               nlr = 0
               DO i = 1,idat
                  IF ( lorr(i) == ilvry .AND. MOD(porg(i),2)            &
                     &  == MOD(iporg,2) .AND. mode(i) == nmod ) THEN
                     nlr = nlr + 1
                     tmp(nlr) = tper(i)
                     key(nlr) = nlr
                     jmap(nlr) = i
                  ENDIF
               ENDDO
               IF ( nlr > 0 ) THEN
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


!#################################################################
      SUBROUTINE GETLIM(Modemx,Idat,Lorr,Mode,Imap,Ilvry)
      IMPLICIT NONE
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

      INTEGER i,Ilvry,im,j,STDOUT,md,n
      PARAMETER (STDOUT=6)

      INTEGER(kind=4) Modemx(2,2),Idat,Lorr(*),Imap(*),Mode(*)
      INTEGER(kind=4) is(100),ie(100)
      DATA is/100*0/,ie/100*0/

      md = 0
      DO i = 1,2
         IF ( Modemx(Ilvry,i) > md ) md = Modemx(Ilvry,i)
      ENDDO
!-----
!     perform linear searches for simplicity
!-----
      DO n = md,1,-1
         DO j = 1,Idat
            IF ( Mode(j) == n .AND. Lorr(j) == Ilvry ) THEN
               im = Imap(j)
               IF ( is(n) == 0 .OR. is(n) > im ) is(n) = im
               IF ( ie(n) == 0 .OR. ie(n) < im ) ie(n) = im
            ENDIF
         ENDDO
      ENDDO
!-----
!     fill out comb
!-----
      DO n = md,2,-1
         IF ( is(n) < is(n-1) ) is(n-1) = is(n)
         IF ( is(n-1) == 0 ) is(n-1) = is(n)
         IF ( ie(n) > ie(n-1) ) ie(n-1) = ie(n)
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
!-----
!     this subroutine takes a sorted list, x(key(i))
!     and determines the unique elements y() in it
!     and returns the unique number ny
!     imap(i) = ny maps original into unique
!-----
      REAL(kind=4) Y(*),X(*)
      INTEGER(kind=4) Key(*),Imap(*),i,j,Nx,Ny

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
         IF ( Y(Ny) < X(j) ) THEN
            Ny = Ny + 1
            Y(Ny) = X(j)
         ENDIF
         Imap(j) = Ny
      ENDDO
      END


!#################################################################
      SUBROUTINE SORT(X,Key,N)
      IMPLICIT NONE

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

      INTEGER N,i,j,ktmp,Key(N)
      REAL tmp, X(N)

      DO i = 1,N
         Key(i) = i
      ENDDO

      DO i = N,1,-1
         DO j = 1,i - 1
            IF ( X(j) > X(j+1) ) THEN
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


!#################################################################
      SUBROUTINE GETBLNK(Instr,Lsep,Ls,Lnobl)
      IMPLICIT NONE
!-----
!     determine first non-blank character
!
!     instr   Ch* Character string to be parsed
!     lsep    I*4 index of last non blank character
!     ls  I*4 length of input string
!     lnobl   I*4 index of first non blank character
!-----

      CHARACTER Instr*(*),tab*1
      INTEGER Lsep,Ls,Lnobl,i,igotit

      tab = CHAR(9)
      Lnobl = Lsep + 1
      igotit = 0
      DO i = Lsep + 1,Ls
         IF ( igotit == 0 ) THEN
            IF ( Instr(i:i) /= ' ' .AND. Instr(i:i) /= tab ) THEN
               Lnobl = i
               igotit = 1
            ENDIF
         ENDIF
      ENDDO
      END
