!*==SRFDIS96.spg  processed by SPAG 6.72Dc at 15:04 on  5 Dec 2017
! mes modifs
! le programme prend ses entrees dans stdin, elles correspondent aux sorties de
! la version modifiee de srfpre96 vers stdout
! sorties vers stdout
! tableau:
!    premiere colonne  : 1 pour LOVE, 2 pour RAYLEIGH, 0 pour un mode manquant (ignorer le reste de la ligne dans ce cas)
!    seconde  colonne  : numero de mode, 0=fondamental,...
!    troisieme et 4eme colonnes : periode en s
!        si la 4eme colonne vaut 0, alors la 3eme colonne est la periode absolue
!        sinon la 3eme colonne est diminuee d un facteur 1-h et la 4eme est augmentee d'un fact. 1+h
!    5eme et 6eme : vitesse de phase correspondant aux colonnes 3 et 4
!----------------------------------------------------------------------c
!                                                                    c
!    COMPUTER PROGRAMS IN SEISMOLOGY                                 c
!    VOLUME IV                                                       c
!                                                                    c
!    PROGRAM: SRFDIS                                                 c
!                                                                    c
!    COPYRIGHT 1986, 1991                                            c
!    D. R. Russell, R. B. Herrmann                                   c
!    Department of Earth and Atmospheric Sciences                    c
!    Saint Louis University                                          c
!    221 North Grand Boulevard                                       c
!    St. Louis, Missouri 63103                                       c
!    U. S. A.                                                        c
!                                                                    c
!----------------------------------------------------------------------c
!
      PROGRAM SRFDIS96
      IMPLICIT NONE
!      real :: start, finish, start1, finish1 ! for cpu timer
!
!   This is a combination of program 'surface80' which search the poles
!   on C-T domain, and the program 'surface81' which search in the F-K
!   domain.  The input data is slightly different with its precessors.
!   -Wang   06/06/83.
!
!   The program calculates the dispersion values for any
!   layered model, any frequency, and any mode.
!
!   This program will accept one liquid layer at the surface.
!   In such case ellipticity of rayleigh wave is that at the
!   top of solid array.  Love wave communications ignore
!   liquid layer.
!
!   Program developed by Robert B Herrmann Saint Louis
!   univ. Nov 1971, and revised by C. Y. Wang on Oct 1981.
!   Modified for use in surface wave inversion, and
!   addition of spherical earth flattening transformation, by
!   David R. Russell, St. Louis University, Jan. 1984.
!
!     Changes
!     28 JAN 2003 - fixed minor but for sphericity correction by
!         saving one parameter in subroutine sphere
!     20 JUL 2004 - removed extraneous line at line 550
!         since dc not defined
!         if(dabs(c1-c2)  <=  dmin1(1.d-6*c1,0.005d+0*dc) )go to 1000
!     28 DEC 2007 - changed the Earth flattening to now use layer
!         midpoint and the Biswas (1972: PAGEOPH 96, 61-74, 1972)
!         density mapping for P-SV  - note a true comparison
!         requires the ability to handle a fluid core for SH and SV
!         Also permit one layer with fluid is base of the velocity is 0.001 km/sec
!-----
!
!     NL  - layers in model
!     NP  - number of unique periods

!-----
!     STDIN  - unit for FORTRAN read from terminal
!     STDOUT  - unit for FORTRAN write to terminal
!     i_fisrst_solid_layer - index of the first solid layer, starting with 1 from top
!     idispl, idispr = number of dispersion points to compute for love and rayleigh
!     nlayer = max number of layers in the mode, used in common with other subroutines
!     ifunc = 1 for love 2 for rayleig
!     vsmin, vsmax = extramal values for S-wave velocity in the model in km/s
!                    in case of a water layer, use P-wave velocity instead of S for vsmin
!     jmn = index of layer where vsmin is encountered, starting from 1, from top to botto
!     jsol = a flag to tell if vsmin corresponds to a liquid layer
!     t1a,t1b = periods at which to compute phase velocity,
!               for group velocity dispersions, these two values are used for differentiation
!               for phase velocity, only t1a is used
      !-----
      INTEGER  STDIN, STDOUT 
      INTEGER i_first_solid_layer,idispl,idispr,nlayer,ifunc,jmn,jsol
      INTEGER i,ie,ierr,ifirst,ift, &
            & igr,iq,iret,is,itst,k,kmax
      INTEGER mode
      INTEGER NP, NL, NLAY

      REAL vsmin,vsmax,t1a,t1b
      REAL cc0,cc1,ddc,h,sone

      PARAMETER (STDIN =5,STDOUT =6)
      PARAMETER (NL=100,NLAY=100)
      PARAMETER (NP=512)

      DOUBLE PRECISION twopi,one,onea
      DOUBLE PRECISION cc,c1,clow,cm,dc,t1
      DOUBLE PRECISION c(NP),cb(NP)
      REAL(kind=4) t(NP), thicknesses(NL),vp_values(NL),vs_values(NL),density_values(NL)
      REAL(kind=4) qbinv(NL),qainv(NL)
      COMMON /MODL  / thicknesses, vp_values, vs_values, density_values
      COMMON /PARA  / nlayer,i_first_solid_layer,twopi

      INTEGER iunit,iiso,idimen,icnvel
      CHARACTER(LEN=80) :: FMT  ! WARNING must be long enough for FMT
!-----
      FMT = "(I2,I2,F11.6,F11.6,F11.6,F11.6)" ! print format for output
      READ(STDIN, *) nlayer
      READ(STDIN, *) thicknesses(1:nlayer-1)
      READ(STDIN, *) vp_values(1:nlayer)
      READ(STDIN, *) vs_values(1:nlayer)
      READ(STDIN, *) density_values(1:nlayer)
      READ(STDIN, *) h
      READ(STDIN, *) idispl
      ddc = 0.005  ! for now

      idispr=idispl ! always true, see srfpre96
      iunit = 0
      iiso = 0
      idimen = 0
      icnvel = 0
      ierr = 0
!-----
!     check for water layer
!-----
      i_first_solid_layer = 1 !first solid layer index
      IF ( vs_values(1) <= 0.0 ) i_first_solid_layer = 2
      twopi = 2.D0*3.141592653589793D0
      one = 1.0D-2
      jmn = 1
      vsmax = -1.E20
      vsmin = 1.E20
!-----
!     find the extremal velocities to assist in starting search
!-----
      DO i = 1, nlayer
         IF ( vs_values(i) >  0.01 .AND. vs_values(i) <  vsmin ) THEN
            vsmin = vs_values(i)
            jmn = i
            jsol = 1
         ELSEIF ( vs_values(i) <= 0.01 .AND. vp_values(i) <  vsmin ) THEN
            vsmin = vp_values(i)
            jmn = i
            jsol = 0
         ENDIF
         IF ( vs_values(i) >  vsmax ) vsmax = vs_values(i)
      ENDDO

      sone = 0.
      DO ifunc = 1, 2 ! 1 = love 2 = rayleigh
          IF   ((ifunc == 1 .AND. idispl >  0) &
            .OR.(ifunc == 2 .AND. idispr >  0)) THEN

               READ (STDIN  ,*) kmax, mode, igr !, h
               READ (STDIN  ,*) (t(i),i=1,kmax)

               IF ( sone <  0.01 ) sone = 2.0
               onea = DBLE(sone)
               !-----
               !     get starting value for phase velocity,
               !         which will correspond to the
               !     VP/VS ratio
               !-----
               IF ( jsol == 0 ) THEN ! water layer
                  cc1 = vsmin
               ELSE ! solid layer solve halfspace period equation
                  CALL GTSOLH(vp_values(jmn),vs_values(jmn),cc1)
               ENDIF
               !-----
               !     back off a bit to get a starting value at a lower phase velocity
               !-----
               cc1 = .95*cc1
               cc1 = .90*cc1
               cc = DBLE(cc1)
               dc = DBLE(ddc)
               dc = DABS(dc)
               c1 = cc
               cm = cc
               DO i = 1, kmax
                  cb(i) = 0.0D0
                  c(i) = 0.0D0
               ENDDO
               ift = 999

               DO iq = 1, mode
!                  call cpu_time(start1)

                  READ (STDIN  ,*) is, ie
                  itst = ifunc
                  DO k = is, ie
                     IF ( k >= ift ) GOTO 5
                     t1 = DBLE(t(k)) ! period at which to compute phase velo
                     IF ( igr >  0 ) THEN
                        t1a = t1/(1.+h) ! for group, shift the period a little bit for differential
                        t1b = t1/(1.-h) ! for group, shift the period a little bit for differential
                        t1 = DBLE(t1a)
                     ELSE
                        t1a = SNGL(t1)
                     ENDIF
                    !-----
                    !     get initial phase velocity estimate to begin search
                    !
                    !     in the notation here, c() is an array of phase velocities
                    !     c(k-1) is the velocity estimate of the present mode
                    !     at the k-1 period, while c(k) is the phase velocity of the
                    !     previous mode at the k period. Since there must be no mode
                    !     crossing, we make use of these values. The only complexity
                    !     is that the dispersion may be reversed.
                    !
                    !     The subroutine getsol determines the zero crossing and refines
                    !     the root.
                    !-----
                     IF ( k == is .AND. iq == 1 ) THEN
                        c1 = cc
                        clow = cc
                        ifirst = 1
                     ELSEIF ( k == is .AND. iq >  1 ) THEN
                        c1 = c(is) + one*dc
                        clow = c1
                        ifirst = 1
                     ELSEIF ( k >  is .AND. iq >  1 ) THEN
                        ifirst = 0
                    ! clow = c(k) + one*dc
                    ! c1 = c(k-1) -onea*dc
                        clow = c(k) + one*dc
                        c1 = c(k-1)
                        IF ( c1 <  clow ) c1 = clow
                     ELSEIF ( k >  is .AND. iq == 1 ) THEN
                        ifirst = 0
                        c1 = c(k-1) - onea*dc
                        clow = cm
                     ENDIF
                     !-----
                     !     bracket root and refine it
                     !-----
!                     !WRITE(*,*) "VERBOSE : av1",t1,c1,clow,dc,cm,vsmax,iret,ifunc,ifirst
                     CALL GETSOL(t1,c1,clow,dc,cm,vsmax,iret,ifunc,ifirst)
!                     !WRITE(*,*) "VERBOSE : ap1",t1,c1,clow,dc,cm,vsmax,iret,ifunc,ifirst
                     IF ( iret == -1 ) GOTO 5 ! if root not found
                     c(k) = c1
                     !-----
                     !     for group velocities compute near above solution
                     !-----
                     IF ( igr >  0 ) THEN
                        ! for group, need to compute phase velo at periods t1a and t1b instead of just t1
                        t1 = DBLE(t1b)
                        ifirst = 0
                        clow = cb(k) + one*dc
                        c1 = c1 - onea*dc
!                        write(*,*) "av2",t1,c1,clow,dc,cm,vsmax,iret,ifunc,ifirst
                        CALL GETSOL(t1,c1,clow,dc,cm,vsmax,iret,ifunc,ifirst)
!                        write(*,*) "ap2",t1,c1,clow,dc,cm,vsmax,iret,ifunc,ifirst
                        !-----
                        !     test if root not found at slightly larger period
                        !-----
                        IF ( iret == -1 ) c1 = c(k)
                        cb(k) = c1
                     ELSE
                        c1 = 0.0D+00
                     ENDIF
                     cc0 = SNGL(c(k))
                     cc1 = SNGL(c1)

                      WRITE(*,FMT) itst,iq - 1,t1a,t1b,cc0,cc1

                  ENDDO
                  GOTO 10

 5                ift = k
                  itst = 0
 10            ENDDO
          ENDIF
      ENDDO
      END


! ########################################################
      SUBROUTINE GTSOLH(A,B,C)
      IMPLICIT NONE

      REAL A, B, C, fac1, fac2, fr, frp, gamma
      INTEGER i
      REAL(kind=4) kappa, k2, gk2
      !WRITE(*,*) "VERBOSE : ENTER GTSOLH" !VERBOSE!

!-----
!     starting solution
!-----
      C = 0.95*B
      DO i = 1, 5
         gamma = B/A
         kappa = C/B
         k2 = kappa**2
         gk2 = (gamma*kappa)**2
         fac1 = SQRT(1.0-gk2)
         fac2 = SQRT(1.0-k2)
         fr = (2.0-k2)**2 - 4.0*fac1*fac2
         frp = -4.0*(2.0-k2)*kappa + 4.0*fac2*gamma*gamma*kappa/fac1 +  &
             & 4.0*fac1*kappa/fac2
         frp = frp/B
         C = C - fr/frp
      ENDDO
      END


! ########################################################
      SUBROUTINE GETSOL(T1,C1,Clow,Dc,Cm,vsmax,Iret,Ifunc,Ifirst)
      IMPLICIT NONE
      REAL vsmax
      INTEGER Ifirst, Ifunc, Iret
      REAL(kind=8) wvno, omega, twopi
      REAL(kind=8) C1, c2, cn, Cm, Dc, T1, Clow, fdir
      REAL(kind=8) DLTAR, del1, del2, del1st, plmn
      SAVE del1st

      !WRITE(*,*) "VERBOSE : ENTER GETSOL" !VERBOSE!

      ! presumed inputs = T1, Clow, DC, CM, vsmax, Ifunc, Ifirst
      ! presumed outputs = C1, Iret

!-----
!     subroutine to bracket dispersion curve
!     and then refine it
!-----
!     t1  - period
!     c1  - initial guess on low side of mode
!     clow - lowest possible value for present mode in a
!           reversed direction search
!     dc  - phase velocity search increment
!     cm  - minimum possible solution
!     vsmax   - maximum shear velocity
!     iret    - 1 = successful
!         - -1= unsuccessful
!     ifunc   - 1 - Love
!         - 2 - Rayleigh
!     ifirst  - 1 this is first period for a particular mode
!         - 0 this is not the first period
!             (this is to define period equation sign
!              for mode jumping test)
!-----
!     to avoid problems in mode jumping with reversed dispersion
!     we note what the polarity of period equation is for phase
!     velocities just beneath the zero crossing at the
!         first period computed.
!-----
!     bracket solution
!-----

      twopi = 2.D0*3.141592653589793D0
      omega = twopi/T1
      wvno = omega/C1
      del1 = DLTAR(wvno,omega,Ifunc)
      IF ( Ifirst == 1 ) del1st = del1
      plmn = DSIGN(1.0D+00,del1st)*DSIGN(1.0D+00,del1)
      IF ( Ifirst == 1 ) THEN
         fdir = +1.0
      ELSEIF ( Ifirst /= 1 .AND. plmn >= 0.0D+00 ) THEN
         fdir = +1.0
      ELSEIF ( Ifirst /= 1 .AND. plmn <  0.0D+00 ) THEN
         fdir = -1.0
      ENDIF
!-----
!     idir indicates the direction of the search for the
!     true phase velocity from the initial estimate.
!     Usually phase velocity increases with period and
!     we always underestimate, so phase velocity should increase
!     (idir = +1). For reversed dispersion, we should look
!     downward from the present estimate. However, we never
!     go below the floor of clow, when the direction is reversed
!-----
!     original syntax
! 100  IF ( idir >  0 ) THEN
!         c2 = C1 + Dc
!      ELSE
!         c2 = C1 - Dc
!      ENDIF
!      IF ( c2 <= Clow ) THEN
!         idir = +1
!         C1 = Clow
!      ENDIF
!      IF ( c2 <= Clow ) GOTO 100

!     syntax 1
100   C2 = C1 + fdir * Dc
      IF ( c2 <= Clow ) THEN
         ! wrong start
         ! we are below Clow, force research direction to upward
         ! put C1 to Clow and restart the loop from beginning
         fdir = +1.0
         C1 = Clow
         GOTO 100
      ENDIF

      omega = twopi/T1
      wvno = omega/c2
      del2 = DLTAR(wvno,omega,Ifunc)

      IF ( DSIGN(1.0D+00,del1) /= DSIGN(1.0D+00,del2) ) THEN
        !-----
        !     root bracketed, refine it
        !-----
         CALL NEVILL(T1,C1,c2,del1,del2,Ifunc,cn)
         C1 = cn
         IF ( C1 >  (vsmax) ) THEN
            Iret = -1  ! failure signal
!            GOTO 99999 ! quit subroutine
         ELSE
            Iret = 1  !            GOTO 99999
!            RETURN
         ENDIF
         GOTO 99999
      ENDIF
      C1 = c2
      del1 = del2
!   check that c1 is in region of solutions
      IF ( C1 <  Cm ) THEN
         Iret = -1 ! failure signal
      ELSE
         IF ( C1 <  (vsmax+Dc) ) GOTO 100
         Iret = -1 ! failure signal
      ENDIF


99999 END


! ########################################################
      SUBROUTINE NEVILL(T,C1,C2,Del1,Del2,Ifunc,Cc)
!-----
!   hybrid method for refining root once it has been bracketted
!   between c1 and c2.  interval halving is used where other schemes
!   would be inefficient.  once suitable region is found neville s
!   iteration method is used to find root.
!   the procedure alternates between the interval halving and neville
!   techniques using whichever is most efficient
!-----
!     the control integer nev means the following:
!
!     nev = 0 force interval halving
!     nev = 1 permit neville iteration if conditions are proper
!     nev = 2 neville iteration is being used
!-----

      IMPLICIT NONE
!*--NEVILL525
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION C1, C2, c3, Cc, Del1, Del2, del3, denom, &
                     & DLTAR, omega, s1, s13, s2, s32, ss1, ss2 ,&
                     & T, twopi, wvno, x
      DOUBLE PRECISION y
      INTEGER Ifunc, j, kk, i_first_solid_layer, m, nlayer, nctrl, nev, NL


      PARAMETER (NL=100)
      REAL(kind=4) thicknesses(NL), vp_values(NL), &
              &vs_values(NL), density_values(NL)
      DIMENSION x(20), y(20)
      COMMON /MODL  / thicknesses, vp_values, vs_values, density_values
      COMMON /PARA  / nlayer, i_first_solid_layer, twopi
!-----
!     initial guess
!-----
      omega = twopi/T
      CALL HALF(C1,C2,c3,del3,omega,Ifunc)
      nev = 1
      nctrl = 1
 100  nctrl = nctrl + 1
      IF ( nctrl >= 100 ) THEN
         Cc = c3
      ELSE
!-----
!     make sure new estimate is inside the previous values. If not
!     perform interval halving
!-----
         IF ( c3 <  DMIN1(C1,C2) .OR. c3 >  DMAX1(C1,C2) ) THEN
            nev = 0
            CALL HALF(C1,C2,c3,del3,omega,Ifunc)
         ENDIF
         s13 = Del1 - del3
         s32 = del3 - Del2
!-----
!     define new bounds according to the sign of the period equation
!-----
         IF ( DSIGN(1.D+00,del3)*DSIGN(1.D+00,Del1) <  0.0D+00 ) THEN
            C2 = c3
            Del2 = del3
         ELSE
            C1 = c3
            Del1 = del3
         ENDIF
!-----
!     check for convergence. vp_values relative error criteria is used
!-----
         IF ( DABS(C1-C2) <= 1.D-6*C1 ) THEN
            Cc = c3
         ELSE
!-----
!     if the slopes are not the same between c1, c3 and c3
!     do not use neville iteration
!-----
            IF ( DSIGN(1.0D+00,s13) /= DSIGN(1.0D+00,s32) ) nev = 0
!-----
!     if the period equation differs by more than a factor of 10
!     use interval halving to avoid poor behavior of polynomial fit
!-----
            ss1 = DABS(Del1)
            s1 = 0.01*ss1
            ss2 = DABS(Del2)
            s2 = 0.01*ss2
            IF ( s1 >  ss2 .OR. s2 >  ss1 .OR. nev == 0 ) THEN
               CALL HALF(C1,C2,c3,del3,omega,Ifunc)
               nev = 1
               m = 1
            ELSE
               IF ( nev == 2 ) THEN
                  x(m+1) = c3
                  y(m+1) = del3
               ELSE
                  x(1) = C1
                  y(1) = Del1
                  x(2) = C2
                  y(2) = Del2
                  m = 1
               ENDIF
!-----
!     perform Neville iteration. Note instead of generating y(x)
!     we interchange the x and y of formula to solve for x(y) when
!     y = 0
!-----
               DO kk = 1, m
                  j = m - kk + 1
                  denom = y(m+1) - y(j)
                  IF ( DABS(denom) <  1.0D-10*ABS(y(m+1)) ) GOTO 110
                  x(j) = (-y(j)*x(j+1)+y(m+1)*x(j))/denom
               ENDDO
               c3 = x(1)
               wvno = omega/c3
               del3 = DLTAR(wvno,omega,Ifunc)
               nev = 2
               m = m + 1
               IF ( m >  10 ) m = 10
               GOTO 100
 110           CALL HALF(C1,C2,c3,del3,omega,Ifunc)
               nev = 1
               m = 1
            ENDIF
            GOTO 100
         ENDIF
      ENDIF
      END


! ########################################################
      SUBROUTINE HALF(C1,C2,C3,Del3,Omega,Ifunc)
      IMPLICIT NONE

      DOUBLE PRECISION C1, C2, C3, Del3, DLTAR, Omega, wvno
      INTEGER Ifunc

      C3 = 0.5*(C1+C2)
      wvno = Omega/C3
      Del3 = DLTAR(wvno,Omega,Ifunc)
      END


! ########################################################
      FUNCTION DLTAR(Wvno,Omega,Ifunc)
!   control the way to P-SV or SH.
!   -> call DLTAR1 for LOVE waves
!      and DLTAR4 for RAYLEIGH waves
      IMPLICIT NONE

      DOUBLE PRECISION DLTAR, DLTAR1, DLTAR4, Omega, Wvno
      INTEGER Ifunc
      !WRITE(*,*) "VERBOSE : ENTER DLTAR" !VERBOSE!

      IF ( Ifunc == 1 ) THEN
         !love wave period equation
         DLTAR = DLTAR1(Wvno,Omega)
      ELSEIF ( Ifunc == 2 ) THEN
         !rayleigh wave period equation
         DLTAR = DLTAR4(Wvno,Omega)
      ENDIF
      !WRITE(*,*) "VERBOSE : EXIT DLTAR", DLTAR !VERBOSE!
      END


! ########################################################
      FUNCTION DLTAR1(Wvno,Omega)
!
!   find SH dispersion values.
!
      IMPLICIT NONE
      DOUBLE PRECISION beta1, cosq, DLTAR1, e1, e10, e2, e20,    &
                     & fac, Omega, q, rb, rho1, sinq, twopi,     &
                     & Wvno, wvnom, wvnop, xkb, xmu, xnor
      DOUBLE PRECISION y, ynor, z
      INTEGER i_first_solid_layer, m, nlayer, mmm1, NL
      PARAMETER (NL=100)
      REAL(kind=4) thicknesses(NL), vp_values(NL), vs_values(NL), density_values(NL) !, RTP(NL), DTP(NL),  BTP(NL)
      COMMON /MODL  / thicknesses, vp_values, vs_values, density_values !, RTP, DTP, BTP
      COMMON /PARA  / nlayer, i_first_solid_layer, twopi

      !WRITE(*,*) "VERBOSE : ENTER DLTAR1" !VERBOSE!

!
!   Haskell-Thompson love wave formulation from halfspace
!   to surface.
!
      beta1 = DBLE(vs_values(nlayer))
      rho1 = DBLE(density_values(nlayer))
      xkb = Omega/beta1
      wvnop = Wvno + xkb
      wvnom = DABS(Wvno-xkb)
      rb = DSQRT(wvnop*wvnom)
      e1 = rho1*rb
      e2 = 1.D+00/(beta1*beta1)
      mmm1 = nlayer - 1
      DO m = mmm1, i_first_solid_layer, -1
         beta1 = DBLE(vs_values(m))
         rho1 = DBLE(density_values(m))
         xmu = rho1*beta1*beta1
         xkb = Omega/beta1
         wvnop = Wvno + xkb
         wvnom = DABS(Wvno-xkb)
         rb = DSQRT(wvnop*wvnom)
         q = DBLE(thicknesses(m))*rb
         IF ( Wvno <  xkb ) THEN
            sinq = DSIN(q)
            y = sinq/rb
            z = -rb*sinq
            cosq = DCOS(q)
         ELSEIF ( Wvno == xkb ) THEN
            cosq = 1.0D+00
            y = DBLE(thicknesses(m))
            z = 0.0D+00
         ELSE
            fac = 0.0D+00
            IF ( q <  16 ) fac = DEXP(-2.0D+0*q)
            cosq = (1.0D+00+fac)*0.5D+00
            sinq = (1.0D+00-fac)*0.5D+00
            y = sinq/rb
            z = rb*sinq
         ENDIF
         e10 = e1*cosq + e2*xmu*z
         e20 = e1*y/xmu + e2*cosq
         xnor = DABS(e10)
         ynor = DABS(e20)
         IF ( ynor >  xnor ) xnor = ynor
         IF ( xnor <  1.D-40 ) xnor = 1.0D+00
         e1 = e10/xnor
         e2 = e20/xnor
      ENDDO
      DLTAR1 = e1
      END


! ########################################################!
      FUNCTION DLTAR4(Wvno,Omga)
!   find P-SV dispersion values.
      IMPLICIT NONE
      DOUBLE PRECISION A0, beta, ca, cosp, CPCq, CPY, CPZ, CQW, &
                     & CQX, cr, DLTAR4, dpth, e, ee, exa, gam,  &
                     & gamm1, gammk, omega, Omga
      DOUBLE PRECISION p, q, ra, rb, rho1, t, twopi, w, w0,    &
                     & Wvno, wvno2, wvnom, wvnop, WY, WZ, xka,   &
                     & xkb, XY, XZ, znul
      INTEGER i, j, i_first_solid_layer, m, nlayer, mmm1, NL
      PARAMETER (NL=100)
      DIMENSION e(5), ee(5), ca(5,5)
      REAL(kind=4) thicknesses(NL), vp_values(NL), &
              &vs_values(NL), density_values(NL)
      COMMON /MODL  / thicknesses, vp_values, &
              &vs_values, density_values
      COMMON /PARA  / nlayer, i_first_solid_layer, twopi
      COMMON /OVRFLW/ A0, CPCq, CPY, CPZ, CQW, CQX, XY, XZ, WY, WZ
!
      !WRITE(*,*) "VERBOSE : ENTER DLTAR4" !VERBOSE!
      omega = Omga
      IF ( omega <  1.0D-4 ) omega = 1.0D-4
      wvno2 = Wvno*Wvno
      xka = omega/DBLE(vp_values(nlayer))
      xkb = omega/DBLE(vs_values(nlayer))
      wvnop = Wvno + xka
      wvnom = DABS(Wvno-xka)
      ra = DSQRT(wvnop*wvnom)
      wvnop = Wvno + xkb
      wvnom = DABS(Wvno-xkb)
      rb = DSQRT(wvnop*wvnom)
      t = DBLE(vs_values(nlayer))/omega
!-----
!   E matrix for the bottom half-space.
!-----
      gammk = 2.D+00*t*t
      gam = gammk*wvno2
      gamm1 = gam - 1.D+00
      rho1 = DBLE(density_values(nlayer))
      e(1) = rho1*rho1*(gamm1*gamm1-gam*gammk*ra*rb)
      e(2) = -rho1*ra
      e(3) = rho1*(gamm1-gammk*ra*rb)
      e(4) = rho1*rb
      e(5) = wvno2 - ra*rb
!-----
!   matrix multiplication from bottom layer upward
!-----
      mmm1 = nlayer - 1
      DO m = mmm1, i_first_solid_layer, -1
         xka = omega/DBLE(vp_values(m))
         xkb = omega/DBLE(vs_values(m))
         t = DBLE(vs_values(m))/omega
         gammk = 2.D+00*t*t
         gam = gammk*wvno2
         wvnop = Wvno + xka
         wvnom = DABS(Wvno-xka)
         ra = DSQRT(wvnop*wvnom)
         wvnop = Wvno + xkb
         wvnom = DABS(Wvno-xkb)
         rb = DSQRT(wvnop*wvnom)
         dpth = DBLE(thicknesses(m))
         rho1 = DBLE(density_values(m))
         p = ra*dpth
         q = rb*dpth
         beta = DBLE(vs_values(m))
!-----
!   evaluate cosP, cosQ,.... in var.
!   evaluate Dunkin's matrix in dnka.
!-----
         CALL VAR(p,q,ra,rb,Wvno,xka,xkb,dpth,w,cosp,exa)
         CALL DNKA(ca,wvno2,gam,gammk,rho1)
         DO i = 1, 5
            cr = 0.0D+00
            DO j = 1, 5
               cr = cr + e(j)*ca(j,i)
            ENDDO
            ee(i) = cr
         ENDDO
         CALL NORMC(ee,exa)
         DO i = 1, 5
            e(i) = ee(i)
         ENDDO
      ENDDO
      IF ( i_first_solid_layer /= 1 ) THEN
!-----
!   include water layer.
!-----
         xka = omega/DBLE(vp_values(1))
         wvnop = Wvno + xka
         wvnom = DABS(Wvno-xka)
         ra = DSQRT(wvnop*wvnom)
         dpth = DBLE(thicknesses(1))
         rho1 = DBLE(density_values(1))
         p = ra*dpth
         beta = DBLE(vs_values(1))
         znul = 1.0D-05
         CALL VAR(p,znul,ra,znul,Wvno,xka,znul,dpth,w,cosp,exa)
         w0 = -rho1*w
         DLTAR4 = cosp*e(1) + w0*e(2)
      ELSE
         DLTAR4 = e(1)
      ENDIF
      END


! ########################################################
      SUBROUTINE VAR(P,Q,Ra,Rb,Wvno,Xka,Xkb,Dpth,W,Cosp,Exa)
!-----
!   find variables cosP, cosQ, sinP, sinQ, etc.
!   as well as cross products required for compound matrix
!-----
!   To handle the hyperbolic functions correctly for large
!   arguments, we use an extended precision procedure,
!   keeping in mind that the maximum precision in double
!   precision is on the order of 16 decimal places.
!
!   So  cosp = 0.5 ( exp(+p) + exp(-p))
!            = exp(p) * 0.5 * ( 1.0 + exp(-2p) )
!   becomes
!       cosp = 0.5 * (1.0 + exp(-2p) ) with an exponent p
!   In performing matrix multiplication, we multiply the modified
!   cosp terms and add the exponents. At the last step
!   when it is necessary to obtain a true amplitude,
!   we then form exp(p). For normalized amplitudes at any depth,
!   we carry an exponent for the numerator and the denominator, and
!   scale the resulting ratio by exp(NUMexp - DENexp)
!
!   The propagator matrices have three basic terms
!
!   HSKA        cosp  cosq
!   DUNKIN      cosp*cosq     1.0
!
!   When the extended floating point is used, we use the
!   largest exponent for each, which is  the following:
!
!   Let pex = p exponent > 0 for evanescent waves = 0 otherwise
!   Let sex = s exponent > 0 for evanescent waves = 0 otherwise
!   Let exa = pex + sex
!
!   Then the modified matrix elements are as follow:
!
!   Haskell:  cosp -> 0.5 ( 1 + exp(-2p) ) exponent = pex
!             cosq -> 0.5 ( 1 + exp(-2q) ) * exp(q-p)
!                                          exponent = pex
!          (this is because we are normalizing all elements in the
!           Haskell matrix )
!    Compound:
!            cosp * cosq -> normalized cosp * cosq exponent = pex + qex
!             1.0  ->    exp(-exa)
!-----
      IMPLICIT NONE
!*--VAR893
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION A0, Cosp, cosq, CPCq, CPY, CPZ, CQW, CQX ,&
                     & Dpth, Exa, fac, P, pex, Q, qmp, Ra, Rb, &
                     & sex, sinp, sinq
      DOUBLE PRECISION W, Wvno, WY, WZ, x, Xka, Xkb, XY, XZ,   &
                     & y, z
!*** End of declarations inserted by SPAG
      COMMON /OVRFLW/ A0, CPCq, CPY, CPZ, CQW, CQX, XY, XZ, WY ,&
                    & WZ

      !WRITE(*,*) "VERBOSE : ENTER VAR" !VERBOSE!

      Exa = 0.0D+00
      A0 = 1.0D+00
!-----
!   examine P-wave eigenfunctions
!      checking whether c> vp c=vp or c < vp
!-----
      pex = 0.0D+00
      sex = 0.0D+00
      IF ( Wvno <  Xka ) THEN
         sinp = DSIN(P)
         W = sinp/Ra
         x = -Ra*sinp
         Cosp = DCOS(P)
      ELSEIF ( Wvno == Xka ) THEN
         Cosp = 1.0D+00
         W = Dpth
         x = 0.0D+00
      ELSEIF ( Wvno >  Xka ) THEN
         pex = P
         fac = 0.0D+00
         IF ( P <  16 ) fac = DEXP(-2.0D+00*P)
         Cosp = (1.0D+00+fac)*0.5D+00
         sinp = (1.0D+00-fac)*0.5D+00
         W = sinp/Ra
         x = Ra*sinp
      ENDIF
!-----
!   examine S-wave eigenfunctions
!      checking whether c > vs, c = vs, c < vs
!-----
      IF ( Wvno <  Xkb ) THEN
         sinq = DSIN(Q)
         y = sinq/Rb
         z = -Rb*sinq
         cosq = DCOS(Q)
      ELSEIF ( Wvno == Xkb ) THEN
         cosq = 1.0D+00
         y = Dpth
         z = 0.0D+00
      ELSEIF ( Wvno >  Xkb ) THEN
         sex = Q
         fac = 0.0D+00
         IF ( Q <  16 ) fac = DEXP(-2.0D+0*Q)
         cosq = (1.0D+00+fac)*0.5D+00
         sinq = (1.0D+00-fac)*0.5D+00
         y = sinq/Rb
         z = Rb*sinq
      ENDIF
!-----
!   form eigenfunction products for use with compound matrices
!-----
      Exa = pex + sex
      A0 = 0.0D+00
      IF ( Exa <  60.0D+00 ) A0 = DEXP(-Exa)
      CPCq = Cosp*cosq
      CPY = Cosp*y
      CPZ = Cosp*z
      CQW = cosq*W
      CQX = cosq*x
      XY = x*y
      XZ = x*z
      WY = W*y
      WZ = W*z
      qmp = sex - pex
      fac = 0.0D+00
      IF ( qmp >  -40.0D+00 ) fac = DEXP(qmp)
      cosq = cosq*fac
      y = fac*y
      z = fac*z
      END


! ########################################################
      SUBROUTINE NORMC(Ee,Ex)
!   This routine is an important step to control over- or
!   underflow.
!   The Haskell or Dunkin vectors are normalized before
!   the layer matrix stacking.
!   Note that some precision will be lost during normalization.
!
      IMPLICIT NONE
      DOUBLE PRECISION Ee, Ex, t1, t2
      INTEGER i
      DIMENSION Ee(5)

      !WRITE(*,*) "VERBOSE : ENTER NORMC"

      Ex = 0.0D+00
      t1 = 0.0D+00
      DO i = 1, 5
         IF ( DABS(Ee(i)) >  t1 ) t1 = DABS(Ee(i))
      ENDDO
      IF ( t1 <  1.D-40 ) t1 = 1.D+00
      DO i = 1, 5
         t2 = Ee(i)
         t2 = t2/t1
         Ee(i) = t2
      ENDDO
!-----
!   store the normalization factor in exponential form.
!-----
      Ex = DLOG(t1)
      END


! ########################################################
      SUBROUTINE DNKA(Ca,Wvno2,Gam,Gammk,Rho)
!     Dunkin's matrix.
      IMPLICIT NONE
      DOUBLE PRECISION A0, a0pq, Ca, CPCq, CPY, CPZ, CQW, CQX,  &
                     & Gam, gamm1, Gammk, gm1sq, gmgm1, gmgmk,    &
                     & one, Rho, rho2, t, twgm1, two
      DOUBLE PRECISION Wvno2, WY, WZ, XY, XZ
      DIMENSION Ca(5,5)
      COMMON /OVRFLW/ A0, CPCq, CPY, CPZ, CQW, CQX, XY, XZ, WY ,&
                    & WZ
      DATA one, two/1.D+00, 2.D+00/

      !WRITE(*,*) "VERBOSE : ENTER DNKA"
      gamm1 = Gam - one
      twgm1 = Gam + gamm1
      gmgmk = Gam*Gammk
      gmgm1 = Gam*gamm1
      gm1sq = gamm1*gamm1
      rho2 = Rho*Rho
      a0pq = A0 - CPCq
      Ca(1,1) = CPCq - two*gmgm1*a0pq - gmgmk*XZ - Wvno2*gm1sq*WY
      Ca(1,2) = (Wvno2*CPY-CQX)/Rho
      Ca(1,3) = -(twgm1*a0pq+Gammk*XZ+Wvno2*gamm1*WY)/Rho
      Ca(1,4) = (CPZ-Wvno2*CQW)/Rho
      Ca(1,5) = -(two*Wvno2*a0pq+XZ+Wvno2*Wvno2*WY)/rho2
      Ca(2,1) = (gmgmk*CPZ-gm1sq*CQW)*Rho
      Ca(2,2) = CPCq
      Ca(2,3) = Gammk*CPZ - gamm1*CQW
      Ca(2,4) = -WZ
      Ca(2,5) = Ca(1,4)
      Ca(4,1) = (gm1sq*CPY-gmgmk*CQX)*Rho
      Ca(4,2) = -XY
      Ca(4,3) = gamm1*CPY - Gammk*CQX
      Ca(4,4) = Ca(2,2)
      Ca(4,5) = Ca(1,2)
      Ca(5,1) = -(two*gmgmk*gm1sq*a0pq+gmgmk*gmgmk*XZ+gm1sq*gm1sq*WY)   &
              & *rho2
      Ca(5,2) = Ca(4,1)
      Ca(5,3) = -(Gammk*gamm1*twgm1*a0pq+Gam*Gammk*Gammk*XZ+gamm1*gm1sq*&
              & WY)*Rho
      Ca(5,4) = Ca(2,1)
      Ca(5,5) = Ca(1,1)
      t = -two*Wvno2
      Ca(3,1) = t*Ca(5,3)
      Ca(3,2) = t*Ca(4,3)
      Ca(3,3) = A0 + two*(CPCq-Ca(1,1))
      Ca(3,4) = t*Ca(2,3)
      Ca(3,5) = t*Ca(1,3)
      END

