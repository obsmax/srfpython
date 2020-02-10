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

!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
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

      PROGRAM SRFPRE96
      IMPLICIT NONE

      INTEGER STDIN,STDOUT
      INTEGER NLC,NLU,NRC,NRU
      PARAMETER (STDIN=5,STDOUT=6)

      READ (STDIN,*) NLC,NLU,NRC,NRU
      CALL GETDSP(NLC,NLU,NRC,NRU)
      END


