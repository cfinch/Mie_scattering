
      SUBROUTINE MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,
     &                  NMOM, IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC,
     &                  PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK,
     &                  SPIKE )

c       NoPMOM version:  saves storage by deleting everything
c       having to do with PMOM calculation:

c       * deletes modules LPCoef, LPCO1T, LPCO2T

c       * NMOM,IPOLZN,PMOM,MOMDIM,CalcMo,NPQuan deleted from argument
c         lists and everywhere else they occur internally, although the
c         first four are retained in the MIEV0 argument list for
c         compatibility with complete version of MIEV0

c       * MIEV0:   shrinks size of LitA,LitB arrays to 2;  deletes call
c                  to LPCoef;  saving of An,Bn in LitA,LitB arrays
c                  deleted

c       * MiPrnt:  PMOM print deleted

c       * CkInMi:  NMOM test changed to allow only zero value (nonzero
c                  value would have no impact but indicates user
c                  conceptual error or negligence)

c       * TestMi:  saving and restoring of NMOM,IPOLZN deleted;
c                  30-loop deleted

c    Computes Mie scattering and extinction efficiencies; asymmetry
c    factor;  forward- and backscatter amplitude;  scattering
c    amplitudes vs. scattering angle for incident polarization parallel
c    and perpendicular to the plane of scattering;
c    some quantities needed in polarized radiative transfer;  and
c    information about whether or not a resonance has been hit.

c    Input and output arguments are described in file MIEV.doc;
c    internal variables are described in MIEV0.f file.
c    Many statements are accompanied by comments referring to 
c    references in MIEV.doc, notably the NCAR Mie report which is now
c    available electronically and which is referred to using the
c    shorthand (Rn), meaning Eq. (n) of the report.

c    CALLING TREE:

c        MIEV0
c            TESTMI
c                TSTBAD
c                MIPRNT
c                ERRMSG
c            CKINMI
c                WRTBAD
c                WRTDIM
c                ERRMSG
c            SMALL1
c            SMALL2
c            ERRMSG
c            BIGA
c                CONFRA
c                    ERRMSG
c            MIPRNT
c ----------------------------------------------------------------------


      IMPLICIT  NONE

c ----------------------------------------------------------------------
c --------  I / O SPECIFICATIONS FOR SUBROUTINE MIEV0  -----------------
c ----------------------------------------------------------------------
      LOGICAL  ANYANG, PERFCT, PRNT(*)
      INTEGER  IPOLZN, MOMDIM, NUMANG, NMOM
      REAL     GQSC, MIMCUT, PMOM( 0:MOMDIM, * ), QEXT, QSCA, SPIKE,
     &         XMU(*), XX
      COMPLEX  CREFIN, SFORW, SBACK, S1(*), S2(*), TFORW(*), TBACK(*)
c ----------------------------------------------------------------------

c                                  ** NOTE --  MaxTrm = 10100  is neces-
c                                  ** sary to do some of the test probs,
c                                  ** but 1100 is sufficient for most
c                                  ** conceivable applications
c     .. Parameters ..

      INTEGER   MAXANG, MXANG2
      PARAMETER ( MAXANG = 501, MXANG2 = MAXANG / 2 + 1 )
      INTEGER   MAXTRM
      PARAMETER ( MAXTRM = 10100 )
      REAL      ONETHR
      PARAMETER ( ONETHR = 1. / 3. )
c     ..
c     .. Local Scalars ..

      LOGICAL   NOABS, PASS1, YESANG
      INTEGER   I, J, N, NANGD2, NTRM
      REAL      CHIN, CHINM1, COEFF, DENAN, DENBN, FN, MIM, MM, MRE,
     &          NP1DN, PSIN, PSINM1, RATIO, RIORIV, RN, RTMP, TAUN,
     &          TCOEF, TWONP1, XINV
      COMPLEX   AN, ANM1, ANP, ANPM, BN, BNM1, BNP, BNPM, CDENAN,
     &          CDENBN, CIOR, CIORIV, CSUM1, CSUM2, CTMP, ZET, 
     &          ZETN, ZETNM1
c     ..
c     .. Local Arrays ..

      REAL      PIN( MAXANG ), PINM1( MAXANG ), RBIGA( MAXTRM )
      COMPLEX   CBIGA( MAXTRM ), LITA( 2 ), LITB( 2 ), SM( MAXANG ),
     &          SMS( MXANG2 ), SP( MAXANG ), SPS( MXANG2 )
c     ..
c     .. External Subroutines ..

      EXTERNAL  BIGA, CKINMI, ERRMSG, MIPRNT, SMALL1, SMALL2, TESTMI
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, CMPLX, CONJG, COS, MAX, MIN, REAL, SIN
c     ..
      SAVE      PASS1

c     .. Statement Functions ..

      REAL      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
c     ..
      DATA      PASS1 /.TRUE./


c                    ** Save some input variables and replace them
c                    ** with values needed to do the self-test

      IF( PASS1 ) CALL TESTMI( .FALSE., XX, CREFIN, MIMCUT, PERFCT,
     &                         ANYANG, NUMANG, XMU, QEXT, QSCA, GQSC,
     &                         SFORW, SBACK, S1, S2, TFORW, TBACK )

   10 CONTINUE
c                                        ** Check input and calculate
c                                        ** certain variables from input

      CALL CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, NMOM, ANYANG,
     &             XMU )


      IF( PERFCT .AND. XX.LE.0.1 ) THEN
c                                            ** Use totally-reflecting
c                                            ** small-particle limit

         CALL SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK,
     &                S1, S2, TFORW, TBACK, LITA, LITB )

         NTRM = 2
         GO TO  100

      END IF


      NOABS = .TRUE.

      IF( .NOT.PERFCT ) THEN

         CIOR = CREFIN
         IF( AIMAG( CIOR ).GT.0.0 ) CIOR   = CONJG( CIOR )

         MRE    =   REAL( CIOR )
         MIM    = - AIMAG( CIOR )
         NOABS  = MIM.LE.MIMCUT
         CIORIV = 1.0 / CIOR
         RIORIV = 1.0 / MRE


         IF( XX*MAX( 1.0, ABS(CIOR) ) .LE. 0.1 ) THEN

c                                    ** Use general-refractive-index
c                                    ** small-particle limit

            CALL SMALL2( XX, CIOR, MIM.GT.MIMCUT, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW,
     &                   TBACK, LITA, LITB )

            NTRM = 2
            GO TO  100

         END IF

      END IF


      NANGD2 = ( NUMANG + 1 ) / 2
      YESANG = NUMANG.GT.0

c                              ** Number of terms in Mie series; Eq R50
      IF( XX.LE.8.0 ) THEN

         NTRM   = XX + 4.*XX**ONETHR + 1.

      ELSE IF( XX.LT.4200. ) THEN

         NTRM   = XX + 4.05*XX**ONETHR + 2.

      ELSE

         NTRM   = XX + 4.*XX**ONETHR + 2.
         
      END IF

      IF( NTRM + 1.GT.MAXTRM ) 
     &    CALL ERRMSG('MIEV0--PARAMETER MaxTrm TOO SMALL',.TRUE.)

c                            ** Calculate logarithmic derivatives of
c                            ** J-Bessel-fcn., A-sub-(1 to NTrm)

      IF( .NOT.PERFCT ) CALL BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA,
     &                             CBIGA )

c                            ** Initialize Ricatti-Bessel functions
c                            ** (psi,chi,zeta)-sub-(0,1) for upward
c                            ** recurrence ( Eq. R19 )
      XINV   = 1.0 / XX
      PSINM1 = SIN( XX )
      CHINM1 = COS( XX )
      PSIN   = PSINM1*XINV - CHINM1
      CHIN   = CHINM1*XINV + PSINM1
      ZETNM1 = CMPLX( PSINM1, CHINM1 )
      ZETN   = CMPLX( PSIN, CHIN )
c                                     ** Initialize previous coeffi-
c                                     ** cients for GQSC series
      ANM1 = ( 0.0, 0.0 )
      BNM1 = ( 0.0, 0.0 )
c                             ** Initialize angular function  pi
c                             ** and sums for S+, S- ( Ref. 2, p. 1507 )
      IF( ANYANG ) THEN

         DO 20 J = 1, NUMANG
c                             ** Eq. R39
            PINM1( J ) = 0.0
            PIN( J )   = 1.0
            SP( J )    = ( 0.0, 0.0 )
            SM( J )    = ( 0.0, 0.0 )
   20    CONTINUE

      ELSE

         DO 30 J = 1, NANGD2
c                             ** Eq. R39
            PINM1( J ) = 0.0
            PIN( J )   = 1.0
            SP( J )    = ( 0.0, 0.0 )
            SM( J )    = ( 0.0, 0.0 )
            SPS( J )   = ( 0.0, 0.0 )
            SMS( J )   = ( 0.0, 0.0 )
   30    CONTINUE

      END IF
      
c                       ** Initialize Mie sums for efficiencies, etc.
      QSCA  = 0.0
      GQSC  = 0.0
      SFORW = ( 0., 0.)
      SBACK = ( 0., 0.)
      CSUM1 = ( 0., 0. )
      CSUM2 = ( 0., 0. )


c ---------  LOOP TO SUM MIE SERIES  -----------------------------------

      MM     = +1.0
      SPIKE  = 1.0

      DO 60  N = 1, NTRM
c                           ** Compute various numerical coefficients
         FN     = N
         RN     = 1.0 / FN
         NP1DN  = 1.0 + RN
         TWONP1 = 2*N + 1
         COEFF  = TWONP1 / ( FN*( N + 1 ) )
         TCOEF  = TWONP1*( FN*( N + 1 ) )

c                              ** Calculate Mie series coefficients
         IF( PERFCT ) THEN
c                                 ** Totally-reflecting case; Eq R/A.1,2

            AN     = ( ( FN*XINV )*PSIN - PSINM1 ) /
     &               ( ( FN*XINV )*ZETN - ZETNM1 )
            BN     = PSIN / ZETN

         ELSE IF( NOABS ) THEN
c                                      ** No-absorption case; Eq (R16)

            CDENAN = ( RIORIV*RBIGA(N) + (FN*XINV) ) * ZETN - ZETNM1
            AN =   ( ( RIORIV*RBIGA(N) + (FN*XINV) ) * PSIN - PSINM1 )
     &             / CDENAN
            CDENBN = (  MRE * RBIGA(N) + (FN*XINV) ) * ZETN - ZETNM1
            BN =   ( (  MRE * RBIGA(N) + (FN*XINV) ) * PSIN - PSINM1 )
     &             / CDENBN

         ELSE
c                                       ** Absorptive case; Eq (R16)

            CDENAN = ( CIORIV * CBIGA(N) + (FN*XINV) ) * ZETN - ZETNM1
            CDENBN = (   CIOR * CBIGA(N) + (FN*XINV) ) * ZETN - ZETNM1
            AN =   ( ( CIORIV * CBIGA(N) + (FN*XINV) ) * PSIN - PSINM1 )
     &             / CDENAN
            BN =   ( (   CIOR * CBIGA(N) + (FN*XINV) ) * PSIN - PSINM1 )
     &             / CDENBN
c                                         ** Eq (R7)

            QSCA = QSCA + TWONP1 * ( SQ( AN ) + SQ( BN ) )

         END IF


         IF( .NOT.PERFCT .AND. N.GT.XX ) THEN
c                                               ** Flag resonance spikes
            DENAN  = ABS( CDENAN )
            DENBN  = ABS( CDENBN )
c                                                   ** Eq. R/B.9
            RATIO  = DENAN / DENBN
c                                                   ** Eq. R/B.10
            IF( RATIO.LE.0.2 .OR. RATIO.GE.5.0 ) 
     &          SPIKE  = MIN( SPIKE, DENAN, DENBN )

         END IF
c                                  ** Increment Mie sums for non-angle-
c                                  ** dependent quantities

c                                                   ** Eq. R/B.2
         SFORW = SFORW + TWONP1 * ( AN + BN )
c                                                   ** Eq. R/B.5,6
         CSUM1 = CSUM1 + TCOEF *( AN - BN )
c                                                   ** Eq. R/B.1
         SBACK = SBACK      + ( MM * TWONP1 ) * ( AN - BN )
c                                                   ** Eq. R/B.7,8
         CSUM2 = CSUM2 + ( MM*TCOEF ) *( AN + BN )

c                                         ** Eq (R8)

         GQSC  = GQSC + ( FN - RN ) * REAL( ANM1 * CONJG( AN )
     &                                    + BNM1 * CONJG( BN ) )
     &          + COEFF * REAL( AN * CONJG( BN ) )


         IF( YESANG ) THEN
c                                      ** Put Mie coefficients in form
c                                      ** needed for computing S+, S-
c                                      ** ( Eq R10 )
            ANP = COEFF * ( AN + BN )
            BNP = COEFF * ( AN - BN )
c                                      ** Increment Mie sums for S+, S-
c                                      ** while upward recursing
c                                      ** angular functions pi and tau
            IF( ANYANG ) THEN
c                                         ** Arbitrary angles

c                                              ** vectorizable loop
               DO 40 J = 1, NUMANG
c                                                 ** Eq. (R37b)

                  RTMP = ( XMU( J ) * PIN( J ) ) - PINM1( J )

c                                                 ** Eq. (R38b)
                  TAUN =  FN * RTMP - PINM1( J )

c                                                   ** Eq (R10)

                  SP( J )  = SP( J ) + ANP * ( PIN( J ) + TAUN )
                  SM( J )  = SM( J ) + BNP * ( PIN( J ) - TAUN )
                  PINM1( J ) = PIN( J )
c                                                 ** Eq. R37c

                  PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
   40          CONTINUE

            ELSE
c                                  ** Angles symmetric about 90 degrees
               ANPM   = MM * ANP
               BNPM   = MM * BNP
c                                          ** vectorizable loop
               DO 50 J = 1, NANGD2
c                                                 ** Eq. (R37b)

                  RTMP = ( XMU( J ) * PIN( J ) ) - PINM1( J )

c                                                 ** Eq. (R38b)
                  TAUN =  FN * RTMP - PINM1( J )

c                                                 ** Eq (R10,12)

                  SP ( J ) = SP ( J ) +  ANP * ( PIN( J ) + TAUN )
                  SMS( J ) = SMS( J ) + BNPM * ( PIN( J ) + TAUN )
                  SM ( J ) = SM ( J ) +  BNP * ( PIN( J ) - TAUN )
                  SPS( J ) = SPS( J ) + ANPM * ( PIN( J ) - TAUN )
                  PINM1( J ) = PIN( J )
c                                                 ** Eq. R37c

                  PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
   50          CONTINUE

            END IF

         END IF
c                          ** Update relevant quantities for next
c                          ** pass through loop
         MM     = -MM
         ANM1   = AN
         BNM1   = BN
c                           ** Upward recurrence for Ricatti-Bessel
c                           ** functions ( Eq. R17 )

         ZET    = ( TWONP1*XINV )*ZETN - ZETNM1
         ZETNM1 = ZETN
         ZETN   = ZET
         PSINM1 = PSIN
         PSIN   = REAL( ZETN )

   60 CONTINUE

c ---------- END LOOP TO SUM MIE SERIES --------------------------------


c                                         ** Eq (R6)
      QEXT   = 2./ XX**2 * REAL( SFORW )

      IF( PERFCT .OR. NOABS ) THEN

         QSCA = QEXT

      ELSE

         QSCA = 2./ XX**2 * QSCA
         
      END IF

      GQSC   = 4./ XX**2 * GQSC
      SFORW = 0.5 * SFORW
      SBACK = 0.5 * SBACK
      TFORW( 1 ) =  0.5*SFORW - 0.125*CSUM1
      TFORW( 2 ) =  0.5*SFORW + 0.125*CSUM1
      TBACK( 1 ) = -0.5*SBACK + 0.125*CSUM2
      TBACK( 2 ) =  0.5*SBACK + 0.125*CSUM2


      IF( YESANG ) THEN
c                                ** Recover scattering amplitudes
c                                ** from S+, S- ( Eq (R11) )

         IF( ANYANG ) THEN
c                                         ** vectorizable loop
            DO 70 J = 1, NUMANG
c                                                   ** Eq (R11)
               S1( J ) = 0.5*( SP( J ) + SM( J ) )
               S2( J ) = 0.5*( SP( J ) - SM( J ) )
   70       CONTINUE

         ELSE
c                                         ** vectorizable loop
            DO 80 J = 1, NANGD2
c                                                   ** Eq (R11)
               S1( J ) = 0.5*( SP( J ) + SM( J ) )
               S2( J ) = 0.5*( SP( J ) - SM( J ) )
   80       CONTINUE
c                                         ** vectorizable loop
            DO 90 J = 1, NANGD2
               S1( NUMANG + 1 - J ) = 0.5*( SPS( J ) + SMS( J ) )
               S2( NUMANG + 1 - J ) = 0.5*( SPS( J ) - SMS( J ) )
   90       CONTINUE

         END IF

      END IF


  100 CONTINUE
      IF( AIMAG( CREFIN ).GT.0.0 ) THEN
c                                         ** Take complex conjugates
c                                         ** of scattering amplitudes
         SFORW  = CONJG( SFORW )
         SBACK  = CONJG( SBACK )

         DO 110 I = 1, 2
            TFORW( I ) = CONJG( TFORW( I ) )
            TBACK( I ) = CONJG( TBACK( I ) )
  110    CONTINUE

         DO 120 J = 1, NUMANG
            S1( J ) = CONJG( S1( J ) )
            S2( J ) = CONJG( S2( J ) )
  120    CONTINUE

      END IF


      IF( PASS1 ) THEN
c                           ** Compare test case results with
c                           ** correct answers and abort if bad;
c                           ** otherwise restore user input and proceed

         CALL TESTMI( .TRUE., XX, CREFIN, MIMCUT, PERFCT, ANYANG,
     &                NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK, S1,
     &                S2, TFORW, TBACK )

         PASS1  = .FALSE.
         GO TO  10

      END IF


      IF ( PRNT(1) .OR. PRNT(2) )
     &   CALL  MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                 QSCA, GQSC, SFORW, SBACK, TFORW, TBACK, S1, S2 )

      RETURN

      END

      SUBROUTINE CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, NMOM,
     &                   ANYANG, XMU )

c        Check for bad input to MIEV0
c        (NoPMOM version)

c     Routines called :  ERRMSG, WRTBAD, WRTDIM


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   ANYANG, PERFCT
      INTEGER   MAXANG, NMOM, NUMANG
      REAL      XX
      COMPLEX   CREFIN
c     ..
c     .. Array Arguments ..

      REAL      XMU( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   INPERR
      INTEGER   I
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      EXTERNAL  WRTBAD, WRTDIM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, REAL
c     ..


      INPERR = .FALSE.

      IF( NUMANG.GT.MAXANG ) INPERR = WRTDIM( 'MaxAng', NUMANG )
      IF( NUMANG.LT.0 ) INPERR = WRTBAD( 'NUMANG' )

      IF( XX.LT.0.) INPERR = WRTBAD( 'XX' )

      IF( .NOT.PERFCT .AND. REAL( CREFIN ).LE.0. )
     &    INPERR = WRTBAD( 'CREFIN' )


      IF( NMOM.NE.0 ) THEN

         INPERR = WRTBAD( 'NMOM' )
         WRITE( *, '(A)' )
     &      ' For nonzero NMOM, use complete version of MIEV0'
      END IF


      IF( ANYANG ) THEN
c                                ** Allow for slight imperfections in
c                                ** computation of cosine
         DO 10 I = 1, NUMANG

             IF ( XMU(I).LT.-1.00001 .OR. XMU(I).GT.1.00001 )
     &            INPERR = WRTBAD( 'XMU' )

   10    CONTINUE

      ELSE

         DO 20 I = 1, ( NUMANG + 1 ) / 2

             IF ( XMU(I).LT.-0.00001 .OR. XMU(I).GT.1.00001 )
     &            INPERR = WRTBAD( 'XMU' )

   20    CONTINUE

      END IF


      IF( INPERR ) CALL ERRMSG( 'MIEV0--INPUT ERROR(S).  Aborting...',
     &                          .TRUE.)

      IF( XX.GT.20000.0 .OR. REAL( CREFIN ).GT.10.0 .OR.
     &    ABS( AIMAG( CREFIN ) ).GT.10.0 )
     &    CALL ERRMSG( 'MIEV0--XX or CREFIN outside tested range',
     &    .FALSE. )

      RETURN
      END

      SUBROUTINE BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA )

c        Calculate logarithmic derivatives of J-Bessel-function

c     Input :  CIOR, XX, NTRM, NOABS, YESANG  (defined in MIEV0)

c    Output :  RBIGA or CBIGA  (defined in MIEV0)

c    Routines called :  CONFRA

c    INTERNAL VARIABLES :

c       CONFRA     Value of Lentz continued fraction for cBigA(NTrm),
c                     used to initialize downward recurrence

c       DOWN       = True, use down-recurrence.  False, do not.

c       F1,F2,F3   Arithmetic statement functions used in determining
c                     whether to use up-  or down-recurrence
c                     ( Eqs. R47a,b,c )

c       MRE        Real refractive index
c       MIM        Imaginary refractive index

c       REZINV     1 / ( MRE * XX ); temporary variable for recurrence
c       ZINV       1 / ( CIOR * XX ); temporary variable for recurrence


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   NOABS, YESANG
      INTEGER   NTRM
      REAL      XX
      COMPLEX   CIOR
c     ..
c     .. Array Arguments ..

      REAL      RBIGA( * )
      COMPLEX   CBIGA( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   DOWN
      INTEGER   N
      REAL      MIM, MRE, REZINV, RTMP
      COMPLEX   CTMP, ZINV
c     ..
c     .. External Functions ..

      COMPLEX   CONFRA
      EXTERNAL  CONFRA
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, COS, EXP, REAL, SIN
c     ..
c     .. Statement Functions ..

      REAL      F1, F2, F3
c     ..
c     .. Statement Function definitions ..

c                                                   ** Eq. R47c
      F1( MRE ) = -8.0 + MRE**2*( 26.22 +
     &            MRE*( -0.4474 + MRE**3*( 0.00204 - 0.000175*MRE ) ) )
     
c                                                   ** Eq. R47b
      F2( MRE ) = 3.9 + MRE*( -10.8 + 13.78*MRE )
c                                                   ** Eq. R47a      
      F3( MRE ) = -15.04 + MRE*( 8.42 + 16.35*MRE )
c     ..

c                                  ** Decide whether BigA can be
c                                  ** calculated by up-recurrence
      MRE = REAL( CIOR )
      MIM = ABS( AIMAG( CIOR ) )

      IF( MRE.LT.1.0 .OR. MRE.GT.10.0 .OR. MIM.GT.10.0 ) THEN

         DOWN = .TRUE.

      ELSE IF( YESANG ) THEN

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F2( MRE ) ) DOWN   = .FALSE.

      ELSE

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F1( MRE ) ) DOWN   = .FALSE.

      END IF


      ZINV   = 1.0 / ( CIOR*XX )
      REZINV = 1.0 / ( MRE*XX )

      IF( DOWN ) THEN
c                          ** Compute initial high-order BigA using
c                          ** Lentz method ( Ref. 1, pp. 17-20 )

         CTMP = CONFRA( NTRM, ZINV )

c                                   *** Downward recurrence for BigA
         IF( NOABS ) THEN
c                                        ** No-absorption case; Eq (R23)
            RBIGA( NTRM ) = REAL( CTMP )

            DO 10 N = NTRM, 2, -1
               RBIGA( N - 1 ) = ( N*REZINV ) -
     &                          1.0 / ( ( N*REZINV ) + RBIGA( N ) )
   10       CONTINUE

         ELSE
c                                         ** Absorptive case; Eq (R23)
            CBIGA( NTRM ) = CTMP

            DO 20 N = NTRM, 2, -1
               CBIGA( N-1 ) = (N*ZINV) - 1.0 / ( (N*ZINV) + CBIGA( N ) )
   20       CONTINUE

         END IF


      ELSE
c                              *** Upward recurrence for BigA
         IF( NOABS ) THEN
c                                  ** No-absorption case; Eq (R20,21)
            RTMP   = SIN( MRE*XX )
            RBIGA( 1 ) = -REZINV + RTMP /
     &                   ( RTMP*REZINV - COS( MRE*XX ) )

            DO 30 N = 2, NTRM
               RBIGA( N ) = -( N*REZINV ) +
     &                      1.0 / ( ( N*REZINV ) - RBIGA( N - 1 ) )
   30       CONTINUE

         ELSE
c                                     ** Absorptive case; Eq (R20,22)

            CTMP = EXP( - (0.,2.) * CIOR * XX )
            CBIGA( 1 ) = - ZINV + (1.-CTMP) /
     &                          ( ZINV * (1.-CTMP) - (0.,1.)*(1.+CTMP) )

            DO 40 N = 2, NTRM
               CBIGA( N ) = - (N*ZINV) + 1.0 / ((N*ZINV) - CBIGA( N-1 ))
   40       CONTINUE

         END IF

      END IF


      RETURN
      END

      COMPLEX FUNCTION CONFRA( N, ZINV )

c         Compute Bessel function ratio A-sub-N from its
c         continued fraction using Lentz method

c         ZINV = Reciprocal of argument of A


c    I N T E R N A L    V A R I A B L E S
c    ------------------------------------

c    CAK      Term in continued fraction expansion of A (Eq. R25)

c    CAPT     Factor used in Lentz iteration for A (Eq. R27)

c    CNUMER   Numerator   in capT  ( Eq. R28A )
c    CDENOM   Denominator in capT  ( Eq. R28B )

c    CDTD     Product of two successive denominators of capT factors
c                 ( Eq. R34C )
c    CNTN     Product of two successive numerators of capT factors
c                 ( Eq. R34B )

c    EPS1     Ill-conditioning criterion
c    EPS2     Convergence criterion

c    KK       Subscript k of cAk  ( Eq. R25B )

c    KOUNT    Iteration counter ( used to prevent infinite looping )

c    MAXIT    Max. allowed no. of iterations

c    MM       + 1  and - 1, alternately
c --------------------------------------------------------------------

      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   N
      COMPLEX   ZINV
c     ..
c     .. Local Scalars ..

      INTEGER   KK, KOUNT, MAXIT, MM
      REAL      EPS1, EPS2
      COMPLEX   CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, REAL
c     ..
      DATA      EPS1 / 1.E-2 / , EPS2 / 1.E-8 /
      DATA      MAXIT / 10000 /


c                                 ** Eq. R25a
      CONFRA = ( N + 1 ) * ZINV
      MM     = - 1
      KK     = 2*N + 3
c                                 ** Eq. R25b, k=2
      CAK    = ( MM*KK ) * ZINV
      CDENOM = CAK
      CNUMER = CDENOM + 1.0 / CONFRA
      KOUNT  = 1

   10 CONTINUE
      KOUNT = KOUNT + 1

      IF( KOUNT.GT.MAXIT )
     &    CALL ERRMSG('ConFra--Iteration failed to converge',.TRUE.)

      MM  = - MM
      KK  = KK + 2
c                                 ** Eq. R25b
      CAK = ( MM*KK ) * ZINV
c                                          ** Eq. R32
      IF( ABS( CNUMER / CAK ).LE.EPS1 .OR.
     &    ABS( CDENOM / CAK ).LE.EPS1 ) THEN

c                                  ** Ill-conditioned case -- stride
c                                  ** two terms instead of one

c                                       ** Eq. R34
         CNTN   = CAK * CNUMER + 1.0
         CDTD   = CAK * CDENOM + 1.0
c                                           ** Eq. R33
         CONFRA = ( CNTN / CDTD ) * CONFRA

         MM  = - MM
         KK  = KK + 2
c                                 ** Eq. R25b
         CAK = ( MM*KK ) * ZINV
c                                      ** Eq. R35
         CNUMER = CAK + CNUMER / CNTN
         CDENOM = CAK + CDENOM / CDTD
         KOUNT  = KOUNT + 1
         GO TO  10

      ELSE
c                           *** Well-conditioned case

c                                  ** Eq. R27
         CAPT   = CNUMER / CDENOM
c                                  ** Eq. R26
         CONFRA = CAPT * CONFRA
c                                  ** Check for convergence; Eq. R31

         IF (      ABS( REAL (CAPT) - 1.0 ).GE.EPS2
     &        .OR. ABS( AIMAG(CAPT) )      .GE.EPS2 )  THEN

c                                        ** Eq. R30
            CNUMER = CAK + 1.0 / CNUMER
            CDENOM = CAK + 1.0 / CDENOM
            GO TO  10

         END IF

      END IF


      RETURN

      END

      SUBROUTINE MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, SFORW, SBACK, TFORW, TBACK, S1,
     &                   S2 )

c         Print scattering quantities of a single particle
c         (NoPMOM version)


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   PERFCT
      INTEGER   NUMANG
      REAL      GQSC, QEXT, QSCA, XX
      COMPLEX   CREFIN, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      LOGICAL   PRNT( * )
      REAL      XMU( * )
      COMPLEX   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I
      REAL      I1, I2
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CONJG, REAL
c     ..

      IF( PERFCT ) WRITE( *, '(''1'',10X,A,1P,E11.4)' )
     &    'Perfectly Conducting Case, size parameter =', XX

      IF( .NOT.PERFCT ) WRITE( *, '(''1'',10X,3(A,1P,E11.4))' )
     &    'Refractive Index:  Real ', REAL( CREFIN ), '  Imag ',
     &    AIMAG( CREFIN ), ',   Size Parameter =', XX


      IF( PRNT( 1 ) .AND. NUMANG.GT.0 ) THEN

         WRITE( *, '(/,A)' )
     &      '    cos(angle)  ------- S1 ---------  ------- S2 ---------'
     &      // '  --- S1*conjg(S2) ---   i1=S1**2   i2=S2**2  (i1+i2)/2'
     &      // '  DEG POLZN'

         DO 10 I = 1, NUMANG
            I1     = REAL( S1( I ) )**2 + AIMAG( S1( I ) )**2
            I2     = REAL( S2( I ) )**2 + AIMAG( S2( I ) )**2
            WRITE( *, '( I4, F10.6, 1P,10E11.3 )'   )
     &              I, XMU(I), S1(I), S2(I), S1(I)*CONJG(S2(I)),
     &              I1, I2, 0.5*(I1+I2), (I2-I1)/(I2+I1)
   10    CONTINUE

      END IF


      IF( PRNT( 2 ) ) THEN

         WRITE ( *, '(/,A,9X,A,17X,A,17X,A,/,(0P,F7.2, 1P,6E12.3) )' )
     &           '  Angle', 'S-sub-1', 'T-sub-1', 'T-sub-2',
     &               0.0,     SFORW,    TFORW(1),  TFORW(2),
     &              180.,     SBACK,    TBACK(1),  TBACK(2)
         WRITE ( *, '(/,4(A,1P,E11.4))' )
     &           ' Efficiency Factors,  extinction:', QEXT,
     &                              '   scattering:', QSCA,
     &                              '   absorption:', QEXT-QSCA,
     &                           '   rad. pressure:', QEXT-GQSC

      END IF


      RETURN

      END

      SUBROUTINE SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
     &                   SBACK, S1, S2, TFORW, TBACK, A, B )

c       Small-particle limit of Mie quantities in totally reflecting
c       limit ( Mie series truncated after 2 terms )

c        A,B       First two Mie coefficients, with numerator and
c                  denominator expanded in powers of XX ( a factor
c                  of XX**3 is missing but is restored before return
c                  to calling program )  ( Ref. 2, p. 1508 )


      IMPLICIT  NONE

c     .. Parameters ..

      REAL      TWOTHR, FIVTHR, FIVNIN
      PARAMETER ( TWOTHR = 2. / 3., FIVTHR = 5. / 3., FIVNIN = 5. / 9. )
c     ..
c     .. Scalar Arguments ..

      INTEGER   NUMANG
      REAL      GQSC, QEXT, QSCA, XX
      COMPLEX   SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL      XMU( * )
      COMPLEX   A( * ), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER   J
      REAL      RTMP
      COMPLEX   CTMP
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CMPLX, CONJG, REAL
c     ..
c     .. Statement Functions ..

      REAL      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
c     ..

c                                                       ** Eq. R/A.5
      A( 1 ) = CMPLX( 0., TWOTHR*( 1. - 0.2*XX**2 ) ) /
     &         CMPLX( 1. - 0.5*XX**2, TWOTHR*XX**3 )
c                                                       ** Eq. R/A.6
      B( 1 ) = CMPLX( 0., -( 1. - 0.1*XX**2 ) / 3. ) /
     &         CMPLX( 1. + 0.5*XX**2, -XX**3 / 3. )
c                                                       ** Eq. R/A.7,8
      A( 2 ) = CMPLX( 0.,   XX**2 / 30. )
      B( 2 ) = CMPLX( 0., - XX**2 / 45. )
c                                                       ** Eq. R/A.9
      QSCA = 6.* XX**4 *( SQ( A(1) ) + SQ( B(1) ) +
     &                      FIVTHR*( SQ( A(2) ) + SQ( B(2) ) ) )
      QEXT = QSCA
c                                                       ** Eq. R/A.10
      GQSC = 6.* XX**4 * REAL( A(1)*CONJG( A(2) + B(1) ) +
     &                        ( B(1) + FIVNIN*A(2) )*CONJG( B(2) ) )

      RTMP   = 1.5 * XX**3
      SFORW  = RTMP * ( A(1) + B(1) + FIVTHR * ( A(2) + B(2) ) )
      SBACK  = RTMP * ( A(1) - B(1) - FIVTHR * ( A(2) - B(2) ) )
      TFORW( 1 ) = RTMP*( B(1) + FIVTHR * ( 2.*B(2) - A(2) ) )
      TFORW( 2 ) = RTMP*( A(1) + FIVTHR * ( 2.*A(2) - B(2) ) )
      TBACK( 1 ) = RTMP*( B(1) - FIVTHR * ( 2.*B(2) + A(2) ) )
      TBACK( 2 ) = RTMP*( A(1) - FIVTHR * ( 2.*A(2) + B(2) ) )


      DO 10 J = 1, NUMANG
c                                                    ** Eq. R/A.11,12
         S1( J ) = RTMP*( A(1) + B(1)*XMU( J ) +
     &                    FIVTHR*( A(2)*XMU( J ) + 
     &                             B(2)*( 2.*XMU( J )**2 - 1.) ) )
         S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) +
     &                    FIVTHR*( B(2)*XMU( J ) + 
     &                             A(2)*( 2.*XMU( J )**2 - 1.) ) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = XX**3 * B(2)

      RETURN
      END

      SUBROUTINE SMALL2( XX, CIOR, CALCQE, NUMANG, XMU, QEXT, QSCA,
     &                   GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK, 
     &                   A, B )

c       Small-particle limit of Mie quantities for general refractive
c       index ( Mie series truncated after 2 terms )

c        A,B       First two Mie coefficients, with numerator and
c                  denominator expanded in powers of XX ( a factor
c                  of XX**3 is missing but is restored before return
c                  to calling program )

c        CIORSQ    Square of refractive index


      IMPLICIT  NONE

c     .. Parameters ..

      REAL      TWOTHR, FIVTHR
      PARAMETER  ( TWOTHR = 2./3., FIVTHR = 5./3. )
c     ..
c     .. Scalar Arguments ..

      LOGICAL   CALCQE
      INTEGER   NUMANG
      REAL      GQSC, QEXT, QSCA, XX
      COMPLEX   CIOR, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL      XMU( * )
      COMPLEX   A( * ), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER   J
      REAL      RTMP
      COMPLEX   CIORSQ, CTMP
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CMPLX, CONJG, REAL
c     ..
c     .. Statement Functions ..

      REAL      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
c     ..


      CIORSQ = CIOR**2
      CTMP   = CMPLX( 0., TWOTHR )*( CIORSQ - 1.0 )

c                                           ** Eq. R42a
      A( 1 ) = CTMP*( 1. - 0.1*XX**2 +
     &         ( CIORSQ / 350.+ 1./ 280.)*XX**4 ) /
     &         ( CIORSQ + 2.+ ( 1.- 0.7*CIORSQ )*XX**2 -
     &         ( CIORSQ**2 / 175.- 0.275*CIORSQ + 0.25 )*XX**4 +
     &         XX**3 * CTMP*( 1.- 0.1*XX**2 ) )

c                                           ** Eq. R42b
      B( 1 ) = ( XX**2 / 30.) * CTMP * ( 1.+
     &         ( CIORSQ / 35.- 1./ 14.)*XX**2 ) /
     &         ( 1.- ( CIORSQ / 15.- 1./ 6.)*XX**2 )

c                                           ** Eq. R42c

      A( 2 ) = ( 0.1*XX**2 )*CTMP*( 1.- XX**2 / 14.) /
     &         ( 2.*CIORSQ + 3.- ( CIORSQ / 7.- 0.5 )*XX**2 )

c                                           ** Eq. R40a

      QSCA = (6.*XX**4) * ( SQ( A(1) ) + SQ( B(1) ) +
     &                       FIVTHR*SQ( A(2) ) )

c                                           ** Eq. R40b
      QEXT = QSCA
      IF( CALCQE ) QEXT = 6.*XX * REAL( A(1) + B(1) + FIVTHR*A(2) )

c                                           ** Eq. R40c

      GQSC = (6.*XX**4) * REAL( A(1)*CONJG( A(2) + B(1) ) )

      RTMP   = 1.5 * XX**3
      SFORW  = RTMP*( A(1) + B(1) + FIVTHR*A(2) )
      SBACK  = RTMP*( A(1) - B(1) - FIVTHR*A(2) )
      TFORW( 1 ) = RTMP*( B(1) - FIVTHR*A(2) )
      TFORW( 2 ) = RTMP*( A(1) + 2.*FIVTHR*A(2) )
      TBACK( 1 ) = TFORW(1)
      TBACK( 2 ) = RTMP*( A(1) - 2.*FIVTHR*A(2) )


      DO 10 J = 1, NUMANG
c                                      ** Eq. R40d,e

         S1( J ) = RTMP*( A(1) + ( B(1) + FIVTHR*A(2) )*XMU( J ) )
         S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) +
     &                    FIVTHR * A(2)*( 2.*XMU( J )**2 - 1.) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = ( 0., 0.)

      RETURN
      END

      SUBROUTINE TESTMI( COMPAR, XX, CREFIN, MIMCUT, PERFCT, ANYANG,
     &                   NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK,
     &                   S1, S2, TFORW, TBACK )

c         Set up to run test case when  COMPAR = False;  when  = True,
c         compare Mie code test case results with correct answers
c         and abort if even one result is inaccurate.
c         (NoPMOM version)

c         The test case is :  Mie size parameter = 10
c                             refractive index   = 1.5 - 0.1 i
c                             scattering angle = 140 degrees
c                             1 Sekera moment

c         Results for this case may be found among the test cases
c         at the end of reference (1).

c         *** NOTE *** When running on some computers, esp. in single
c         precision, the Accur criterion below may have to be relaxed.
c         However, if Accur must be set larger than 10**-3 for some
c         size parameters, your computer is probably not accurate
c         enough to do Mie computations for those size parameters.

c     Routines called :  ERRMSG, MIPRNT, TSTBAD


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   ANYANG, COMPAR, PERFCT
      INTEGER   NUMANG
      REAL      GQSC, MIMCUT, QEXT, QSCA, XX
      COMPLEX   CREFIN, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL      XMU( * )
      COMPLEX   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   ANYSAV, OK, PERSAV
      INTEGER   N, NUMSAV
      REAL      ACCUR, CALC, EXACT, MIMSAV, TESTGQ, TESTQE, TESTQS,
     &          XMUSAV, XXSAV
      COMPLEX   CRESAV, TESTS1, TESTS2, TESTSB, TESTSF
c     ..
c     .. Local Arrays ..

      LOGICAL   PRNT( 2 )
      COMPLEX   TESTTB( 2 ), TESTTF( 2 )
c     ..
c     .. External Functions ..

      LOGICAL   TSTBAD
      EXTERNAL  TSTBAD
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, MIPRNT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, REAL
c     ..
c     .. Statement Functions ..

      LOGICAL   WRONG
c     ..
      SAVE      XXSAV, CRESAV, MIMSAV, PERSAV, ANYSAV, NUMSAV, XMUSAV

      DATA   TESTQE / 2.459791 /,
     &       TESTQS / 1.235144 /,
     &       TESTGQ / 1.139235 /,
     &       TESTSF / ( 61.49476, -3.177994 ) /,
     &       TESTSB / ( 1.493434, 0.2963657 ) /,
     &       TESTS1 / ( -0.1548380, -1.128972) /,
     &       TESTS2 / ( 0.05669755, 0.5425681) /,
     &       TESTTF / ( 12.95238, -136.6436 ), ( 48.54238, 133.4656 ) /,
     &       TESTTB / ( 41.88414, -15.57833 ), ( 43.37758, -15.28196 )/
     
      DATA      ACCUR / 1.E-4 /
c     ..
c     .. Statement Function definitions ..

      WRONG( CALC, EXACT ) = ABS( ( CALC - EXACT ) / EXACT ).GT.ACCUR
c     ..


      IF( .NOT.COMPAR ) THEN
c                                   ** Save certain user input values
         XXSAV  = XX
         CRESAV = CREFIN
         MIMSAV = MIMCUT
         PERSAV = PERFCT
         ANYSAV = ANYANG
         NUMSAV = NUMANG
         XMUSAV = XMU( 1 )
c                                   ** Reset input values for test case
         XX     = 10.0
         CREFIN = ( 1.5, -0.1 )
         MIMCUT = 0.0
         PERFCT = .FALSE.
         ANYANG = .TRUE.
         NUMANG = 1
         XMU( 1 ) = -0.7660444

      ELSE
c                                    ** Compare test case results with
c                                    ** correct answers and abort if bad
         OK     = .TRUE.

        IF ( WRONG( QEXT,TESTQE ) )
     &       OK =  TSTBAD( 'QEXT', ABS((QEXT - TESTQE) / TESTQE) )
     
        IF ( WRONG( QSCA,TESTQS ) )
     &       OK =  TSTBAD( 'QSCA', ABS((QSCA - TESTQS) / TESTQS) )
     
        IF ( WRONG( GQSC,TESTGQ ) )
     &       OK =  TSTBAD( 'GQSC', ABS((GQSC - TESTGQ) / TESTGQ) )

         IF( WRONG( REAL( SFORW ),REAL( TESTSF ) ) .OR.
     &       WRONG( AIMAG( SFORW ),AIMAG( TESTSF ) ) )
     &       OK = TSTBAD( 'SFORW', ABS( ( SFORW - TESTSF ) / TESTSF ) )

         IF( WRONG( REAL( SBACK ),REAL( TESTSB ) ) .OR.
     &       WRONG( AIMAG( SBACK ),AIMAG( TESTSB ) ) )
     &       OK = TSTBAD( 'SBACK', ABS( ( SBACK - TESTSB ) / TESTSB ) )

         IF( WRONG( REAL( S1( 1 ) ),REAL( TESTS1 ) ) .OR.
     &       WRONG( AIMAG( S1( 1 ) ),AIMAG( TESTS1 ) ) )
     &       OK = TSTBAD( 'S1', ABS( ( S1( 1 ) - TESTS1 ) / TESTS1 ) )

         IF( WRONG( REAL( S2( 1 ) ),REAL( TESTS2 ) ) .OR.
     &       WRONG( AIMAG( S2( 1 ) ),AIMAG( TESTS2 ) ) )
     &       OK = TSTBAD( 'S2', ABS( ( S2( 1 ) - TESTS2 ) / TESTS2 ) )


         DO 10 N = 1, 2

           IF ( WRONG(  REAL(TFORW(N)),  REAL(TESTTF(N)) ) .OR.
     &          WRONG( AIMAG(TFORW(N)), AIMAG(TESTTF(N)) ) )
     &          OK =  TSTBAD( 'TFORW', ABS( (TFORW(N) - TESTTF(N)) /
     &                                       TESTTF(N) ) )
           IF ( WRONG(  REAL(TBACK(N)),  REAL(TESTTB(N)) ) .OR.
     &          WRONG( AIMAG(TBACK(N)), AIMAG(TESTTB(N)) ) )
     &          OK =  TSTBAD( 'TBACK', ABS( (TBACK(N) - TESTTB(N)) /
     &                                        TESTTB(N) ) )

   10    CONTINUE


         IF( .NOT.OK ) THEN

            PRNT( 1 ) = .TRUE.
            PRNT( 2 ) = .TRUE.

            CALL  MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                    QSCA, GQSC, SFORW, SBACK, TFORW, TBACK, S1,S2)

            CALL ERRMSG( 'MIEV0 -- Self-test failed', .TRUE.)

         END IF
c                                       ** Restore user input values
         XX     = XXSAV
         CREFIN = CRESAV
         MIMCUT = MIMSAV
         PERFCT = PERSAV
         ANYANG = ANYSAV
         NUMANG = NUMSAV
         XMU( 1 ) = XMUSAV

      END IF

      RETURN
      END

