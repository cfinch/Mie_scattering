
      SUBROUTINE MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,
     &                  NMOM, IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC,
     &                  PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK,
     &                  SPIKE )

c    Computes Mie scattering and extinction efficiencies; asymmetry
c    factor;  forward- and backscatter amplitude;  scattering
c    amplitudes vs. scattering angle for incident polarization parallel
c    and perpendicular to the plane of scattering;
c    coefficients in the Legendre polynomial expansions of either the
c    unpolarized phase function or the polarized phase matrix;
c    some quantities needed in polarized radiative transfer;  and
c    information about whether or not a resonance has been hit.

c    Input and output variables are described in file MIEV.doc. 
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
c            LPCOEF
c                LPCO1T
c                LPCO2T
c                ERRMSG
c            MIPRNT


c      I N T E R N A L   V A R I A B L E S
c      -----------------------------------

c  AN,BN           Mie coefficients a-sub-n, b-sub-n ( Ref. 1, Eq. 16 )

c  ANM1,BNM1       Mie coefficients  a-sub-(n-1),
c                     b-sub-(n-1);  used in GQSC sum

c  ANP             Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
c  BNP             Coeffs. in S- expansion ( Ref. 2, p. 1507 )
c  ANPM            Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
c                     when  MU  is replaced by  - MU
c  BNPM            Coeffs. in S- expansion ( Ref. 2, p. 1507 )
c                     when  MU  is replaced by  - MU

c  CALCMO(K)       TRUE, calculate moments for K-th phase quantity
c                     (derived from IPOLZN)

c  CBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
c                     ( COMPLEX version )

c  CDENAN,         (COMPLEX) denominators of An,Bn
c   CDENBN

c  CIOR            Complex index of refraction with negative
c                     imaginary part (Van de Hulst convention)
c  CIORIV          1 / cIoR

c  COEFF           ( 2N + 1 ) / ( N ( N + 1 ) )

c  CSUM1,2         temporary sum variables for TFORW, TBACK

c  FN              Floating point version of loop index for
c                     Mie series summation

c  LITA,LITB(N)    Mie coefficients An, Bn, saved in arrays for
c                     use in calculating Legendre moments PMOM

c  MAXTRM          Max. possible no. of terms in Mie series

c  MM              (-1)^(n+1), where n is Mie series sum index 

c  MIM             Magnitude of imaginary refractive index
c  MRE             Real part of refractive index

c  MAXANG          Max. possible value of input variable NUMANG
c  NANGD2          (NUMANG+1)/2 ( no. of angles in 0-90 deg; ANYANG=F )

c  NOABS           TRUE, sphere non-absorbing (determined by MIMCUT)

c  NP1DN           ( N + 1 ) / N

c  NPQUAN          Highest-numbered phase quantity for which moments are
c                     to be calculated (the largest digit in IPOLZN
c                     if  IPOLZN .NE. 0)

c  NTRM            No. of terms in Mie series

c  PASS1           TRUE on first entry, FALSE thereafter; for self-test

c  PIN(J)          Angular function pi-sub-n ( Ref. 2, Eq. 3 )
c                     at J-th angle
c  PINM1(J)        pi-sub-(n-1) ( see PIn ) at J-th angle

c  PSINM1          Ricatti-Bessel function psi-sub-(n-1), argument XX
c  PSIN            Ricatti-Bessel function psi-sub-n of argument XX
c                     ( Ref. 1, p. 11 ff. )

c  RBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
c                     ( REAL version, for when imag refrac index = 0 )

c  RIORIV          1 / Mre

c  RN              1 / N

c  RTMP            (REAL) temporary variable

c  SP(J)           S+  for J-th angle  ( Ref. 2, p. 1507 )
c  SM(J)           S-  for J-TH angle  ( Ref. 2, p. 1507 )
c  SPS(J)          S+  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )
c  SMS(J)          S-  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )

c  TAUN            Angular function tau-sub-n ( Ref. 2, Eq. 4 )
c                     at J-th angle

c  TCOEF           N ( N+1 ) ( 2N+1 ) (for summing TFORW,TBACK series)

c  TWONP1          2N + 1

c  YESANG          TRUE if scattering amplitudes are to be calculated

c  ZETNM1          Ricatti-Bessel function  zeta-sub-(n-1) of argument
c                     XX  ( Ref. 2, Eq. 17 )
c  ZETN            Ricatti-Bessel function  zeta-sub-n of argument XX
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

c                                  ** NOTE --  MAXTRM = 10100  is neces-
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
      INTEGER   I, J, N, NANGD2, NPQUAN, NTRM
      REAL      CHIN, CHINM1, COEFF, DENAN, DENBN, FN, MIM, MM, MRE,
     &          NP1DN, PSIN, PSINM1, RATIO, RIORIV, RN, RTMP, TAUN,
     &          TCOEF, TWONP1, XINV
      COMPLEX   AN, ANM1, ANP, ANPM, BN, BNM1, BNP, BNPM, CDENAN,
     &          CDENBN, CIOR, CIORIV, CSUM1, CSUM2, CTMP, ZET, 
     &          ZETN, ZETNM1
c     ..
c     .. Local Arrays ..

      LOGICAL   CALCMO( 4 )
      REAL      PIN( MAXANG ), PINM1( MAXANG ), RBIGA( MAXTRM )
      COMPLEX   CBIGA( MAXTRM ), LITA( MAXTRM ), LITB( MAXTRM ),
     &          SM( MAXANG ), SMS( MXANG2 ), SP( MAXANG ), SPS( MXANG2 )
c     ..
c     .. External Subroutines ..

      EXTERNAL  BIGA, CKINMI, ERRMSG, LPCOEF, MIPRNT, SMALL1, SMALL2,
     &          TESTMI
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
      DATA      PASS1 / .TRUE. /


c                    ** Save some input variables and replace them
c                    ** with values needed to do the self-test

      IF( PASS1 ) CALL TESTMI( .FALSE., XX, CREFIN, MIMCUT, PERFCT,
     &                         ANYANG, NMOM, IPOLZN, NUMANG, XMU, QEXT,
     &                         QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW,
     &                         TBACK, PMOM, MOMDIM )

   10 CONTINUE
c                                        ** Check input and calculate
c                                        ** certain variables from input

      CALL CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM, NMOM,
     &             IPOLZN, ANYANG, XMU, CALCMO, NPQUAN )


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

         IF( AIMAG(CIOR).GT.0.0 ) CIOR = CONJG( CIOR )

         MRE    = REAL( CIOR )
         MIM    = -AIMAG( CIOR )
         NOABS  = MIM.LE.MIMCUT
         CIORIV = 1.0 / CIOR
         RIORIV = 1.0 / MRE

         IF( XX*MAX( 1.0, ABS(CIOR) ).LE.0.1 ) THEN

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

         NTRM = XX + 4.*XX**ONETHR + 1.

      ELSE IF( XX.LT.4200. ) THEN

         NTRM = XX + 4.05*XX**ONETHR + 2.

      ELSE

         NTRM = XX + 4.*XX**ONETHR + 2.

      END IF

      IF( NTRM+1 .GT. MAXTRM )
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
            PIN( J ) = 1.0

            SP( J ) = ( 0.0, 0.0 )
            SM( J ) = ( 0.0, 0.0 )
   20    CONTINUE

      ELSE

         DO 30 J = 1, NANGD2
c                             ** Eq. R39
            PINM1( J ) = 0.0
            PIN( J ) = 1.0

            SP( J ) = ( 0.0, 0.0 )
            SM( J ) = ( 0.0, 0.0 )
            SPS( J ) = ( 0.0, 0.0 )
            SMS( J ) = ( 0.0, 0.0 )
   30    CONTINUE

      END IF

c                       ** Initialize Mie sums for efficiencies, etc.
      QSCA  = 0.0
      GQSC  = 0.0
      SFORW = ( 0., 0. )
      SBACK = ( 0., 0. )
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
         COEFF  = TWONP1 / ( FN * ( N + 1 ) )
         TCOEF  = TWONP1 * ( FN * ( N + 1 ) )

c                           ** Calculate Mie series coefficients
         IF( PERFCT ) THEN
c                                 ** Totally-reflecting case; Eq R/A.1,2

            AN = ( ( FN*XINV )*PSIN - PSINM1 ) /
     &           ( ( FN*XINV )*ZETN - ZETNM1 )
            BN = PSIN / ZETN

         ELSE IF( NOABS ) THEN
c                                      ** No-absorption case; Eq (R16)

            CDENAN = ( RIORIV*RBIGA(N) + ( FN*XINV ) ) * ZETN - ZETNM1
            AN   = ( ( RIORIV*RBIGA(N) + ( FN*XINV ) ) * PSIN - PSINM1 )
     &             / CDENAN
            CDENBN = ( MRE*RBIGA(N) + ( FN*XINV ) ) * ZETN - ZETNM1
            BN   = ( ( MRE*RBIGA(N) + ( FN*XINV ) ) * PSIN - PSINM1 )
     &             / CDENBN

         ELSE
c                                       ** Absorptive case; Eq (R16)

            CDENAN = ( CIORIV*CBIGA( N ) + ( FN*XINV ) )*ZETN - ZETNM1
            CDENBN =   ( CIOR*CBIGA( N ) + ( FN*XINV ) )*ZETN - ZETNM1
            AN   = ( ( CIORIV*CBIGA( N ) + ( FN*XINV ) )*PSIN - PSINM1 )
     &             / CDENAN
            BN     = ( ( CIOR*CBIGA( N ) + ( FN*XINV ) )*PSIN - PSINM1 )
     &             / CDENBN
c                                         ** Eq (R7)

            QSCA   = QSCA + TWONP1*( SQ( AN ) + SQ( BN ) )

         END IF
c                       ** Save Mie coefficients for PMOM calculation

         LITA( N ) = AN
         LITB( N ) = BN


         IF( .NOT.PERFCT .AND. N.GT.XX ) THEN
c                                               ** Flag resonance spikes
            DENAN  = ABS( CDENAN )
            DENBN  = ABS( CDENBN )
c                                                   ** Eq. R/B.9
            RATIO  = DENAN / DENBN
c                                                   ** Eq. R/B.10
            IF( RATIO.LE.0.2 .OR. RATIO.GE.5.0 )
     &          SPIKE = MIN( SPIKE, DENAN, DENBN )

         END IF
c                                  ** Increment Mie sums for non-angle-
c                                  ** dependent quantities

c                                                   ** Eq. R/B.2
         SFORW = SFORW + TWONP1*( AN + BN )
c                                                   ** Eq. R/B.5,6
         CSUM1 = CSUM1 + TCOEF *( AN - BN )
c                                                   ** Eq. R/B.1
         SBACK = SBACK + ( MM*TWONP1 )*( AN - BN )
c                                                   ** Eq. R/B.7,8
         CSUM2 = CSUM2 + ( MM*TCOEF ) *( AN + BN )

c                                         ** Eq (R8)

         GQSC  = GQSC  + ( FN - RN ) * REAL( ANM1 * CONJG( AN ) +
     &                                       BNM1 * CONJG( BN ) )
     &           + COEFF * REAL( AN * CONJG( BN ) )


         IF( YESANG ) THEN
c                                      ** Put Mie coefficients in form
c                                      ** needed for computing S+, S-
c                                      ** ( Eq R10 )
            ANP = COEFF*( AN + BN )
            BNP = COEFF*( AN - BN )

c                                      ** Increment Mie sums for S+, S-
c                                      ** while upward recursing
c                                      ** angular functions pi and tau
            IF( ANYANG ) THEN
c                                         ** Arbitrary angles

c                                              ** vectorizable loop
               DO 40 J = 1, NUMANG
c                                                 ** Eq. (R37b)

                  RTMP = ( XMU(J) * PIN(J) ) - PINM1( J )

c                                                 ** Eq. (R38b)
                  TAUN   = FN * RTMP - PINM1( J )

c                                                   ** Eq (R10)

                  SP( J ) = SP( J ) + ANP * ( PIN( J ) + TAUN )
                  SM( J ) = SM( J ) + BNP * ( PIN( J ) - TAUN )

                  PINM1( J ) = PIN( J )
c                                                 ** Eq. R37c

                  PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
   40          CONTINUE

            ELSE
c                                  ** Angles symmetric about 90 degrees
               ANPM = MM*ANP
               BNPM = MM*BNP
c                                          ** vectorizable loop
               DO 50 J = 1, NANGD2
c                                                 ** Eq. (R37b)

                  RTMP = ( XMU(J) * PIN(J) ) - PINM1( J )

c                                                 ** Eq. (R38b)
                  TAUN = FN * RTMP - PINM1( J )

c                                                 ** Eq (R10,12)

                  SP ( J ) = SP ( J ) + ANP * ( PIN( J ) + TAUN )
                  SMS( J ) = SMS( J ) + BNPM *( PIN( J ) + TAUN )
                  SM ( J ) = SM ( J ) + BNP * ( PIN( J ) - TAUN )
                  SPS( J ) = SPS( J ) + ANPM *( PIN( J ) - TAUN )

                  PINM1( J ) = PIN( J )
c                                                 ** Eq. R37c

                  PIN( J ) = ( XMU(J) * PIN(J) ) + NP1DN * RTMP
   50          CONTINUE

            END IF

         END IF
c                          ** Update relevant quantities for next
c                          ** pass through loop
         MM   = - MM
         ANM1 = AN
         BNM1 = BN
c                           ** Upward recurrence for Ricatti-Bessel
c                           ** functions ( Eq. R17 )

         ZET    = ( TWONP1*XINV ) * ZETN - ZETNM1
         ZETNM1 = ZETN
         ZETN   = ZET
         PSINM1 = PSIN
         PSIN   = REAL( ZETN )

   60 CONTINUE

c ---------- END LOOP TO SUM MIE SERIES --------------------------------


c                                         ** Eq (R6)
      QEXT = 2. / XX**2*REAL( SFORW )

      IF( PERFCT .OR. NOABS ) THEN

         QSCA = QEXT

      ELSE

         QSCA = 2./ XX**2 * QSCA

      END IF

      GQSC   = 4./ XX**2 * GQSC
      SFORW  = 0.5*SFORW
      SBACK  = 0.5*SBACK
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
c                                 ** Calculate Legendre moments

  100 CONTINUE
      IF( NMOM.GT.0 ) CALL LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM, CALCMO,
     &                             NPQUAN, LITA, LITB, PMOM )


      IF( AIMAG( CREFIN ).GT.0.0 ) THEN
c                                         ** Take complex conjugates
c                                         ** of scattering amplitudes

         SFORW = CONJG( SFORW )
         SBACK = CONJG( SBACK )

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

         CALL TESTMI( .TRUE., XX, CREFIN, MIMCUT, PERFCT, ANYANG, NMOM,
     &                IPOLZN, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
     &                SBACK, S1, S2, TFORW, TBACK, PMOM, MOMDIM )

         PASS1  = .FALSE.
         GO TO  10

      END IF


      IF( PRNT( 1 ) .OR. PRNT( 2 ) ) 
     &  CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &               QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &               SFORW, SBACK, TFORW, TBACK, S1, S2 )

      RETURN

      END

      SUBROUTINE CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM,
     &                   NMOM, IPOLZN, ANYANG, XMU, CALCMO, NPQUAN )

c        Check for bad input to MIEV0 and calculate CALCMO, NPQUAN

c     Routines called :  ERRMSG, WRTBAD, WRTDIM


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   ANYANG, PERFCT
      INTEGER   IPOLZN, MAXANG, MOMDIM, NMOM, NPQUAN, NUMANG
      REAL      XX
      COMPLEX   CREFIN
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      REAL      XMU( * )
c     ..
c     .. Local Scalars ..

      CHARACTER STRING*4
      LOGICAL   INPERR
      INTEGER   I, IP, J, L
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      EXTERNAL  WRTBAD, WRTDIM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, ICHAR, MAX, REAL
c     ..


      INPERR = .FALSE.

      IF( NUMANG.GT.MAXANG ) INPERR = WRTDIM( 'MaxAng', NUMANG )
      IF( NUMANG.LT.0 ) INPERR = WRTBAD( 'NUMANG' )

      IF( XX.LT.0. ) INPERR = WRTBAD( 'XX' )

      IF( .NOT.PERFCT .AND. REAL( CREFIN ).LE.0. )
     &    INPERR = WRTBAD( 'CREFIN' )

      IF( MOMDIM.LT.0 ) INPERR = WRTBAD( 'MOMDIM' )


      IF( NMOM.NE.0 ) THEN

         IF( NMOM.LT.0 .OR. NMOM.GT.MOMDIM ) INPERR = WRTBAD( 'NMOM' )

         IF( ABS( IPOLZN ).GT.4444 ) INPERR = WRTBAD( 'IPOLZN' )

         NPQUAN = 0

         DO 10 L = 1, 4
            CALCMO( L ) = .FALSE.
   10    CONTINUE

         IF( IPOLZN.NE.0 ) THEN
c                                 ** Parse out IPOLZN into its digits
c                                 ** to find which phase quantities are
c                                 ** to have their moments calculated

            WRITE( STRING, '(I4)' ) ABS( IPOLZN )

            DO 20 J = 1, 4
               IP = ICHAR( STRING( J:J ) ) - ICHAR( '0' )

               IF( IP.GE.1 .AND. IP.LE.4 ) CALCMO( IP ) = .TRUE.

               IF( IP.EQ.0 .OR. ( IP.GE.5 .AND. IP.LE.9 ) )
     &             INPERR = WRTBAD( 'IPOLZN' )

               NPQUAN = MAX( NPQUAN, IP )
   20       CONTINUE

         END IF

      END IF


      IF( ANYANG ) THEN
c                                ** Allow for slight imperfections in
c                                ** computation of cosine
         DO 30 I = 1, NUMANG

            IF( XMU( I ).LT.-1.00001 .OR. XMU( I ).GT.1.00001 )
     &          INPERR = WRTBAD( 'XMU' )

   30    CONTINUE

      ELSE

         DO 40 I = 1, ( NUMANG + 1 ) / 2

            IF( XMU( I ).LT.-0.00001 .OR. XMU( I ).GT.1.00001 )
     &          INPERR = WRTBAD( 'XMU' )

   40    CONTINUE

      END IF


      IF( INPERR ) CALL ERRMSG( 'MIEV0--Input error(S).  Aborting...',
     &                          .TRUE. )

      IF( XX.GT.20000.0 .OR. REAL( CREFIN ).GT.10.0 .OR.
     &    ABS( AIMAG( CREFIN ) ).GT.10.0 )
     &    CALL ERRMSG( 'MIEV0--XX or CREFIN outside tested range',
     &    .FALSE.)

      RETURN
      END

      SUBROUTINE LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM, CALCMO, NPQUAN, A,
     &                   B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities ( Ref. 5 formulation )

c     INPUT:  NTRM                    Number terms in Mie series
c             NMOM, IPOLZN, MOMDIM    MIEV0 arguments
c             CALCMO                  Flags calculated from IPOLZN
c             NPQUAN                  Defined in MIEV0
c             A, B                    Mie series coefficients

c     OUTPUT: PMOM                   Legendre moments (MIEV0 argument)

c     Routines called :  ERRMSG, LPCO1T, LPCO2T

c     *** NOTES ***

c         (1)  Eqs. 2-5 are in error in Dave, Appl. Opt. 9,
c         1888 (1970).  Eq. 2 refers to M1, not M2;  eq. 3 refers to
c         M2, not M1.  In eqs. 4 and 5, the subscripts on the second
c         term in square brackets should be interchanged.

c         (2)  The general-case logic in this subroutine works correctly
c         in the two-term Mie series case, but subroutine LPCO2T
c         is called instead, for speed.

c         (3)  Subroutine  LPCO1T, to do the one-term case, is never
c         called within the context of MIEV0, but is included for
c         complete generality.

c         (4)  Some improvement in speed is obtainable by combining the
c         310- and 410-loops, if moments for both the third and fourth
c         phase quantities are desired, because the third phase quantity
c         is the real part of a complex series, while the fourth phase
c         quantity is the imaginary part of that very same series.  But
c         most users are not interested in the fourth phase quantity,
c         which is related to circular polarization, so the present
c         scheme is usually more efficient.


c           ** Definitions of local variables ***

c      AM(M)       Numerical coefficients  a-sub-m-super-l
c                     in Dave, Eqs. 1-15, as simplified in Ref. 5.

c      BI(I)       Numerical coefficients  b-sub-i-super-l
c                     in Dave, Eqs. 1-15, as simplified in Ref. 5.

c      BIDEL(I)    1/2 Bi(I) times factor capital-del in Dave

c      CM,DM()     Arrays C and D in Dave, Eqs. 16-17 (Mueller form),
c                     calculated using recurrence derived in Ref. 5

c      CS,DS()     Arrays C and D in Ref. 4, Eqs. A5-A6 (Sekera form),
c                     calculated using recurrence derived in Ref. 5

c      C,D()       Either CM,DM or CS,DS, depending on IPOLZN

c      EVENL       True for even-numbered moments;  false otherwise

c      IDEL        1 + little-del  in Dave

c      MAXTRM      Max. no. of terms in Mie series

c      MAXMOM      Max. no. of non-zero moments

c      NUMMOM      Number of non-zero moments

c      RECIP(K)    1 / K


      IMPLICIT  NONE

c     .. Parameters ..

      INTEGER   MAXTRM, MAXMOM, MXMOM2, MAXRCP
      PARAMETER ( MAXTRM = 1102, MAXMOM = 2*MAXTRM, MXMOM2 = MAXMOM / 2,
     &            MAXRCP = 4*MAXTRM + 2 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   IPOLZN, MOMDIM, NMOM, NPQUAN, NTRM
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      REAL      PMOM( 0:MOMDIM, * )
      COMPLEX   A( * ), B( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   EVENL, PASS1
      INTEGER   I, IDEL, IMAX, J, K, L, LD2, M, MMAX, NUMMOM
      REAL      SUM
c     ..
c     .. Local Arrays ..

      REAL      AM( 0:MAXTRM ), BI( 0:MXMOM2 ), BIDEL( 0:MXMOM2 ),
     &          RECIP( MAXRCP )
      COMPLEX   C( MAXTRM ), CM( MAXTRM ), CS( MAXTRM ), D( MAXTRM ),
     &          DM( MAXTRM ), DS( MAXTRM )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, LPCO1T, LPCO2T
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CONJG, MAX, MIN, MOD, REAL
c     ..
c     .. Equivalences ..

      EQUIVALENCE ( C, CM ), ( D, DM )
c     ..
      SAVE      PASS1, RECIP

      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN

         DO 10 K = 1, MAXRCP
            RECIP( K ) = 1.0 / K
   10    CONTINUE

         PASS1  = .FALSE.

      END IF


      DO 30 J = 1, MAX( 1, NPQUAN )

         DO 20 L = 0, NMOM
            PMOM( L, J ) = 0.0
   20    CONTINUE

   30 CONTINUE


      IF( NTRM.EQ.1 ) THEN

         CALL LPCO1T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

         RETURN

      ELSE IF( NTRM.EQ.2 ) THEN

         CALL LPCO2T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

         RETURN

      END IF


      IF( NTRM + 2.GT.MAXTRM )
     &    CALL ERRMSG('LPCoef--PARAMETER MaxTrm too small',.TRUE.)

c                                     ** Calculate Mueller C, D arrays
      CM( NTRM + 2 ) = ( 0., 0. )
      DM( NTRM + 2 ) = ( 0., 0. )
      CM( NTRM + 1 ) = ( 1. - RECIP( NTRM+1 ) ) * B( NTRM )
      DM( NTRM + 1 ) = ( 1. - RECIP( NTRM+1 ) ) * A( NTRM )
      CM( NTRM ) = ( RECIP( NTRM ) + RECIP( NTRM+1 ) ) * A( NTRM ) +
     &             ( 1. - RECIP( NTRM ) )*B( NTRM-1 )
      DM( NTRM ) = ( RECIP( NTRM ) + RECIP( NTRM+1 ) ) * B( NTRM ) +
     &             ( 1. - RECIP( NTRM ) )*A( NTRM-1 )

      DO 40 K = NTRM-1, 2, -1
         CM( K ) = CM( K+2 ) - ( 1. + RECIP(K+1) ) * B( K+1 )
     &                       + ( RECIP(K) + RECIP(K+1) ) * A( K )
     &                       + ( 1. - RECIP(K) ) * B( K-1 )
         DM( K ) = DM( K+2 ) - ( 1. + RECIP(K+1) ) * A( K+1 )
     &                       + ( RECIP(K) + RECIP(K+1) ) * B( K )
     &                       + ( 1. - RECIP(K) ) * A( K-1 )
   40 CONTINUE

      CM( 1 ) = CM( 3 ) + 1.5 * ( A( 1 ) - B( 2 ) )
      DM( 1 ) = DM( 3 ) + 1.5 * ( B( 1 ) - A( 2 ) )


      IF( IPOLZN.GE.0 ) THEN

         DO 50 K = 1, NTRM + 2
            C( K ) = ( 2*K - 1 ) * CM( K )
            D( K ) = ( 2*K - 1 ) * DM( K )
   50    CONTINUE

      ELSE
c                                    ** Compute Sekera C and D arrays
         CS( NTRM + 2 ) = ( 0., 0. )
         DS( NTRM + 2 ) = ( 0., 0. )
         CS( NTRM + 1 ) = ( 0., 0. )
         DS( NTRM + 1 ) = ( 0., 0. )

         DO 60 K = NTRM, 1, -1
            CS( K ) = CS( K+2 ) + ( 2*K + 1 ) * ( CM( K+1 ) - B( K ) )
            DS( K ) = DS( K+2 ) + ( 2*K + 1 ) * ( DM( K+1 ) - A( K ) )
   60    CONTINUE

         DO 70 K = 1, NTRM + 2
            C( K ) = ( 2*K - 1 ) * CS( K )
            D( K ) = ( 2*K - 1 ) * DS( K )
   70    CONTINUE

      END IF


      IF( IPOLZN.LT.0 ) NUMMOM = MIN( NMOM, 2*NTRM - 2 )
      IF( IPOLZN.GE.0 ) NUMMOM = MIN( NMOM, 2*NTRM )

      IF( NUMMOM.GT.MAXMOM )
     &    CALL ERRMSG('LPCoef--PARAMETER MaxTrm too small',.TRUE.)


c                          ** Loop over moments

      DO 240 L = 0, NUMMOM

         LD2 = L / 2
         EVENL  = MOD( L, 2 ).EQ.0
c                                    ** Calculate numerical coefficients
c                                    ** a-sub-m and b-sub-i in Dave
c                                    ** double-sums for moments
         IF( L.EQ.0 ) THEN

            IDEL = 1

            DO 80 M = 0, NTRM
               AM( M ) = 2.0 * RECIP( 2*M + 1 )
   80       CONTINUE

            BI( 0 ) = 1.0

         ELSE IF( EVENL ) THEN

            IDEL = 1

            DO 90 M = LD2, NTRM
               AM( M ) = ( 1. + RECIP( 2*M - L + 1 ) ) * AM( M )
   90       CONTINUE

            DO 100 I = 0, LD2 - 1
               BI( I ) = ( 1. - RECIP( L - 2*I ) ) * BI( I )
  100       CONTINUE

            BI( LD2 ) = ( 2. - RECIP( L ) ) * BI( LD2 - 1 )

         ELSE

            IDEL = 2

            DO 110 M = LD2, NTRM
               AM( M ) = ( 1. - RECIP( 2*M + L + 2 ) ) * AM( M )
  110       CONTINUE

            DO 120 I = 0, LD2
               BI( I ) = ( 1. - RECIP( L + 2*I + 1 ) ) * BI( I )
  120       CONTINUE

         END IF
c                                     ** Establish upper limits for sums
c                                     ** and incorporate factor capital-
c                                     ** del into b-sub-i
         MMAX = NTRM - IDEL
         IF( IPOLZN.GE.0 ) MMAX = MMAX + 1
         IMAX = MIN( LD2, MMAX - LD2 )

         IF( IMAX.LT.0 ) GO TO  250

         DO 130 I = 0, IMAX
            BIDEL( I ) = BI( I )
  130    CONTINUE

         IF( EVENL ) BIDEL( 0 ) = 0.5*BIDEL( 0 )

c                                    ** Perform double sums just for
c                                    ** phase quantities desired by user
         IF( IPOLZN.EQ.0 ) THEN

            DO 150 I = 0, IMAX
c                                           ** vectorizable loop

               SUM = 0.0

               DO 140 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                      ( REAL( C(M-I+1) * CONJG( C(M+I+IDEL) ) )
     &                      + REAL( D(M-I+1) * CONJG( D(M+I+IDEL) ) ) )
  140          CONTINUE

               PMOM( L, 1 ) = PMOM( L, 1 ) + BIDEL( I ) * SUM

  150       CONTINUE

            PMOM( L, 1 ) = 0.5*PMOM( L, 1 )
            GO TO  240

         END IF


         IF( CALCMO( 1 ) ) THEN

            DO 170 I = 0, IMAX

               SUM = 0.0
c                                           ** vectorizable loop
               DO 160 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                        REAL( C(M-I+1) * CONJG( C(M+I+IDEL) ) )
  160          CONTINUE

               PMOM( L, 1 ) = PMOM( L, 1 ) + BIDEL( I ) * SUM

  170       CONTINUE

         END IF


         IF( CALCMO( 2 ) ) THEN

            DO 190 I = 0, IMAX

               SUM = 0.0
c                                           ** vectorizable loop
               DO 180 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                        REAL( D(M-I+1) * CONJG( D(M+I+IDEL) ) )
  180          CONTINUE

               PMOM( L, 2 ) = PMOM( L, 2 ) + BIDEL( I ) * SUM

  190       CONTINUE

         END IF


         IF( CALCMO( 3 ) ) THEN

            DO 210 I = 0, IMAX

               SUM = 0.0
c                                           ** vectorizable loop
               DO 200 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                      ( REAL( C(M-I+1) * CONJG( D(M+I+IDEL) ) )
     &                      + REAL( C(M+I+IDEL) * CONJG( D(M-I+1) ) ) )
  200          CONTINUE

               PMOM( L, 3 ) = PMOM( L, 3 ) + BIDEL( I ) * SUM

  210       CONTINUE

            PMOM( L, 3 ) = 0.5*PMOM( L, 3 )
         END IF


         IF( CALCMO( 4 ) ) THEN

            DO 230 I = 0, IMAX

               SUM= 0.0
c                                           ** vectorizable loop
               DO 220 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                      ( AIMAG( C(M-I+1) * CONJG( D(M+I+IDEL) ) )
     &                      + AIMAG( C(M+I+IDEL) * CONJG( D(M-I+1) ) ) )
  220          CONTINUE

               PMOM( L, 4 ) = PMOM( L, 4 ) + BIDEL( I ) * SUM

  230       CONTINUE

            PMOM( L, 4 ) = - 0.5 * PMOM( L, 4 )

         END IF

  240 CONTINUE


  250 CONTINUE

      RETURN
      END

      SUBROUTINE LPCO1T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities in special case where
c         no. terms in Mie series = 1

c        INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
c                CALCMO                   Flags calculated from IPOLZN
c                A(1), B(1)               Mie series coefficients

c        OUTPUT: PMOM                     Legendre moments


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   IPOLZN, MOMDIM, NMOM
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      REAL      PMOM( 0:MOMDIM, * )
      COMPLEX   A( * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   L, NUMMOM
      REAL      A1SQ, B1SQ
      COMPLEX   A1B1C, CTMP
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CONJG, MIN, REAL
c     ..
c     .. Statement Functions ..

      REAL      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
c     ..


      A1SQ   = SQ( A( 1 ) )
      B1SQ   = SQ( B( 1 ) )
      A1B1C  = A( 1 ) * CONJG( B( 1 ) )


      IF( IPOLZN.LT.0 ) THEN

         IF( CALCMO( 1 ) ) PMOM( 0, 1 ) = 2.25*B1SQ

         IF( CALCMO( 2 ) ) PMOM( 0, 2 ) = 2.25*A1SQ

         IF( CALCMO( 3 ) ) PMOM( 0, 3 ) = 2.25*REAL( A1B1C )

         IF( CALCMO( 4 ) ) PMOM( 0, 4 ) = 2.25*AIMAG( A1B1C )

      ELSE

         NUMMOM = MIN( NMOM, 2 )

c                             ** Loop over moments
         DO 10  L = 0, NUMMOM

            IF( IPOLZN.EQ.0 ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 1.5*( A1SQ + B1SQ )

               IF( L.EQ.1 ) PMOM( L, 1 ) = 1.5*REAL( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = 0.15*( A1SQ + B1SQ )

               GO TO  10

            END IF


            IF( CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 2.25*( A1SQ + B1SQ / 3.)

               IF( L.EQ.1 ) PMOM( L, 1 ) = 1.5*REAL( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = 0.3*B1SQ

            END IF


            IF( CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 2 ) = 2.25*( B1SQ + A1SQ / 3. )

               IF( L.EQ.1 ) PMOM( L, 2 ) = 1.5*REAL( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 2 ) = 0.3*A1SQ

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 3.0*REAL( A1B1C )

               IF( L.EQ.1 ) PMOM( L, 3 ) = 0.75*( A1SQ + B1SQ )

               IF( L.EQ.2 ) PMOM( L, 3 ) = 0.3*REAL( A1B1C )

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = -1.5*AIMAG( A1B1C )

               IF( L.EQ.1 ) PMOM( L, 4 ) = 0.0

               IF( L.EQ.2 ) PMOM( L, 4 ) = 0.3*AIMAG( A1B1C )

            END IF


   10    CONTINUE

      END IF

      RETURN
      END

      SUBROUTINE LPCO2T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities in special case where
c         no. terms in Mie series = 2

c        INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
c                CALCMO                   Flags calculated from IPOLZN
c                A(1-2), B(1-2)           Mie series coefficients

c        OUTPUT: PMOM                     Legendre moments


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   IPOLZN, MOMDIM, NMOM
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      REAL      PMOM( 0:MOMDIM, * )
      COMPLEX   A( * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   L, NUMMOM
      REAL      A2SQ, B2SQ, PM1, PM2
      COMPLEX   A2C, B2C, CA, CAC, CAT, CB, CBC, CBT, CG, CH, CTMP
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CONJG, MIN, REAL
c     ..
c     .. Statement Functions ..

      REAL      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
c     ..


      CA   = 3.*A( 1 ) - 5.*B( 2 )
      CAT  = 3.*B( 1 ) - 5.*A( 2 )
      CAC  = CONJG( CA )
      A2SQ = SQ( A( 2 ) )
      B2SQ = SQ( B( 2 ) )
      A2C  = CONJG( A( 2 ) )
      B2C  = CONJG( B( 2 ) )


      IF( IPOLZN.LT.0 ) THEN

c                                   ** Loop over Sekera moments
         NUMMOM = MIN( NMOM, 2 )

         DO 10 L = 0, NUMMOM

            IF( CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 0.25 * ( SQ( CAT )
     &                                      + (100./3.)* B2SQ )

               IF( L.EQ.1 ) PMOM( L, 1 ) = (5./3.)*REAL( CAT*B2C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = (10./3.)*B2SQ

            END IF


            IF( CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 2 ) = 0.25 * ( SQ( CA )
     &                                      + (100./3.) * A2SQ )

               IF( L.EQ.1 ) PMOM( L, 2 ) = (5./3.)*REAL( CA*A2C )

               IF( L.EQ.2 ) PMOM( L, 2 ) = (10./3.)*A2SQ

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 0.25 * REAL( CAT * CAC
     &                                      + (100./3.) * B(2) * A2C )

               IF( L.EQ.1 ) PMOM( L, 3 ) = 5./6.*
     &                                     REAL( B(2)*CAC + CAT*A2C )

               IF( L.EQ.2 ) PMOM( L, 3 ) = 10./3.* REAL( B(2)*A2C )

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = -0.25 * AIMAG( CAT * CAC
     &                                      + (100./3.)* B(2) * A2C )

               IF( L.EQ.1 ) PMOM( L, 4 ) = -5./ 6.*
     &                                     AIMAG( B(2)*CAC + CAT*A2C )

               IF( L.EQ.2 ) PMOM( L, 4 ) = -10./ 3.* AIMAG( B(2)*A2C )

            END IF


   10    CONTINUE


      ELSE

         CB  = 3.*B( 1 ) + 5.*A( 2 )
         CBT = 3.*A( 1 ) + 5.*B( 2 )
         CBC = CONJG( CB )
         CG  = ( CBC*CBT + 10.*( CAC*A( 2 ) + B2C*CAT ) ) / 3.
         CH  = 2.*( CBC*A( 2 ) + B2C*CBT )

c                               ** Loop over Mueller moments
         NUMMOM = MIN( NMOM, 4 )

         DO 20 L = 0, NUMMOM


            IF( IPOLZN.EQ.0 .OR. CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PM1 = 0.25*SQ( CA ) + SQ( CB ) / 12.
     &                            + (5./3.)*REAL( CA*B2C ) + 5.*B2SQ

               IF( L.EQ.1 ) PM1 = REAL( CB * ( CAC / 6.+ B2C ) )

               IF( L.EQ.2 ) PM1 = SQ( CB ) / 30.+ (20./7.)*B2SQ
     &                            + (2./3.)*REAL( CA*B2C )

               IF( L.EQ.3 ) PM1 = (2./7.) * REAL( CB*B2C )

               IF( L.EQ.4 ) PM1 = (40./63.) * B2SQ

               IF( CALCMO( 1 ) ) PMOM( L, 1 ) = PM1

            END IF


            IF( IPOLZN.EQ.0 .OR. CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PM2 = 0.25*SQ( CAT ) + SQ( CBT ) / 12.
     &                           + ( 5./ 3.) * REAL( CAT*A2C )
     &                           + 5.*A2SQ

               IF( L.EQ.1 ) PM2 = REAL( CBT *
     &                                 ( CONJG( CAT ) / 6.+ A2C ) )

               IF( L.EQ.2 ) PM2 = SQ( CBT ) / 30.
     &                            + ( 20./7.) * A2SQ
     &                            + ( 2./3.) * REAL( CAT*A2C )

               IF( L.EQ.3 ) PM2 = (2./7.) * REAL( CBT*A2C )

               IF( L.EQ.4 ) PM2 = (40./63.) * A2SQ

               IF( CALCMO( 2 ) ) PMOM( L, 2 ) = PM2

            END IF


            IF( IPOLZN.EQ.0 ) THEN

               PMOM( L, 1 ) = 0.5*( PM1 + PM2 )
               GO TO  20

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 0.25 * REAL( CAC*CAT + CG
     &                                         + 20.* B2C * A(2) )

               IF( L.EQ.1 ) PMOM( L, 3 ) = REAL( CAC*CBT + CBC*CAT
     &                                          + 3.*CH ) / 12.

               IF( L.EQ.2 ) PMOM( L, 3 ) = 0.1 * REAL( CG
     &                                      + (200./7.) * B2C * A(2) )

               IF( L.EQ.3 ) PMOM( L, 3 ) = REAL( CH ) / 14.

               IF( L.EQ.4 ) PMOM( L, 3 ) = 40./63.* REAL( B2C*A(2) )

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = 0.25 * AIMAG( CAC*CAT + CG
     &                                      + 20.* B2C * A(2) )

               IF( L.EQ.1 ) PMOM( L, 4 ) = AIMAG( CAC*CBT + CBC*CAT
     &                                           + 3.*CH ) / 12.

               IF( L.EQ.2 ) PMOM( L, 4 ) = 0.1 * AIMAG( CG
     &                                     + (200./7.) * B2C * A(2) )

               IF( L.EQ.3 ) PMOM( L, 4 ) = AIMAG( CH ) / 14.

               IF( L.EQ.4 ) PMOM( L, 4 ) = 40./63.* AIMAG( B2C*A(2) )

            END IF


   20    CONTINUE

      END IF

      RETURN
      END

      SUBROUTINE BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA )

c        Calculate logarithmic derivatives of J-Bessel-function

c     Input :  CIOR, XX, NTRM, NOABS, YESANG  (defined in MIEV0)

c    Output :  RBIGA or CBIGA  (defined in MIEV0)

c    Routines called :  CONFRA


c    INTERNAL VARIABLES :

c       CONFRA     Value of Lentz continued fraction for CBIGA(NTRM),
c                     used to initialize downward recurrence

c       DOWN       = True, use down-recurrence.  False, do not.

c       F1,F2,F3   Arithmetic statement functions used in determining
c                     whether to use up-  or down-recurrence
c                     ( Ref. 2, Eqs. 6-8 )

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
         IF( MIM*XX .LT. F2( MRE ) ) DOWN = .FALSE.

      ELSE

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F1( MRE ) ) DOWN = .FALSE.

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
c                            *** Upward recurrence for BigA
         IF( NOABS ) THEN
c                                  ** No-absorption case; Eq (R20,21)
            RTMP = SIN( MRE*XX )
            RBIGA( 1 ) = - REZINV + RTMP /
     &                   ( RTMP*REZINV - COS( MRE*XX ) )

            DO 30 N = 2, NTRM
               RBIGA( N ) = -( N*REZINV ) +
     &                      1.0 / ( ( N*REZINV ) - RBIGA( N - 1 ) )
   30       CONTINUE

         ELSE
c                                     ** Absorptive case; Eq (R20,22)

            CTMP = EXP( - (0.,2.)*CIOR*XX )
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
     &                   QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &                   SFORW, SBACK, TFORW, TBACK, S1, S2 )

c         Print scattering quantities of a single particle


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   PERFCT
      INTEGER   IPOLZN, MOMDIM, NMOM, NUMANG
      REAL      GQSC, QEXT, QSCA, XX
      COMPLEX   CREFIN, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * ), PRNT( * )
      REAL      PMOM( 0:MOMDIM, * ), XMU( * )
      COMPLEX   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      CHARACTER FMAT*22
      INTEGER   I, J, M
      REAL      FNORM, I1, I2
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
            I1 = REAL( S1( I ) )**2 + AIMAG( S1( I ) )**2
            I2 = REAL( S2( I ) )**2 + AIMAG( S2( I ) )**2
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

         IF( NMOM.GT.0 ) THEN

            WRITE( *, '(/,A)' ) ' Normalized moments of :'

            IF( IPOLZN.EQ.0 ) WRITE( *, '(''+'',27X,A)' )
     &          'Phase Fcn'

            IF( IPOLZN.GT.0 ) WRITE( *, '(''+'',33X,A)' )
     &          'M1           M2          S21          D21'

            IF( IPOLZN.LT.0 ) WRITE( *, '(''+'',33X,A)' )
     &          'R1           R2           R3           R4'

            FNORM = 4./ ( XX**2 * QSCA )

            DO 30  M = 0, NMOM

               WRITE( *, '(A,I4)' ) '      Moment no.', M

               DO 20 J = 1, 4

                  IF( CALCMO( J ) ) THEN
                     WRITE( FMAT, '(A,I2,A)' )
     &                      '( ''+'', T', 24+(J-1)*13, ', 1P,E13.4 )'
                     WRITE( *, FMAT ) FNORM * PMOM( M, J )
                  END IF

   20          CONTINUE
   30       CONTINUE

         END IF

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
      PARAMETER ( TWOTHR = 2./3., FIVTHR = 5./3., FIVNIN = 5./9. )
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
      B( 1 ) = CMPLX( 0., - ( 1. - 0.1*XX**2 ) / 3.) /
     &         CMPLX( 1. + 0.5*XX**2, - XX**3 / 3.)
c                                                       ** Eq. R/A.7,8
      A( 2 ) = CMPLX( 0.,   XX**2 / 30.)
      B( 2 ) = CMPLX( 0., - XX**2 / 45.)
c                                                       ** Eq. R/A.9
      QSCA = 6.* XX**4 *( SQ( A(1) ) + SQ( B(1) ) +
     &           FIVTHR*( SQ( A(2) ) + SQ( B(2) ) ) )
      QEXT = QSCA
c                                                       ** Eq. R/A.10
      GQSC = 6.* XX**4 *REAL( A(1)*CONJG( A(2) + B(1) ) +
     &         ( B(1) + FIVNIN*A(2) )*CONJG( B(2) ) )

      RTMP   = 1.5 * XX**3
      SFORW  = RTMP*( A(1) + B(1) + FIVTHR*( A(2) + B(2) ) )
      SBACK  = RTMP*( A(1) - B(1) - FIVTHR*( A(2) - B(2) ) )
      TFORW( 1 ) = RTMP*( B(1) + FIVTHR*( 2.*B(2) - A(2) ) )
      TFORW( 2 ) = RTMP*( A(1) + FIVTHR*( 2.*A(2) - B(2) ) )
      TBACK( 1 ) = RTMP*( B(1) - FIVTHR*( 2.*B(2) + A(2) ) )
      TBACK( 2 ) = RTMP*( A(1) - FIVTHR*( 2.*A(2) + B(2) ) )


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
      PARAMETER ( TWOTHR = 2./3., FIVTHR = 5./3.)
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
      A( 1 ) = CTMP*( 1.- 0.1*XX**2 +
     &         ( CIORSQ / 350. + 1./280.)*XX**4 ) /
     &         ( CIORSQ + 2.+ ( 1.- 0.7*CIORSQ )*XX**2 -
     &         ( CIORSQ**2 / 175.- 0.275*CIORSQ + 0.25 )*XX**4 +
     &         XX**3 * CTMP * ( 1.- 0.1*XX**2 ) )

c                                           ** Eq. R42b
      B( 1 ) = ( XX**2 / 30. )*CTMP*( 1.+
     &         ( CIORSQ / 35. - 1./ 14.)*XX**2 ) /
     &         ( 1.- ( CIORSQ / 15. - 1./6.)*XX**2 )

c                                           ** Eq. R42c

      A( 2 ) = ( 0.1*XX**2 )*CTMP*( 1.- XX**2 / 14. ) /
     &         ( 2.*CIORSQ + 3.- ( CIORSQ / 7.- 0.5 ) * XX**2 )

c                                           ** Eq. R40a

      QSCA = (6.*XX**4) * ( SQ( A(1) ) + SQ( B(1) ) +
     &                     FIVTHR * SQ( A(2) ) )

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
     &                    FIVTHR*A(2)*( 2.*XMU( J )**2 - 1.) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = ( 0., 0.)

      RETURN
      END

      SUBROUTINE TESTMI( COMPAR, XX, CREFIN, MIMCUT, PERFCT, ANYANG,
     &                   NMOM, IPOLZN, NUMANG, XMU, QEXT, QSCA, GQSC,
     &                   SFORW, SBACK, S1, S2, TFORW, TBACK, PMOM,
     &                   MOMDIM )

c         Set up to run test case when  COMPAR = False;  when  = True,
c         compare Mie code test case results with correct answers
c         and abort if even one result is inaccurate.

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
      INTEGER   IPOLZN, MOMDIM, NMOM, NUMANG
      REAL      GQSC, MIMCUT, QEXT, QSCA, XX
      COMPLEX   CREFIN, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL      PMOM( 0:MOMDIM, * ), XMU( * )
      COMPLEX   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   ANYSAV, OK, PERSAV
      INTEGER   IPOSAV, M, N, NMOSAV, NUMSAV
      REAL      ACCUR, CALC, EXACT, MIMSAV, TESTGQ, TESTQE, TESTQS,
     &          XMUSAV, XXSAV
      COMPLEX   CRESAV, TESTS1, TESTS2, TESTSB, TESTSF
c     ..
c     .. Local Arrays ..

      LOGICAL   CALCMO( 4 ), PRNT( 2 )
      REAL      TESTPM( 0:1 )
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
      SAVE      XXSAV, CRESAV, MIMSAV, PERSAV, ANYSAV, NMOSAV, IPOSAV,
     &          NUMSAV, XMUSAV

      DATA      TESTQE / 2.459791 /,
     &          TESTQS / 1.235144 /,
     &          TESTGQ / 1.139235 /,
     &          TESTSF / ( 61.49476, -3.177994 ) /,
     &          TESTSB / ( 1.493434, 0.2963657 ) /,
     &          TESTS1 / ( -0.1548380, -1.128972 ) /,
     &          TESTS2 / ( 0.05669755, 0.5425681 ) /,
     &          TESTTF / ( 12.95238, -136.6436 ),
     &                   ( 48.54238, 133.4656 ) /,
     &          TESTTB / ( 41.88414, -15.57833 ),
     &                   ( 43.37758, -15.28196 ) /,
     &          TESTPM / 227.1975, 183.6898 /

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
         NMOSAV = NMOM
         IPOSAV = IPOLZN
         NUMSAV = NUMANG
         XMUSAV = XMU( 1 )
c                                   ** Reset input values for test case
         XX     = 10.0
         CREFIN = ( 1.5, -0.1 )
         MIMCUT = 0.0
         PERFCT = .FALSE.
         ANYANG = .TRUE.
         NMOM   = 1
         IPOLZN = -1
         NUMANG = 1
         XMU( 1 ) = -0.7660444

      ELSE
c                                    ** Compare test case results with
c                                    ** correct answers and abort if bad
         OK = .TRUE.

         IF( WRONG( QEXT,TESTQE ) )
     &       OK = TSTBAD( 'QEXT', ABS( ( QEXT - TESTQE ) / TESTQE ) )

         IF( WRONG( QSCA,TESTQS ) )
     &       OK = TSTBAD( 'QSCA', ABS( ( QSCA - TESTQS ) / TESTQS ) )

         IF( WRONG( GQSC,TESTGQ ) )
     &       OK = TSTBAD( 'GQSC', ABS( ( GQSC - TESTGQ ) / TESTGQ ) )

         IF( WRONG( REAL( SFORW ),REAL( TESTSF ) ) .OR.
     &       WRONG( AIMAG( SFORW ),AIMAG( TESTSF ) ) )
     &       OK = TSTBAD( 'SFORW', ABS( ( SFORW - TESTSF ) / TESTSF ) )

         IF( WRONG( REAL( SBACK ),REAL( TESTSB ) ) .OR.
     &       WRONG( AIMAG( SBACK ),AIMAG( TESTSB ) ) )
     &       OK = TSTBAD( 'SBACK', ABS( ( SBACK - TESTSB ) / TESTSB ) )

         IF( WRONG( REAL( S1(1) ),REAL( TESTS1 ) ) .OR.
     &       WRONG( AIMAG( S1(1) ),AIMAG( TESTS1 ) ) )
     &       OK = TSTBAD( 'S1', ABS( ( S1(1) - TESTS1 ) / TESTS1 ) )

         IF( WRONG( REAL( S2(1) ),REAL( TESTS2 ) ) .OR.
     &       WRONG( AIMAG( S2(1) ),AIMAG( TESTS2 ) ) )
     &       OK = TSTBAD( 'S2', ABS( ( S2(1) - TESTS2 ) / TESTS2 ) )


         DO 10  N = 1, 2

            IF( WRONG( REAL( TFORW(N) ),REAL( TESTTF(N) ) ) .OR.
     &          WRONG( AIMAG( TFORW(N) ),
     &          AIMAG( TESTTF(N) ) ) ) OK = TSTBAD( 'TFORW',
     &          ABS( ( TFORW(N) - TESTTF(N) ) / TESTTF(N) ) )

            IF( WRONG( REAL( TBACK(N) ),REAL( TESTTB(N) ) ) .OR.
     &          WRONG( AIMAG( TBACK(N) ),
     &          AIMAG( TESTTB(N) ) ) ) OK = TSTBAD( 'TBACK',
     &          ABS( ( TBACK(N) - TESTTB(N) ) / TESTTB(N) ) )

   10    CONTINUE


         DO 20 M = 0, 1

            IF ( WRONG( PMOM(M,1), TESTPM(M) ) )
     &           OK =  TSTBAD( 'PMOM', ABS( (PMOM(M,1)-TESTPM(M)) /
     &                                      TESTPM(M) ) )

   20    CONTINUE


         IF( .NOT.OK ) THEN

            PRNT( 1 ) = .TRUE.
            PRNT( 2 ) = .TRUE.
            CALCMO( 1 ) = .TRUE.
            CALCMO( 2 ) = .FALSE.
            CALCMO( 3 ) = .FALSE.
            CALCMO( 4 ) = .FALSE.

            CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &                   SFORW, SBACK, TFORW, TBACK, S1, S2 )

            CALL ERRMSG( 'MIEV0 -- Self-test failed',.TRUE.)

         END IF
c                                       ** Restore user input values
         XX     = XXSAV
         CREFIN = CRESAV
         MIMCUT = MIMSAV
         PERFCT = PERSAV
         ANYANG = ANYSAV
         NMOM   = NMOSAV
         IPOLZN = IPOSAV
         NUMANG = NUMSAV
         XMU( 1 ) = XMUSAV

      END IF

      RETURN
      END

