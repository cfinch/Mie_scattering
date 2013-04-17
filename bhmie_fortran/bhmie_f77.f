      SUBROUTINE BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
C Declare parameters:
C Note: important that MXNANG be consistent with dimension of S1 and S2
C       in calling routine!
      INTEGER MXNANG,NMXX
      PARAMETER(MXNANG=1000,NMXX=15000)
C Arguments:
      INTEGER NANG
      REAL GSCA,QBACK,QEXT,QSCA,X
      COMPLEX REFREL
      COMPLEX S1(2*MXNANG-1),S2(2*MXNANG-1)
C Local variables:
      INTEGER J,JJ,N,NSTOP,NMX,NN
      REAL APSI,APSI1,CHI,CHI0,CHI1,DANG,FN,P,PII,
     &     RN,THETA,XSTOP,YMOD
      REAL AMU(MXNANG),PI(MXNANG),PI0(MXNANG),PI1(MXNANG),TAU(MXNANG)
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
      COMPLEX AN,AN1,BN,BN1,XI,XI1,Y
      COMPLEX D(NMXX)
C***********************************************************************
C Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
C    to calculate scattering and absorption by a homogenous isotropic
C    sphere.
C Given:
C    X = 2*pi*a/lambda
C    REFREL = (complex refr. index of sphere)/(real index of medium)
C    NANG = number of angles between 0 and 90 degrees
C           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
C           if called with NANG<2, will set NANG=2 and will compute
C           scattering for theta=0,90,180.
C Returns:
C    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
C                                scatt. E perp. to scatt. plane)
C    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
C                                scatt. E parr. to scatt. plane)
C    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
C    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
C    QBACK = (dC_sca/domega)/pi*a**2
C          = backscattering efficiency
C    GSCA = <cos(theta)> for scattering
C
C Original program taken from Bohren and Huffman (1983), Appendix A
C Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
C in order to compute <cos(theta)>
C 91/05/07 (BTD): Modified to allow NANG=1
C 91/08/15 (BTD): Corrected error (failure to initialize P)
C 91/08/15 (BTD): Modified to enhance vectorizability.
C 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
C 91/08/15 (BTD): Changed definition of QBACK.
C 92/01/08 (BTD): Note that this version has been superceded by
C                 fully double precision version = bhmie.f which,
C                 unfortunately, is not standard f77.
C                 However, retain this in case standard f77 version
C                 is required for porting to some other system.
C***********************************************************************
C*** Safety checks
      IF(NANG.GT.MXNANG)STOP'***Error: NANG > MXNANG in bhmie'
      IF(NANG.LT.2)NANG=2
C*** Obtain pi:
      PII=4.E0*ATAN(1.E0)
      DX=X
      Y=X*REFREL
      YMOD=ABS(Y)
C
C*** Series expansion terminated after NSTOP terms
C    Logarithmic derivatives calculated from NMX on down
      XSTOP=X+4.E0*X**0.3333+2.0
C*** Original code:
C      NMX=AMAX1(XSTOP,YMOD)+15
C      NSTOP=XSTOP
C*** Experimental code:
      NMX=1.0*AMAX1(XSTOP,YMOD)+15
      NSTOP=1.0*XSTOP
C
      IF(NMX.GT.NMXX)THEN
          WRITE(0,*)'Error: NMX > NMXX=',NMXX,' for |m|x=',YMOD
          STOP
      ENDIF
C*** Require NANG.GE.1 in order to calculate scattering intensities
      DANG=0.
      IF(NANG.GT.1)DANG=.5E0*PII/FLOAT(NANG-1)
      DO 1000 J=1,NANG
          THETA=FLOAT(J-1)*DANG
          AMU(J)=COS(THETA)
 1000 CONTINUE
      DO 1100 J=1,NANG
          PI0(J)=0.E0
          PI1(J)=1.E0
 1100 CONTINUE
      NN=2*NANG-1
      DO 1200 J=1,NN
          S1(J)=(0.E0,0.E0)
          S2(J)=(0.E0,0.E0)
 1200 CONTINUE
C
C*** Logarithmic derivative D(J) calculated by downward recurrence
C    beginning with initial value (0.,0.) at J=NMX
C
      D(NMX)=(0.E0,0.E0)
      NN=NMX-1
      DO 2000 N=1,NN
          RN=NMX-N+1
          D(NMX-N)=(RN/Y)-(1.E0/(D(NMX-N+1)+RN/Y))
 2000 CONTINUE
C
C*** Riccati-Bessel functions with real argument X
C    calculated by upward recurrence
C
      PSI0=DCOS(DX)
      PSI1=DSIN(DX)
      CHI0=-SIN(X)
      CHI1=COS(X)
C APSI0 never used, so this line removed from program:
C      APSI0=PSI0
      APSI1=PSI1
C XI0 never used, so this line removed from program:
C      XI0=CMPLX(APSI0,-CHI0)
      XI1=CMPLX(APSI1,-CHI1)
      QSCA=0.E0
      GSCA=0.E0
      P=-1.
      DO 3000 N=1,NSTOP
          DN=N
          RN=N
          FN=(2.E0*RN+1.E0)/(RN*(RN+1.E0))
          PSI=(2.E0*DN-1.E0)*PSI1/DX-PSI0
          APSI=PSI
          CHI=(2.E0*RN-1.E0)*CHI1/X-CHI0
          XI=CMPLX(APSI,-CHI)
C
C*** Store previous values of AN and BN for use
C    in computation of g=<cos(theta)>
          IF(N.GT.1)THEN
              AN1=AN
              BN1=BN
          ENDIF
C
C*** Compute AN and BN:
          AN=(D(N)/REFREL+RN/X)*APSI-APSI1
          AN=AN/((D(N)/REFREL+RN/X)*XI-XI1)
          BN=(REFREL*D(N)+RN/X)*APSI-APSI1
          BN=BN/((REFREL*D(N)+RN/X)*XI-XI1)
C
C*** Augment sums for Qsca and g=<cos(theta)>
          QSCA=QSCA+(2.*RN+1.)*(CABS(AN)**2+CABS(BN)**2)
          GSCA=GSCA+((2.*RN+1.)/(RN*(RN+1.)))*
     &         (REAL(AN)*REAL(BN)+AIMAG(AN)*AIMAG(BN))
          IF(N.GT.1)THEN
              GSCA=GSCA+((RN-1.)*(RN+1.)/RN)*
     &        (REAL(AN1)*REAL(AN)+AIMAG(AN1)*AIMAG(AN)+
     &         REAL(BN1)*REAL(BN)+AIMAG(BN1)*AIMAG(BN))
          ENDIF
C
C*** Now calculate scattering intensity pattern
C    First do angles from 0 to 90
          DO 2500 J=1,NANG
              JJ=2*NANG-J
              PI(J)=PI1(J)
              TAU(J)=RN*AMU(J)*PI(J)-(RN+1.E0)*PI0(J)
              S1(J)=S1(J)+FN*(AN*PI(J)+BN*TAU(J))
              S2(J)=S2(J)+FN*(AN*TAU(J)+BN*PI(J))
 2500     CONTINUE
C
C*** Now do angles greater than 90 using PI and TAU from
C    angles less than 90.
C    P=1 for N=1,3,...; P=-1 for N=2,4,...
          P=-P
          DO 2600 J=1,NANG-1
              JJ=2*NANG-J
              S1(JJ)=S1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
              S2(JJ)=S2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
 2600     CONTINUE
          PSI0=PSI1
          PSI1=PSI
          APSI1=PSI1
          CHI0=CHI1
          CHI1=CHI
          XI1=CMPLX(APSI1,-CHI1)
C
C*** Compute pi_n for next value of n
C    For each angle J, compute pi_n+1
C    from PI = pi_n , PI0 = pi_n-1
          DO 2800 J=1,NANG
              PI1(J)=((2.*RN+1.)*AMU(J)*PI(J)-(RN+1.)*PI0(J))/RN
              PI0(J)=PI(J)
 2800     CONTINUE
 3000 CONTINUE
C
C*** Have summed sufficient terms.
C    Now compute QSCA,QEXT,QBACK,and GSCA
      GSCA=2.*GSCA/QSCA
      QSCA=(2.E0/(X*X))*QSCA
      QEXT=(4.E0/(X*X))*REAL(S1(1))
      QBACK=CABS(S1(2*NANG-1))*CABS(S1(2*NANG-1))/(PII*X*X)
      RETURN
      END
