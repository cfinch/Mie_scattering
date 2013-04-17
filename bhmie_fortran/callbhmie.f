      PROGRAM CALLBHMIE
C Parameters:
      INTEGER MXNANG
      PARAMETER(MXNANG=1000)
C Variables:
      INTEGER IREADEP,J,NAN,NANG,NANG0
      REAL AJ,ANG,DANG,GSCA,PI,POL,
     &     QABS,QBACK,QEXT,QSCA,RAD,REFMED,
     &     S11,S12,S33,S34,WAVEL,X
      COMPLEX REFREL,CXEPS,S1(2*MXNANG-1),S2(2*MXNANG-1)
C***********************************************************************
C Program to interactively call Bohren-Huffman Mie theory program
C
C CALLBHMIE will interactively prompt for:
C 1. refractive index of surrounding medium
C 2. either refractive index or dielectric constant of sphere
C 3. radius of sphere
C 4. wavelength (in vacuo)
C 5. number of angles at which to calculate scattering intensities
C
C CALLBHMIE will return:
C 1. Q_ext, Q_abs, Q_sca, g, Q_back
C 2. If NANG>0, then will also return scattering matrix elements
C    S_11, S_33, S_34, and POL
C
C Adapted by B.T.Draine, Princeton Univ. Obs.
C***********************************************************************
      PI=4.E0*ATAN(1.E0)
      OPEN(UNIT=7,FILE='callbhmie.out',STATUS='UNKNOWN')
      WRITE(*,*)' Enter (real) refractive index of surrounding medium'
      READ(*,*)REFMED
      WRITE(*,*)' Wish to enter refr.index or epsilon? (0 or 1)'
      READ(*,*)IREADEP
 1000 IF(IREADEP.LE.0)THEN
          WRITE(*,*)' Enter complex refractive index of sphere'
          READ(*,*)REFREL
      ELSE
          WRITE(*,*)' Enter complex epsilon of sphere'
          READ(*,*)CXEPS
          REFREL=SQRT(CXEPS)
      ENDIF
      REFREL=REFREL/REFMED
      WRITE(*,6012)REFREL
      WRITE(7,6012)REFREL
 6012 FORMAT(' Complex refractive index=',1P2E10.3)
      WRITE(*,*)' Enter radius'
      READ(*,*)RAD
      WRITE(*,*)' Enter wavelength'
      READ(*,*)WAVEL
      WRITE(*,*)' Enter NANG = number of angles between 0 and 90'
      READ(*,*)NANG0
      IF(NANG0.GT.MXNANG)STOP'***Error: NANG > MXNANG'
      NANG=NANG0
      IF(NANG0.LT.2)NANG=2
      X=2.E0*PI*RAD*REFMED/WAVEL
      WRITE(7,6013)RAD,WAVEL,X
      WRITE(*,6013)RAD,WAVEL,X
C**********
C NANG=number of angles between 0 and 90 degrees (incl. 0 and 90)
C Scattering matrix elements are calculated for 2*NANG-1 angles
C including 0, 90, and 180 degrees.
C**********
      IF(NANG.GT.1)DANG=0.5E0*PI/FLOAT(NANG-1)
      CALL BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
      QABS=QEXT-QSCA
      WRITE(*,6065)QEXT,QSCA,QABS,GSCA,QBACK
      WRITE(7,6065)QEXT,QSCA,QABS,GSCA,QBACK
C POL=degree of polarization (for incident upolarized light)
      IF(NANG0.GT.1)THEN
          NAN=2*NANG-1
          WRITE(*,6017)
          WRITE(7,6017)
          DO 355 J=1,NAN
          AJ=J
          S11=0.5E0*CABS(S2(J))*CABS(S2(J))
          S11=S11+0.5E0*CABS(S1(J))*CABS(S1(J))
          S12=0.5E0*CABS(S2(J))*CABS(S2(J))
          S12=S12-0.5E0*CABS(S1(J))*CABS(S1(J))
          POL=-S12/S11
          S33=REAL(S2(J)*CONJG(S1(J)))
          S34=AIMAG(S2(J)*CONJG(S1(J)))
          ANG=DANG*(AJ-1.E0)*180.E0/PI
          WRITE(7,6075)ANG,S11,POL,S33,S34
          WRITE(*,6075)ANG,S11,POL,S33,S34
  355     CONTINUE
      ENDIF
 6013 FORMAT(' radius=',1PE11.4,' lambda=',E11.4,' x=',E11.4)
 6017 FORMAT(2X,'theta',7X,'S11',11X,'POL',11X,'S33',11X,'S34')
 6065 FORMAT(/,'Qext=',1PE11.4,' Qsca=',E11.4,' Qabs=',E11.4,
     &' <cos>=',E11.4,/,17X,'Qbk =',E11.4)
 6075 FORMAT(1X,F6.2,2X,1PE12.5,2X,E12.5,2X,E12.5,2X,E12.5)
      STOP
      END
