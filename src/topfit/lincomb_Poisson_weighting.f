      SUBROUTINE LINCOMB_POISSON_WEIGHTING(BASF,MF,NF,DATA,WGHT,ND,OK,
     &                                                           DEBUG)
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
      DIMENSION BASF(MF,ND), DATA(ND), WGHT(ND)
      LOGICAL OK, DEBUG, OKI
C 
C     This subroutine fits DATA(n), n=1,...,ND, using the following linear
C     combination of basic functions
C 
C                 NF
C     MODEL(n) = SUM BASF(i,n)*P(i)
C                i=1
C 
C     where P(i), are the fitting parameters, and the array BASF(i,n)
C     contains the values of NF linear independent BASis Functions,
C     calculated for ND values of (an) independent variable(s).
C     The first dimension MF of the array BASF(MF,ND) must be equal to or
C     greater than NF.
C 
C     The MODEL(n) is linear in all the fitting parameters P(i).
C 
C     Assuming that DATA(n) obey Poissonian statistics, the weights are
C     calculated as 
C 
C     WGHT(n) = 1/MODEL(n),
C 
C     which results in a weighted least squares problem with variable
C    (self-consistent) weights.
C 
C     The sum that must be minimized is defined as
C 
C          ND  [DATA(n)-MODEL(n)]^2
C     S = SUM ----------------------
C         n=1       MODEL(n)
C 
C     The condition of the minimum: the first derivatives must equal zero:
C 
C      dS      ND {     [DATA(n)]^2  }
C     ----- = SUM {1 - --------------} * BASF(i,n) = 0
C     dP(i)   n=1 {     [MODEL(n)]^2 }
C 
C     Newton–Raphson method is used to find the P(i) values that result in
C     zero first derivatives; this method requires the calculation of the
C     second derivatives:
C
C        d^2S         ND  [DATA(n)]^2
C     ---------- = 2*SUM -------------- * BASF(i,n) * BASF(j,n)
C     dP(i)dP(j)     n=1  [MODEL(n)]^3
C_______________________________________________________________________
C 
C                   POISSON ARTIFACT CORRECTION
C 
C     Poisson distribution
C 
C               exp(-a) * a^n
C     p(n|a) = ---------------
C                     n! 
C 
C     inf
C     SUM p(n|a) = 1
C     n=0
C 
C     inf
C     SUM n * p(n|a) = a
C     n=0
C 
C     inf
C     SUM n^2 * p(n|a) = a^2+a
C     n=0
C 
C            inf  (n-f)^2              a^2+a
C     S(f) = SUM --------- * p(n|a) = ------- - 2*a + f
C            n=0     f                   f
C 
C     dS(f)      a^2+a
C     ----- = - ------- + 1 = 0    ==>  a^2+a = f^2
C      df         f^2
C 
C     Solution:  a = sqrt(f^2-0.25)-0.5;
C 
C     WGHT(n) = 1/(SQRT(MODEL(n)**2+0.25)-0.5)
C 
C     ------------------------------------------------------------------
C     Copyright 2025 Dmitri Toptygin
C 
C     Permission is hereby granted, free of charge, to any person obtaining a
C     copy of this software and associated documentation files (the "Software"),
C     to deal in the Software without restriction, including without limitation
C     the rights to use, copy, modify, merge, publish, distribute, sublicense,
C     and/or sell copies of the Software, and to permit persons to whom the
C     Software is furnished to do so, subject to the following conditions:
C 
C     The above copyright notice and this permission notice shall be included
C     in all copies or substantial portions of the Software.
C 
C     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
C     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
C     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
C     THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
C     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
C     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
C     DEALINGS IN THE SOFTWARE.
C     ------------------------------------------------------------------
      DIMENSION FM(:), P(:), Q(:), H(:,:), V(:), U(:)
      ALLOCATABLE FM, P, Q, H, V, U
      PARAMETER(SHIFT1=0.0625D0, THRESH1=0.03125D0, NABORT1=256)
      PARAMETER(ACC2=1.0D-10, ACCF=1.074607828D0, NABORT2=256)
C 
      OK=.FALSE.
      IF(NF.GT.MF)GOTO 420
      NNZ=0
      S=0.0D0
      DO 10,N=1,ND
      D=DATA(N)
      IF(D.GT.0.0D0)THEN
        NNZ=NNZ+1
        WGHT(N)=1.0D0/D
        S=S+D
      ELSE
        WGHT(N)=0.0D0
      ENDIF
  10  CONTINUE
      IF(NF.GE.NNZ)THEN
        IF(DEBUG)WRITE(*,'("NF=",I4," >/= NNZ=",I4)')NF,NNZ
        GOTO 420
      ENDIF
C 
      ALLOCATE(FM(ND), P(NF), Q(NF), H(NF,NF), V(NF), U(NF), STAT=IERR)
      IF(IERR.NE.0)GOTO 601
      RNDF=DBLE(ND-NF)
C 
C-------------------- PART 1: Linear least squares ---------------------
C 
      ITER=0
      SHIFT=SHIFT1
      NMIN=0
      YMIN=0.0D0
 101  CONTINUE
      IF(DEBUG)WRITE(*,102)ITER,SHIFT,S/RNDF,NMIN,YMIN
 102  FORMAT('LSE Iter =',I3.2,', Shift =',F8.4,', χ² =',F9.3,
     &       ', N_min =',I6,', F_min=',F8.4)
C 
      DO 110,I=1,NF
      P(I)=0.0D0
      DO 111,J=1,I
      H(J,I)=0.0D0
 111  CONTINUE
 110  CONTINUE
      DO 100,N=1,ND
      D=MAX(DATA(N),SHIFT)
      W=1.0D0/D
      DO 120,I=1,NF
      WB=W*BASF(I,N)
      P(I)=P(I)+WB*D
      DO 121,J=1,I
      H(J,I)=H(J,I)+WB*BASF(J,N)
 121  CONTINUE
 120  CONTINUE
 100  CONTINUE
      DO 130,I=1,NF
      DO 131,J=1,I-1
      H(I,J)=H(J,I)
 131  CONTINUE
 130  CONTINUE
C 
      CALL LINEQSYS(H,P,NF,NF,OKI)
      IF(.NOT.OKI)THEN
        IF(DEBUG)WRITE(*,'("Singular Hessian Matrix 1")')
        GOTO 410
      ENDIF
C 
      YMIN=1.0D308
      S=0.0D0
      DO 160,N=1,ND
      Y=0.0D0
      DO 161,I=1,NF
      Y=Y+BASF(I,N)*P(I)
 161  CONTINUE
      FM(N)=Y
      IF(Y.LT.YMIN)THEN
        YMIN=Y
        NMIN=N
      ENDIF
      D=MAX(DATA(N),SHIFT)
      Z=D-Y
      S=S+Z*Z/D
 160  CONTINUE
      IF(S-S.NE.0.0D0)THEN
        IF(DEBUG)WRITE(*,'("NaN #1")')
        GOTO 410
      ENDIF
      IF(YMIN.LE.THRESH1)THEN
        ITER=ITER+1
        IF(ITER.GT.NABORT1)THEN
          IF(DEBUG)WRITE(*,'("Iter > NABORT1")')
          GOTO 410
        ENDIF
        SHIFT=SHIFT-YMIN+SHIFT1
        GOTO 101
      ENDIF
C 
C-------------------- PART 2: Newton–Raphson method --------------------
C 
      ITER=0
      ACC=ACC2
      S=1.0D+308
 201  CONTINUE
      ITER=ITER+1
      S0=S
C 
      DO 210,I=1,NF
      V(I)=0.0D0
      DO 211,J=1,I
      H(J,I)=0.0D0
 211  CONTINUE
 210  CONTINUE
      DO 200,N=1,ND
      D=DATA(N)
      IF(D.GT.0.0D0)THEN
        F=FM(N)
        R=D/F
        R2=R*R
        Y=R2-1.0D0
        Z=2.0D0*R2/F
        DO 220,I=1,NF
        BNI=BASF(I,N)
        V(I)=V(I)+Y*BNI
        ZBNI=Z*BNI
        DO 221,J=1,I
        H(J,I)=H(J,I)+ZBNI*BASF(J,N)
 221    CONTINUE
 220    CONTINUE
      ELSE
        DO 230,I=1,NF
        V(I)=V(I)-BASF(I,N)
 230    CONTINUE
      ENDIF
 200  CONTINUE
      DO 240,I=1,NF
      U(I)=V(I)
      DO 241,J=1,I-1
      H(I,J)=H(J,I)
 241  CONTINUE
 240  CONTINUE
C 
      CALL LINEQSYS(H,V,NF,NF,OKI)
      IF(.NOT.OKI)THEN
        IF(DEBUG)WRITE(*,'("Singular Hessian Matrix 2")')
        GOTO 410
      ENDIF
      UV=0.0D0
      DO 245,I=1,NF
      UV=UV+U(I)*V(I)
 245  CONTINUE
      UV=0.5D0*UV
C 
      ALAMBD=1.0D0
 250  CONTINUE
      DO 251,I=1,NF
      Q(I)=P(I)+ALAMBD*V(I)
 251  CONTINUE
C 
      YMIN=1.0D308
      S=0.0D0
      DO 260,N=1,ND
      Y=0.0D0
      DO 261,I=1,NF
      Y=Y+BASF(I,N)*Q(I)
 261  CONTINUE
      IF(Y.LE.0.0D0)THEN
        ALAMBD=0.5D0*ALAMBD
        IF(ALAMBD.LT.1.0D-6)THEN
          IF(DEBUG)WRITE(*,'("Lambda < 1.0D-6  #2")')
          GOTO 410
        ENDIF
        GOTO 250
      ENDIF
      FM(N)=Y
      Z=DATA(N)-Y
      S=S+Z*Z/Y
 260  CONTINUE
C 
      IF(S-S.NE.0.0D0)THEN
        IF(DEBUG)WRITE(*,'("NaN #2")')
        GOTO 410
      ENDIF
      IF(S.GT.S0*1.00000001D0)THEN
        ALAMBD=0.5D0*ALAMBD
        IF(ALAMBD.LT.1.0D-6)THEN
          IF(DEBUG)WRITE(*,'("Lambda < 1.0D-6")')
          GOTO 410
        ENDIF
        GOTO 250
      ENDIF
C 
      IF(DEBUG)WRITE(*,271)ITER,S/RNDF,UV/RNDF,ALAMBD
 271  FORMAT('Newton Iter =',I3.2,', χ² =',F13.6,' Expect Δ =',F13.6,
     &' λ =',F9.6)
C 
      DO 270,I=1,NF
      P(I)=Q(I)
 270  CONTINUE
C 
      THR=ACC*RNDF
      IF(S0-S.LE.THR.AND.UV.LE.THR.AND.ALAMBD.GE.1.0D0)THEN
        CONTINUE
      ELSE
        IF(ITER.GE.NABORT2)THEN
          IF(DEBUG)WRITE(*,'("Iter > NABORT2")')
          GOTO 410
        ENDIF
        ACC=ACC*ACCF
        GOTO 201
      ENDIF
C 
C------ PART 3: Weight calculation & Poisson artifact correction -------
C 
      DO 310,N=1,ND
      F=FM(N)
      WGHT(N)=1.0D0/(SQRT(F*F+0.25D0)-0.5D0)
 310  CONTINUE
      OK=.TRUE.
C 
 410  DEALLOCATE(FM, P, Q, H, V, U, STAT=IERR)
      IF(IERR.NE.0)GOTO 601
 420  RETURN
C 
 601  STOP 'Memory allocation error in SUBROUTINE LINCOMB_POISSON_WEIGHT
     &ING.  '
      END

