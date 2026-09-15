      SUBROUTINE TOPFIT(MODEL,PAR,STD,FLO,NP,MAXLTRAIN,GLFIRST,DOMAIN,
     &ARGM,DATA,PRED,WGHT,ND,NABORT,ACCUR,ALPHA0,CHI2,OK,SWRITE,DWRITE,
     &SAVEINVHESSM)
C     ------------------------------------------------------------------
C 
C     TOPFIT() is a Weighted Nonlinear Least Squares Estimator written
C     by Dmitri Toptygin. TOP are the first letters ofthe author's name,
C     FIT are the first letters of the word fitting like in data fitting
C     by a nonlinear model function of multiple unknown parameters.
C 
C     The first version of subroutines TOPFIT() was written in
C     FORTRAN 66 by Dmitri Toptygin in 1984. It used REAL(4) arithmetic,
C     fixed-dimension arrays, and was based on Gauss-Newton algorithm.
C 
C     Different versions of subroutine TOPFIT() have been written. The
C     difference is in the algorithm (Gauss-Newton or Marquardt or
C     Hessian matrix eigendecomposition algorithm), whether the internal
C     arrays are fixed size or allocatable, whether the subroutine is
C     optimized for use in global analysis, etc.
C     
C     Subroutines TOPFIT() have identical lists of formal arguments,
C     which makes them interchangeable. Different subroutines TOPFIT()
C     are saved under different file names. 
C 
C     "topfit_marquardt.f" is a REAL(8) subroutine with Allocatable
C     arrays and a Global Accaleration feature. It is based on a
C     modification of the Marquardt algorithm.
C
C     To take advantage of the Global Acceleration, all global parameters
C     must either preceed local parameters (GLFIRST=.TRUE.) or succeed
C     local parameters (GLFIRST=.FALSE.). If global and local parameters
C     are mixed, then there will be no benefit in computation speed.
C     If the program is not used for global analysis, then it makes no
C     difference whether GLFIRST is .TRUE. or .FALSE..
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
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
      EXTERNAL MODEL, DOMAIN
      LOGICAL FLO, GLFIRST, DOMAIN, OK, SWRITE, DWRITE, SAVEINVHESSM
      DIMENSION PAR(NP),STD(NP),FLO(NP)
      DIMENSION ARGM(*)
      DIMENSION DATA(ND),PRED(ND),WGHT(ND)
      ALLOCATABLE TRY(:),DELTA(:),DCOPY(:),DELTA0(:),DERIV(:,:)
      ALLOCATABLE HESSIAN(:,:),HSCOPY(:,:)
      LOGICAL OKEEP, OCM, ONDOM, ONAN
C 
      PARAMETER(ALPHAMAX=1.0D0,ALPHAMIN=1.776D-15)
      PARAMETER(ALPHAFINCR=4.0D0,ALPHAFDECR=0.250D0)
      PARAMETER(NFAILMAX=24)
C 
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,VRYLRG=1.78D+308)
 4    FORMAT(   /,   A,' CHI SQUARE MINIMIZATION ABORTED.')
 5    FORMAT(   /,   A,' THE DOMAIN OF TOLERATED VALUES IN PARAMETER SPA
     &CE. CHI SQUARE MINIMIZATION ABORTED.')
C6    FORMAT(   /,   A,' IS CURRENTLY SET TO ',I4)
C 
      CHI2=ZERO
      IF(NP.LE.0)THEN
        WRITE(*,4)'ZERO OR NEGATIVE NUMBER OF PARAMETERS.'
        GOTO 801
      ENDIF
C 
      ONDOM=.NOT.DOMAIN(PAR,NP)
      IF(ONDOM)THEN
        WRITE(*,5)'INITIAL GUESSES OUTSIDE'
        GOTO 801
      ENDIF
C 
      IF(MAXLTRAIN.LE.0)THEN
        WRITE(*,4)'ZERO OR NEGATIVE MAXLTRAIN.'
        GOTO 801
      ENDIF
C 
      IF(ND.LE.0)THEN
        WRITE(*,4)'ZERO OR NEGATIVE NUMBER OF DATA POINTS.'
        GOTO 801
      ENDIF
C 
      ALPHA=ALPHA0
C 
      IF(GLFIRST)THEN
      LIM1=NP
      LIM2=1
      INCR=-1
      ELSE
      LIM1=1
      LIM2=NP
      INCR=+1
      ENDIF
C 
      MAXNP=NP
      ALLOCATE(TRY(MAXNP),DELTA(MAXNP),DCOPY(MAXNP),DELTA0(MAXNP),
     &DERIV(MAXNP,MAXLTRAIN),
     &HESSIAN(MAXNP,MAXNP),HSCOPY(MAXNP,MAXNP),STAT=IERR)
      IF(IERR.NE.0)THEN
        WRITE(*,4)'INSUFFICIENT MEMORY FOR HESSIAN AND OTHER MATRIXES.'
        GOTO 801
      ENDIF
C 
      OK=.TRUE.
      IF(SWRITE)WRITE(*,51)
      IF(DWRITE)THEN
        OPEN(792,FILE='chi2list.txt',STATUS='UNKNOWN')
        WRITE(792,51)
      ENDIF
 51   FORMAT('Iteration',10X,'Chi square ',2X,' Lambda ',2X,'  Alpha  ')
 52   FORMAT(   I7.4,2X ,           F21.9 ,2X,   F8.6   ,2X,   1PE9.3  )
 53   FORMAT(   I7.4,2X ,      1PD21.14   ,2X, 0PF8.6   ,2X,   1PE9.3  )
C 
      OCM=.FALSE.
      DO 90,I=1,NP
      DELTA(I)=ZERO
 90   CONTINUE
      ALAMBD=ZERO
      CHIOLD=VRYLRG
      CHIMIN=VRYLRG
      CHIDEC=ZERO
      ITER=-1
      NFAIL=0
 100  ITER=ITER+1
      IF(ITER.GT.NABORT)THEN
        WRITE(*,101)NABORT
 101    FORMAT(   /,' AFTER',I6,' ITERATIONS THE MINIMUM IS NOT FOUND. C
     &HI SQUARE MINIMIZATION ABORTED.')
        OK=.FALSE.
        GOTO 500
      ENDIF
 110  DO 111,I=1,NP
      TRY(I)=PAR(I)+ALAMBD*DELTA(I)
 111  CONTINUE
      ONDOM=.NOT.DOMAIN(TRY,NP)
      IF(ONDOM)THEN
        ITER=ITER-1
        ONAN=.FALSE.
        GOTO 320
      ENDIF
      CHISUM=ZERO
      NCHI=0
      DO 120,I=1,NP
      DELTA(I)=ZERO
      DO 121,J=1,I
      HESSIAN(J,I)=ZERO
 121  CONTINUE
 120  CONTINUE
C 
      M=0
      DO 200,N=1,ND
      IF(N.GT.M)THEN
        L=0
        CALL MODEL(TRY,NP,MAXNP,ARGM,N,PRED(N),DERIV,LTRAIN)
        IF(LTRAIN.LE.0.OR.LTRAIN.GT.MAXLTRAIN)THEN
          WRITE(*,4)'LTRAIN VALUE ERROR.'
          GOTO 802
        ENDIF
        M=M+LTRAIN
      ENDIF
      L=L+1
      W=WGHT(N)
      IF(W.GT.ZERO)THEN
        Y=PRED(N)
        ONAN=(Y-Y).NE.ZERO
        IF(ONAN)GOTO 320
        Y=DATA(N)-Y
        CHISUM=CHISUM+W*Y*Y
        NCHI=NCHI+1
        DO 210,I=1,NP
        IF(FLO(I))THEN
          X=DERIV(I,L)
          ONAN=(X-X).NE.ZERO
          IF(ONAN)GOTO 320
          IF(X.NE.ZERO)THEN
            X=W*X
            DELTA(I)=DELTA(I)+X*Y
            DO 211,J=1,I
            IF(FLO(J))THEN
              Z=DERIV(J,L)
              IF(Z.NE.ZERO)HESSIAN(J,I)=HESSIAN(J,I)+X*Z
            ENDIF
 211        CONTINUE
          ENDIF
        ENDIF
 210    CONTINUE
      ELSEIF(W.LT.ZERO)THEN
        WRITE(*,4)'NEGATIVE WEIGHT ENCOUNTERED.'
        GOTO 802
      ENDIF
 200  CONTINUE
      DO 290,I=1,NP
      DO 291,J=1,I-1
      HESSIAN(I,J)=HESSIAN(J,I)
 291  CONTINUE
      IF(FLO(I).AND.HESSIAN(I,I).GT.ZERO)NCHI=NCHI-1
 290  CONTINUE
      IF(NCHI.LT.1)NCHI=1
      RNCHI=DBLE(NCHI)
      CHICUR=CHISUM/RNCHI
C 
      IF(SWRITE)THEN
        IF(CHICUR.LE.99999999999.0D0)THEN
          WRITE(*,52)ITER,CHICUR,ALAMBD,ALPHA
        ELSE
          WRITE(*,53)ITER,CHICUR,ALAMBD,ALPHA
        ENDIF
      ENDIF
      IF(DWRITE)THEN
        WRITE(792,52)ITER,CHICUR,ALAMBD,ALPHA
        OPEN(793,FILE='curr_par.txt',STATUS='UNKNOWN')
        WRITE(793,301)(TRY(I),I=1,NP)
 301    FORMAT(G20.12)
        WRITE(793,302)ITER,CHICUR,ALAMBD,ALPHA
 302    FORMAT('Iteration',I6,'   Chisquare =',F16.6,'   Lambda =',F9.6,
     &'   Alpha =',1PE10.3)
        CLOSE(793)
      ENDIF
C 
      CHIMIN=MIN(CHIMIN,CHISUM)
      TOLDIF=ACCUR*MAX(CHIMIN,RNCHI)
      OKEEP=(CHISUM-CHIMIN).LE.TOLDIF
      IF(OKEEP)THEN
        CHI2=CHICUR
        DO 310,I=1,NP
        PAR(I)=TRY(I)
        DO 311,J=1,NP
        HSCOPY(I,J)=HESSIAN(I,J)
 311    CONTINUE
        DCOPY(I)=DELTA(I)
 310    CONTINUE
        OCM=.TRUE.
        IF(ALAMBD.GE.ONE.AND.NFAIL.EQ.0)
     &  ALPHA=MAX(ALPHAMIN,ALPHA*ALPHAFDECR)
        GOTO 340
      ENDIF
 320  OKEEP=.FALSE.
      NFAIL=NFAIL+1
      IF(NFAIL.GT.NFAILMAX)THEN
        IF(ONDOM)THEN
          WRITE(*,5)'CHI SQUARE GRADIENT LEADS OUTSIDE'
        ELSEIF(ONAN)THEN
          WRITE(*,4)'MODEL FUNCTION RETURNED NaN (NOT A NUMBER).'
        ELSEIF(.NOT.OCM)THEN
          WRITE(*,4)'SINGULAR HESSIAN MATRIX.'
        ELSE
          WRITE(*,4)'CHI SQUARE GRADIENT LEADS TO A REGION IN PARAMETER 
     &SPACE WHERE THE HESSIAN MATRIX IS SINGULAR.'
        ENDIF
        OK=.FALSE.
        GOTO 500
      ENDIF
C     ALPHA=MAX(ALPHA*ALPHAFINCR,ALPHA0)
      ALPHA=ALPHA*ALPHAFINCR
      IF(ALPHA.GT.ALPHAMAX)THEN
        ALAMBD=0.5D0*ALAMBD
        ALPHA=ALPHAMAX
      ENDIF
      DO 330,I=1,NP
      DO 331,J=1,NP
      HESSIAN(I,J)=HSCOPY(I,J)
 331  CONTINUE
      DELTA(I)=DCOPY(I)
      TRY(I)=PAR(I)
 330  CONTINUE
 340  CONTINUE
C 
      DO 350,I=1,NP
      IF(.NOT.FLO(I))THEN
        DO 351,J=1,NP
        HESSIAN(I,J)=ZERO
        HESSIAN(J,I)=ZERO
 351    CONTINUE
        HESSIAN(I,I)=ONE
        DELTA(I)=ZERO
      ENDIF
 350  CONTINUE
C 
      X=ONE+ALPHA
      DO 360,I=1,NP
      IF(FLO(I))THEN
        Y=X*HESSIAN(I,I)
        IF(Y.GT.ZERO)THEN
          HESSIAN(I,I)=Y
        ELSE
          HESSIAN(I,I)=ONE
        ENDIF
      ENDIF
      DELTA0(I)=DELTA(I)
 360  CONTINUE
C 
      DO 400,I=LIM1,LIM2,INCR
      X=ONE/HESSIAN(I,I)
      IF((X-X).NE.ZERO)GOTO 320
      DO 410,J=I+INCR,LIM2,INCR
      HESSIAN(I,J)=X*HESSIAN(I,J)
 410  CONTINUE
      DELTA(I)=X*DELTA(I)
      DO 420,K=LIM1,LIM2,INCR
      IF(K.NE.I)THEN
        X=HESSIAN(K,I)
        IF(X.NE.ZERO)THEN
          DO 421,J=I+INCR,LIM2,INCR
          HESSIAN(K,J)=HESSIAN(K,J)-X*HESSIAN(I,J)
 421      CONTINUE
          DELTA(K)=DELTA(K)-X*DELTA(I)
        ENDIF
      ENDIF
 420  CONTINUE
 400  CONTINUE
      DO 430,I=1,NP
      X=DELTA(I)
      IF((X-X).NE.ZERO)GOTO 320
 430  CONTINUE
C 
      IF(OKEEP)THEN
C 
      CHIDEC=ZERO
      DO 490,I=1,NP
      CHIDEC=CHIDEC+DELTA(I)*DELTA0(I)
 490  CONTINUE
      IF(ABS(CHISUM-CHIOLD).LE.TOLDIF.AND.CHIDEC.LE.TOLDIF.AND.
     &         ALAMBD.GE.ONE.AND.ALPHA.LE.MAX(ALPHA0,ALPHAMIN))GOTO 500
      CHIOLD=CHISUM
      NFAIL=0
      ALAMBD=ONE
C 
      ENDIF
C 
      GOTO 100
C 
 500  IF(DWRITE)CLOSE(792)
      X=ONE+1.776D-15
      DO 510,I=1,NP
      IF(OCM.AND.FLO(I))THEN
        DO 511,J=1,NP
        IF(OCM.AND.FLO(J))THEN
          HESSIAN(I,J)=HSCOPY(I,J)
        ELSE
          HESSIAN(I,J)=ZERO
        ENDIF
 511    CONTINUE
        HESSIAN(I,I)=MAX(X*HSCOPY(I,I),1.0D-200)
      ELSE
        DO 512,J=1,NP
        HESSIAN(I,J)=ZERO
 512    CONTINUE
        HESSIAN(I,I)=ONE
      ENDIF
 510  CONTINUE
C 
      DO 600,I=LIM1,LIM2,INCR
      X=ONE/HESSIAN(I,I)
      HESSIAN(I,I)=ONE
      DO 610,J=LIM1,LIM2,INCR
      HESSIAN(I,J)=X*HESSIAN(I,J)
 610  CONTINUE
      DO 620,K=LIM1,LIM2,INCR
      IF(K.NE.I)THEN
        X=HESSIAN(K,I)
        HESSIAN(K,I)=ZERO
        IF(X.NE.ZERO)THEN
          DO 621,J=LIM1,LIM2,INCR
          HESSIAN(K,J)=HESSIAN(K,J)-X*HESSIAN(I,J)
 621      CONTINUE
        ENDIF
      ENDIF
 620  CONTINUE
 600  CONTINUE
C 
      DO 630,I=LIM1,LIM2,INCR
      IF(.NOT.FLO(I))THEN
        DO 631,J=LIM1,LIM2,INCR
        HESSIAN(J,I)=ZERO
        HESSIAN(I,J)=ZERO
 631    CONTINUE
      ENDIF
 630  CONTINUE
C 
      IF(SAVEINVHESSM)THEN
        OPEN(794,FILE='invhessm.bin',FORM='UNFORMATTED',ACCESS='DIRECT',
     &  RECL=8,STATUS='UNKNOWN',ERR=681)
        IREC=0
        DO 670,J=1,NP
        DO 671,I=1,NP
        IREC=IREC+1
        WRITE(794,REC=IREC,ERR=681)HESSIAN(I,J)
 671    CONTINUE
 670    CONTINUE
 681    CLOSE(794,STATUS='KEEP',ERR=682)
 682    CONTINUE
      ENDIF
C 
      X=MAX(CHI2,ONE)
      DO 690,I=1,NP
      IF(OCM.AND.FLO(I))THEN
        STD(I)=MIN(DSQRT(X*HESSIAN(I,I)),1.0D+99)
      ELSEIF(FLO(I))THEN
        STD(I)=1.0D+99
      ELSE
        STD(I)=ZERO
      ENDIF
 690  CONTINUE
      DEALLOCATE(TRY,DELTA,DCOPY,DELTA0,DERIV,HESSIAN,HSCOPY)
      GOTO 990
C 
 802  IF(DWRITE)CLOSE(792)
      DEALLOCATE(TRY,DELTA,DCOPY,DELTA0,DERIV,HESSIAN,HSCOPY)
 801  OK=.FALSE.
      DO 810,I=1,NP
      IF(FLO(I))THEN
        STD(I)=1.0D+99
      ELSE
        STD(I)=ZERO
      ENDIF
 810  CONTINUE
C 
 990  RETURN
      END

