C     The three functions in this library calculate the following:
C        1. Chi squared Distribution Function,
C        2. Cumulative Chi squared Distribution Function,
C        3. Inverse Cumulative Chi squared Distribution Function.
C     Function 3 depends on functions 2 and 1, therefore all three functions
C     are combined in one source file.
C 
C     Originally written by Dmitri Toptygin in 2020.
C 
C     ------------------------------------------------------------------
C     Copyright 2022 Dmitri Toptygin
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
      FUNCTION CHI2_DF(X,NDF)
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
C 
C     Chi squared Distribution Function
C 
C     NDF = Number of Degrees of Freedom
C 
C     CHI2_DF = the density of probability that chi^2 equals X
C 
      IF(X.GT.0.0D0.AND.NDF.GE.1)THEN
      HNDF=0.5D0*DBLE(NDF)
      HX=0.5D0*X
      A=(HNDF-1.0D0)*LOG(HX)-HX-LOG_GAMMA(HNDF)
      CHI2_DF=0.5D0*EXP(A)
      ELSE
      CHI2_DF=0.0D0
      ENDIF
      RETURN
      END
C     ------------------------------------------------------------------
      FUNCTION CHI2_CDF(X,NDF)
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
C 
C     Chi squared Cumulative Distribution Function
C 
C     NDF = Number of Degrees of Freedom
C 
C     CHI2_CDF = the probability that chi^2 is less than X
C 
      IF(X.GT.0.0D0.AND.NDF.GE.1)THEN
      HNDF=0.5D0*DBLE(NDF)
      HX=0.5D0*X
      M=MAX(0,NINT(HX-HNDF))
      SUMA=0.0D0
      IF(M.GT.0)THEN
      R=1.0D0
      DO 10,I=M,1,-1
      R=R*(DBLE(I)+HNDF)/HX
      SOLD=SUMA
      SUMA=SUMA+R
      IF(SUMA.EQ.SOLD)GOTO 11
  10  CONTINUE
  11  CONTINUE
      ENDIF
      SUMB=0.0D0
      R=1.0D0
      J=M
  20  J=J+1
      R=R*HX/(DBLE(J)+HNDF)
      SOLD=SUMB
      SUMB=SUMB+R
      IF(SUMB.NE.SOLD)GOTO 20
      S=SUMA+1.0D0+SUMB
      U=DBLE(M)+HNDF
      A=U*LOG(HX)-HX-LOG_GAMMA(U+1.0D0)
      CHI2_CDF=MIN(S*EXP(A),1.0D0)
      ELSEIF(X.GT.0.0D0)THEN
      CHI2_CDF=1.0D0
      ELSE
      CHI2_CDF=0.0D0
      ENDIF
      RETURN
      END
C     ------------------------------------------------------------------
      FUNCTION CHI2_ICDF(P,NDF)
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
C 
C     Chi squared Inverse Cumulative Distribution Function
C 
C     NDF = Number of Degrees of Freedom
C 
C     P is the probability that chi^2 is less than CHI2_ICDF
C 
      IF(P.GT.0.0D0.AND.P.LT.1.0D0.AND.NDF.GE.1)THEN
      ACCUR=1.0D-15*DBLE(NDF+100)
      X=MAX(1.0D-6,DBLE(NDF-2))
  10  DELTA=P-CHI2_CDF(X,NDF)
      X=X+DELTA/CHI2_DF(X,NDF)
      IF(ABS(DELTA).GT.ACCUR)GOTO 10
      CHI2_ICDF=X
      ELSEIF(NDF.GE.1)THEN
        IF(P.LE.0.0D0)THEN
          CHI2_ICDF=0.0D0
        ELSE
          CHI2_ICDF=1.0D38
        ENDIF
      ELSE
      CHI2_ICDF=0.0D0
      ENDIF
      RETURN
      END

