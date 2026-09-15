      SUBROUTINE LINEQSYS(H,V,N,IDIM,OK)
C 
C     The subroutine solves a system of N linear algebraic equations
C     with N unknowns and at the same time inverts N*N matrix H.
C 
C     SUM H(I,J)*Vout(J)=Vinp(I)
C     Hout=THE INVERSE OF Hin
C     GOOD FOR POSITIVE DEFINITE H ONLY
C 
C     Originally written in FORTRAN 66 by Dmitri Toptygin in ~ 1980.
C 
C     In ~ 1990 upgraded from REAL(4) to REAL(8) arithmetic.
C 
C     In ~ 2018 slightly modified to avoid compilation-time warnings
C     issued by modern fortran compilers like gfortran, specifically,
C     "Shared DO termination label".
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
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
      DIMENSION H(IDIM,IDIM),V(IDIM)
      LOGICAL OK
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      OK=.FALSE.
      DO 100,I=1,N
      X=H(I,I)
      H(I,I)=ONE
      IF(X.LE.ZERO)RETURN
      DO 101,J=1,N
      H(I,J)=H(I,J)/X
101   CONTINUE
      V(I)=V(I)/X
      DO 110,J=1,N
      IF(J.NE.I)THEN
      X=H(J,I)
      H(J,I)=ZERO
      DO 111,K=1,N
      H(J,K)=H(J,K)-H(I,K)*X
111   CONTINUE
      V(J)=V(J)-V(I)*X
      ENDIF
110   CONTINUE
100   CONTINUE
      OK=.TRUE.
      RETURN
      END

