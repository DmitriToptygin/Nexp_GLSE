C     This subroutine is for use with global analysis programs;
C     it simplifies transitions between global and local parameters.
C 
C     GMAPDEFINE(L,N) defines/redefines the global parameter map dimensions
C     L (number of Local parameters) and N (Number of locations)
C 
C     GMAPWRITE(I,J) stores global numbers I(1),I(2),...,I(L) for the
C     L local parameters corresponding to location J for future use
C 
C     GMAPREAD(I,J) returns global numbers I(1),I(2),...,I(L) for the
C     L local parameters corresponding to location J, where J=1,2,...,N
C 
C     GMAPSHOWDIM(L,N) returns current map dimensions L and N
C 
C     Written by Dmitri Toptygin in 2020.
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
      SUBROUTINE GMAPDEFINE(LOCPAR,NLOC)
      IMPLICIT INTEGER(4) (A-Z)
      DIMENSION I(*), MAP(:,:)
      ALLOCATABLE MAP
      SAVE L, N, MAP
      DATA N /0/
C 
      IF(N.NE.NLOC.OR.L.NE.LOCPAR)THEN
      L=LOCPAR
      N=NLOC
      IF(ALLOCATED(MAP))DEALLOCATE(MAP)
      IF(N.GT.0)ALLOCATE(MAP(L,N))
      ENDIF
      RETURN
C 
      ENTRY GMAPWRITE(I,J)
      IF(J.GE.1.AND.J.LE.N)THEN
      DO 10,K=1,L
      MAP(K,J)=I(K)
  10  CONTINUE
      ENDIF
      RETURN
C 
      ENTRY GMAPREAD(I,J)
      IF(J.GE.1.AND.J.LE.N)THEN
      DO 20,K=1,L
      I(K)=MAP(K,J)
  20  CONTINUE
      ELSE
      DO 30,K=1,L
      I(K)=0
  30  CONTINUE
      ENDIF
      RETURN
C 
      ENTRY GMAPSHOWDIM(LO,NL)
      LO=L
      NL=N
      RETURN
      END

