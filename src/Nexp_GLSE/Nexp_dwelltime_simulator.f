      PROGRAM Nexp_dwelltime_simulator
C 
C     This program generates a set of random dwell times.
C     The number of dwell times is determined by the user.
C 
C     The theoretical survival function is a linear combination of Nexp
C     exponential functions; the number Nexp is determined by the user.
C     The amplitudes and lifetimes for all exponential terms are
C     determined by the user.
C 
C     For details read the documentation file "Nexp_GLSE_manual.pdf"
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
      ALLOCATABLE ALP, TAU, ISEED
      DIMENSION ALP(:), TAU(:), ISEED(:)
      CHARACTER(128)FILENAME
   1  FORMAT(A)
   3  FORMAT(1X)
   4  FORMAT(A)
   5  FORMAT('Α(',I1.1,'),  τ(',I1.1,') ? ')
   6  FORMAT('Α(',I2.2,'),  τ(',I2.2,') ? ')
   7  FORMAT('Α(',I3.3,'),  τ(',I3.3,') ? ')
   8  FORMAT('Α(',I4.4,'),  τ(',I4.4,') ? ')
   9  FORMAT('Random number generator seed (',I2,' integer numbers)? ')
C 
      CALL RANDOM_SEED(SIZE=ISEEDDIM)
      ALLOCATE(ISEED(ISEEDDIM), STAT=IERR)
      IF(IERR.NE.0)GOTO 600
      CALL RANDOM_INIT(.FALSE.,.TRUE.)
      CALL RANDOM_SEED(GET=ISEED)
C 
      WRITE(*,3,ERR=620)
      WRITE(*,4,ADVANCE='NO',ERR=620)'Number of exponentials? '
      READ(*,*,ERR=601,END=601)NEXP
      WRITE(*,*,ERR=620)NEXP
      IF(NEXP.LT.1.OR.NEXP.GT.9999)GOTO 607
      ALLOCATE(ALP(NEXP), TAU(NEXP), STAT=IERR)
      IF(IERR.NE.0)GOTO 600
C 
      SUMA=0.0D0
      DO 50,I=1,NEXP
      IF(NEXP.LE.9)THEN
        WRITE(*,5,ADVANCE='NO',ERR=620)I,I
      ELSEIF(NEXP.LE.99)THEN
        WRITE(*,6,ADVANCE='NO',ERR=620)I,I
      ELSEIF(NEXP.LE.999)THEN
        WRITE(*,7,ADVANCE='NO',ERR=620)I,I
      ELSE
        WRITE(*,8,ADVANCE='NO',ERR=620)I,I
      ENDIF
      READ(*,*,ERR=601,END=601)A,T
      WRITE(*,*,ERR=620)A,T
      IF(T.LE.0.0D0)GOTO 607
      ALP(I)=A
      TAU(I)=T
      SUMA=SUMA+A
  50  CONTINUE
      IF(SUMA.LE.0.0D0)GOTO 607
      DO 51,I=1,NEXP
      ALP(I)=ALP(I)/SUMA
  51  CONTINUE
C 
      WRITE(*,4,ADVANCE='NO',ERR=620)'Number of dwell times? '
      READ(*,*,ERR=601,END=601)ND
      WRITE(*,*,ERR=620)ND
      IF(ND.LT.1)GOTO 607
C 
      WRITE(*,9,ADVANCE='NO',ERR=620)ISEEDDIM
      READ(*,*,ERR=70,END=70)(ISEED(I),I=1,ISEEDDIM)
  70  CONTINUE
      WRITE(*,*,ERR=620)(ISEED(I),I=1,ISEEDDIM)
      CALL RANDOM_SEED(PUT=ISEED)
C 
      WRITE(*,4,ADVANCE='NO',ERR=620)'Output file name? '
      READ(*,1,ERR=601,END=601)FILENAME
      LFN=LEN_TRIM(FILENAME)
      WRITE(*,1,ERR=620)FILENAME(1:LFN)
      OPEN(2,FILE=FILENAME(1:LFN),ERR=621)
C 
      DO 100,N=1,ND
      CALL RANDOM_NUMBER(Y)
      X=0.0D0
 110  CONTINUE
      F=0.0D0
      D=0.0D0
      DO 111,I=1,NEXP
      T= TAU(I)
      AE=ALP(I)*EXP(-X/T)
      F=F+AE
      D=D-AE/T
 111  CONTINUE
      DIFF=Y-F
      X=X+DIFF/D
      IF(ABS(DIFF).GT.1.0D-14)GOTO 110
      WRITE(2,*,ERR=621)X
 100  CONTINUE
C 
      CLOSE(2,ERR=621)
C 
      DEALLOCATE(ISEED, ALP, TAU, STAT=IERR)
      IF(IERR.NE.0)GOTO 600
C 
      STOP 'OK.  '
 600  STOP 'ARRAY ALLOCATION ERROR.  '
 601  STOP 'ERROR while reading from the standard input file.  '
 607  STOP 'ERROR: number out of range.  '
 620  STOP 'ERROR while writing to the standard output.  '
 621  STOP 'ERROR while writing to the output file.  '
      END

