C     The two subroutines in this library are useful for timing numeric
C     computations and displaying the results in common time units.
C 
C     Originally written in FORTRAN 77 by Dmitri Toptygin in 1993.
C 
C     In 2003 upgraded from Lahey FORTRAN 77 compiler to
C     Intel FORTRAN 95 compiler.
C 
C     No changes were necessary for modern fortran compilers like gfortran.
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
      SUBROUTINE TIMER(ITIME)
      IMPLICIT INTEGER(4)(I-N)
      CHARACTER(8)CDATE
      CHARACTER(10)CTIME
      CHARACTER(5)CZONE
      DIMENSION IARR(8)
      CALL DATE_AND_TIME(CDATE,CTIME,CZONE,IARR)
      ITIME=IARR(5)*3600000+IARR(6)*60000+IARR(7)*1000+IARR(8)
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE SPELLTIME(ITIME,HRMNSCMS)
      IMPLICIT INTEGER(4)(I-N)
      CHARACTER(17)HRMNSCMS
 1    FORMAT(I2.2,'h ',I2.2,'m ',I2.2,'s ',I3.3,'ms')
C 
C     Lahey FORTRAN 77 JTIME=ITIME*10
C     Intel FORTRAN 95 JTIME=ITIME
C 
      JTIME=ITIME
      IF(JTIME.LT.0)JTIME=JTIME+86400000
      IHR=JTIME/3600000
      JTIME=JTIME-IHR*3600000
      IMIN=JTIME/60000
      JTIME=JTIME-IMIN*60000
      ISEC=JTIME/1000
      JTIME=JTIME-ISEC*1000
      IMSEC=JTIME
      WRITE(HRMNSCMS,1)IHR,IMIN,ISEC,IMSEC
      RETURN
      END

