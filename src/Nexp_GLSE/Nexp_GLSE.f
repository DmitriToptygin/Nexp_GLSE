      PROGRAM Nexp_GLSE
C 
C     A Global Least Square Estimator for fitting multiple lifetime or
C     dwelltime distributions by a linear combination of N exponentials.
C 
C     This program assumes that the Impulse Response Function equals
C     Dirac  delta function, therefore no convolution is performed.
C 
C     The weighting of the squared residuals is based on the assumption
C     of Poissonian distribution of the counts in the bins.
C 
C     The bins may have equal or unequal widths.
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
      DIMENSION PAR(:), STD(:), FLO(:), IDP(:), MGP(:)
      LOGICAL FLO, OK
      DIMENSION ARGM(:,:), DATA(:), PRED(:), WGHT(:)
      DIMENSION HESSIAN(:,:)
      DIMENSION JD(:), KD(:), LD(:)
      ALLOCATABLE  PAR,STD,FLO, IDP,MGP, ARGM, DATA,PRED,WGHT, HESSIAN
      ALLOCATABLE JD, KD, LD
      PARAMETER(MAXNEXP=32, NARG=2)
      DIMENSION LOCMAP(2,MAXNEXP)
      EXTERNAL NEXP_GMODEL, NEXP_GDOMAIN
      LOGICAL GLFIRST, SWRITE, DWRITE, SAVEINVHESSM
      DATA GLFIRST /.TRUE./, MAXLTRAIN /1/,
     & NABORT /1024/, ACCUR/1.0D-8/, AMARQ/1.0D-8/,
     & SWRITE /.TRUE./, DWRITE /.FALSE./, SAVEINVHESSM /.TRUE./
      CHARACTER(1) ANSWER, C1LOGIC
      CHARACTER(256) FNAME, OUTFILE
      CHARACTER(22) EXCLFILE
      CHARACTER(17)HRMNSCMS
      LOGICAL IDENTBINDEF, TWOCOL, DEBUG
   1  FORMAT(A)
   3  FORMAT(1X)
   4  FORMAT(A)
   5  FORMAT(A,'_',I4.4,A1,'.txt')
C 
      INQUIRE(FILE='DEBUG',EXIST=DEBUG)
C 
      WRITE(*,3,ERR=620)
      WRITE(*,4,ADVANCE='NO',ERR=620)
     &'Number of data sets to be analyzed? '
      READ(*,*,ERR=601,END=601)NSET
      WRITE(*,*,ERR=620)NSET
      IF(NSET.LT.1.OR.NSET.GT.9999)GOTO 607
C 
      IF(NSET.EQ.1)THEN
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &  'Number of exponentials? '
      ELSE
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &'Number of exponentials for one data set (same N for all sets)? '
      ENDIF
      READ(*,*,ERR=601,END=601)NEXP
      WRITE(*,*,ERR=620)NEXP
      IF(NEXP.LT.1.OR.NEXP.GT.MAXNEXP)GOTO 607
C 
      CALL GMAPDEFINE(2*NEXP,NSET)
C 
      OPEN(2,FILE='results.txt',ERR=621)
      IF(NSET.GT.1)THEN
      WRITE(2,'("Global",I3,"-exponential fitting of lifetime or dwell t
     &ime data,",I5," data sets.")',ERR=621)NEXP, NSET
      ELSE
      WRITE(2,'(I2,"-exponential fitting of lifetime or dwell time data.
     &")',ERR=621)NEXP
      ENDIF
      WRITE(2,'("Variances are calculated assuming Poissonian error stat
     &istics.")',ERR=621)
C 
      IF(NSET.EQ.1)THEN
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &  'Total number of model parameters? '
      ELSE
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &  'Total number of model parameters (global+local)? '
      ENDIF
      READ(*,*,ERR=601,END=601)NP
      WRITE(*,*,ERR=620)NP
      IF(NP.LT.2.OR.NP.GT.8192)GOTO 607
      WRITE(2,'(I4," model parameters.")',ERR=621)NP
C 
      IF(NSET.EQ.1)THEN
        IDENTBINDEF=.TRUE.
      ELSE
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &  'Are bins defined identically for all data sets (Y/N)? '
        READ(*,1,ERR=601,END=601)ANSWER
        WRITE(*,4,ERR=620)ANSWER
        IF    (ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')THEN
          IDENTBINDEF=.TRUE.
        ELSEIF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
          IDENTBINDEF=.FALSE.
        ELSE 
          GOTO 606
        ENDIF
      ENDIF
C 
      IF(IDENTBINDEF)THEN
        IF(NSET.EQ.1)THEN
          WRITE(*,4,ADVANCE='NO',ERR=620) 'Number of bins? '
        ELSE
          WRITE(2,'("Identical bin definitions for all data sets.")',
     &    ERR=621)
          WRITE(*,4,ADVANCE='NO',ERR=620)
     &   'Number of bins per one data set? '
        ENDIF
        READ(*,*,ERR=601,END=601)NBIN
        WRITE(*,*,ERR=620)NBIN
        IF(NBIN.LT.1.OR.NBIN.GT.16777216)GOTO 607
        NDMAX=NSET*NBIN
      ELSE
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &  'Estimated maximum total number of bins in all data sets? '
        READ(*,*,ERR=601,END=601)NDMAX
        WRITE(*,*,ERR=620)NDMAX
      ENDIF
      IF(NDMAX.LT.1.OR.NDMAX.GT.16777216)GOTO 607
C 
      ALLOCATE(PAR(NP),STD(NP),FLO(NP),IDP(NP),MGP(NP),ARGM(NARG,NDMAX),
     & DATA(NDMAX),PRED(NDMAX),WGHT(NDMAX), HESSIAN(NP,NP),
     & JD(NSET), KD(NSET), LD(NSET), STAT=IERR)
      IF(IERR.NE.0)GOTO 608
      DO 10,I=1,NP
      PAR(I)=0.0D0
      FLO(I)=.FALSE.
      IDP(I)=0
      MGP(I)=0
  10  CONTINUE
C 
      IF(IDENTBINDEF)THEN
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &  'Create equally spaced bins or read bin times from File (C/F)? '
        READ(*,1,ERR=601,END=601)ANSWER
        WRITE(*,4,ERR=620)ANSWER
        IF    (ANSWER.EQ.'C'.OR.ANSWER.EQ.'c')THEN
          WRITE(*,4,ADVANCE='NO',ERR=620)'Bin width? '
          READ(*,*,ERR=601,END=601)WBIN
          WRITE(*,*,ERR=620)WBIN
          IF(WBIN.LE.0.0D0)GOTO 607
          DO 21,N=1,NBIN
          ARGM(2,N)=DBLE(N)*WBIN
  21      CONTINUE
          WRITE(2,'(I8," bins; bin width =", 1PE13.6)',ERR=621)NBIN,WBIN
        ELSEIF(ANSWER.EQ.'F'.OR.ANSWER.EQ.'f')THEN
          WRITE(*,4,ADVANCE='NO',ERR=620)
     &    'Name of the file with bin end times? '
          READ(*,1,ERR=601,END=601)FNAME
          LFN=LEN_TRIM(FNAME)
          WRITE(*,4,ERR=620)FNAME(1:LFN)
          OPEN(1,FILE=FNAME(1:LFN),ERR=602)
          T0=0.0D0
          DO 22,N=1,NBIN
          READ(1,*,ERR=603,END=604)T
          IF(T.LE.T0)GOTO 605
          T0=T
          ARGM(2,N)=T
  22      CONTINUE
          CLOSE(1)
          WRITE(2,'(I8," bins; bin end times from file ",A)',ERR=621)
     &    NBIN,FNAME(1:LFN)
        ELSE
          GOTO 606
        ENDIF
      ENDIF
C 
      ND=0
      NDNZW=0
      DO 100,ISET=1,NSET
      RSET=DBLE(ISET)
      ND1=ND+1
      WRITE(*,3,ERR=620)
      WRITE(2,3,ERR=621)
      IF(NSET.GT.1)THEN
        WRITE(*,101,ERR=620)ISET
        WRITE(2,101,ERR=621)ISET
      ENDIF
 101  FORMAT('**** Data set ',I4.4,' ****')
      IF(IDENTBINDEF)THEN
        DO 110,N=1,NBIN
        ND=ND+1
        IF(ND.GT.NDMAX)GOTO 611
        ARGM(1,ND)=RSET
        ARGM(2,ND)=ARGM(2,N)
        DATA(ND)=0.0D0
 110    CONTINUE
      ELSE
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &  'Number of bins for this data set? '
        READ(*,*,ERR=601,END=601)NBIN
        WRITE(*,*,ERR=620)NBIN
        IF(NBIN.LT.1.OR.NBIN.GT.16777216)GOTO 607
        WRITE(*,4,ADVANCE='NO',ERR=620)
     &  'Create equally spaced bins or read bin times from File (C/F)? '
        READ(*,1,ERR=601,END=601)ANSWER
        WRITE(*,4,ERR=620)ANSWER
        IF    (ANSWER.EQ.'C'.OR.ANSWER.EQ.'c')THEN
          WRITE(*,4,ADVANCE='NO',ERR=620)'Bin width? '
          READ(*,*,ERR=601,END=601)WBIN
          WRITE(*,*,ERR=620)WBIN
          IF(WBIN.LE.0.0D0)GOTO 607
          DO 121,N=1,NBIN
          ND=ND+1
          IF(ND.GT.NDMAX)GOTO 611
          ARGM(1,ND)=RSET
          ARGM(2,ND)=DBLE(N)*WBIN
          DATA(ND)=0.0D0
 121      CONTINUE
          WRITE(2,'(I8," bins; bin width =", 1PE13.6)',ERR=621)NBIN,WBIN
        ELSEIF(ANSWER.EQ.'F'.OR.ANSWER.EQ.'f')THEN
          WRITE(*,4,ADVANCE='NO',ERR=620)
     &    'Name of the file with bin end times? '
          READ(*,1,ERR=601,END=601)FNAME
          LFN=LEN_TRIM(FNAME)
          WRITE(*,4,ERR=620)FNAME(1:LFN)
          OPEN(1,FILE=FNAME(1:LFN),ERR=602)
          T0=0.0D0
          DO 122,N=1,NBIN
          READ(1,*,ERR=603,END=604)T
          IF(T.LE.T0)GOTO 605
          T0=T
          ND=ND+1
          IF(ND.GT.NDMAX)GOTO 611
          ARGM(1,ND)=RSET
          ARGM(2,ND)=T
          DATA(ND)=0.0D0
 122      CONTINUE
          CLOSE(1)
          WRITE(2,'(I8," bins; bin end times from file ",A)',ERR=621)
     &    NBIN,FNAME(1:LFN)
        ELSE
          GOTO 606
        ENDIF
      ENDIF
C 
      WRITE(*,4,ADVANCE='NO',ERR=620)
     &'Lifetimes / Dwell times file name? '
      READ(*,1,ERR=601,END=601)FNAME
      LFN=LEN_TRIM(FNAME)
      WRITE(*,4,ERR=620)FNAME(1:LFN)
      OPEN(1,FILE=FNAME(1:LFN),ERR=602)
C 
      WRITE(*,4,ADVANCE='NO',ERR=620)
     &          'Is this a List (1-col) or Binned data (2-col)? '
      READ(*,1,ERR=601,END=601)ANSWER
      WRITE(*,4,ERR=620)ANSWER
      IF    (ANSWER.EQ.'1'.OR.ANSWER.EQ.'L'.OR.ANSWER.EQ.'l')THEN
        TWOCOL=.FALSE.
        NCOL=1
      ELSEIF(ANSWER.EQ.'2'.OR.ANSWER.EQ.'B'.OR.ANSWER.EQ.'b')THEN
        TWOCOL=.TRUE.
        NCOL=2
      ELSE 
        GOTO 606
      ENDIF
C 
      NREAD=0
      NOFFRANGE=0
 130  CONTINUE
      IF(TWOCOL)THEN
        READ(1,*,ERR=603,END=140)T,RNCOUNT
        NCOUNT=NINT(RNCOUNT)
        IF(NCOUNT.LT.0)GOTO 615
      ELSE
        READ(1,*,ERR=603,END=140)T
        NCOUNT=1
      ENDIF
      IF(T.LE.0.0D0)GOTO 609
      NREAD=NREAD+NCOUNT
      TR=0.99999999999999D0*T
      DO 131,N=ND1,ND
      IF(TR.LE.ARGM(2,N))THEN
        DATA(N)=DATA(N)+DBLE(NCOUNT)
        GOTO 130
      ENDIF
 131  CONTINUE
      NOFFRANGE=NOFFRANGE+NCOUNT
      GOTO 130
 140  CONTINUE
      CLOSE(1)
      WRITE(2,'("Data file: ",A," (",I1.1,"-column).")',ERR=621)
     &FNAME(1:LFN),NCOL
      IF(NREAD.LT.1)GOTO 604
        WRITE(*,141,ERR=620)NREAD
        WRITE(2,141,ERR=621)NREAD
      IF(NOFFRANGE.GT.0)THEN
        WRITE(*,142,ERR=620)NOFFRANGE
        WRITE(2,142,ERR=621)NOFFRANGE
      ENDIF
 141  FORMAT(I11,' lifetimes / dwell times read from file.')
 142  FORMAT(I11,' lifetimes / dwell times are out of range; it is recom
     &mended to increase the number of bins.')
C 
      NDD=ND-ND1+1
      CALL EXPSERIES_WEIGHTING(ARGM(1,ND1),DATA(ND1),WGHT(ND1),NDD,OK,
     & DEBUG)
      IF(.NOT.OK)THEN
        WRITE(2,4,ERR=621)'ERROR: Weight estimation failed. There is som
     &ething wrong with the data.  '
        GOTO 612
      ENDIF
C 
      L=0
      DO 150,N=ND1,ND
      IF(WGHT(N).GT.0.0D0)L=L+1
 150  CONTINUE
      JD(ISET)=ND1
      KD(ISET)=ND
      LD(ISET)=L+ND1-1
      NDNZW=NDNZW+L
      IF(L.LT.NBIN)THEN
        WRITE(*,151,ERR=620)L
        WRITE(2,151,ERR=621)L
      ELSE
        WRITE(*,152,ERR=620)L
        WRITE(2,152,ERR=621)L
      ENDIF
 151  FORMAT('Variances were reliably estimated for the first',I8,
     &' bins;',/,18X,'the data in remaining bins will be ignored.')
 152  FORMAT('Variances were reliably estimated for all',I8,' bins.')
C 
      LEX=0
      WRITE(EXCLFILE,'("excluded_bins_",I4.4,".txt")')ISET
      OPEN(1,FILE=EXCLFILE,STATUS='OLD',ERR=161)
      GOTO 162
 161  CONTINUE
      EXCLFILE='excluded_bins.txt'
      OPEN(1,FILE='excluded_bins.txt',STATUS='OLD',ERR=165)
 162  CONTINUE
      READ(1,*,END=163,ERR=616)IEXCLB
      IF(IEXCLB.LE.0) GOTO 616
      NEXCLB=ND1-1+IEXCLB
      IF(NEXCLB.GT.ND)GOTO 162
      W=WGHT(NEXCLB)
      WGHT(NEXCLB)=0.0D0
      IF(W.GT.0.0D0)LEX=LEX+1
      GOTO 162
 163  CLOSE(1)
      IF(LEX.GT.0)THEN
        LTEF=LEN_TRIM(EXCLFILE)
        WRITE(*,164,ERR=620)LEX,EXCLFILE(1:LTEF)
        WRITE(2,164,ERR=621)LEX,EXCLFILE(1:LTEF)
        NDNZW=NDNZW-LEX
      ENDIF
 164  FORMAT(I11,' bin(s) will be excluded, according to the list in "',
     &A,'".')
 165  CONTINUE
C 
      IF(ISET.EQ.1)THEN
        WRITE(*,3,ERR=620)
        WRITE(*,'(A,I4.4)',ERR=620)
     &  'All parameters are numbered sequentially from 0001 to ',NP
      ENDIF
      IF(NSET.GT.1.AND.ISET.EQ.1)THEN
        WRITE(*,4,ERR=620)'Global parameters are those that correspond t
     &o multiple data sets.'
        WRITE(*,4,ERR=620)'For faster operation global parameters should
     & have smaller numbers.'
      ENDIF
C 
      DO 170,IEXP=1,NEXP
      WRITE(*,171,ADVANCE='NO',ERR=620)IEXP
 171  FORMAT('Exponential ',I2.2,': parameter numbers for Α and τ? ')
      READ(*,*,ERR=601,END=601)LALP,LTAU
      WRITE(*,*,ERR=620)LALP,LTAU
      IF(LALP.LT.1.OR.LALP.GT.NP)GOTO 607
      IF(IDP(LALP).EQ.0)THEN
        IDP(LALP)=1
      ELSE
        IF(IDP(LALP).NE.1)GOTO 610
      ENDIF
      MGP(LALP)=MGP(LALP)+1
      LOCMAP(1,IEXP)=LALP
      IF(LTAU.LT.1.OR.LTAU.GT.NP)GOTO 607
      IF(IDP(LTAU).EQ.0)THEN
        IDP(LTAU)=2
      ELSE
        IF(IDP(LTAU).NE.2)GOTO 610
      ENDIF
      MGP(LTAU)=MGP(LTAU)+1
      LOCMAP(2,IEXP)=LTAU
      IF(IEXP.EQ.1)THEN
        WRITE(2,172,ERR=621)IEXP,LALP,IEXP,LTAU
      ELSE
        WRITE(2,173,ERR=621)IEXP,LALP,IEXP,LTAU
      ENDIF
 172  FORMAT('Model parameters: ','Α(',I2.2,') = PAR(',I4.4,'),  τ(',
     &I2.2,') = PAR(',I4.4,')')
 173  FORMAT('                  ','Α(',I2.2,') = PAR(',I4.4,'),  τ(',
     &I2.2,') = PAR(',I4.4,')')
 170  CONTINUE
C 
      CALL GMAPWRITE(LOCMAP,ISET)
 100  CONTINUE
C 
      WRITE(*,3,ERR=620)
      WRITE(*,4,ERR=620)'Enter initial guesses for parameters followed b
     &y space and T (tunable) or F (fixed).'
      WRITE(2,3,ERR=621)
      WRITE(2,4,ERR=621)'Initial guesses:'
      NFPAR=0
      DO 180,I=1,NP
      IF(IDP(I).NE.0)THEN
        WRITE(*,181,ADVANCE='NO',ERR=620)I
        READ(*,*,ERR=601,END=601)PAR(I),FLO(I)
        WRITE(*,182,ERR=620)PAR(I),C1LOGIC(FLO(I))
        WRITE(2,183,ERR=621)I,PAR(I),C1LOGIC(FLO(I))
        IF(FLO(I))NFPAR=NFPAR+1
      ENDIF
 180  CONTINUE
      NDEGF=NDNZW-NFPAR
 181  FORMAT('    PAR(',I4.4,') ? ')
 182  FORMAT(1PE13.6,2X,A1)
 183  FORMAT('PAR(',I4.4,') =',1PE13.6,2X,A1)
C 
      WRITE(*,3,ERR=620)
      WRITE(2,3,ERR=621)
      CALL TIMER(ITIME1)
C 
      CALL TOPFIT(NEXP_GMODEL,PAR,STD,FLO,NP,MAXLTRAIN,GLFIRST,
     &NEXP_GDOMAIN,ARGM,DATA,PRED,WGHT,ND,NABORT,ACCUR,AMARQ,CHI2,OK,
     &SWRITE,DWRITE,SAVEINVHESSM)
C 
      CALL TIMER(ITIME2)
      CALL SPELLTIME(ITIME2-ITIME1,HRMNSCMS)
      WRITE(*,3,ERR=620)
      IF(OK)THEN
        WRITE(*,'("The χ² minimization took ",A)',ERR=620)HRMNSCMS
        WRITE(2,'("The χ² minimization took ",A)',ERR=621)HRMNSCMS
      ELSE
        WRITE(*,201,ERR=620)
        WRITE(2,201,ERR=621)
      ENDIF
      WRITE(*,3,ERR=620)
      WRITE(2,3,ERR=621)
      RNDF=DBLE(NDEGF)
      C950=CHI2_ICDF(0.950D0,NDEGF)/RNDF
      C995=CHI2_ICDF(0.995D0,NDEGF)/RNDF
      WSSQ=CHI2*RNDF
      PCHI=1.0D0-CHI2_CDF(WSSQ,NDEGF)
      AIC=WSSQ+REAL(2*NFPAR-NDNZW)
      WRITE(*,202,ERR=620)CHI2
      WRITE(*,203,ERR=620)NDNZW
      WRITE(*,204,ERR=620)NFPAR
      WRITE(*,205,ERR=620)NDEGF
      WRITE(*,206,ERR=620)PCHI*100.0D0
      WRITE(*,207,ERR=620)95.0D0,C950
      WRITE(*,207,ERR=620)99.5D0,C995
      WRITE(*,208,ERR=620)WSSQ
      WRITE(*,209,ERR=620)AIC
      WRITE(2,202,ERR=621)CHI2
      WRITE(2,203,ERR=621)NDNZW
      WRITE(2,204,ERR=621)NFPAR
      WRITE(2,205,ERR=621)NDEGF
      WRITE(2,206,ERR=621)PCHI*100.0D0
      WRITE(2,207,ERR=621)95.0D0,C950
      WRITE(2,207,ERR=621)99.5D0,C995
      WRITE(2,208,ERR=621)WSSQ
      WRITE(2,209,ERR=621)AIC
 201  FORMAT('WARNING: χ² minimization failed. Try different ',
     &'initial guesses or try fixing some parameters.')
 202  FORMAT('Reduced χ² =',F13.3)
 203  FORMAT('Number of bins with nonzero weights NB =',I8)
 204  FORMAT('Number of  free  model  parameters  NP =',I8)
 205  FORMAT('Number of  degrees  of  freedom  NB-NP =',I8)
 206  FORMAT('The probability of the above χ² or greater is'
     &,F7.3,'%')
 207  FORMAT('With the probability of',F5.1,
     &'% reduced χ² must not exceed',F7.3)
 208  FORMAT('Absolute χ² =',F14.2)
 209  FORMAT('Akaike I.C. =',F14.2)
C 
      OPEN(3,FILE='chisquare.txt',ERR=623)
      WRITE(3,*,ERR=623)CHI2
      WRITE(3,*,ERR=623)NDNZW
      WRITE(3,*,ERR=623)NFPAR
      WRITE(3,*,ERR=623)NDEGF
      WRITE(3,*,ERR=623)PCHI
      WRITE(3,*,ERR=623)C950
      WRITE(3,*,ERR=623)C995
      WRITE(3,*,ERR=623)WSSQ
      WRITE(3,*,ERR=623)AIC
      WRITE(3,1,ERR=623)C1LOGIC(OK)
      CLOSE(3,ERR=623)
C 
      WRITE(*,3,ERR=620)
      WRITE(*,211,ERR=620)
      WRITE(2,3,ERR=621)
      WRITE(2,211,ERR=621)
      DO 210,I=1,NP
      IF(IDP(I).NE.0)THEN
        IF(FLO(I))THEN
          WRITE(*,212,ERR=620)I,PAR(I),STD(I)
          WRITE(2,212,ERR=621)I,PAR(I),STD(I)
        ELSE
          WRITE(*,213,ERR=620)I,PAR(I)
          WRITE(2,213,ERR=621)I,PAR(I)
        ENDIF
      ENDIF
 210  CONTINUE
 211  FORMAT('Best-fit parameter values with standard deviations:')
 212  FORMAT('PAR(',I4.4,') =',1PE13.6,' ±',1PE13.6)
 213  FORMAT('PAR(',I4.4,') =',1PE13.6,' known a priori')
C 
      OPEN(1,FILE='invhessm.bin',FORM='UNFORMATTED',ACCESS='DIRECT',
     &RECL=8,STATUS='UNKNOWN')
      IREC=0
      DO 221,J=1,NP
      DO 222,I=1,NP
      IREC=IREC+1
      READ(1,REC=IREC)HESSIAN(I,J)
 222  CONTINUE
 221  CONTINUE
      CLOSE(1)
C 
      DO 300,ISET=1,NSET
      CALL GMAPREAD(LOCMAP,ISET)
      WRITE(*,3,ERR=620)
      WRITE(2,3,ERR=621)
      IF(NSET.GT.1)THEN
        WRITE(*,301,ERR=620)ISET
        WRITE(2,301,ERR=621)ISET
      ENDIF
 301  FORMAT('Best-fit parameter values for data set ',I4.4)
      SUMA=0.0D0
      SUMB=0.0D0
      SUMH=0.0D0
      SUMQ=0.0D0
      DO 310,IEXP=1,NEXP
      LMAI=LOCMAP(1,IEXP)
      LMTI=LOCMAP(2,IEXP)
      AI=PAR(LMAI)
      RI=1.0D0/PAR(LMTI)
      BI=AI*RI
      SUMA=SUMA+AI
      SUMB=SUMB+BI
      DO 311,JEXP=1,NEXP
      LMAJ=LOCMAP(1,JEXP)
      LMTJ=LOCMAP(2,JEXP)
      AJ=PAR(LMAJ)
      RJ=1.0D0/PAR(LMTJ)
      BJ=AJ*RJ
      SUMH=SUMH+ HESSIAN(LMAI,LMAJ)
      SUMQ=SUMQ+(HESSIAN(LMAI,LMAJ)
     &          -HESSIAN(LMTI,LMAJ)*BI
     &          -HESSIAN(LMAI,LMTJ)*BJ
     &          +HESSIAN(LMTI,LMTJ)*BI*BJ)*RI*RJ
 311  CONTINUE
 310  CONTINUE
      DO 320,IEXP=1,NEXP
      LMAI=LOCMAP(1,IEXP)
      LMTI=LOCMAP(2,IEXP)
      AI=PAR(LMAI)
      TI=PAR(LMTI)
      RI=1.0D0/TI
      BI=AI*RI
      SUMG=0.0D0
      SUMP=0.0D0
      DO 321,JEXP=1,NEXP
      LMAJ=LOCMAP(1,JEXP)
      LMTJ=LOCMAP(2,JEXP)
      AJ=PAR(LMAJ)
      RJ=1.0D0/PAR(LMTJ)
      BJ=AJ*RJ
      SUMG=SUMG+ (HESSIAN(LMAI,LMAJ)+HESSIAN(LMAJ,LMAI))
      SUMP=SUMP+((HESSIAN(LMAI,LMAJ)+HESSIAN(LMAJ,LMAI))
     &          -(HESSIAN(LMTI,LMAJ)+HESSIAN(LMAJ,LMTI))*BI
     &          -(HESSIAN(LMAI,LMTJ)+HESSIAN(LMTJ,LMAI))*BJ
     &          +(HESSIAN(LMTI,LMTJ)+HESSIAN(LMTJ,LMTI))*BI*BJ)*RI*RJ
 321  CONTINUE
      FA=AI/SUMA
      FB=BI/SUMB
      VB= (HESSIAN(LMAI,LMAI)-BI*(HESSIAN(LMAI,LMTI)+HESSIAN(LMTI,LMAI)
     &                       -BI* HESSIAN(LMTI,LMTI)))*RI*RI
      VFA=(HESSIAN(LMAI,LMAI)-FA*(SUMG-FA*SUMH))/(SUMA*SUMA)
      VFB=(VB                -FB*(SUMP-FB*SUMQ))/(SUMB*SUMB)
      IF(VB.GT.0.0D0)THEN
        SDB=SQRT(VB)
      ELSE
        SDB=0.0D0
      ENDIF
      IF(VFA.GT.0.0D0)THEN
       SDFA=SQRT(VFA)
      ELSE
       SDFA=0.0D0
      ENDIF
      IF(VFB.GT.0.0D0)THEN
       SDFB=SQRT(VFB)
      ELSE
       SDFA=0.0D0
      ENDIF
      WRITE(*,302,ERR=620)IEXP,PAR(LMAI),STD(LMAI)
      WRITE(*,303,ERR=620)IEXP,FA,SDFA
      WRITE(*,304,ERR=620)IEXP,TI,STD(LMTI)
      WRITE(*,305,ERR=620)IEXP,RI,STD(LMTI)*RI*RI
      WRITE(*,306,ERR=620)IEXP,BI,SDB
      WRITE(*,307,ERR=620)IEXP,FB,SDFB
      WRITE(2,302,ERR=621)IEXP,PAR(LMAI),STD(LMAI)
      WRITE(2,303,ERR=621)IEXP,FA,SDFA
      WRITE(2,304,ERR=621)IEXP,TI,STD(LMTI)
      WRITE(2,305,ERR=621)IEXP,RI,STD(LMTI)*RI*RI
      WRITE(2,306,ERR=620)IEXP,BI,SDB
      WRITE(2,307,ERR=620)IEXP,FB,SDFB
 320  CONTINUE
 300  CONTINUE
 302  FORMAT('absolute Α(',I2.2,') =',1PE13.6,' ±',1PE13.6)
 303  FORMAT('fraction Α(',I2.2,') =',1PE13.6,' ±',1PE13.6)
 304  FORMAT('         τ(',I2.2,') =',1PE13.6,' ±',1PE13.6)
 305  FORMAT('         k(',I2.2,') =',1PE13.6,' ±',1PE13.6)
 306  FORMAT('absolute α(',I2.2,') =',1PE13.6,' ±',1PE13.6)
 307  FORMAT('fraction α(',I2.2,') =',1PE13.6,' ±',1PE13.6)
C 
      CLOSE(2,ERR=621)
C 
      OPEN(2,FILE='parameters.txt',ERR=622)
      DO 350,I=1,NP
      WRITE(2,*,ERR=622)PAR(I),STD(I)
 350  CONTINUE
      CLOSE(2,ERR=622)
C 
      IF(NSET.GT.1)OPEN(3,FILE='local_chisquare.txt',ERR=624)
C 
      DO 400,ISET=1,NSET
      WSSQLOC=0.0D0
      DEGFLOC=0.0D0
C 
      WRITE(OUTFILE,5)'binned',ISET,'B'
      LOF=LEN_TRIM(OUTFILE)
      OPEN(21,FILE=OUTFILE(1:LOF),ERR=613)
      WRITE(OUTFILE,5)'binned',ISET,'D'
      LOF=LEN_TRIM(OUTFILE)
      OPEN(22,FILE=OUTFILE(1:LOF),ERR=613)
      T0=0.0D0
      WRITE(21,*,ERR=614)T0,0
      DO 410,N=JD(ISET),KD(ISET)
      T=ARGM(2,N)
      TH=0.5D0*(T+T0)
      X=DATA(N)
      IX=NINT(X)
      WRITE(21,*,ERR=614)T0,IX
      WRITE(21,*,ERR=614)T, IX
      WRITE(22,*,ERR=614)TH,IX
      T0=T
 410  CONTINUE
      WRITE(21,*,ERR=614)T0,0
      CLOSE(21,ERR=614)
      CLOSE(22,ERR=614)
C 
      WRITE(OUTFILE,5)'expser',ISET,'B'
      LOF=LEN_TRIM(OUTFILE)
      OPEN(23,FILE=OUTFILE(1:LOF),ERR=613)
      WRITE(OUTFILE,5)'expser',ISET,'D'
      LOF=LEN_TRIM(OUTFILE)
      OPEN(24,FILE=OUTFILE(1:LOF),ERR=613)
      WRITE(OUTFILE,5)'Nexpfit',ISET,'B'
      LOF=LEN_TRIM(OUTFILE)
      OPEN(25,FILE=OUTFILE(1:LOF),ERR=613)
      WRITE(OUTFILE,5)'Nexpfit',ISET,'D'
      LOF=LEN_TRIM(OUTFILE)
      OPEN(26,FILE=OUTFILE(1:LOF),ERR=613)
      WRITE(OUTFILE,5)'WghResid',ISET,'B'
      LOF=LEN_TRIM(OUTFILE)
      OPEN(27,FILE=OUTFILE(1:LOF),ERR=613)
      WRITE(OUTFILE,5)'WghResid',ISET,'D'
      LOF=LEN_TRIM(OUTFILE)
      OPEN(28,FILE=OUTFILE(1:LOF),ERR=613)
      T0=0.0D0
      WRITE(23,*,ERR=614)T0,0.0D0
      WRITE(25,*,ERR=614)T0,0.0D0
      WRITE(27,*,ERR=614)T0,0.0D0
      DO 420,N=JD(ISET),LD(ISET)
      T=ARGM(2,N)
      TH=0.5D0*(T+T0)
      W=WGHT(N)
      IF(W.GT.0.0D0)THEN
        Y=1.0D0/W
      ELSE
        Y=0.0D0
      ENDIF
      Z=PRED(N)
      DMZ=DATA(N)-Z
      WD=SQRT(W)*DMZ
      WSSQLOC=WSSQLOC+W*DMZ*DMZ
      IF(W.GT.0.0D0)DEGFLOC=DEGFLOC+1.0D0
      WRITE(23,*,ERR=614)T0,Y
      WRITE(23,*,ERR=614)T, Y
      WRITE(24,*,ERR=614)TH,Y
      WRITE(25,*,ERR=614)T0,Z
      WRITE(25,*,ERR=614)T, Z
      WRITE(26,*,ERR=614)TH,Z
      WRITE(27,*,ERR=614)T0,WD
      WRITE(27,*,ERR=614)T, WD
      WRITE(28,*,ERR=614)TH,WD
      T0=T
 420  CONTINUE
      WRITE(23,*,ERR=614)T0,0.0D0
      WRITE(25,*,ERR=614)T0,0.0D0
      WRITE(27,*,ERR=614)T0,0.0D0
      CLOSE(23,ERR=614)
      CLOSE(24,ERR=614)
      CLOSE(25,ERR=614)
      CLOSE(26,ERR=614)
      CLOSE(27,ERR=614)
      CLOSE(28,ERR=614)
C 
      CALL GMAPREAD(LOCMAP,ISET)
      DO 430,IEXP=1,NEXP
      DO 431,J=1,2
      L=LOCMAP(J,IEXP)
      IF(FLO(L))DEGFLOC=DEGFLOC-1.0D0/DBLE(MGP(L))
 431  CONTINUE
 430  CONTINUE
C 
      NDEGFLOC=NINT(DEGFLOC)
      REDCHI2LOC=WSSQLOC/DEGFLOC
      IF(NSET.GT.1)WRITE(3,441,ERR=624)ISET,REDCHI2LOC,
     &(1.0D0-CHI2_CDF(REDCHI2LOC*DBLE(NDEGFLOC),NDEGFLOC))*100.0D0
 441  FORMAT(I4,F13.3,1X,F7.3,'%')
C 
 400  CONTINUE
C 
      IF(NSET.GT.1)CLOSE(3,ERR=624)
C 
      DEALLOCATE(PAR,STD,FLO,IDP,MGP,ARGM,DATA,PRED,WGHT,HESSIAN,
     & JD, KD, LD, STAT=IERR)
      IF(IERR.NE.0)GOTO 608
C 
      CALL GMAPDEFINE(0,0)
C 
      STOP 'OK.  '
 601  STOP 'ERROR while reading from the standard input file.  '
 602  STOP 'ERROR while opening the above file as an existing file.  '
 603  STOP 'ERROR while reading from the above file.  '
 604  STOP 'ERROR: the above file is too short.  '
 605  STOP 'ERROR: bin end time must be greater than bin start time.  '
 606  STOP 'ERROR: unexpected answer.  '
 607  STOP 'ERROR: number out of range.  '
 608  STOP
     $'MEMORY ALLOCATION ERROR: too many parameters and rip curves.  '
 609  WRITE(*,*)'Input file line No',NREAD+1,', lifetime/dwelltime=',T
      STOP 'ERROR: zero or negative lifetime/dwelltime encountered.  '
 610  STOP 'ERROR: two different physical quantities are equated to the 
     &same numbered parameter.  '
 611  STOP 'ERROR: Estimated maximum total number of bins in all data se
     &ts has been exceeded.  '
 612  STOP 'ERROR: Weight estimation failed. There is something wrong wi
     &th the data.  '
 613  WRITE(*,*,ERR=620)OUTFILE(1:LOF)
      STOP 'ERROR while opening the above file.  '
 614  STOP 'ERROR while writing to an output file.  '
 615  STOP 'ERROR: a negative count in a bin encountered.  '
 616  STOP 'ERROR while reading from a file "excluded_bins*.txt".  '
 620  STOP 'ERROR while writing to the standard output.  '
 621  STOP 'ERROR while writing to file "results.txt".  '
 622  STOP 'ERROR while writing to file "parameters.txt".  '
 623  STOP 'ERROR while writing to file "chisquare.txt".  '
 624  STOP 'ERROR while writing to file "local_chisquare.txt".  '
      END
C     ------------------------------------------------------------------
      SUBROUTINE EXPSERIES_WEIGHTING(ARGM,DATA,WGHT,ND,OK,DEBUG)
C 
C     Copyright 2025 (C) Dmitri Toptygin, GPL license version 3.0.
C 
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
      DIMENSION ARGM(2,ND), DATA(ND), WGHT(ND)
      LOGICAL OK, DEBUG, OVERWLIMIT
      DIMENSION BASF(:,:)
      ALLOCATABLE BASF
      PARAMETER (TRATIO=1.8D0, EXTENSION=1.05D0, TFACTOR=2.3D0)
      PARAMETER (WLIMIT=10.0D0)
C 
      OK=.FALSE.
      LASTNZ=0
      DO 10,N=1,ND
      IF(DATA(N).GT.0.0D0)LASTNZ=N
  10  CONTINUE
      IF(LASTNZ.LE.2)GOTO 200
      TDCUT=EXTENSION*ARGM(2,LASTNZ)
      NDCUT=0
      DO 20,N=1,ND
      IF(ARGM(2,N).LE.TDCUT)NDCUT=N
  20  CONTINUE
C 
      TRL=LOG(TRATIO)
      TB1=ARGM(2,1)
      TBL=ARGM(2,LASTNZ)
      NFMAX=NINT(LOG(TBL/TB1)/TRL)+1
      ALLOCATE(BASF(NFMAX,NDCUT), STAT=IERR)
      IF(IERR.NE.0)GOTO 601
C 
      DO 50,I=1,NFMAX
      R=EXP(-TRL*DBLE(I-1))/TB1
      F0=1.0D0
      DO 51,N=1,NDCUT
      F=EXP(-R*ARGM(2,N))
      BASF(I,N)=F0-F
      F0=F
  51  CONTINUE
  50  CONTINUE
C 
      NF=NFMAX
 100  CONTINUE
 101  FORMAT('Exponential Series N =',I3.2,' τmin =',1PE10.3,' τmax =',
     &1PE10.3,' τn+1/τn =',0PF6.3,' Tmax =',1PE10.3)
C 
      TMAX=TFACTOR*EXP(TRL*DBLE(NF-1))*TB1
      NLAST=0
      DO 120,N=1,NDCUT
      IF(ARGM(2,N).LE.TMAX)NLAST=N
 120  CONTINUE
      IF(DEBUG)
     &WRITE(*,101)NF,TB1,TB1*TRATIO**(NF-1),TRATIO,ARGM(2,NLAST)
C 
      CALL LINCOMB_POISSON_WEIGHTING(BASF,NFMAX,NF,DATA,WGHT,NLAST,OK,
     &                                                          DEBUG)
      IF(.NOT.OK)THEN
        NF=NF-1
        IF(NF.GE.1)GOTO 100
      ENDIF
C 
      OVERWLIMIT=.FALSE.
      DO 131,N=1,NLAST
      OVERWLIMIT=OVERWLIMIT.OR.WGHT(N).GT.WLIMIT
      IF(OVERWLIMIT)WGHT(N)=0.0D0
 131  CONTINUE
      DO 132,N=NLAST+1,ND
      WGHT(N)=0.0D0
 132  CONTINUE
C 
      DEALLOCATE(BASF, STAT=IERR)
      IF(IERR.NE.0)GOTO 601
 200  RETURN
C 
 601  STOP 'Memory allocation error in EXPSERIES_WEIGHTING.  '
      END
C     ------------------------------------------------------------------
      FUNCTION NEXP_GDOMAIN(PAR,NP)
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
      LOGICAL NEXP_GDOMAIN, OK
      DIMENSION PAR(NP)
      PARAMETER(MAXNEXP=32)
      DIMENSION LOCMAP(2,MAXNEXP)
C 
      CALL GMAPSHOWDIM(LOCPAR,NSET)
      NEXP=LOCPAR/2
C 
      OK=.TRUE.
      DO 100,I=1,NSET
      CALL GMAPREAD(LOCMAP,I)
      DO 110,J=1,NEXP
      OK=OK.AND.PAR(LOCMAP(2,J)).GT.0.0D0
 110  CONTINUE
 100  CONTINUE
C 
      NEXP_GDOMAIN=OK
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE NEXP_GMODEL(PAR,NP,MAXNP,ARGM,N,PRED,DERIV,LTRAIN)
      IMPLICIT REAL(8)(A-H,O-Z), INTEGER(4)(I-N)
      DIMENSION PAR(NP), DERIV(NP)
      PARAMETER(MAXNEXP=32, NARG=2)
      DIMENSION ARGM(NARG,*)
      DIMENSION LOCMAP(2,MAXNEXP), TWOSINH(MAXNEXP), TWOCOSH(MAXNEXP)
      SAVE ISETOLD, NEXP, LOCMAP, TBEG, DELTAT, TWOSINH, TWOCOSH
      LOGICAL NEWDELTA
C 
      LTRAIN=1
      DO 10,I=1,NP
      DERIV(I)=0.0D0
  10  CONTINUE
C 
      IF(N.EQ.1)THEN
        CALL GMAPSHOWDIM(NLOCPAR,NSET)
        NEXP=NLOCPAR/2
      ENDIF
C 
      ISET=NINT(ARGM(1,N))
      IF(N.EQ.1.OR.ISET.NE.ISETOLD)THEN
        CALL GMAPREAD(LOCMAP,ISET)
        TBEG=0.0D0
        DELTAT=0.0D0
        DO 20,I=1,NEXP
        TWOSINH(I)=0.0D0
        TWOCOSH(I)=2.0D0
  20    CONTINUE
      ENDIF
      ISETOLD=ISET
C 
      TEND=ARGM(2,N)
      DELTANEW=TEND-TBEG
      NEWDELTA=ABS(DELTANEW-DELTAT).GE.2.0D-14*TEND
      IF(NEWDELTA)DELTAT=DELTANEW
      PRED=0.0D0
      DO 100,I=1,NEXP
      J=LOCMAP(1,I)
      K=LOCMAP(2,I)
      ALPHA=PAR(J)
      TAU  =PAR(K)
      TWOTAU=TAU+TAU
      HDT=DELTAT/TWOTAU
      IF(HDT.LT.18.420680744D0)THEN
        RTM=(TBEG+TEND)/(TWOTAU)
        ETM=EXP(-RTM)
        IF(NEWDELTA)THEN
          IF(HDT.LT.1.0D0)THEN
            HDT2=HDT*HDT
            X=2.0D0
            Y=HDT2/3.0D0
            Z=Y
 110        CONTINUE
            X=X+2.0D0
            Y=HDT2*Y/(X*(X+1.0D0))
            Z=Z+Y
            IF(Y.GT.8.4D-15)GOTO 110
            TSH=HDT*(2.0D0+Z)
          ELSE
            TSH=2.0D0*SINH(HDT)
          ENDIF
          TCH=SQRT(4.0D0+TSH*TSH)
          TWOSINH(I)=TSH
          TWOCOSH(I)=TCH
        ELSE
          TSH=TWOSINH(I)
          TCH=TWOCOSH(I)
        ENDIF
        EXPDIFF=ETM*TSH
        DEXDIFF=ETM*(RTM*TSH-HDT*TCH)/TAU
      ELSE
        RTB=TBEG/TAU
        ETB=EXP(-RTB)
        EXPDIFF=ETB
        DEXDIFF=ETB*RTB/TAU
      ENDIF
      PRED=PRED+        ALPHA*EXPDIFF
      DERIV(J)=DERIV(J)+      EXPDIFF
      DERIV(K)=DERIV(K)+ALPHA*DEXDIFF
 100  CONTINUE
      TBEG=TEND
C 
      RETURN
      END
C     ------------------------------------------------------------------
      FUNCTION C1LOGIC(B)
      CHARACTER(1) C1LOGIC
      LOGICAL B
      IF(B)THEN
      C1LOGIC='T'
      ELSE
      C1LOGIC='F'
      ENDIF
      RETURN
      END

