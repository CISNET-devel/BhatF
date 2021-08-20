      SUBROUTINE  MIGRAD(NDIM,NPAR,NFCN,X,AMIN)

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(NMAX=100,SENS=2.,XINF=32.)
C     CHARACTER*28 STAMP,STATUS
      CHARACTER*28 STATUS
      CHARACTER*10 LABELS(NMAX)
      CHARACTER    RNF(NMAX)
      CHARACTER*10 ANSWER
      CHARACTER*17 STR1
      CHARACTER*8  STR2
      CHARACTER*4  STR3
      CHARACTER(len=256) WRTBUFFER(16) ! WARNING: max 16 lines output at a time

CC        PERFORMS A LOCAL FUNCTION MINIMIZATION USING BASICALLY THE
CC        METHOD OF DAVIDON-FLETCHER-POWELL AS MODIFIED BY FLETCHER
CC        REF. -- FLETCHER, COMP.J. 13,317 (1970)   "SWITCHING METHOD"

      DIMENSION GS(250),R(NMAX),XXS(NMAX),FLNU(NMAX),VG(NMAX),
     $          VII(NMAX),IVN(NMAX)
      DIMENSION X(NMAX),DIRIN(NMAX),DSTEPS(NMAX)
      DIMENSION G(NMAX),G2(NMAX),V(NMAX,NMAX),IDN(32)
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)
      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /IFLAGS/INDEX,IMCMC,IBOOT,ISEED,IERR  
      DATA SLAMIN,SLAMAX,TLAMIN,TLAMAX/0.2D0, 3.0D0, 0.05D0, 6.0D0/

      interface
         subroutine logmessage(msg, nlines)
         character(len=*), intent(in) :: msg(:)
         integer, intent(in), optional :: nlines
         end subroutine logmessage
      end interface
      
      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,AMIN); NFCN = NFCN + 1

      ISW2=0
      NRETRY=0
      NFCNMX=10000
      VTEST=0.001D0
      APSI =.05d0	!or .01 for better convergence?
      UP=1.D0

      IF (NPAR .LE. 0)  RETURN
      NPFN = NFCN
      PARN=NPAR
      NPFN = NFCN
      PARN=NPAR

      RHO2 = 10.D0*APSI
      ROSTOP = 1.0D-5 * APSI

      TRACE=1.D0
      FS = AMIN
      wrtbuffer(:)=' '
      write (wrtbuffer,450)
      call logmessage(wrtbuffer)
      IF(IBOOT.EQ.0) then
         WRITE(WRTBUFFER, 470)  ROSTOP, APSI, VTEST
         call logmessage(wrtbuffer)
      endif
         
      GO TO 2
     
    1 WRITE(WRTBUFFER, 520)
      call logmessage(wrtbuffer)

    2 NPARD = NPAR

! determine step sizes via DSTEP routine

      CALL DQSTEP(NFCN,NDIM,NPAR,X,SENS,XINF,DSTEPS)

      DO  3 I= 1, NPAR
         
         D=DSTEPS(I)

      IF (ISW2 .GE. 1)  THEN
         D = 0.002D0 * DSQRT(DABS(V(I,I))*UP)
      ENDIF

    3 DIRIN(I) = D
C                                        . . . . . . STARTING GRADIENT
      NTRY = 0
    4 NEGG2 = 0

      DO 10  ID= 1, NPARD
      I = ID + NPAR - NPARD
      D = DIRIN(I)
      XTF = X(IVN(I))
      X(IVN(I)) = XTF + D
      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,FS1); NFCN = NFCN + 1
      X(IVN(I)) = XTF - D
      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,FS2); NFCN = NFCN + 1
      X(IVN(I)) = XTF
      GS(I) = (FS1-FS2)/(2. * D)
      G2(I) = (FS1 + FS2 - 2.*AMIN) / D**2
      IF (G2(I) .GT. 0.D0)  GO TO 10	
C                                        . . . SEARCH IF G2 .LE. 0. . .
      WRITE(WRTBUFFER, 520)
      call logmessage(wrtbuffer)
      NEGG2 = NEGG2 + 1
      NTRY = NTRY + 1
      IF (NTRY .GT. 6)  GO TO 230               !EGL: increased NTRY to 6
      D = 5.*DABS(DIRIN(I))  !EGL 50. -> 5.
      XBEG = XTF
      IF (GS(I) .LT. 0.)  DIRIN(I) = -DIRIN(I)
      KG = 0
      NF = 0
      NS = 0
    5 X(IVN(I)) = XTF + D

      IF(DABS(X(IVN(I))).GT.XINF) THEN
         WRITE(WRTBUFFER,*)'parameter',ivn(i),' near boundary'
         call logmessage(wrtbuffer,1)
         X(IVN(I))=SIGN(1.D0,X(IVN(I)))*XINF
      ENDIF

      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,F); NFCN = NFCN + 1

      IF (F.EQ.inf .or. F.EQ.nan) GO TO 9 
      IF (F.EQ.-inf .or. F.EQ.-nan) GO TO 9 

      IF (F .LE. AMIN. AND. DABS(X(IVN(I))).LE.XINF)  GO TO 6

C         FAILURE
 9    IF (KG .EQ. 1)  GO TO 8
      KG = -1
      NF = NF + 1
      D = -0.4*D
      IF (NF .LT. 10)  GO TO 5
      D = 1000.*D
      GO TO 7

C         SUCCESS
    6 XTF = X(IVN(I))
      D = 3.0*D
      AMIN = F
      KG = 1
      NS = NS + 1
      IF (NS .LT. 10)  GO TO 5
      IF (AMIN .LT. FS)  GO TO 8
      D = 0.001*D
    7 XTF = XBEG
      G2(I) = 1.0
      NEGG2 = NEGG2 - 1
    8 X(IVN(I)) = XTF
      DIRIN(I) = 0.1D0*D
      FS = AMIN
   10 CONTINUE

      IF (NEGG2 .GE. 1)  GO TO 4
      NTRY = 0
      MATGD = 1
C                                        . . . . . . DIAGONAL MATRIX
      IF (ISW2 .GT. 1)  GO TO 15
   11 NTRY = 1
      MATGD = 0
      DO 13 I= 1, NPAR
      DO 12 J= 1, NPAR
   12 V(I,J) = 0.D0
   13 V(I,I) = 2.D0/G2(I)
C                                        . . . GET SIGMA AND SET UP LOOP
 15   SIGMA = 0.D0
      DO 18 I= 1, NPAR
      IF (V(I,I) .LE. 0.D0)  GO TO 11
      RI = 0.D0
      DO 17 J= 1, NPAR
      XXS(I) = X(IVN(I))
   17 RI= RI+ V(I,J) * GS(J)
   18 SIGMA = SIGMA + GS(I) *RI * 0.5D0
      IF (SIGMA .GE. 0.D0)  GO TO 20
      WRITE(WRTBUFFER, 520)
      call logmessage(wrtbuffer)
      IF (NTRY.EQ.0)  GO TO 11
      ISW2 = 0
      GO TO 230
   20 ISW2 = 1
      ITER = 0

      IF(IBOOT.EQ.0) THEN
        WRITE(WRTBUFFER,*) 'GRADIENT SEARCH INTERMEDIATE OUTPUT:'
        call logmessage(wrtbuffer,1)
        WRITE(WRTBUFFER,9010)
        call logmessage(wrtbuffer)
        IV1=IVN(1)
        WRITE(WRTBUFFER,9020)
     &       ITER,AMIN,LABELS(IV1),SU(IV1),GS(1),G2(1),NFCN
        call logmessage(wrtbuffer)
      DO I=2,NPAR
        IVI=IVN(I)
        WRITE(WRTBUFFER,9030) LABELS(IVI),SU(IVI),GS(I),G2(I)
        call logmessage(wrtbuffer)
      ENDDO
      ENDIF
C                                        . . . . .  START MAIN LOOP

   24 CONTINUE
      GDEL = 0.D0
      DO 30  I=1,NPAR
      IVI=IVN(I)
      RI = 0.D0
      DO 25 J=1,NPAR
   25 RI = RI + V(I,J) *GS(J)
      DIRIN(I) = -0.5D0*RI
      GDEL = GDEL + DIRIN(I)*GS(I)
C                                        		.LINEAR SEARCH ALONG -VG  . . .
      X(IVI) =XXS(I) + DIRIN(I)
 30   IF(DABS(X(IVI)).GT.XINF) X(IVI)=SIGN(1.D0,X(IVI))*XINF

      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,F);NFCN=NFCN+1

C                                        . QUADR INTERP USING SLOPE GDEL
      DENOM = 2.D0*(F-AMIN-GDEL)
      IF (DENOM .LE. 0.D0)  GO TO 35
      SLAM = -GDEL/DENOM
      IF (SLAM .GT. SLAMAX)  GO TO 35
      IF (SLAM .LT. SLAMIN)  SLAM=SLAMIN
      GO TO 40
   35 SLAM = SLAMAX
   40 IF (DABS(SLAM-1.D0) .LT. 0.1D0)  GO TO 70
      DO 45 I= 1, NPAR
      IVI=IVN(I)
      X(IVI) =XXS(I) + SLAM*DIRIN(I) 
 45   IF(DABS(X(IVI)).GT.XINF) X(IVI)=SIGN(1.D0,X(IVI))*XINF

      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,F2); NFCN = NFCN + 1

C                                        . QUADR INTERP USING 3 POINTS
      AA = FS/SLAM
      BB = F/(1.D0-SLAM)
      CC = F2/ (SLAM*(SLAM-1.D0))
      DENOM = 2.D0*(AA+BB+CC)
      IF (DENOM .LE. 0.D0)  GO TO 48
      TLAM = (AA*(SLAM+1.D0) + BB*SLAM + CC)/DENOM
      IF (TLAM .GT. TLAMAX)  GO TO 48
      IF (TLAM .LT. TLAMIN)  TLAM=TLAMIN
      GO TO 50
   48 TLAM = TLAMAX
   50 CONTINUE
      DO 51 I= 1, NPAR
      IVI=IVN(I)
      X(IVI) = XXS(I)+TLAM*DIRIN(I)
 51   IF(DABS(X(IVI)).GT.XINF) X(IVI)=SIGN(1.D0,X(IVI))*XINF

      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,F3); NFCN = NFCN + 1

      IF (F.GE.AMIN .AND. F2.GE.AMIN .AND. F3.GE.AMIN) GO TO 200
      IF (F .LT. F2 .AND. F .LT. F3)  GO TO 61
      IF (F2 .LT. F3)  GO TO 58
   55 F = F3
      SLAM = TLAM
      GO TO 65
   58 F = F2
      GO TO 65
   61 SLAM = 1.D0
   65 DO 67 I= 1, NPAR
      IVI=IVN(I)
      DIRIN(I) = DIRIN(I)*SLAM
      X(IVI) = XXS(I) + DIRIN(I)
 67   IF(DABS(X(IVI)).GT.XINF) X(IVI)=SIGN(1.D0,X(IVI))*XINF

   70 AMIN = F
      ISW2 = 2
      IF (SIGMA+FS-AMIN .LT. ROSTOP)  THEN
         WRITE(WRTBUFFER,*) 'convergence type 1'
         call logmessage(wrtbuffer,1)
         GO TO 170
      ENDIF
      IF (SIGMA+RHO2+FS-AMIN .GT. APSI)  GO TO 75
      IF (TRACE .LT. VTEST)  THEN
         WRITE(WRTBUFFER,*) 'convergnece type 2'
         call logmessage(wrtbuffer,1)
         GO TO 170
      ENDIF

 75   IF (NFCN-NPFN .GE. NFCNMX)  GO TO 190
      ITER = ITER + 1

C                                        . . . GET GRADIENT AND SIGMA .

C ---  COMPUTE FIRST AND SECOND (DIAGONAL) DERIVATIVES 

      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,AMIN); NFCN = NFCN + 1

! determine step sizes via DSTEP routine
      
      IF(ISW2.EQ.0) THEN
      CALL DQSTEP(NFCN,NDIM,NPAR,X,SENS,XINF,DSTEPS)
      ELSE
!      WRITE(WRTBUFFER,*) 'using V() for computation of D'
      ENDIF

      DO I= 1, NPAR

      IF (ISW2 .GE. 1)  THEN
         D = 0.002D0 * DSQRT(DABS(V(I,I))*UP)
      ELSE
      D=DSTEPS(I)
      ENDIF

      XTF = X(IVN(I))
      X(IVN(I)) = XTF + D
      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,FS1); NFCN = NFCN + 1
      X(IVN(I)) = XTF - D
      CALL BTRAFO(NDIM,X,SU)
      CALL FUNC(SU,NDIM,FS2); NFCN = NFCN + 1
      X(IVN(I)) = XTF
      G(I) = (FS1-FS2)/(2.D0 * D)
      G2(I) = (FS1 + FS2 - 2.D0*AMIN) / D**2
      ENDDO

! ---------------------------------------------------

      RHO2 = SIGMA
      SIGMA = 0.D0
      GVG = 0.D0
      DELGAM = 0.D0
      DO 100 I= 1, NPAR
      RI = 0.D0
      VGI = 0.D0
      DO 90 J= 1, NPAR
      VGI = VGI + V(I,J)*(G(J)-GS(J))
   90 RI = RI + V(I,J) *G (J)
      R(I) = RI * 0.5D0
      VG(I) = VGI*0.5D0
      GAMI = G(I) - GS(I)
      GVG = GVG + GAMI*VG(I)
      DELGAM = DELGAM + DIRIN(I)*GAMI
  100 SIGMA = SIGMA + G(I)*R(I)
      IF (SIGMA .LT. 0.D0)  GO TO 1
      IF (GVG .LE. 0.D0)  GO TO 105
      IF (DELGAM .LE. 0.D0)  GO TO 105
      GO TO 107
  105 IF (SIGMA .LT. 0.1D0*ROSTOP)  THEN
         WRITE(WRTBUFFER,*) 'convergence crit. 3'
         call logmessage(wrtbuffer,1)
         GO TO 170
      ENDIF 
      GO TO 1
  107 CONTINUE
C                                        .  UPDATE COVARIANCE MATRIX
      TRACE=0.D0
      DO 120 I= 1, NPAR
      VII(I) = V(I,I)
      DO  120  J=1,NPAR
      D = DIRIN(I)*DIRIN(J)/DELGAM - VG(I)*VG(J)/GVG
  120 V(I,J) = V(I,J) + 2.D0*D
      IF (DELGAM .LE. GVG)  GO TO 135
      DO 125 I= 1, NPAR
  125 FLNU(I) = DIRIN(I)/DELGAM - VG(I)/GVG
      DO 130 I= 1, NPAR
      DO 130 J= 1, NPAR
  130 V(I,J) = V(I,J) + 2.D0*GVG*FLNU(I)*FLNU(J)
  135 CONTINUE
      DO 140 I= 1, NPAR
      XXS(I) = X(IVN(I))
      GS(I) = G(I)
  140 TRACE = TRACE + ((V(I,I)-VII(I))/(V(I,I)+VII(I)))**2
      TRACE = DSQRT(TRACE/PARN)

      FS = F
      GO TO 24
C                                        . . . . .  END MAIN LOOP

  170 write(WRTBUFFER, 500)   !EGL: convergence
      call logmessage(wrtbuffer)

      ISW2 = 3
      IF (MATGD .GT. 0)  THEN
         STATUS='converged'
         GO TO 435
      ENDIF

      NPARGD = NPAR*(NPAR+5)/2
      IF (NFCN-NPFN .GE. NPARGD)  THEN
         STATUS='converged'
         GO TO 435
      ENDIF

      write(WRTBUFFER, 180)
      call logmessage(wrtbuffer)
  180 FORMAT ('COVARIANCE MATRIX INACCURATE.  BHAT WILL ',
     1'TRY TO RECALCULATE')
C      CALL HESSE
C      CALL MPRINT(1,AMIN)
C      WRITE(WRTBUFFER,*) 'DO YOU WANT COVARIANCEES? (Y/N)'
C      READ (*,460) ICHAR
C      IF(ICHAR.EQ.'Y') CALL MATOUT(0.0D0, 1)

      IF (ISW2 .GE. 2)  THEN
         ISW2 = 3
         STATUS='not converged'
         GO TO 435
      ENDIF

 190  ISW1 = 1          !EGL: non-convergence
      GO TO 230

 200  write(WRTBUFFER, 650)
      call logmessage(wrtbuffer)
      DO 210 I= 1, NPAR
 210     X(IVN(I)) = XXS(I)
         ISW2 = 1
         WRITE(WRTBUFFER,*) 'SIGMA: ',SIGMA
         call logmessage(wrtbuffer,1)
      IF (SIGMA .LT. ROSTOP)  THEN
         WRITE(WRTBUFFER,*) 'convergence type 4'
         call logmessage(wrtbuffer,1)
         GO TO 170
      ENDIF
      IF (MATGD .GT. 0)  GO TO 2

 230  write(WRTBUFFER, 510)         !EGL: retry
      call logmessage(wrtbuffer)

      IF (ISW2 .LE. 1)  THEN    !EGL: non-convergence after retry
         STATUS='not converged'
         GO TO 435
      ENDIF
 
      write(WRTBUFFER, 511)
      call logmessage(wrtbuffer)
         STATUS='not converged'

  435 CALL BTRAFO(NDIM,X,SU) 
      CALL FUNC(SU,NDIM,AMIN); NFCN=NFCN+1

c       GRADIENT SEARCH FINAL OUTPUT:

C      CALL BTIME(STAMP)
C        WRITE(WRTBUFFER,9009) STAMP,STATUS,AMIN
        WRITE(WRTBUFFER,9009) STATUS,AMIN
        call logmessage(wrtbuffer)
        WRITE(WRTBUFFER,9010)
        call logmessage(wrtbuffer)
        J=1
        IF(IVN(J).EQ.1) THEN
          WRITE(WRTBUFFER,9020) ITER,AMIN,LABELS(1),SU(1),GS(1),G2(1),NFCN
          call logmessage(wrtbuffer)
          J=J+1
        ELSE
          WRITE(WRTBUFFER,9021) ITER,AMIN,LABELS(1),SU(1),NFCN
          call logmessage(wrtbuffer)
        ENDIF

      DO I=2,NDIM
        IF(IVN(J).EQ.I) THEN
          WRITE(WRTBUFFER,9030) LABELS(I),SU(I),GS(J),G2(J)
          call logmessage(wrtbuffer)
          J=J+1
        ELSE
          WRITE(WRTBUFFER,9031) LABELS(I),SU(I)
          call logmessage(wrtbuffer)
        ENDIF
      ENDDO
      
        OPEN(23,file='bhat.final',status='unknown')

C        WRITE(23,9009) STAMP,STATUS,AMIN
        WRITE(23,9009) STATUS,AMIN
        WRITE(23,9010)
        IV1=IVN(1)
        WRITE(23,9020) ITER,AMIN,LABELS(IV1),SU(IV1),GS(1),G2(1),NFCN
      DO I=2,NPAR
        IVI=IVN(I)
        WRITE(23,9030) LABELS(IVI),SU(IVI),GS(I),G2(I)
      ENDDO

        CLOSE(23)

      IF(IBOOT.EQ.2) THEN
      WRITE(9,6000) AMIN,(SU(IVN(I)),I=1,NPAR)
      RETURN	!GO BACK TO BOOTSTRAPLOOP
      ENDIF

      
      RETURN
  450 FORMAT('START DAVIDON-FLETCHER-POWELL ALGORITHM',/)
C  460 FORMAT(A1)
  470 FORMAT ('CONVERGENCE CRITERIA: ',
     +    ' --  ESTIMATED DISTANCE TO MINIMUM (EDM) .LT.',E9.2,/,
     +22X,' --                               OR EDM .LT.',E9.2,/,
     +22X,' --  AND FRAC. CHANGE IN VARIANCE MATRIX .LT.',E9.2,/)
  500 FORMAT ('GRADIENT SEARCH COMPLETED',/)
  510 FORMAT ('GRADIENT SEARCH NONCONVERGENT',/)
  511 FORMAT ('GRADIENT SEARCH NONCONVERGENT - STOP',/)
  520 FORMAT ('COVARIANCE MATRIX IS NOT POSITIVE-DEFINITE',/)
  650 FORMAT ('GRADIENT SEARCH FAILS TO FIND IMPROVEMENT',/)
C 5000  FORMAT(I6,30E18.8)
6000  FORMAT(F15.5,30E18.8,/)
C 7000  FORMAT(/,'DO YOU WANT GRAPHICAL INFORMATION (Y/N) ')
C 7010  FORMAT(A10)
C 8000  FORMAT('plot "tmp_" using 1:3 title "',a5,'" w l ',
C     1 10(a17,i5,a8,a5,'"',a4))
C 8001  FORMAT('"tmp_" using 1:',i5,' w l')
C 8002  FORMAT('pause -1')
C 8010  FORMAT(/,'GRAPHING  LABELS:   ',2A10)
C 8020  FORMAT(i5,5x,a5)
C 9009  FORMAT(/,'BHAT RUN: ',A28,T42,'STATUS: ',A15,T67,'VALUE:',E16.8)
9009  FORMAT('BHAT RUN STATUS: ',A15,T67,'VALUE:',E16.8,/)
9010  FORMAT('IT',T11,'VALUE',T18,'LBL',T26,
     & 'ESTIMATES',T41,
     & 'DERIVATIVES',T57,
     & 'CURVE',T67,'CALLS',/)
 9020 FORMAT(I2,T4,E12.5,T18,A5,E12.5,T40,F12.2,E10.2,I10,/)
 9021 FORMAT(I2,T4,E12.5,T18,A5,E12.5,T40,
     $       '  .... PARAMETER FIXED',T62,I10,/)
 9030 FORMAT(T18,A5,E12.5,T40,F12.2,E10.2,/)
 9031 FORMAT(T18,A5,E12.5,T40,'  .... PARAMETER FIXED',/)
      END


c$$$      SUBROUTINE BTIME(STAMP)
c$$$       CHARACTER*28 STAMP
c$$$
c$$$      IDUM=SYSTEM('date > bhat.stamp')
c$$$      OPEN(20,file='bhat.stamp',status='unknown')
c$$$      READ(20,202) STAMP
c$$$ 202  FORMAT(A28)
c$$$      CLOSE(20)
c$$$
c$$$      END SUBROUTINE BTIME
