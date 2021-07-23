      SUBROUTINE DQSTEP(NFCN,NDIM,NPAR,X0,SENS,XINF,DSTEPS)

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=100)

      CHARACTER*10 LABELS(NMAX)
      CHARACTER RNF(NMAX)

      DIMENSION X(NMAX),X0(NMAX),DSTEPS(NMAX),IVN(NMAX)
      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)

      STEP=.01D0; DSTEPS=.01D0
      
! OBTAIN REFERENCE POINT VALUE

      CALL BTRAFO(NDIM,X0,SU)
      CALL FUNC(SU,NDIM,F0); NFCN = NFCN + 1

      DO I=1,NPAR
         
         STEPI=STEP
 2       X(1:NDIM)=X0(1:NDIM)

         X(IVN(I))=X0(IVN(I))-STEPI; X1=-STEPI
         CALL BTRAFO(NDIM,X,SU)
         CALL FUNC(SU,NDIM,F1); NFCN = NFCN + 1

         X(IVN(I))=X0(IVN(I))+STEPI; X2= STEPI
         CALL BTRAFO(NDIM,X,SU)
         CALL FUNC(SU,NDIM,F2); NFCN = NFCN + 1

         IF(F2.EQ. nan .OR. F2.EQ.-nan) THEN
           PRINT*,'NaNs: reducing step size ',I
            STEPI=STEPI/10.
            GOTO 2
         ENDIF

         IF(DABS(F2-F0).GT.5*SENS.OR.DABS(F1-F0).GT.5*SENS) THEN
c           PRINT*,'reducing step size ',I
            STEPI=STEPI/10.
            GOTO 2
         ENDIF

         B=((F1-F0)*X2-(F2-F0)*X1)/(X1*X1*X2-X2*X2*X1)
         A=(F1-F0)/X1-B*X1

c --- roots
         XS1=0.5*(-A-DSQRT(A*A+4*B*SENS))/B
         XS2=0.5*(-A+DSQRT(A*A+4*B*SENS))/B

         IF(DABS(XS1).LE.DABS(XS2)) THEN
            XS=XS1
         ELSE
            XS=XS2
         ENDIF
         
c --- see where we end up
         X(IVN(I))=X0(IVN(I))+XS
         CALL BTRAFO(NDIM,X,SU)
         CALL FUNC(SU,NDIM,F2); NFCN = NFCN + 1

         IF(F2.EQ. nan .OR. F2.EQ.-nan) THEN
            PRINT*,'oops: unable to find stepsize, use default ',I
            CYCLE
         ELSE

            DSTEPS(I)=XS
            PRINT*,'DSTEP:',IVN(I),XS,(F2-F0)
         ENDIF

       ENDDO
       RETURN
       END


