      SUBROUTINE DSTEP(NFCN,NDIM,NPAR,X0,SENS,XINF,DSTEPS)

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=100)

      CHARACTER*10 LABELS(NMAX)
      CHARACTER RNF(NMAX)

      DIMENSION X(NMAX),X0(NMAX),DSTEPS(NMAX),IVN(NMAX)
      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)

      STEP=1.; JMAX=10
      
! OBTAIN REFERENCE POINT VALUE

      CALL BTRAFO(NDIM,X0,SU)
      CALL FUNC(SU,NDIM,F0); NFCN = NFCN + 1

      DO I=1,NPAR
         
         STEPI=STEP
         X(1:NDIM) =X0(1:NDIM)

         XL=X0(IVN(I)); XU=XL+STEP; X(IVN(I))=XU
 2       CALL BTRAFO(NDIM,X,SU)
         CALL FUNC(SU,NDIM,F); NFCN = NFCN + 1

         IF(F.EQ. nan .OR. F.EQ.-nan) THEN
            PRINT*,'oops: stepsize too large'
            STEPI=STEPI/10.
            XU=XL+STEPI; X(IVN(I))=XU
            GOTO 2
         ENDIF

         DF=DABS(F-F0)
         IF(DF.LT.SENS) GOTO 1

         X(IVN(I))=0.5*(XL+XU)
         J=0
         DO WHILE(DABS(DF-SENS).GT..05) !within a few % should be fine

         CALL BTRAFO(NDIM,X,SU)
         CALL FUNC(SU,NDIM,F); NFCN = NFCN + 1
         DF=DABS(F-F0)

         IF(DF-SENS.GT.0.) THEN
            XU=X(IVN(I)); X(IVN(I))=0.5*(XL+X(IVN(I)))
         ELSE
            XL=X(IVN(I)); X(IVN(I))=0.5*(XU+X(IVN(I)))
         ENDIF

         J=J+1
         IF(J.GT.JMAX) THEN
            PRINT*,'could not find dstep for parameter',IVN(I)
            GOTO 1
         ENDIF

         ENDDO

 1       DSTEPS(I)=X(IVN(I))-X0(IVN(I))
         PRINT*,'DSTEP:',IVN(I),DSTEPS(I),DF

       ENDDO
       RETURN
       END


