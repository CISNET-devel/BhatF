        SUBROUTINE HESSIAN(NPAR,NDIM,IVN,X0,Y,DF,DDF,W,Z,NCT)

C       NDIM:           TOTAL NO. OF PARAMETER IN PROBLEM
C       NPAR:           SUBSET OF PARAMETERS OF INTEREST
C       NCT:            COUNTS FUNCTION CALLS
C !! IVN should be passed through common

        IMPLICIT REAL*8(A-H,O-Z)
        PARAMETER(NMAX=100)
        INTEGER NRC,NM,IWORK(NMAX)

        DIMENSION X0(NMAX),Y(NMAX),XN(NMAX,2*NMAX*NMAX)
        DIMENSION DF(NMAX),DDF(NMAX,NMAX),FDD(NMAX,NMAX),F(2*NMAX*NMAX) 
        DIMENSION W(NMAX),Z(NMAX,NMAX),IVN(NMAX)
        DIMENSION FWORK1(NMAX), FWORK2(NMAX)

        COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)
        COMMON /DELTA/DEL(NMAX)

        NLM=2*NPAR*NPAR
        NP2=2*NPAR

        CALL BTRAFO(NDIM,X0,SU)
        CALL FUNC(SU,NDIM,F0); NCT=NCT+1
        PRINT*,' '
        PRINT*,'INTERNAL HESSIAN NOW COMPUTED AT FCN VALUE: ',F0

        DO 20 IN=1,NPAR
        DO 20 JN=1,NLM
20      XN(IN,JN)=X0(IVN(IN))

C ***   ON AXES

        DO 30 I=1,NPAR

                J=2*I
                K=J-1
                XN(I,J)=XN(I,J)-DEL(I)  !POINT IM
30              XN(I,K)=XN(I,K)+DEL(I)  !POINT IP

C ***   OFF AXES
        IF(NPAR.GE.2) THEN
        MC=NP2+1
        DO 40 I=2,NPAR
                DO 40 J=1,I-1
                XN(I,MC)=XN(I,MC)+DEL(I)        !POINT MC1
                XN(J,MC)=XN(J,MC)+DEL(J)
              MC=MC+1
                XN(I,MC)=XN(I,MC)-DEL(I)        !POINT MC2
                XN(J,MC)=XN(J,MC)+DEL(J)
              MC=MC+1
                XN(I,MC)=XN(I,MC)-DEL(I)        !POINT MC3
                XN(J,MC)=XN(J,MC)-DEL(J)
              MC=MC+1
                XN(I,MC)=XN(I,MC)+DEL(I)        !POINT MC4
                XN(J,MC)=XN(J,MC)-DEL(J)
40            MC=MC+1
        ENDIF

        Y=X0
        DO  I=1,NLM
        Y(IVN(1:NPAR))=XN(1:NPAR,I)
        CALL BTRAFO(NDIM,Y,SU)
        CALL FUNC(SU,NDIM,F(I)); NCT=NCT+1
        ENDDO

C        GEL2=2.D0*GEL
C        GELQ=GEL*GEL                                          

C ***   FIRST AND DIAGONAL SECOND DERIVATIVES 
C
        DO 70 I=1,NPAR
        I2=I*2
        I1=I2-1
        DF(I)=(F(I2)-F(I1))/2.D0/DEL(I)
70      DDF(I,I)=(F(I2)+F(I1)-2.D0*F0)/DEL(I)/DEL(I)

C ***   SECOND OFFDIAGONAL DERIVATIVES

        IF(NPAR.GE.2) THEN
        MC=NP2+1
        DO 80 I=2,NPAR
                DO 80 J=1,I-1
                MC1=MC
                MC2=MC1+1
                MC3=MC2+1
                MC4=MC3+1
           DDF(I,J)=(F(MC1)+F(MC3)-F(MC2)-F(MC4))/DEL(I)/DEL(J)/4.D0
                DDF(J,I)=DDF(I,J)
80              MC=MC4+1
        ENDIF

!       DO I=1,NPAR
!       IF(DDF(I,I).GT.0.D0) THEN
!       DEL(I)=.002D0/DSQRT(DDF(I,I))
!       ELSE
!       PRINT*,'problem in asigning a new delta(i)',DDF(I,I)
!       ENDIF
!       ENDDO

C ***   CHECK EIGENVALUES

        NRC=NPAR
        NM=NMAX

        FDD=DDF
        CALL RS(NM,NRC,FDD,W,0,Z,fwork1,fwork2,ie)

        PRINT 1000
        IF(W(1).LE.0.)
     $       PRINT*,'Hessian not positive definite'//
     $              ' or near singular'
c       PRINT 9000,(W(I),I=1,NPAR)
c       PRINT 9050,((Z(J,I),J=1,NPAR),I=1,NPAR)
        RETURN
1000    FORMAT(/)
9000    FORMAT(' INTERNAL EIGENVALUES:',/,30E12.4,/)
9050    FORMAT(30E12.4)
        END

