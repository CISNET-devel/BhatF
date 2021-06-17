        SUBROUTINE SIMPLEX(NCT,X,IV,P,Y,MP,NP,NPAR,NDIM,FTOL,ITER)
C
C ----  SIMPLE DOWNHILL SIMPLEX METHOD BY NELDER AND MEAD
C
        IMPLICIT REAL*8(A-H,O-Z)
        PARAMETER(NMAX=100,NMAX1=NMAX+1)
c	PARAMETER(ALPHA=1.0,BETA=.25,GAMMA=2.0,ITMAX=400)
C	PARAMETER(ALPHA=1.0,BETA=.50,GAMMA=2.0,ITMAX=800)
	PARAMETER(ALPHA=1.0,BETA=.50,GAMMA=1.5,ITMAX=1000)
        DIMENSION P(NMAX1,NMAX),Y(NMAX1),PR(NMAX),PRR(NMAX),PBAR(NMAX)
        DIMENSION IV(NMAX),X(NMAX)
        COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)

        MPTS=MP
        ITER=0

1       ILO=1
        IF(Y(1).GT.Y(2)) THEN
                IHI=1
                INHI=2
        ELSE
                IHI=2
                INHI=1
        ENDIF

        DO 11 I=1,MPTS
                IF(Y(I).LT.Y(ILO)) ILO=I
                IF(Y(I).GT.Y(IHI)) THEN
                        INHI=IHI
                        IHI=I
                ELSE IF(Y(I).GT.Y(INHI)) THEN
                        IF(I.NE.IHI) INHI=I
                ENDIF
11      CONTINUE

!!        RTOL=2.*DABS(Y(IHI)-Y(ILO))/(DABS(Y(IHI))+DABS(Y(ILO)))
        RTOL=DABS(Y(IHI)-Y(ILO))

        IF(RTOL.LT.FTOL) RETURN
        IF(ITER.EQ.ITMAX) THEN
           PRINT*,'exceeding ',itmax,' iterations'
           RETURN
        ENDIF
        ITER=ITER+1

        DO 12 J=1,NP
        PBAR(J)=0.
12      CONTINUE

        DO 14 I=1,MPTS
        IF(I.NE.IHI) THEN
                DO 13 J=1,NP
13              PBAR(J)=PBAR(J)+P(I,J)
        ENDIF
14      CONTINUE

        DO 15 J=1,NP
        PBAR(J)=PBAR(J)/NP
        PR(J)=(1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15      CONTINUE

        X(IV(1:NP))=PR(1:NP)
        CALL BTRAFO(NDIM,X,SU)
        CALL FUNC(SU,NDIM,YPR); NCT=NCT+1

        IF(YPR.LE.Y(ILO)) THEN

        PRR(1:NP)=GAMMA*PR(1:NP)+(1.-GAMMA)*PBAR(1:NP)

        X(IV(1:NP))=PRR(1:NP)
        CALL BTRAFO(NDIM,X,SU)
        CALL FUNC(SU,NDIM,YPRR); NCT=NCT+1

        IF(YPRR.LT.Y(ILO)) THEN

                DO 17 J=1,NP
                        P(IHI,J)=PRR(J)
17              CONTINUE
                Y(IHI)=YPRR
        ELSE
                DO 18 J=1,NP
                        P(IHI,J)=PR(J)
18              CONTINUE
                Y(IHI)=YPR
        ENDIF

        ELSE IF(YPR.GE.Y(INHI)) THEN
                IF(YPR.LT.Y(IHI)) THEN
                DO 19 J=1,NP
                        P(IHI,J)=PR(J)
19              CONTINUE
                Y(IHI)=YPR
        ENDIF

        DO 21 J=1,NP
                PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
21      CONTINUE

        X(IV(1:NP))=PRR(1:NP)
        CALL BTRAFO(NDIM,X,SU)
        CALL FUNC(SU,NDIM,YPRR); NCT=NCT+1

        IF(YPRR.LT.Y(IHI)) THEN
                DO 22 J=1,NP
                        P(IHI,J)=PRR(J)
22              CONTINUE
                Y(IHI)=YPRR
        ELSE
                DO 24 I=1,MPTS
                IF(I.NE.ILO) THEN
                        DO 23 J=1,NP
                                PR(J)=0.5*(P(I,J)+P(ILO,J))
                                P(I,J)=PR(J)
23                      CONTINUE
c                        CALL SELECT(NP,IV,PR,X)
                        X(IV(1:NP))=PR(1:NP)
                        CALL BTRAFO(NDIM,X,SU)
                        CALL FUNC(SU,NDIM,Y(I)); NCT=NCT+1
                ENDIF
24              CONTINUE
        ENDIF
        
        ELSE
        
                DO 25 J=1,NP
                        P(IHI,J)=PR(J)
25              CONTINUE
                Y(IHI)=YPR
        ENDIF
        
        GOTO 1
        END                

