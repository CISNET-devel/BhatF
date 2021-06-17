C ********************************************************************** 
C                                           
C                                           
C
C       PROGRAM:     PLKHCI.F
C                    MODIFIED NEWTON-RAPHSON OPTIMIZER
C                    TO OBTAIN PROFILE-LIKELIHOOD-BASED CONFIDENCE INTERVALS
C
C       METHOD:      D.J. VENZON and S.H. MOOLGAVKAR
C                    " A METHOD FOR COMPUTING PROFILE-LIKELIHOOD-BASED
C                      CONFIDENCE INTERVALS "
C                    in THE JOURNAL OF THE ROYAL STATISTICAL SOCIETY SERIES C
C                       VOLUME 37, NO.1, 1988, PP. 87-94
C       
C       DATE:        FHCRC/PHS OCT.17TH 1988
C
C       USAGE:       Warning!!!! this is a rough version
C                    NOT EDITED YET!
C                    
C                    NB identifies which parameter (among NPAR) is to be used for profile
C                    likelihood calc. 
C     
C                    NDIM:	total number of parameters in model
C                    NPAR:      subset of parameters (NDIM>NPAR)
C
C
C *********************** SPECIFICATIONS *************************************


      SUBROUTINE PLKHCI(X0_MLE,F_MLE,NPAR,NDIM,NCT)

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(NMAX=100,NMAX2=NMAX*NMAX)
      PARAMETER (SENS=0.01,XINF=32.)

      CHARACTER*10 LABELS(NMAX),FLBL
      CHARACTER*14 SUU_FLAG,SUL_FLAG
      CHARACTER*10 ANSWER
      CHARACTER RNF(NMAX)

      DIMENSION Y(NMAX),DF(NMAX),DDF(NMAX,NMAX),Z(NMAX,NMAX)
     &,FDD(NMAX,NMAX),DISC(NMAX),WAREA(NMAX),W(NMAX),X0(NMAX),IVN(NMAX)
     &,BIDY(NMAX,NMAX),X0_MLE(NMAX),DDL(NMAX,NMAX),DBO(NMAX),DBL(NMAX)
     &,RHS(NMAX),SU_MLE(NMAX),DISC0(NMAX),X00(NMAX)

      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /IFLAGS/INDEX,IMCMC,IBOOT,ISEED,IERR 
      COMMON /DELTA/DEL(NMAX)
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)

      DATA DF/NMAX*0.D0/ DDF/NMAX2*1.D3/ BIDY/NMAX2*0.d0/

      ND2=2*NPAR
      ND1=NPAR-1
      NLM=(NPAR*NPAR+3*NPAR)/2

C --------------------------------------------------------------------------
      WRITE(6, FMT='(A54)',ADVANCE='YES') 'Bhat will attempt to compute profile likelihood bounds'
 24   WRITE(6, FMT='(A7)', ADVANCE='NO') 'Label: '
      READ(5, FMT='(A10)', ADVANCE='YES', ERR=25) FLBL

      IF(FLBL.EQ.'') GOTO 25
           DO II=1,NPAR
              IF(FLBL.EQ.LABELS(IVN(II))) THEN
                 WRITE(6,FMT='(A15,2X,A10)') 'using parameter',FLBL
                 inb=ii; exit
              ENDIF
           ENDDO

        PRINT*,' assuming max. log. likelihood: ',F_MLE
	
        SUU_FLAG='non-converged'; SUL_FLAG='non-converged'

	CALL BTRAFO(NDIM,X0_MLE,SU_MLE)

        CALL DQSTEP(NCT,NDIM,NPAR,X0_MLE,SENS,XINF,DEL)

 	CALL HESSIAN(NPAR,NDIM,IVN,X0_MLE,Y,DF,DDF,W,Z,NCT)

c *** shoot along tangent to initialize first step

c ****  construct DDL

	IE=0
	DO 82 I=1,NPAR
	I1=I-IE

		IF(I.EQ.INB) THEN
		IE=1
		GOTO 82
		ENDIF

		JE=0
		DO 84 J=1,NPAR
		J1=J-JE

			IF(J.EQ.INB) THEN
			JE=1
			GOTO 84
			ENDIF
	
	DDL(I1,J1)=DDF(I,J)
 84             CONTINUE  
 82   CONTINUE  

c ****  construct DBO

	J=1
	DO I=1,NPAR
	IF(I.NE.INB) THEN
	DBO(J)=DDF(INB,I)
	J=J+1
	ENDIF
	ENDDO

C ****  construct DBB

	DBB=DDF(INB,INB)

C ***   matrix inversion: DDL inverse

        BIDY=0.D0
        DO I=1,ND1
	BIDY(I,I)=1.D0
        ENDDO
        FDD=DDL

	CALL AXEB(FDD,ND1,NMAX,BIDY,ND1,NMAX,0,WAREA,IE) !check warea

        SUM_I=0.D0; DBL=0.D0
        DO I=1,ND1
         DO J=1,ND1
          DBL(I)=DBL(I)+BIDY(I,J)*DBO(J)
         ENDDO
         SUM_I=SUM_I+DBO(I)*DBL(I)
        ENDDO

c      PRINT*,'value for stop criterium (recommend .0001)?'
c      READ*,EPS
      eps=.001
c      PRINT*,'step reduction (recommend .2)?'
c      READ*,RNEPS
        PRINT*,'No. of iterations? (10 recommended)'
        READ*,N_max


        DO 98 IDIR=1,-1,-2
           
        X0=X0_MLE   
        X0(IVN(INB))=X0(IVN(INB))+IDIR*DSQRT(3.84D0/(DBB-SUM_I))/2.D0
	j=1
	do i=1,npar
	if(i.ne.inb) then
        X0(IVN(I))=X0(IVN(I))-IDIR*DBL(J)*DSQRT(3.84D0/(DBB-SUM_I))/2.D0
	j=j+1
        endif
	enddo

	CALL BTRAFO(NDIM,X0,SU)
        CALL FUNC(SU,NDIM,F0); NCT=NCT+1

c         PRINT*,'      modified    mle       f0: ',F0
c	DO I=1,NPAR
c         write(6,8071) IVN(I),SU(IVN(I)),SU_MLE(IVN(I))
c	ENDDO

c *** search via Newton-Raphson

 99     PRINT*,' '
        IF(IDIR==-1) THEN
           PRINT*,'begin Newton-Raphson search for lower confidence bound'
        ELSE
           PRINT*,'begin Newton-Raphson search for upper confidence bound'
        ENDIF

      PRINT*,' '
  
        DO N=1,N_max

 	CALL HESSIAN(NPAR,NDIM,IVN,X0,Y,DF,DDF,W,Z,NCT)

c *** compute DDL 

	IE=0
	DO 92 I=1,NPAR
	I1=I-IE

		IF(I.EQ.INB) THEN
		IE=1
		GOTO 92
		ENDIF

		JE=0
		DO 94 J=1,NPAR
		J1=J-JE

			IF(J.EQ.INB) THEN
			JE=1
			GOTO 94
			ENDIF
	
	DDL(I1,J1)=DDF(I,J)
 94             CONTINUE  
 92   CONTINUE 
 
c ****  construct DBO

	J=1
	DO I=1,NPAR
	IF(I.NE.INB) THEN
	DBO(J)=DDF(INB,I)
	J=J+1
	ENDIF
	ENDDO

c ***  construct matrix FDD

        BIDY=0.D0
        DO I=1,NPAR
           BIDY(I,I)=1.D0
        ENDDO

        FDD(1,1)=-DF(INB)
        RHS(1) =(F_MLE-F0+1.92D0) 
        J=2
        DO I=1,NPAR               !first row
           IF(I.NE.INB) THEN
              FDD(1,J)=-DF(I)
              RHS(J)=DF(I)
              J=J+1
           ENDIF
        ENDDO
           
        DO 202 I=2,NPAR              !lower left block
        DO 202 J=2,NPAR
 202    FDD(I,J)=DDL(I-1,J-1)
      
        DO I=2,NPAR
           FDD(I,1)=DBO(I-1)
        ENDDO

	CALL AXEB(FDD,NPAR,NMAX,BIDY,NPAR,NMAX,0,WAREA,IE) !check warea

c *** shoot ahead

        RNEPS=1.; IF(N.LE.2) RNEPS=.2
        RHS0=DABS(RHS(1)); NEPSC=0; DISC=0.d0; X00=X0
        DO J=1,NPAR
         DO K=1,NPAR
          DISC(J)=DISC(J)+BIDY(J,K)*RHS(K)
         ENDDO
        ENDDO
        DISC0=DISC
        RHS1_MIN=1.D6

 100    DISC=RNEPS*DISC0

C       PRINT 8050
C       PRINT 8060,DISC(J),DUM

        NZ=0
        X0(IVN(INB))=X00(IVN(INB))+DISC(1)
        IF(DABS(DISC(1)).LE.EPS) NZ=NZ+1

        J=2
        DO I=1,NPAR
        IF(I.NE.INB) THEN   
        X0(IVN(I))=X00(IVN(I))+DISC(J)
        IF(DABS(DISC(J)).LE.EPS) NZ=NZ+1
        J=J+1
        ENDIF
        ENDDO

C *** update f0

	CALL BTRAFO(NDIM,X0,SU)
        CALL FUNC(SU,NDIM,F0); NCT=NCT+1
        RHS1=DABS(F_MLE-F0+1.92D0)

        IF(RHS1.GT.3.84 .AND. NEPSC.LT.11) THEN
           RNEPS=.75*RNEPS
           PRINT*,'reducing step size:',RNEPS,RHS1 
           NEPSC=NEPSC+1
           IF(RHS1 <= RHS1_MIN) THEN
              RNEPS_MIN=RNEPS; RHS1_MIN=RHS1
           ENDIF
           GOTO 100
c        ELSEIF(RHS0.GE.RHS1. AND. NEPSC.EQ.0) THEN
c           RNEPS=RNEPS**.5
         ELSEIF(RHS1.GT.3.84 .AND. NEPSC.EQ.11) THEN
            PRINT*,'problem closing in, will continue with min distance'  
            DISC=RNEPS_MIN*DISC0
            NZ=0
            X0(IVN(INB))=X00(IVN(INB))+DISC(1)
            IF(DABS(DISC(1)).LE.EPS) NZ=NZ+1

            J=2
            DO I=1,NPAR
               IF(I.NE.INB) THEN   
                  X0(IVN(I))=X00(IVN(I))+DISC(J)
                  IF(DABS(DISC(J)).LE.EPS) NZ=NZ+1
                  J=J+1
               ENDIF
            ENDDO

            CALL BTRAFO(NDIM,X0,SU)
            CALL FUNC(SU,NDIM,F0); NCT=NCT+1

        ENDIF

C ***   STOP-CRITERION FOR ITERATION
        IF(NZ.EQ.NPAR .AND. NEPSC.EQ.0) THEN

c -----------------------------------------------------------------------
        PRINT*,'CONVERGENCE: ',SU(IVN(INB))

        WRITE(6,9010)
        WRITE(6,9020) N,(F0-F_MLE),LABELS(IVN(1)),SU(IVN(1)),X0(IVN(1)),
     &  DF(1),DDF(1,1),NCT

	DO I=2,NPAR
           IVI=IVN(I)
      write(6,9030) LABELS(IVI),SU(IVI),X0(IVI),DF(I),DDF(I,I)
	ENDDO
        IF(IDIR== 1) THEN
           SUU=SU(IVN(INB)); SUU_flag='converged'
        ELSEIF(IDIR==-1) THEN
           SUL=SU(IVN(INB)); SUL_flag='converged'
        ENDIF

        GOTO 98

        END IF
c -----------------------------------------------------------------------

        WRITE(6,9010)
        WRITE(6,9020) N,(F0-F_MLE),LABELS(IVN(1)),SU(IVN(1)),X0(IVN(1)),
     &  DF(1),DDF(1,1),NCT

	DO I=2,NPAR
           IVI=IVN(I)
      write(6,9030) LABELS(IVI),SU(IVI),X0(IVI),DF(I),DDF(I,I)
	ENDDO
        PRINT*,' '
        PRINT*,'Confidence Bound Estimate: ', SU(IVN(INB))
        PRINT*,' '

        ENDDO !N=

        PRINT*,'another Newton-Raphson run?'
	READ(5,7010) ANSWER

	IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y') THEN
           GOTO 99
        ELSE
           IF(IDIR== 1) THEN
              SUU=SU(IVN(INB)); SUU_flag='non-converged'
           ELSEIF(IDIR==-1) THEN
              SUL=SU(IVN(INB)); SUL_flag='non-converged'
           ENDIF
           GOTO 98
        ENDIF

 98     CONTINUE
        WRITE(6, FMT='(/,A37,A10)',ADVANCE='YES') 'ESTIMATED PROFILE LIKELIHOOD BOUNDS: ',FLBL 
        WRITE(6, FMT='(A20,E16.8,2x,A13)',ADVANCE='YES') 'lower plkh bound: ',sul,sul_flag
        WRITE(6, FMT='(A20,E16.8,2x,A13)',ADVANCE='YES') 'upper plkh bound: ',suu,suu_flag

        GOTO 24

 25     RETURN
7010	FORMAT(A10)
8071    FORMAT(I5,2F12.6)
9010  FORMAT(/,' IT',T11,'VALUE',T18,'LBL',T26,'ESTIMATES',T39,
     & 'LOG-ESTIMATES',T61,'DERIVATIVES',T77,'CURVE',T87,'CALLS',/)
9020  FORMAT(I3,T5,F11.5,T18,A5,E12.5,T40,E12.5,T60,F12.2,E10.2,I10)
9030  FORMAT(T18,A5,E12.5,T40,E12.5,T60,F12.2,E10.2)
        END



	
