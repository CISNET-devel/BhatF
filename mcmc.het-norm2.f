c *** MARKOV-CHAIN MONTE-CARLO METHOD FOR SAMPLING LIKELIHHOODS WITH 
c *** HETEROGENEITY PARAMETERS

c     uses mvn Q sampler on original scale
c     passes nobs and vector lambda to FUNC
c     needs generalization to more than one component

c     imcmc=0: do mcmc (assuming homogeneity)
c     imcmc=1: do heterogeneity problem via mcmc
c     imcmc=2: do heterogeneity problem via mcmc

      SUBROUTINE MCMC2(NDIM,NPAR,NCT,X0,F0,MEANV,COVM,X)
	
      IMPLICIT REAL*8 (A-H,O-Z)
      
      PARAMETER(NMAX=100, XINF=16.)

      CHARACTER*10 LABELS(NMAX),LABELS_ORG(NMAX),ANSWER
      CHARACTER L_PRIOR, RNF(NMAX)
      INTEGER NPAR, NCT, NRC, IWORK(NMAX), iseed1, iseed2
      INTEGER LH,IVN(NMAX),IRE(NPAR),NOBS,ID1(10),ID2(10)
      REAL*8, DIMENSION(200) :: LKH,LKH0
      REAL*8, DIMENSION(NMAX,200) :: LAMBDA,LAMBDA0,TLAM,TLAM0
      REAL MEANV(NPAR),COVM(NPAR,NPAR),X(NPAR),WORK(NPAR)
      REAL PARM(NPAR * (NPAR+3)/2+1), genunf, gengam
      EXTERNAL genunf,gengam

      REAL*8, DIMENSION(NMAX) :: Y,X0,XP,SU0_MC,SM
      DIMENSION SQ(NMAX),ISQ(NMAX),JSQ(NMAX),SU0(NMAX),SU00(NMAX)
      DIMENSION W(NMAX),S(NMAX,NMAX),FWORK1(NMAX),FWORK2(NMAX)
      DIMENSION E1(NMAX),E2(NMAX),EU0(NMAX),EU(NMAX)
      DIMENSION C_SUM(NMAX,NMAX),theta(NMAX),ESU0(NMAX)

      DIMENSION DDF(NMAX,NMAX),DF(NMAX)
      DIMENSION WAREA(NMAX),BIDY(NMAX,NMAX),X00(NMAX)

      REAL, DIMENSION(NMAX) :: XMU,SIG

      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /IFLAGS/INDEX,IMCMC,IBOOT,ISEED,IERR  
      COMMON /DELTA/DEL(NMAX)
      COMMON /HETERO/NOBS,ID1,ID2,LAMBDA,LKH
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)
      COMMON /PMCMC/X_PRIOR(NMAX),S_PRIOR(NMAX),L_PRIOR(NMAX)

      DATA SM/NMAX*0.D0/

      INCLUDE '/usr/pgi/linux86/include/lib3f.h'

 	INTERFACE
 	 EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_INI(NPAR_)
          INTEGER, INTENT(IN) :: NPAR_
 	 END SUBROUTINE GNUPLOT_INI
 	END INTERFACE

 	INTERFACE
 	 EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT()
 	 END SUBROUTINE GNUPLOT
 	END INTERFACE
c -----------------------------------------------------------------

      PRINT*,'length of cyclic M chain (variable MC)?'
      READ(5,*) MC
C *** CYCLE THROUGH PARAMETER SPACE SEQUENTIALLY

c      PRINT*,'scale that defines width of proposal distr.?'
c      READ(5,*) DSC

      PRINT*,'training period I (recommend 1000 cycles)?'
      READ(5,*) MCT0

      PRINT*,'training period II (recommend 2000 cycles)?'
      READ(5,*) MC1

      PRINT*,'record every mth cycle:'
      READ(5,*) mth
	
      PRINT*,'iseed1, iseed2:'
      READ(5,*) iseed1,iseed2

!      PRINT*,'scale up prior means:'
!      READ*,SCL

!        OPEN(7,FILE="mcmc.inp",status="unknown")
!        i=1
! 2      READ (7, *, end=55) x_prior(i),s_prior(i),l_prior(i)
!        IF(S_PRIOR(I).GT.0.) THEN
!           IMCMC=2
!           x_prior(i)=x_prior(i)*scl
!        ENDIF
!        i=i+1
!        GOTO 2

! 55   CLOSE(7)

      DO I=1,NDIM
         IF(RNF(I).EQ.'N'. OR. RNF(I).EQ.'G') THEN
            IMCMC=2
            X_PRIOR(I)=SU(I)
            S_PRIOR(I)=SE(I)
          ENDIF
      ENDDO
        
      NCMC=0; MCT=0
      F_MIN=1.D10
      ACC_P=0.D0
        
      print*,'setall'
      CALL SETALL(ISEED1,ISEED2)

      CALL BTRAFO(NDIM,X0,SU0)
      CALL FUNC(SU0,NDIM,F00); NCT=NCT+1
      F0=F00

      CALL GNUPLOT_INI(NPAR)

      N_ACC=0; N_ACC_LAMBDA=0

c --- compute variances

	DEL(1:NPAR)=.001d0  !perhaps use dstep
        MEANV(1:NPAR)=0.

 	CALL HESSIAN(NPAR,NDIM,IVN,X0,Y,DF,DDF,W,Z,NCT)

c        PRINT 9001,(DF(I),I=1,NPAR)
c        PRINT 9002,(DDF(I,I),I=1,NPAR)
c	 PRINT 1000

C ---  compute covariance 

        BIDY=0.D0
        DO I=1,NPAR
           BIDY(I,I)=1.D0
        ENDDO

	IERR=0
	CALL AXEB(DDF,NPAR,NMAX,BIDY,NPAR,NMAX,0,WAREA,IERR)

        COVM=.5*BIDY  !use a progressive algorithm for adjusting scale

c --- initialize multivariate normal sampler

        CALL SETGMN(MEANV,COVM,NPAR,PARM)

C ***   transform back to original scale (use Jacobian)
C
C	DO I=1,NPAR
C	  DYDXI=1.D0/(SU0(I)-S1(I))+1.D0/(S2(I)-SU0(I))
C	    DO J=1,NPAR	
C	      BIDY(I,J)=0.d0
C	      DYDXJ=1.D0/(SU0(J)-S1(J))+1.D0/(S2(J)-SU0(J))
C	      DDF(I,J)=DYDXI*DDF(I,J)*DYDXJ
C	    ENDDO
C	ENDDO

c --- initialize covariance matrix for first run

          SM=0.D0; C_SUM=0.D0

C --- start 1. training M chain -------------------

        DO N=1,MCT0

c --- use multivariate normal now 
c     compute width of random uniform proposal kernel

        ACC_P=0.                !by default reject
        CALL GENMN(PARM,X,WORK)

c --- compute proposal x' (=y). If unacceptable, goto 4  

        X00=X0; Y=X0
	DO I=1,NPAR
	Y(IVN(I))=X0(IVN(I))+1.*X(I)
	IF(ABS(Y(IVN(I))).GT.XINF) THEN
           print*,'mcmc close to boundary, move rejected!',ivn(i)
           GOTO 4
        ENDIF  
	ENDDO

c ---   get likelihood or energy 

	IERR=0
        CALL BTRAFO(NDIM,Y,SU)
	CALL FUNC(SU,NDIM,FP); NCT=NCT+1

c --- boundary checks from within func

        IF(IERR.EQ.1) GOTO 4 !checks on positivity problems

C ---	R ratio 

        ACC_P=MIN(1.D0,DEXP((-FP+F0)))
 4	IF(ACC_P.EQ.1.D0) THEN
	 X0=Y; F0=FP; n_acc=n_acc+1
	ELSE
	 UD01=GENUNF(0.,1.) !RAN2(ISEED)
	 IF(UD01.LE.ACC_P) THEN
	   X0=Y; F0=FP; n_acc=n_acc+1
	 ENDIF
	ENDIF

        CALL BTRAFO(NDIM,X0,SU0)

c --- update the best likelihood estimate (ble)

        IF(ACC_P.EQ.1.) THEN
         IF(FP.LT.F_MIN) THEN
            SU0_MC=SU0; F_MIN=FP
         ENDIF
        ENDIF

c ---------------------------------------------

        IF(MOD(N,100).EQ.0) THEN
         WRITE(6,1007) N,N_ACC,N_ACC_LAMBDA,F0,(SU0(IVN(I)),I=1,NPAR)
	 N_ACC=0  !resetting the acceptance counter to zero
        ENDIF
	
        CALL BTRAFO(NDIM,X00,SU00)
        CALL FTRAFO(NDIM,X0,SU0)

         DO I=1,NPAR
	  SM(I)=SM(I)+(SU0(IVN(I))-SU00(IVN(I)))
          THETA(I)=SM(I)/FLOAT(N)  !new mean
         ENDDO

         DO I=1,NPAR
          IVI=IVN(I)
          DO J=1,NPAR
           IVJ=IVN(J)
           C_SUM(I,J)=C_SUM(I,J)+((SU0(IVI)-SU00(IVI))-THETA(I))*
     &                          ((SU0(IVJ)-SU00(IVJ))-THETA(J))
	  ENDDO
         ENDDO

c ---------------------------------------

         ENDDO !n     end of 1. training run
         
	 PRINT 999

         SM=0.D0; COVM=C_SUM/FLOAT(MCT0) 

         DO I=1,NPAR
	  PRINT 9050,(COVM(I,J)/SQRT(COVM(I,I))/
     &   SQRT(COVM(J,J)),J=1,NPAR)
         ENDDO
          
C ---	get eigenvalues and eigenvectors of covariance matrix

	NRC=NPAR; IERR=0
        CALL RS(NMAX,NRC,C_SUM,W,0,S,fwork1,fwork2,ie)
	PRINT*,'RS error: ',ie

	PRINT 999
	PRINT 9040,(W(I),I=1,NPAR)
	PRINT 999
        DO I=1,NPAR
	PRINT 9050,(S(I,J),J=1,NPAR)
        ENDDO
	PRINT 999

c --- initialize multivariate normal sampler (on original scale)

        CALL SETGMN(MEANV,COVM,NPAR,PARM)

c ---  procede with Q sampling from multivariate normal with COVM

	MC_Smax=2; MC_S=0; N_ACC=0; N_ACC_LAMBDA=0

	print 999
	print*,'second pilot chain: includes update for heterog. parameters'
	print 999

c ---  check for heterog. parameters

        LH=0
        DO I=1,NDIM
         IF(RNF(I).EQ.'H') THEN    !note, if RNF==H then random eff. model  
            LH=LH+1; IRE(LH)=I; IMCMC=1
            PRINT*,'treating ',I,' as heterogeneous parameter'
         ENDIF
        ENDDO

c ---  initialization of sa and sb / here only for 1 parameter / generalize later

        IF(IMCMC.EQ.1) THEN
c --- on logit scale
        SIG(1:LH)=SE(IRE(1:LH))

c ---  generate lambda_i's, if L>0 
        DO I=1,LH
           DO K=1,NOBS
           TLAM(I,K)=GENNOR(0.,SIG(I))
           LAMBDA(I,K)=SU0(IRE(I))+TLAM(I,K)
           ENDDO
        ENDDO

	CALL FUNC(SU0,NDIM,F0); NCT=NCT+1 !fills lkh's
        print*,'first imcmc=1 function call: ',F0

        ELSEIF(IMCMC.EQ.2) THEN

	CALL FUNC(SU0,NDIM,F0); NCT=NCT+1
        print*,'first imcmc=2 function call: ',F0
           G0=0.
           DO I=1,NDIM
                 IF(RNF(I).EQ.'G') THEN
                 a_prior=(x_prior(i)/s_prior(i))**2
                 b_prior=x_prior(i)/s_prior(i)/s_prior(i)
                 G0=G0+su0(i)*b_prior-(a_prior-1.)*log(su0(i))
                 ELSEIF(RNF(I).EQ.'N') THEN
                 G0=G0-((su0(i)-x_prior(i))**2)/2./s_prior(i)/s_prior(i)
                 ENDIF
           ENDDO

        ELSE

	CALL FUNC(SU0,NDIM,F0); NCT=NCT+1
        print*,'first imcmc=0 function call: ',f0

        ENDIF

        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
        OPEN(9,FILE='mcmc.het.out',STATUS='unknown',POSITION='append')

        IF(IMCMC.EQ.1) THEN
	WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
        WRITE(9,1001) (SU0(IRE(I)),I=1,LH), (SIG(I),I=1,LH)
        ELSE
	WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
        ENDIF

        close(8); close(9)

c ---  reinitialize COV stuff

        print*,'scale for m.v. sampler (recommend 2.):'
        read*,scl

        CALL GNUPLOT()

 66     SM=0.D0; C_SUM=0.D0

c =============================================================================	
C --- start revised chain ------------------------------------------------------ 
c =============================================================================	

        DO N=1,MC1

        ACC_P=0.D0; IERR=0
 	CALL GENMN(PARM,X,WORK)

        IF(IMCMC.EQ.1) lkh0=lkh

c --- step 1 (theta's) / always carried out

c --- propose new parameter vector (Y) from mvn distr. on orig. scale
c     all parameters that are 'not stopped' are used
c     compute proposal x' (=y). If unacceptable, goto 14  

        SU00=SU0; Y=SU0
	DO I=1,NPAR
	Y(IVN(I))=SU0(IVN(I))+scl*X(I)

	 IF(Y(IVN(I)).LT.S1(IVN(I))) THEN
           PRINT*,'down move rejected!',IVN(i)
           GOTO 14
         ELSEIF(Y(IVN(I)).GT.S2(IVN(I))) THEN
           PRINT*,'up move rejected!',IVN(i)
           GOTO 14
         ENDIF

	ENDDO
        
c --- compute mean+random deviation (mu+lambda) using proposals on means

        IF(IMCMC.EQ.1) THEN

        DO I=1,LH    !on logit scale
           LAMBDA(I,1:NOBS)=Y(IRE(I))+TLAM(I,1:NOBS)
        ENDDO

	CALL FUNC(Y,NDIM,FP); NCT=NCT+1  ! uses lambda

        G0=0.; GP=0.
            
        ELSEIF(IMCMC.EQ.2) THEN

	   CALL FUNC(Y,NDIM,FP); NCT=NCT+1
           GP=0.
           DO I=1,NDIM
                 IF(RNF(I).EQ.'G') THEN
                 a_prior=(x_prior(i)/s_prior(i))**2
                 b_prior=x_prior(i)/s_prior(i)/s_prior(i)
                 GP=GP+y(ivn(i))*b_prior-(a_prior-1.)*log(y(ivn(i)))
                 ELSEIF(RNF(I).EQ.'N') THEN
                 GP=GP-((y(ivn(i))-x_prior(i))**2)/2./s_prior(i)/s_prior(i)
                 ENDIF
           ENDDO
           
        ELSE         
	   CALL FUNC(Y,NDIM,FP); NCT=NCT+1
           G0=0.; GP=0.
        ENDIF
         
        Fnew=FP+GP; Fold=F0+G0   !check sign

c --- boundary checks from within subroutine func()
        IF(IERR.EQ.1) GOTO 14 !checks on positivity problems etc.

C ---	R ratio ---------------------------------------------------

        ACC_P=MIN(1.D0,DEXP((-Fnew+Fold)))

 14	IF(ACC_P.EQ.1.D0) THEN
	 SU0=Y; F0=FP; G0=GP; lkh0=lkh; n_acc=n_acc+1
	ELSE
	 IF(genunf(0.,1.).LE.ACC_P) THEN
	   SU0=Y; F0=FP; G0=GP; lkh0=lkh; n_acc=n_acc+1
	 ENDIF
	ENDIF
c -----------------------------------------------------------------

        IF(IMCMC.EQ.1) THEN
            
c ---  step 2: generate new lambda_i's | theta, mu and sigma

        ACC_P=0.; IERR=0; LAMBDA0=LAMBDA; TLAM0=TLAM

c ---  generate lambda_i's, if L>0 
        DO I=1,LH
           DO K=1,NOBS
           TLAM(I,K)=GENNOR(0.,SIG(I))
           LAMBDA(I,K)=SU0(IRE(I))+TLAM(I,K)
           IF(LAMBDA(I,K).LE.S1(IRE(I))) THEN
!              print*,'lambda down move rejected',i,k
              LAMBDA(I,K)=LAMBDA0(I,K)
           ENDIF
           IF(LAMBDA(I,K).GE.S2(IRE(I))) THEN
!              print*,'lambda up move rejected',i,k
              LAMBDA(I,K)=LAMBDA0(I,K)
           ENDIF
           ENDDO
        ENDDO

c --- because we sample on logit scale the generated lambdas stay within their boundaries.

c --- compute remaining factors:
c --- get new lkh's for current point (su0) and with these proposed lambda_i's
c     for lambda's that donnot change this may be simplified

	CALL FUNC(SU0,NDIM,FDUM); NCT=NCT+1 

c --- use Metropolis-Hastings here to accept or reject the lambda_i's 
c     individual by individual, for the parameter block (H parameters)
        N_ACC_LAMBDA=0  !needs to be adjusted for rejected lambdas
        DO K=1,NOBS

         ACC_P=MIN(1.D0,DEXP((-LKH(K)+LKH0(K))))

         IF(ACC_P.EQ.1.D0) THEN
          LAMBDA0(1:LH,K)=LAMBDA(1:LH,K); TLAM0(1:LH,K)=TLAM(1:LH,K)
          LKH0(K)=LKH(K)
          N_ACC_LAMBDA=N_ACC_LAMBDA+1
	 ELSE

	  IF(GENUNF(0.,1.).LE.ACC_P) THEN
           LAMBDA0(1:LH,K)=LAMBDA(1:LH,K);TLAM0(1:LH,K)=TLAM(1:LH,K) 
           LKH0(K)=LKH(K)
           N_ACC_LAMBDA=N_ACC_LAMBDA+1
	  ENDIF
	 ENDIF
        ENDDO
        LAMBDA=LAMBDA0; TLAM=TLAM0; LKH=LKH0

        N_ACC_LAMBDA=100*(N_ACC_LAMBDA/FLOAT(NOBS)) !average acceptance per cycle
        F0=SUM(LKH(1:NOBS))

c ---  step3: sample eta(=1/2/sig^2) from gamma distr.

        DO I=1,LH  !for each heterogeneous parameter
            TLAM2_SUM=SUM(TLAM(I,1:NOBS)**2)
            eta=gengam(real(tlam2_sum),real(0.5*nobs+1.))
            sig(i)=1./sqrt(2.*eta)
        ENDDO

        ENDIF                   !imcmc
! -----------------------------------------------------------------

        ESU0=ESU0+SU0
        NCMC=NCMC+1

c --- pursue minimum neg. log likelihood (approx. mle's)

        IF(ACC_P.EQ.1.) THEN
         IF(F0.LT.F_MIN) THEN
           SU0_MC=SU0; F_MIN=F0
         ENDIF
        ENDIF

c ---------------------------------------------

        if(mod(n,100).eq.0) then

c	DO I=1,NPAR
c	if(a(i).lt.25. and. n.le.500) w(i)=.25*w(i) !balancing acceptance rates
c	if(a(i).gt.75. and. n.le.500) w(i)= 4.*w(i) !balancing acceptance rates
c	ENDDO

        WRITE(6,1001) (SU0(IRE(I)),I=1,LH)
        WRITE(6,1001) (SIG(I),I=1,LH)
        WRITE(6,1007) N,N_ACC,N_ACC_LAMBDA,F0,(SU0(IVN(I)),I=1,NPAR)
	N_ACC=0; N_ACC_LAMBDA=0  !resetting the acceptance counter to zero

        endif
	
c --- write

	if(mod(n,mth).eq.0) then
        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
        OPEN(9,FILE='mcmc.het.out',STATUS='unknown',POSITION='append')

        IF(IMCMC.EQ.1) THEN

	WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
        WRITE(9,1001) (SU0(IRE(I)),I=1,LH), (SIG(I),I=1,LH)
        ELSE
	WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
        ENDIF
        close(8); close(9)
	endif

c --- monitor covariance of 'increment' distribution

        IF(MC_S.LT.MC_Smax) THEN

         DO I=1,NPAR
	  SM(I)=SM(I)+(SU0(IVN(I))-SU00(IVN(I)))
          THETA(I)=SM(I)/FLOAT(N)  !new mean
         ENDDO

         DO I=1,NPAR
          IVI=IVN(I)
          DO J=1,NPAR
           IVJ=IVN(J)
          C_SUM(I,J)=C_SUM(I,J)+((SU0(IVI)-SU00(IVI))-theta(I))*
     &                          ((SU0(IVJ)-SU00(IVJ))-theta(J))
	  ENDDO
         ENDDO

        ENDIF
c ---------------------------------------

	ENDDO !n

	MC_S=MC_S+1

        IF(MC_S.LT.MC_Smax) THEN

	 PRINT 999
         SM=0.d0; COVM=C_SUM/FLOAT(MC1) 

         DO I=1,NPAR
	  PRINT 9050,(COVM(I,J)/SQRT(COVM(I,I))/
     &   SQRT(COVM(J,J)),J=1,NPAR)
         ENDDO
          
c --- re-initialize multivariate normal sampler

        CALL SETGMN(MEANV,COVM,NPAR,PARM)

        ENDIF
	
	print 999
	print*,'chain: ',MC_S
	print 999

	MC1=MC

c --- continue or exit

	IF(MC_S.LT.MC_Smax) THEN

c --- allow for change in scale
        print*,'scale for m.v. sampler (recommend 2.):'
        read*,scl
        GOTO 66

        ELSE
        print*,'continue (y/n):'
        read(5,67) answer
        ENDIF
 67     format(a1)

        if(answer.eq.'y'.or.answer.eq.'Y') then
        print*,'number of cycles:'
        read*,mc1
        goto 66
        endif

c ---  final analysis - return approx. mle

        PRINT*,' '
	
        CALL FTRAFO(NDIM,X0,SU0_MC)

        PRINT 1003
	DO I=1,NDIM
        SU(I)=SU0_MC(I)
	PRINT 1002,LABELS(I),ESU0(I)/FLOAT(NCMC),SU(I)
	ENDDO
        
        F0=F_MIN
        
        IDUM=SYSTEM('killall gnuplot')

	RETURN
 999    FORMAT(/)
 1000	FORMAT(I5,F12.4,30E12.4)
 1001	FORMAT(30E16.8)
 1002   FORMAT(A10,2E12.5)
 1003  FORMAT('pram',t15,' mc mean',t24,' appr. mode')
 1007  FORMAT(3I5,F12.4,30E12.4)
c 1008  FORMAT('accepted: ',30I4)
 9040  FORMAT(/,' eigenvalues:  ',/,30E12.4,/)
 9050  FORMAT(30E14.6)
	END

