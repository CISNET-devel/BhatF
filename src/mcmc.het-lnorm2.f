c *** MARKOV-CHAIN MONTE-CARLO METHOD FOR SAMPLING LIKELIHHOODS WITH 
c *** HETEROGENEITY PARAMETERS

c     uses mvn Q sampler on original scale
c     passes heterogeneity vector lambda to FUNC

c     imcmc=0: do mcmc (assuming homogeneity)
c     imcmc=1: do heterogeneity problem via mcmc
c     note: uses log-normal hetero distributions

c     W: Hessian eigen values
c     Z: Hessian eigen vectors
 
      SUBROUTINE MCMC_HET_LNORM(NDIM,NPAR,LH,NCT,X0,F0,MEANV,COVM,X)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      
      PARAMETER(NMAX=100, XINF=16., SENS=0.1)

      CHARACTER*10 LABELS(NMAX),LABELS_ORG(NMAX),ANSWER
      CHARACTER L_PRIOR, RNF(NMAX)
      INTEGER NPAR, NCT, NRC, IWORK(NMAX), iseed1, iseed2
      INTEGER LH,IVN(NMAX),IRE(NPAR),NOBS,ID1(10),ID2(10)
      INTEGER N_ACC_LAMBDA(LH)
      REAL*8, DIMENSION(200) :: LKH,LKH0
      REAL*8, DIMENSION(10,200) :: LAMBDA,LAMBDA0
      REAL MEANV(NPAR),COVM(NPAR,NPAR),X(NPAR),WORK(NPAR)
      REAL PARM(NPAR * (NPAR+3)/2+1), genunf, gengam
      EXTERNAL genunf,gengam

      REAL*8, DIMENSION(NMAX) :: Y,X0,SU0_MC,SM
      DIMENSION SQ(NMAX),ISQ(NMAX),JSQ(NMAX),SU0(NMAX),SU00(NMAX)
      DIMENSION W(NMAX),Z(NMAX,NMAX)
      DIMENSION S(NMAX,NMAX),FWORK1(NMAX),FWORK2(NMAX)
      DIMENSION E1(NMAX),E2(NMAX),EU0(NMAX),EU(NMAX)
      DIMENSION C_SUM(NMAX,NMAX),theta(NMAX),ESU0(NMAX)

      DIMENSION DDF(NMAX,NMAX),DF(NMAX)
      DIMENSION WAREA(NMAX),BIDY(NMAX,NMAX),X00(NMAX)

      REAL*4, DIMENSION(NMAX) :: SIG,TSIG

      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /IFLAGS/INDEX,IMCMC,IBOOT,ISEED,IERR  
      COMMON /DELTA/DEL(NMAX)
      COMMON /HETERO/NOBS,ID1,ID2,LAMBDA,LKH
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)
      COMMON /PMCMC/X_PRIOR(NMAX),S_PRIOR(NMAX),L_PRIOR(NMAX)

      DATA SM/NMAX*0.D0/

!!      INTERFACE 
!!         INTEGER FUNCTION SYSTEM (COMMANDA) 
!!         CHARACTER(LEN=*) COMMANDA 
!!         END FUNCTION SYSTEM 
!!      END INTERFACE

c -----------------------------------------------------------------
!!      INCLUDE '/usr/pgi/linux86/include/lib3f.h'

!! 	INTERFACE
!! 	 EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_INI(NPAR_)
!!          INTEGER, INTENT(IN) :: NPAR_
!! 	 END SUBROUTINE GNUPLOT_INI
!! 	END INTERFACE

!! 	INTERFACE
!! 	 EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT()
!! 	 END SUBROUTINE GNUPLOT
!! 	END INTERFACE
c -----------------------------------------------------------------

C *** CYCLE THROUGH PARAMETER SPACE SEQUENTIALLY

c      PRINT*,'scale that defines width of proposal distr.?'
c      READ(5,*) DSC

      PRINT*,'training period I (recommend 1000 cycles)?'
      READ(5,*) MCT0

      PRINT*,'training period II (recommend 2000 cycles)?'
      READ(5,*) MC1

      PRINT*,'record every mth cycle:'
      READ(5,*) mth
      
      PRINT*,'ID, iseed1, iseed2:'
      READ(5,*) mc_id,iseed1,iseed2

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
      
      IF(COVM(1,1) .NE. 0.) THEN
         PRINT*,'already have COVM, '//
     $        'will directly procede to pilot run 2'
         PRINT 999
         sm=0.
         GOTO 75
      ENDIF

c --- compute variances

      CALL DQSTEP(NCT,NDIM,NPAR,X0,SENS,XINF,DEL)
        print*,'del',del(1:npar)
        MEANV(1:NPAR)=0.

 	CALL HESSIAN(NPAR,NDIM,IVN,X0,Y,DF,DDF,W,Z,NCT)

c ---  if DDF not pos definite, retain diagonal elements if positive

        IF(W(1).LE.0.) THEN
           PRINT*,'logit-Hessian not positive definite'
           PRINT*,'eigenvalues:',w(1:npar)
c           DO I=1,NPAR
c              IF(DDF(I,I).le.0.) THEN
c                 PRINT*,'curvature negative: use default 10e-6'
c                 DDF(I,I)=1.d-6
c              ENDIF
c           ENDDO
c           PRINT*,'continue to use diagonal elements only'
c           DO I=2,NPAR
c              DO J=1,I-1
c                 DDF(I,J)=0; DDF(J,I)=0
c              ENDDO
c           ENDDO
           PRINT*,'construct a proposal COVM using dqstep'
           CALL DQSTEP(NCT,NDIM,NPAR,X0,2.d0,16.d0,DEL)
           COVM=0
           DO I=1,NPAR
             COVM(I,I)=DEL(I)*DEL(I)
           ENDDO

        ELSE
                 
c ---  compute covariance using Hessian 

        BIDY=0.D0
        DO I=1,NPAR
           BIDY(I,I)=1.D0
        ENDDO
      IERR=0
      CALL AXEB(DDF,NPAR,NMAX,BIDY,NPAR,NMAX,0,WAREA,IERR)
        COVM=.5*BIDY  !use a 'bootstrapping' algorithm for adjusting scale
        
        ENDIF
c --- initialize multivariate normal sampler
        pause

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

      print 999
      print*,'first pilot chain: using MLEs'//
     $     ' and log-transformed covariance'
      print 999

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
         WRITE(6,1007) N,N_ACC,(N_ACC_LAMBDA(I),I=1,LH)
         WRITE(6,1006) F0,(SU0(IVN(I)),I=1,NPAR)
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

c -----------------------------------------------------------------------

         ENDDO !n     end of 1. training run
         
       PRINT 999

         SM=0.D0; COVM=C_SUM/FLOAT(MCT0) 
          
C ---	get eigenvalues and eigenvectors of covariance matrix

      NRC=NPAR; IERR=0
        CALL RS(NMAX,NRC,C_SUM,W,0,S,fwork1,fwork2,ie)
      PRINT*,'RS error: ',ie

      PRINT 999
      PRINT 9040,(W(I),I=1,NPAR)
      PRINT 999
c        DO I=1,NPAR
c 	 PRINT 9050,(S(I,J),J=1,NPAR)
c        ENDDO
c	 PRINT 999

c --- initialize multivariate normal sampler (on original scale)

 75     CONTINUE
c         DO I=1,NPAR
c	  PRINT 9050,(COVM(I,J)/SQRT(COVM(I,I))/
c     &   SQRT(COVM(J,J)),J=1,NPAR)
c         ENDDO
        CALL SETGMN(MEANV,COVM,NPAR,PARM)

c ---  procede with Q sampling from multivariate normal with COVM

      MC_Smax=2; MC_S=0; N_ACC=0; N_ACC_LAMBDA=0

      print 999
      print*,'second pilot chain: includes update '//
     $     'for heterog. parameters'
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

        SIG(1:LH)=SE(IRE(1:LH))

c ---  generate lambda_i's from log-normal distribution
        DO I=1,LH
           DO K=1,NOBS
           LAMBDA(I,K)=SU0(IRE(I))*DEXP(SIG(I)*GENNOR(0.,1.))
           ENDDO
        ENDDO

      CALL FUNC(SU0,NDIM,F0); NCT=NCT+1 !fills lkh's
        print*,'first (imcmc=1) function call: ',F0,sum(lkh(1:nobs))

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

        IF(IMCMC.EQ.1) THEN
        DO I=1,LH
           WSIG=EXP(SIG(I)*SIG(I))
           TSIG(I)=SU0(IRE(I))*DSQRT(WSIG*WSIG-WSIG)
        ENDDO
      WRITE(8,1001) mc_id,F0,(SU0(IVN(I)),I=1,NPAR),(TSIG(I),I=1,LH)
        ELSE
      WRITE(8,1001)  mc_id,F0,(SU0(IVN(I)),I=1,NPAR)
        ENDIF

        close(8)

c ---  reinitialize COV stuff

        print*,'scale for m.v. sampler (recommend 2.):'
        print*,'(if via global recommend 0.25)'
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
c     All parameters SU0 labeled 'H' are replaced by their respective means
c     which are Gibbs sampled.

        SU00=SU0; Y=SU0
      DO I=1,NPAR
      Y(IVN(I))=SU0(IVN(I))+scl*X(I)

        IF(RNF(IVN(I)).NE.'H') THEN
        IF(Y(IVN(I)).LT.S1(IVN(I))) THEN
           PRINT*,'down move rejected!',IVN(i)
           GOTO 14
          ELSEIF(Y(IVN(I)).GT.S2(IVN(I))) THEN
           PRINT*,'up move rejected!',IVN(i)
           GOTO 14
          ENDIF
        ENDIF

      ENDDO
        
      CALL FUNC(Y,NDIM,FP); NCT=NCT+1

        IF(IMCMC.EQ.2) THEN

           GP=0.
           DO I=1,NDIM
                 IF(RNF(I).EQ.'G') THEN
                 a_prior=(x_prior(i)/s_prior(i))**2
                 b_prior=x_prior(i)/s_prior(i)/s_prior(i)
                 GP=GP+y(ivn(i))*b_prior-(a_prior-1.)*log(y(ivn(i)))
                 ELSEIF(RNF(I).EQ.'N') THEN
                 GP=GP-((y(ivn(i))-x_prior(i))**2)/
     $                   2./s_prior(i)/s_prior(i)
                 ENDIF
           ENDDO
           
        ELSE         
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
! -----------------------------------------------------------------

        IF(IMCMC.EQ.1) THEN

c ---  step1: sample eta (log mean) | sig, lambda using Gibbs

           DO I=1,LH  !for each heterogeneous parameter
            expect_log_lambda=sum(dlog(lambda(i,1:nobs)))/nobs
            eta=gennor(real(expect_log_lambda),
     $                 real(sig(i)/sqrt(float(nobs)))) 
            su0(ire(i))=dexp(eta)

c ---  step2: sample 1/2/sig^2 from gamma distr.

           ETA=DLOG(SU0(IRE(I))) !log median, current
           RHOI=0.
           DO K=1,NOBS
              YLAMIK=DLOG(LAMBDA(I,K))
              RHOI=RHOI+(YLAMIK-ETA)**2
           ENDDO
            ETAI=gengam(real(RHOI),real(0.5*NOBS+1.))
            SIG(I)=1./SQRT(2.*ETAI)

c ---  step3: sample lambda_i's from log-normal distribution

        ACC_P=0.; IERR=0; LAMBDA0=LAMBDA

        DO K=1,NOBS
           LAMBDA(I,K)=SU0(IRE(I))*DEXP(SIG(I)*GENNOR(0.,1.))
        ENDDO

c --- compute remaining factors:
c --- get new lkh's for current point (su0) and with these proposed lambda_i's
c     for lambda's that do not change this may be simplified (???)

      CALL FUNC(SU0,NDIM,FDUM); NCT=NCT+1 

c --- use Metropolis-Hastings here to accept or reject the lambda_i's 
c     individual by individual, for a specific 'H' parameter

        N_ACC_LAMBDA(I)=0
        DO K=1,NOBS

         ACC_P=MIN(1.D0,DEXP((-LKH(K)+LKH0(K))))

         IF(ACC_P.EQ.1.D0) THEN
          LAMBDA0(I,K)=LAMBDA(I,K)
          LKH0(K)=LKH(K)
          N_ACC_LAMBDA(I)=N_ACC_LAMBDA(I)+1
       ELSE

        IF(GENUNF(0.,1.).LE.ACC_P) THEN
           LAMBDA0(I,K)=LAMBDA(I,K)
           LKH0(K)=LKH(K)
           N_ACC_LAMBDA(I)=N_ACC_LAMBDA(I)+1
        ENDIF
       ENDIF
        ENDDO
        LAMBDA=LAMBDA0; LKH=LKH0
        N_ACC_LAMBDA(I)=100*(N_ACC_LAMBDA(I)/FLOAT(NOBS)) !average acceptance per cycle
        ENDDO !I=1,LH

        F0=SUM(LKH(1:NOBS))

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

        DO I=1,LH
           WSIG=EXP(SIG(I)*SIG(I))
           PRINT*,SU0(IRE(I))*DSQRT(WSIG*WSIG-WSIG),SIG(I)
        ENDDO
        WRITE(6,1007) N,N_ACC,(N_ACC_LAMBDA(I),I=1,LH)
        WRITE(6,1006) F0,(SU0(IVN(I)),I=1,NPAR)
      N_ACC=0; N_ACC_LAMBDA=0  !resetting the acceptance counter to zero

        endif
      
c --- write

      if(mod(n,mth).eq.0) then
        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')

        IF(IMCMC.EQ.1) THEN

        DO I=1,LH
           WSIG=EXP(SIG(I)*SIG(I))
           TSIG(I)=SU0(IRE(I))*DSQRT(WSIG*WSIG-WSIG)
        ENDDO
      WRITE(8,1001) mc_id,F0,(SU0(IVN(I)),I=1,NPAR),(TSIG(I),I=1,LH)
        ELSE
      WRITE(8,1001) mc_id,F0,(SU0(IVN(I)),I=1,NPAR)
        ENDIF
        close(8)
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
     &             SQRT(COVM(J,J)),J=1,NPAR)
           ENDDO
          
c --- re-initialize multivariate normal sampler

           CALL SETGMN(MEANV,COVM,NPAR,PARM)
c --- allow for change in scale
           print*,'scale for m.v. sampler (recommend 2.):'
           read*,scl
           GOTO 66
        ELSE
           print*,'continue (y/n):'
           read(5,67) answer
 67        format(a1)

           if(answer.eq.'y'.or.answer.eq.'Y') then
              print*,'number of cycles:'
              read*,mc1
              goto 66
           endif
        ENDIF
      
      print 999
      print*,'chain: ',MC_S
      print 999

c ---  exit -- final analysis - return approx. mle

        PRINT*,' '
      
        CALL FTRAFO(NDIM,X0,SU0_MC)

        PRINT 1003
      DO I=1,NDIM
        SU(I)=SU0_MC(I)
      PRINT 1002,LABELS(I),ESU0(I)/FLOAT(NCMC),SU(I)
      ENDDO
        
        F0=F_MIN
        
        CALL SYSTEM('killall gnuplot')

      RETURN
 999    FORMAT(/)
 1000	FORMAT(I5,F12.4,30E12.4)
 1001	FORMAT(I5,30E16.8)
 1002   FORMAT(A10,2E12.5)
 1003  FORMAT('pram',t15,' mc mean',t24,' appr. mode')
 1007  FORMAT('cycle:',I6,' MVN accept:',I4,' (lambda) accept:',30I4)
 1006  FORMAT(F12.4,30E12.4)
c 1008  FORMAT('accepted: ',30I4)
 9040  FORMAT(/,' eigenvalues:  ',/,30E12.4,/)
 9050  FORMAT(30E14.6)
      END

