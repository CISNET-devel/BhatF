c *** MARKOV-CHAIN MONTE-CARLO METHOD FOR SAMPLING LIKELIHHOODS WITH 
c *** HETEROGENEITY PARAMETERS

c     uses mvn Q sampler on original scale
c     passes nobs and vector lambda to FUNC
c     needs generalization to more than one component

c     imcmc=0: do mcmc (assuming homogeneity)
c     imcmc=1: do heterogeneity problem via mcmc
c     imcmc=2: do heterogeneity problem via mcmc

c     note alngam requires real*8 arguments

      SUBROUTINE MCMC2(NDIM,NPAR,NCT,X0,F0,MEANV,COVM,X)
	
      IMPLICIT REAL*8 (A-H,O-Z)
      
      PARAMETER(NMAX=100, XINF=16.)

      CHARACTER*10 LABELS(NMAX),ANSWER
      CHARACTER L_PRIOR, RNF(NMAX)
      INTEGER NPAR, NCT, NRC, IWORK(NMAX), iseed1, iseed2
      INTEGER LH,IRE(NPAR),NOBS,ID1(10),ID2(10),IREJECT(200)
      REAL*8, DIMENSION(200) :: LKH,LKH0
      REAL*8, DIMENSION(10,200) :: LAMBDA,LAMBDA0
      REAL MEANV(NPAR),COVM(NPAR,NPAR),X(NPAR),WORK(NPAR)
      REAL PARM(NPAR * (NPAR+3)/2+1), genunf, gengam
      EXTERNAL genunf,gengam,alngam

      REAL*8, DIMENSION(NMAX) :: Y, X0, SU0_MC, SM, XM
      INTEGER, DIMENSION(NMAX) :: IVN, JSQ, ISQ
      DIMENSION SQ(NMAX),SU0(NMAX),SU00(NMAX)
      DIMENSION W(NMAX),S(NMAX,NMAX),FWORK1(NMAX),FWORK2(NMAX)
      DIMENSION E1(NMAX),E2(NMAX),EU0(NMAX),EU(NMAX)
      DIMENSION C_SUM(NMAX,NMAX),theta(NMAX),ESU0(NMAX)

      DIMENSION DDF(NMAX,NMAX),DF(NMAX)
      DIMENSION WAREA(NMAX),BIDY(NMAX,NMAX),X00(NMAX)

      REAL*8, DIMENSION(NMAX) :: SR,SR_,SA,SA0

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

c      PRINT*,'length of cyclic M chain (variable MC)?'
c      READ(5,*) MC
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
      F0=F00; X00=X0

      CALL GNUPLOT_INI(NPAR)

      N_ACC=0; N_ACC_LAMBDA=0; N_ACC_SR=0; N_ACC_SA=0

      IF(COVM(1,1) .NE. 0.) THEN
         PRINT*,'already have COVM, will directly procede to pilot run 2'
         PRINT 999
         SM=0.
         GOTO 75
      ENDIF

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

        COVM=BIDY  !use cov at max. likelihood

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
c -----------------------------------------------------------------------------
c --- initialize indep. chains

        DO M=1,MCHAINS

c --- use multivariate normal

 3      CALL GENMN(PARM,X,WORK)

c --- compute proposal x' (=y) given x (=x0). If unacceptable, goto 4  

        X0=X00 !X00 is fixed, do not change
	DO I=1,NPAR
	X0(IVN(I))=X00(IVN(I))+X(I)
	IF(ABS(X0(IVN(I))).GT.XINF) THEN
           print*,'mcmc close to virtual boundary, move rejected!',ivn(i)
           goto 3    !resample
        ENDIF  
	ENDDO

c ---   get likelihood 

	IERR=0
        CALL BTRAFO(NDIM,X0,SU0)
	CALL FUNC(SU0,NDIM,F0); NCT=NCT+1

c --- boundary checks from within func

        IF(IERR.EQ.1) GOTO 3 !checks on errors within func

c --- save initial points

        Y0(M,1:NPAR)=X0(IVN(1:NPAR)); FY0(M)=F0

        ENDDO
c -----------------------------------------------------------------------------


c --- start indiv. chains -----------------------------------------------------

        DO N=1,MCT0

c --- use multivariate normal

        ACC_P=0.                !by default reject
        CALL GENMN(PARM,X,WORK)

        X0=Y; CALL BTRAFO(NDIM,Y,SU00)  !Y is current, X0 is proposal
	DO I=1,NPAR
	X0(IVN(I))=Y(IVN(I))+.1*X(I)
	IF(ABS(X0(IVN(I))).GT.XINF) THEN
           print*,'mcmc close to virtual boundary, move rejected!',ivn(i)
           goto 4
        ENDIF  
	ENDDO

c ---   get likelihood 

	IERR=0
        CALL BTRAFO(NDIM,X0,SU)
	CALL FUNC(SU,NDIM,FP); NCT=NCT+1

c --- boundary checks from within func

        IF(IERR.EQ.1) GOTO 4 !checks on errors within func

C ---	R ratio 

        ACC_P=MIN(1.D0,DEXP((-FP+F0)))
 4	IF(ACC_P.EQ.1.D0) THEN
	 Y=X0; F0=FP; n_acc=n_acc+1
	ELSE
	 UD01=GENUNF(0.,1.) !RAN2(ISEED)
	 IF(UD01.LE.ACC_P) THEN
	   Y=X0; F0=FP; n_acc=n_acc+1
	 ENDIF
	ENDIF

        CALL BTRAFO(NDIM,Y,SU0)

        IF(MOD(N,100).EQ.0) THEN
         WRITE(6,1007) N,N_ACC,N_ACC_LAMBDA,N_ACC_SR,N_ACC_SA,F0,(SU0(IVN(I)),I=1,NPAR)
	 N_ACC=0  !resetting the acceptance counter to zero
        ENDIF
	
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
====================================================================================
         ENDDO !n     end of 1. training run
         
	 PRINT 999

         SM=0.D0; COVM=C_SUM/FLOAT(MCT0) 

 75      DO I=1,NPAR
	  PRINT 9050,(COVM(I,J)/SQRT(COVM(I,I))/
     &   SQRT(COVM(J,J)),J=1,NPAR)
         ENDDO
          
C ---	get eigenvalues and eigenvectors of covariance matrix

	NRC=NPAR; IE=0
        CALL RS(NMAX,NRC,C_SUM,W,0,S,fwork1,fwork2,ie)
	PRINT*,'RSM ierr: ',ierr

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

	MC_Smax=3; MC_S=0; N_ACC=0; N_ACC_LAMBDA=0

	print 999
	print*,'second pilot chain: includes update for heterog. parameters'
	print 999

c -------------------------------------------------------------------------
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

c        PRINT*,'choose shape parameter a:'
c        READ*,SR00

        SR(1:LH)=(SU0(IRE(1:LH))/SE(IRE(1:LH)))**2   !shape a: arbitrary starting value
        SA(1:LH)=SR(1:LH)/SU0(IRE(1:LH)) !scale b: choose large for narrow distribution 

c ---  generate lambda_i's, if L>0 / initially all lambda's equal or close
        DO I=1,LH
         LAMBDA(I,1:NOBS)=SU0(IRE(I)) !gengam(real(sa),real(sr))
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
	write(8,1001) mc_id,f0,(su0(ivn(i)),i=1,npar)
        write(9,1001) mc_id,(sr(i)/sa(i),sqrt(sr(i))/sa(i),i=1,lh)
        ELSE
	write(8,1001) mc_id,f0,(su0(ivn(i)),i=1,npar)
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
c     All parameters SU0 labeled 'H' are replaced by their respective shape 
c     parameters.

c --- note that if nuissance parameter out of bounds then reject

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
                 GP=GP-((y(ivn(i))-x_prior(i))**2)/2./s_prior(i)/s_prior(i)
                 ENDIF
           ENDDO
           
        ELSE         
           G0=0.; GP=0.
        ENDIF
         
        IF(IMCMC.EQ.1) THEN
           FH0=0.; FHP=0.; sr_=sr
        IF(MC_S.GE.1) THEN !use mvn sampler to sample a (sr)
c --- step 1: include shape parameters in block of homogeneous parameters 
c             given theta's,lambda's and b's

           DO I=1,LH  !for each heterogeneous parameter

           sll=sum(dlog(lambda(i,1:nobs)*sa(i)))
            
c --- from mvn sampler for all heterog.  parameters
            sr_(i)=y(ire(i))

c --- if chosen a is less than lower bound, reject move       
            if(sr_(i).le.1.d-6) then
              print*,'down move of shape p. rejected!',sr_(i)
              goto 14
            endif
            
            FH0=FH0+sr(i)*sll-nobs*alngam(sr(i)) !check prec.
            FHP=FHP+sr_(i)*sll-nobs*alngam(sr_(i)) !check prec.
        ENDDO

        ENDIF
        ENDIF
        Fnew=FP-FHP-GP; Fold=F0-FH0-G0

c --- boundary checks from within subroutine func()
        IF(IERR.EQ.1) GOTO 14 !checks on positivity problems etc.

C ---	R ratio 

        ACC_P=MIN(1.D0,DEXP((-Fnew+Fold)))

 14	IF(ACC_P.EQ.1.D0) THEN
	 SU0=Y; sr(1:lh)=sr_(1:lh); F0=FP; G0=GP; lkh0=lkh; n_acc=n_acc+1
	ELSE
	 IF(genunf(0.,1.).LE.ACC_P) THEN
	   SU0=Y; sr(1:lh)=sr_(1:lh); F0=FP; G0=GP; lkh0=lkh; n_acc=n_acc+1
	 ENDIF
	ENDIF
! -----------------------------------------------------------------

        IF(IMCMC.EQ.1) THEN
        IF(MC_S.LE.0) THEN !use MH sampler for a in first pilot run

           su0(ire(1:lh))=sr(1:lh)

c ---  step1: generate sr (a) | sa (b),lambda using MH sampler

           DO I=1,LH  !for each heterogeneous parameter
            ACC_G=0.
            sll=sum(dlog(lambda(i,1:nobs)*sa(i)))
            
c --- try exponential hyperprior for sr (a)
c --- try runf proposal with width dep. on current point
            sr_(i)=genunf(real(.95*sr(i)),real(1.05*sr(i)))
c           print*,sr(i),sr_(i),sll

c --- if chosen a is less than 1.d-6 (lower bound) reject       
            if(sr_(i).le.1.d-6) goto 12
            
            G0=sr(i)*sll-nobs*alngam(sr(i)) !check prec.
            GP=sr_(i)*sll-nobs*alngam(sr_(i)) !check prec.

c            print*,'sll:',sll,' sr:',sr(i),' sa:',sa(i)
c            print*,g0,gp,gp-g0
            ACC_G=MIN(1.D0,sr(i)*DEXP(GP-G0)/sr_(i))

 12         IF(ACC_G.EQ.1.D0) THEN
               sr(i)=sr_(i); n_acc_sr=n_acc_sr+1
            ELSE
               IF(genunf(0.,1.).LE.ACC_G) THEN
                  sr(i)=sr_(i); n_acc_sr=n_acc_sr+1
               ENDIF
            ENDIF
           ENDDO

           su0(ire(1:lh))=sr(1:lh)

        ENDIF

c ---  step2: generate sa (b) | sr (a),lambda using Gibbs sampler

        DO I=1,LH
         sum_lambda=sum(lambda(i,1:nobs))   
         sa(i)=gengam(real(sum_lambda),real(nobs*sr(i)+1.)) 
         if(mod(n,100).eq.0) then; print*,i,sr(i),sa(i),sqrt(sr(i))/sa(i); endif
        ENDDO

c ---  step 3: generate new lambda_i's | sa(b),sr(a)

        ACC_P=0.; IERR=0

        LAMBDA0=LAMBDA
        IREJECT=0
        DO I=1,LH
        DO K=1,NOBS
!         lkh0(i)   =lkh(i) !we have this already form sampling Y
         lambda(i,k)=gengam(real(sa(i)),real(sr(i)))
c --- the question is whether or not to use a truncated likelihood. If
c     truncated then rejection occurs if lambda_i's out of bounds.
c         if(lambda(i,k).gt.s2(ire(i)).or.lambda(i,k).lt.s1(ire(i))) then
c            print*,'lambda rejected ',i,k,lambda(i,k)
c            lambda(i,k)=lambda0(i,k) !reject
c            ireject(k)=1
c         endif
        ENDDO
        ENDDO

c ---  get new lkh's for current point (su0) and with proposed lambda_i's
c      for lambda's that donnot change this can be simplified

	CALL FUNC(SU0,NDIM,FDUM); NCT=NCT+1 

c --- use Metropolis-Hastings here to accept or reject the lambda_i's
        N_ACC_LAMBDA=0
        DO K=1,NOBS
         IF(IREJECT(K).EQ.1) GOTO 13

         ACC_P=MIN(1.D0,DEXP((-LKH(K)+LKH0(K))))

c         if(k.eq.200) then
c            print*,n,lambda0(1,k),lambda(1,k)
c            print*,lkh0(k),lkh(k),DEXP((-LKH(K)+LKH0(K)))
c         endif

 13      IF(ACC_P.EQ.1.D0) THEN
c --- all het. parameters are updated in a block indiv. by indiv.
          LAMBDA0(1:LH,K)=LAMBDA(1:LH,K); LKH0(K)=LKH(K)
          N_ACC_LAMBDA=N_ACC_LAMBDA+1
	 ELSE

	  IF(GENUNF(0.,1.).LE.ACC_P) THEN
           LAMBDA0(1:LH,K)=LAMBDA(1:LH,K); LKH0(K)=LKH(K)
           N_ACC_LAMBDA=N_ACC_LAMBDA+1
	  ENDIF
	 ENDIF
        ENDDO
        LAMBDA=LAMBDA0 !now current
        LKH=LKH0

        N_ACC_LAMBDA=100*(N_ACC_LAMBDA/FLOAT(NOBS)) !average acceptance per cycle
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

c	DO I=1,NPAR
c	if(a(i).lt.25. and. n.le.500) w(i)=.25*w(i) !balancing acceptance rates
c	if(a(i).gt.75. and. n.le.500) w(i)= 4.*w(i) !balancing acceptance rates
c	ENDDO
        IF(IMCMC.EQ.1) THEN
           XM=SU0
           XM(IRE(1:LH))= SR(1:LH)/SA(1:LH)
        ENDIF

        WRITE(6,1007) N,N_ACC,N_ACC_LAMBDA,N_ACC_SR,N_ACC_SA,F0,(XM(IVN(I)),I=1,NPAR)
	N_ACC=0; N_ACC_LAMBDA=0; N_ACC_SR=0; N_ACC_SA=0  !resetting the acceptance counter to zero

        endif
	
c --- write

	if(mod(n,mth).eq.0) then
        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
        OPEN(9,FILE='mcmc.het.out',STATUS='unknown',POSITION='append')

        IF(IMCMC.EQ.1) THEN
           XM=SU0
           XM(IRE(1:LH))= SR(1:LH)/SA(1:LH)  
           write(8,1001) mc_id,f0,(xm(ivn(i)),i=1,npar)
           write(9,1001) mc_id,(sr(i)/sa(i),sqrt(sr(i))/sa(i),i=1,lh)
        ELSE
	write(8,1001) mc_id,f0,(su0(ivn(i)),i=1,npar)
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

	ENDDO !n =============================================================

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
        print*,'number of cycles:'
        read*,mc1
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
        SU(I)=SU0_MC(I) !!!!!!!!!!!?????????????
	PRINT 1002,LABELS(I),ESU0(I)/FLOAT(NCMC),SU(I)
	ENDDO
        
        F0=F_MIN
        
        IDUM=SYSTEM('killall gnuplot')

	RETURN
 999    FORMAT(/)
 1000	FORMAT(I5,F12.4,30E12.4)
 1001	FORMAT(I5,30E16.8)
 1002   FORMAT(A10,2E12.5)
 1003  FORMAT('pram',t15,' mc mean',t24,' appr. mode')
 1007  FORMAT(5I5,F12.4,30E12.4)
c 1008  FORMAT('accepted: ',30I4)
 9040  FORMAT(/,' eigenvalues:  ',/,30E12.4,/)
 9050  FORMAT(30E14.6)
	END

