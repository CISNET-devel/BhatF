c *** MARKOV-CHAIN MONTE-CARLO METHOD FOR DECORATING INDIV. LIKELIHHOODS WITH 
c *** BERKSON ERROR (LOG-NORMAL) IN DOSE

c     uses mvn Q sampler on original scale
c     passes nobs and vector lambda to FUNC

c     imcmc=0: do mcmc (assuming no uncertainty)
c     imcmc=1: do Berkson via mcmc

      SUBROUTINE BERKSON(NDIM,NPAR,NCT,X0,F0,MEANV,COVM,X)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      
      PARAMETER(NMAX=100, XINF=16.,NSAMPLE=4000,NDOSE=20)

      CHARACTER*10 LABELS(NMAX),ANSWER,ANSWER_SIGMA
      CHARACTER L_PRIOR, RNF(NMAX)
      INTEGER NPAR, NCT, NRC, IWORK(NMAX), iseed1, iseed2
      INTEGER IVN(NMAX),IRE(NPAR),NOBS,ID1(10),ID2(10)
      REAL*8, DIMENSION(NSAMPLE) :: LKH,LKH0
      REAL*8, DIMENSION(NDOSE,NSAMPLE) :: DOSE,LAMBDA,LAMBDA0
      REAL MEANV(NPAR),COVM(NPAR,NPAR),X(NPAR),WORK(NPAR)
      REAL PARM(NPAR * (NPAR+3)/2+1), genunf, gengam
      EXTERNAL genunf,gengam

      REAL*8, DIMENSION(NMAX) :: Y,X0,XP,SU0_MC,SM
      DIMENSION SQ(NMAX),ISQ(NMAX),JSQ(NMAX),SU0(NMAX),SU00(NMAX)
      DIMENSION W(NMAX),S(NMAX,NMAX),FWORK1(NMAX),FWORK2(NMAX)
      DIMENSION C_SUM(NMAX,NMAX),theta(NMAX),ESU0(NMAX)

      DIMENSION DDF(NMAX,NMAX),DF(NMAX)
      DIMENSION WAREA(NMAX),BIDY(NMAX,NMAX),X00(NMAX)

      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /IFLAGS/INDEX,IMCMC,IBOOT,ISEED,IERR  
      COMMON /DELTA/DEL(NMAX)
      COMMON /HETERO/NOBS,ID1,ID2,LAMBDA,LKH,DOSE
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
!!         INTEGER, INTENT(IN) :: NPAR_
!! 	 END SUBROUTINE GNUPLOT_INI
!! 	END INTERFACE

!!	INTERFACE
!! 	 EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT()
!! 	 END SUBROUTINE GNUPLOT
!! 	END INTERFACE
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
      
      PRINT*,'ID, iseed1, iseed2:'
      READ(5,*) id,iseed1,iseed2
        
      NCMC=0; MCT=0
      F_MIN=1.D10
      ACC_P=0.D0
        
      print*,'setall'
      CALL SETALL(ISEED1,ISEED2)

      CALL BTRAFO(NDIM,X0,SU0)
      CALL FUNC(SU0,NDIM,F00); NCT=NCT+1
      F0=F00
      print*,f0

      CALL GNUPLOT_INI(NPAR)

      N_ACC=0; N_ACC_SIGMA=0; N_ACC_LAMBDA=0

      
      if(covm(1,1) .ne. 0.) then
         print*,'already have COVM, '//
     $        'will directly procede to pilot run 2'
         print 999
         sm=0.
         goto 75
      endif

c ---  compute likelihood based covariance from observed information matrix 

      DEL(1:NPAR)=.001d0  !perhaps use dstep
        MEANV(1:NPAR)=0.

 	CALL HESSIAN(NPAR,NDIM,IVN,X0,Y,DF,DDF,W,Z,NCT)

c        PRINT 9001,(DF(I),I=1,NPAR)
c        PRINT 9002,(DDF(I,I),I=1,NPAR)
c	 PRINT 1000

        BIDY=0.D0
        DO I=1,NPAR
           BIDY(I,I)=1.D0
        ENDDO

      IERR=0
      CALL AXEB(DDF,NPAR,NMAX,BIDY,NPAR,NMAX,0,WAREA,IERR)

        COVM=.5*BIDY  !use a progressive algorithm for adjusting scale

c --- initialize multivariate normal sampler

        CALL SETGMN(MEANV,COVM,NPAR,PARM)

c --- initialize covariance matrix for first run

          SM=0.D0; C_SUM=0.D0

C --- start 1. training M chain -------------------

        DO N=1,MCT0

c --- use multivariate normal now 

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

c ---   get likelihood 

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

c ---------------------------------------------

        IF(MOD(N,100).EQ.0) THEN
         WRITE(6,1007) N,N_ACC,N_ACC_SIGMA,N_ACC_LAMBDA,F0,
     $                 (SU0(IVN(I)),I=1,NPAR)
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

         ENDDO !n     end of 1. training run
       PRINT 999

         SM=0.D0; COVM=C_SUM/FLOAT(MCT0) 

 75      DO I=1,NPAR
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

c --- initialize multivariate normal sampler (on original scale)

        CALL SETGMN(MEANV,COVM,NPAR,PARM)

c ---  procede with Q sampling from multivariate normal with COVM

      MC_Smax=2; MC_S=0; N_ACC=0; N_ACC_SIGMA=0; N_ACC_LAMBDA=0

c -------------------------------------------------------------------------------
      print 999
      print*,'second pilot chain (includes Berkson error)'
      print 999

c ---  generate d_i's (called lambda, using recorded doses) 

        print*,'std of lognormal (sigma)' !assume all uncertainty comes from single error model
        read*,sig0
        print*,'fixed or not (y/n)'
        read(5,67) answer_sigma

        IMCMC=1

        LAMBDA=DOSE
        DO K=1,NOBS
           DO L=1,NDOSE
           IF(DOSE(L,K) .NE. 0.) THEN
           LAMBDA(L,K)=DOSE(L,K)*DEXP(-0.5*SIG0*SIG0)*
     $                  DEXP(SIG0*GENNOR(0.,1.))
           ENDIF
           ENDDO
        ENDDO
c        print*,lambda(1:5,1)
   
      CALL FUNC(SU0,NDIM,F0); NCT=NCT+1 !fills lkh's
        print*,'first (imcmc=1) function call, check: ',F0,
     $         sum(lkh(1:nobs))

        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
        WRITE(8,1001) ID,F0,(SU0(IVN(I)),I=1,NPAR),SIG0
        close(8)

c ---  reinitialize COV stuff

        print*,'scale for m.v. sampler (recommend 2.):'
        read*,scl

        CALL GNUPLOT()

 66     SM=0.D0; C_SUM=0.D0

c =============================================================================	
C --- start chain inluding Berkson error  -------------------------------------
c =============================================================================	

        DO N=1,MC1

        ACC_P=0.D0; IERR=0
 	CALL GENMN(PARM,X,WORK)

        lkh0=lkh

c --- step 1 (sample theta's) / always carried out

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
        
      CALL FUNC(Y,NDIM,FP); NCT=NCT+1  ! uses lambda

c --- boundary checks from within subroutine func()
        IF(IERR.EQ.1) GOTO 14 !checks on positivity problems etc.

C ---	R ratio ---------------------------------------------------

        ACC_P=MIN(1.D0,DEXP((-FP+F0)))

 14	IF(ACC_P.EQ.1.D0) THEN
       SU0=Y; F0=FP; lkh0=lkh; n_acc=n_acc+1
      ELSE
       IF(genunf(0.,1.).LE.ACC_P) THEN
         SU0=Y; F0=FP; lkh0=lkh; n_acc=n_acc+1
       ENDIF
      ENDIF

c ---  step 2: generate new lambda_i's (with 'Berkson doses') | theta, sigma

        ACC_P=0.; IERR=0; LAMBDA0=LAMBDA

        DO K=1,NOBS
           DO L=1,NDOSE
           IF(DOSE(L,K) .NE. 0.) THEN
           LAMBDA(L,K)=DOSE(L,K)*DEXP(-0.5*SIG0*SIG0)*
     $                  DEXP(SIG0*GENNOR(0.,1.))
           ENDIF
           ENDDO
        ENDDO

c --- obtain indiv.likelihood contrib. LKH
      CALL FUNC(SU0,NDIM,FDUM); NCT=NCT+1 

c --- use Metropolis-Hastings here to accept or reject the lambda_i's 
c     individual by individual, for the parameter block (model parameters)
        N_ACC_LAMBDA=0
        DO K=1,NOBS

         ACC_P=MIN(1.D0,DEXP((-LKH(K)+LKH0(K))))

         IF(ACC_P.EQ.1.D0) THEN
          LAMBDA0(1:NDOSE,K)=LAMBDA(1:NDOSE,K)
          LKH0(K)=LKH(K)
          N_ACC_LAMBDA=N_ACC_LAMBDA+1
       ELSE

        IF(GENUNF(0.,1.).LE.ACC_P) THEN
           LAMBDA0(1:NDOSE,K)=LAMBDA(1:NDOSE,K)
           LKH0(K)=LKH(K)
           N_ACC_LAMBDA=N_ACC_LAMBDA+1
        ENDIF
       ENDIF
        ENDDO
        LAMBDA=LAMBDA0; LKH=LKH0

        N_ACC_LAMBDA=100*(N_ACC_LAMBDA/FLOAT(NOBS)) !average acceptance per cycle
        F0=SUM(LKH(1:NOBS))

        IF(ANSWER_SIGMA.EQ.'N' .OR. ANSWER_SIGMA.EQ.'n') THEN
c ---  step3: sample eta(=1/2/sig^2) from gamma distr.
        
        ACC_P=0.            

        ERR_SUM=0.  !!!!!!!!!!!!!!!!!!! dose(?,k) !!!!!!!!!!!!!!!!!!!!
            DO K=1,NOBS
               IF(DOSE(1,K).GT.0.) THEN
               R=dlog(lambda(1,k)/dose(1,k)) + 0.5*SIG0*SIG0
               ERR_SUM=ERR_SUM+DABS(R)**2
               ELSE
               PAUSE 'dose=0'
               ENDIF
            ENDDO
        ETA=GENGAM(REAL(ERR_SUM),REAL(0.5*NOBS+1.))
        SIGP=1./SQRT(2.*ETA)

        ACC_P=MIN(1.,EXP(-(SIGP*SIGP-SIG0*SIG0)/8.))
        
        IF(ACC_P.EQ.1.D0) THEN
           SIG0=SIGP; N_ACC_SIGMA=N_ACC_SIGMA+1
       ELSE
        IF(GENUNF(0.,1.).LE.ACC_P) THEN
             SIG0=SIGP; N_ACC_SIGMA=N_ACC_SIGMA+1
          ENDIF
        ENDIF

c -----------------------------------------------------------------
        ENDIF

        NCMC=NCMC+1

        if(mod(n,100).eq.0) then
        WRITE(6,1007) N,N_ACC,N_ACC_SIGMA,N_ACC_LAMBDA,F0,
     $                (SU0(IVN(I)),I=1,NPAR),SIG0
      N_ACC=0; N_ACC_SIGMA=0; N_ACC_LAMBDA=0  !resetting the acceptance counter to zero
        endif
      
c --- write

      if(mod(n,mth).eq.0) then
        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
      WRITE(8,1001) ID,F0,(SU0(IVN(I)),I=1,NPAR),SIG0
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

        CALL SYSTEM('killall gnuplot')

      RETURN
 999    FORMAT(/)
 1000	FORMAT(I5,F12.4,30E12.4)
 1001	FORMAT(I5,30E16.8)
 1002   FORMAT(A10,2E12.5)
 1007   FORMAT(4I5,F12.4,30E12.4)
 9040   FORMAT(/,' eigenvalues:  ',/,30E12.4,/)
 9050   FORMAT(30E14.6)
      END

