c *** MARKOV-CHAIN MONTE-CARLO METHOD SEARCH ALGORITHM

c     uses mvn Q sampler on logit scale

      SUBROUTINE MCMC(NPAR,NCT,X0,F0,MEANV,COVM,X)
	
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(NMAX=100)

      CHARACTER*10 LABELS(NMAX),LABELS_ORG(NMAX),ANSWER
      INTEGER NPAR, NCT, NRC, IWORK(NMAX), iseed1, iseed2
      REAL MEANV(NPAR),COVM(NPAR,NPAR),X(NPAR),WORK(NPAR)
      REAL PARM(100), genunf
      EXTERNAL genunf

      DIMENSION Y(NMAX),X0(NMAX),IVN(NMAX),SU0_MC(NMAX),SUM(NMAX)
      DIMENSION SQ(NMAX),ISQ(NMAX),JSQ(NMAX),SU0(NMAX)
      DIMENSION W(NMAX),S(NMAX,NMAX),FWORK(8*NMAX)
      DIMENSION E1(NMAX),E2(NMAX),EU0(NMAX),EU(NMAX)
      DIMENSION C_SUM(NMAX,NMAX),theta(NMAX),ESU0(NMAX)

      DIMENSION DDF(NMAX,NMAX),DF(NMAX)
      DIMENSION WAREA(NMAX),BIDY(NMAX,NMAX),X00(NMAX)

      COMMON /LABELS/LABELS,IVN
      COMMON /IFLAGS/INDEX,IMCMC,IBOOT,ISEED,IERR  
      COMMON /DELTA/DEL(NMAX)
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)

      DATA SUM/NMAX*0.D0/ BETA/1.d0/

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

c      PRINT*,'inverse temperature for heating likelihood'
c      READ(5,*) BETA

      PRINT*,'record every mth cycle:'
      READ(5,*) mth
	
      PRINT*,'iseed1, iseed2:'
      READ(5,*) iseed1,iseed2

      NCMC=0
      MCT=0
      F_MIN=1.D10
      ACC_P=0.
        
      print*,'setall'
      CALL SETALL(ISEED1,ISEED2)

      CALL BTRAFO(NPAR,X0,SU0)
      CALL FUNC(SU0,NPAR,F00) !current
      NCT=NCT+1
      F0=F00

      CALL GNUPLOT_INI(NPAR)

      N_ACC=0

c --- compute variances

	DO I=1,NPAR
	DEL(I)=.001d0
        MEANV(I)=0.
c	print*,i,ivn(i),del(i)
	ENDDO

 	CALL HESSIAN(NPAR,NPAR,IVN,X0,Y,DF,DDF,W,Z,NCT)

c        PRINT 9001,(DF(I),I=1,NPAR)
c        PRINT 9002,(DDF(I,I),I=1,NPAR)
c	 PRINT 1000

C ---  compute covariance 

	DO I=1,NPAR
	 BIDY(I,I)=1.d0
	ENDDO

	IERR=0
	CALL AXEB(DDF,NPAR,NMAX,BIDY,NPAR,NMAX,0,WAREA,IERR)

        DO I=1,NPAR
         DO J=1,NPAR
         COVM(I,J)=BIDY(I,J)
         ENDDO
        ENDDO

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


c --- initialize covariance matrix

       DO J=1,NPAR
          SUM(J)=0.
       DO I=1,NPAR
          C_SUM(J,I)=0.d0
       ENDDO
      ENDDO


C --- start training M chain -------------------

        DO N=1,MCT0

c --- USE MULTIVARIATE NORMAL NOW !compute width of random uniform proposal kernel

        ACC_P=0.                !by default reject
        CALL GENMN(PARM,X,WORK)

c --- compute proposal x' (=y). If unacceptable, goto 4  

	DO I=1,NPAR
        X00(I)=X0(I)
	Y(I)=X0(I)+1.*X(I)
	if(abs(y(i)).gt.16.) then
           print*,'mcmc close to boundary, move rejected!',i
           goto 4
        endif  
	ENDDO

c	print*,su(1),su(2),su(3)

c ---   get likelihood or energy 

	IERR=0
        CALL BTRAFO(NPAR,Y,SU)
	CALL FUNC(SU,NPAR,FP)
	NCT=NCT+1

c --- boundary checks from within func

        IF(IERR.EQ.1) GOTO 4 !checks on positivity problems

C ---	R ratio 

        ACC_P=MIN(1.D0,DEXP((-FP+F0)*beta))
 4	IF(ACC_P.EQ.1.D0) THEN

	 DO K=1,NPAR
	  X0(K)=Y(K)
	 ENDDO

         F0=FP
	 n_acc=n_acc+1

	ELSE
	 UD01=genunf(0.,1.) !RAN2(ISEED)
	 IF(UD01.LE.ACC_P) THEN

	  DO K=1,NPAR
	   X0(K)=Y(K)
	  ENDDO

          F0=FP
	  n_acc=n_acc+1

	 ENDIF

	ENDIF

        CALL BTRAFO(NPAR,X0,SU0)

c --- update the best likelihood estimate (ble)

        IF(ACC_P.EQ.1.) THEN
         IF(FP.LT.F_MIN) THEN
           DO J=1,NPAR
            SU0_MC(J)=SU0(J)
           ENDDO
           F_MIN=FP
         ENDIF
        ENDIF

c ---------------------------------------------

c	if(mod(n,mth).eq.0) then
c	WRITE(8,1001) F0,(SU0(I),I=1,NPAR)  !,sum1/float(n),
c     &  sum2/float(n),sum3/float(n),dexp(-f0*(1.d0-beta))
c	endif

c ---------------------------------------

        if(mod(n,100).eq.0) then
        WRITE(6,1007) N,N_ACC,F0,(SU0(I),I=1,NPAR)
	N_ACC=0  !resetting the acceptance counter to zero
        endif
	
        CALL FTRAFO(NPAR,X0,SU0)

         DO I=1,NPAR
	  SUM(I)=SUM(I)+(X0(I)-X00(I))
          theta(i)=SUM(I)/float(n)  !new mean
         ENDDO

         DO I=1,NPAR
          DO J=1,NPAR
          C_SUM(i,j)=C_SUM(i,j)+((X0(I)-X00(I))-theta(I))*
     &                          ((X0(J)-X00(J))-theta(J))
	  ENDDO
         ENDDO

c ---------------------------------------

         ENDDO !n     end of (burn-in + training) period

	 PRINT 999
         DO I=1,NPAR
          SUM(I)=0.d0 
          DO J=1,NPAR
           COVM(i,j)=C_SUM(i,j)/FLOAT(MCT0) 
          ENDDO
	 ENDDO

         DO I=1,NPAR
	  PRINT 9050,(COVM(i,j)/SQRT(COVM(I,I))/
     &   SQRT(COVM(J,J)),J=1,NPAR)
         ENDDO
          
C ---	get eigenvalues and eigenvectors of covariance matrix

	NRC=NPAR
	IERR=0
        CALL RSM(NMAX,NRC,C_SUM,W,NRC,S,fwork,iwork,ierr)
	PRINT*,'RSM ierr: ',ierr

	PRINT 999
	PRINT 9040,(W(I),I=1,NPAR)
	PRINT 999
        DO I=1,NPAR
	PRINT 9050,(S(I,J),J=1,NPAR)
        ENDDO
	PRINT 999

c --- initialize multivariate normal sampler

        CALL SETGMN(MEANV,COVM,NPAR,PARM)

c        PAUSE 'ready?'

c ---------------------------------------

c ---  procede with Q sampling from multivariate normal with COVM

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
	WRITE(8,1001) F0,(SU0(I),I=1,NPAR)
        close(8)

	MC_Smax=2
	MC_S=0
 	N_ACC=0

	print 999
	print*,'second pilot chain:'
	print 999

c ---  reinitialize COV stuff

        print*,'scale for m.v. sampler (recommend 2.):'
        read*,scl

        CALL GNUPLOT()

 66	DO I=1,NPAR
	SUM(I)=0.
	DO J=1,NPAR
	C_SUM(I,J)=0.
	ENDDO
	ENDDO
	
C --- start revised chain 


        DO N=1,MC1

        ACC_P=0.                !by default reject
 	CALL GENMN(PARM,X,WORK)

c --- compute proposal x' (=y). If unacceptable, got 14  

	DO I=1,NPAR
        X00(I)=X0(I)
	Y(I)=X0(I)+scl*X(I)
	if(abs(y(i)).gt.16.) then
           print*,'mcmc close to boundary, move rejected!',i
           goto 14
        endif  
	ENDDO

c ---   get neg. log likelihood or energy value

	IERR=0
        CALL BTRAFO(NPAR,Y,SU)
	CALL FUNC(SU,NPAR,FP)
	NCT=NCT+1

c --- boundary checks from within subroutine func()

        IF(IERR.EQ.1) GOTO 14 !checks on positivity problems

C ---	R ratio 

        ACC_P=MIN(1.D0,DEXP((-FP+F0)*beta))

 14	IF(ACC_P.EQ.1.D0) THEN

	 DO K=1,NPAR
	  X0(K)=Y(K)
	 ENDDO

         F0=FP
	 n_acc=n_acc+1

	ELSE
	 UD01=genunf(0.,1.)  !RAN2(ISEED)
	 IF(UD01.LE.ACC_P) THEN

	  DO K=1,NPAR
	   X0(K)=Y(K)
	  ENDDO

          F0=FP
	  n_acc=n_acc+1

	 ENDIF

	ENDIF

        CALL BTRAFO(NPAR,X0,SU0)

        DO J=1,NPAR
           ESU0(J)=ESU0(J)+SU0(J)
        ENDDO
        NCMC=NCMC+1

c --- update minimum in neg. log likelihood (approx. mle's)

        IF(ACC_P.EQ.1.) THEN
         IF(FP.LT.F_MIN) THEN
           DO J=1,NPAR
           SU0_MC(J)=SU0(J)
           ENDDO
           F_MIN=FP
         ENDIF
        ENDIF

c ---------------------------------------------

        if(mod(n,100).eq.0) then

c	DO I=1,NPAR
c	if(a(i).lt.25. and. n.le.500) w(i)=.25*w(i) !balancing acceptance rates
c	if(a(i).gt.75. and. n.le.500) w(i)= 4.*w(i) !balancing acceptance rates
c	ENDDO
        WRITE(6,1007) N,N_ACC,F0,(SU0(I),I=1,NPAR)
	N_ACC=0  !resetting the acceptance counter to zero

        endif
	
c --- WRITE

	if(mod(n,mth).eq.0) then
        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
	WRITE(8,1001) F0,(SU0(I),I=1,NPAR)
        close(8)
	endif

c --- compute covariance for n>500

        IF(MC_S.LT.MC_Smax) THEN

        CALL FTRAFO(NPAR,X0,SU0)

         DO I=1,NPAR
	  SUM(I)=SUM(I)+(X0(I)-X00(I))
          theta(i)=SUM(I)/float(n)  !new mean
         ENDDO

         DO I=1,NPAR
          DO J=1,NPAR
          C_SUM(i,j)=C_SUM(i,j)+((X0(I)-X00(I))-theta(I))*
     &                          ((X0(J)-X00(J))-theta(J))
	  ENDDO
         ENDDO

        ENDIF
c ---------------------------------------

	ENDDO !n

	MC_S=MC_S+1

        IF(MC_S.LT.MC_Smax) THEN

	 PRINT 999
         DO I=1,NPAR
          SUM(I)=0.d0 
          DO J=1,NPAR
           COVM(i,j)=C_SUM(i,j)/FLOAT(MC1) 
          ENDDO
	 ENDDO

         DO I=1,NPAR
	  PRINT 9050,(COVM(i,j)/SQRT(COVM(I,I))/
     &   SQRT(COVM(J,J)),J=1,NPAR)
         ENDDO
          
c --- initialize multivariate normal sampler

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
	
        CALL FTRAFO(NPAR,X0,SU0_MC)

        PRINT 1003
	DO I=1,NPAR
        SU(I)=SU0_MC(I)
	PRINT 1002,LABELS(I),ESU0(I)/FLOAT(NCMC),SU(I)
	ENDDO
        
        F0=F_MIN
        
        CLOSE(8)

	RETURN
 999    FORMAT(/)
 1000	FORMAT(I5,F12.4,30E12.4)
 1001	FORMAT(30E16.8)
 1002   FORMAT(A10,2E12.5)
 1003  FORMAT('pram',t15,' mc mean',t24,' appr. mode')
 1007  FORMAT(2I5,F12.4,30E12.4)
c 1008  FORMAT('accepted: ',30I4)
 9040  FORMAT(/,' eigenvalues:  ',/,30E12.4,/)
 9050  FORMAT(30E14.6)
	END


      EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_INI(NPAR)
      INTEGER NPAR
c --- initialize gnuplot graphics

      INCLUDE '/usr/local/pgi/linux86/include/lib3f.h'
      IDUM=SYSTEM('if [ -f mcmc.out ]; then mv -f mcmc.out mcmc.out.old; fi')

        OPEN(20,FILE='mcmc.gnu',STATUS='UNKNOWN')
        WRITE(20,1004) 
        DO J=2,NPAR+1
        WRITE(20,1005) J-1,J
        ENDDO
        WRITE(20,1006)
        CLOSE(20)

 1004   FORMAT(' set title "MC Markov chain (-log L)"',/,
     1  'plot "mcmc.out" using 1 w l',/,
     2  'pause 10')
 1005   FORMAT(
     3  ' set title "MC Markov chain (parameter ',i2,')"',/,
     4  'plot "mcmc.out" using ',i2,' w l',/,
     5  'pause 10')
 1006   FORMAT('reread')
        END SUBROUTINE GNUPLOT_INI


      EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT()
      INCLUDE '/usr/local/pgi/linux86/include/lib3f.h'
      IDUM=SYSTEM('gnuplot mcmc.gnu &')
      END SUBROUTINE GNUPLOT




