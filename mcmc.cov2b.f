c *** MARKOV-CHAIN MONTE-CARLO METHOD SEARCH ALGORITHM

      SUBROUTINE MCMC(NDIM,NPAR,NCT,X0,F0,MEANV,COVM,X)
	
      IMPLICIT REAL*8 (A-H,O-Z)

      SAVE
      PARAMETER(NMAX=100,BETA=1.,XINF=32.)

      CHARACTER*10 LABELS(NMAX),ANSWER
      CHARACTER RNF(NMAX)
      INTEGER NPAR, NCT, NRC, IWORK(NMAX), iseed1, iseed2
      REAL MEANV(NPAR),COVM(NPAR,NPAR),X(NPAR),WORK(NPAR)
      REAL PARM(NPAR * (NPAR+3)/2+1), genunf
      EXTERNAL genunf

      DIMENSION Y(NMAX),X0(NMAX),IVN(NMAX),SU0_MC(NMAX),SM(NMAX)
      DIMENSION SQ(NMAX),ISQ(NMAX),JSQ(NMAX),SU0(NMAX),SU00(NMAX)
      DIMENSION W(NMAX),S(NMAX,NMAX),FWORK1(NMAX),FWORK2(NMAX)
      DIMENSION E1(NMAX),E2(NMAX),EU0(NMAX),EU(NMAX)
      DIMENSION C_SUM(NMAX,NMAX),theta(NMAX),ESU0(NMAX)

      DIMENSION DDF(NMAX,NMAX),DF(NMAX)
      DIMENSION WAREA(NMAX),BIDY(NMAX,NMAX),X00(NMAX)

      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /IFLAGS/INDEX,IMCMC,IBOOT,ISEED,IERR  
      COMMON /DELTA/DEL(NMAX)
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)

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

c --- interactive input

      PRINT*,'length of cyclic M chain (variable MC)?'
      READ(5,*) MC

      PRINT*,'training period I (recommend 1000 cycles)?'
      READ(5,*) MCT0

      PRINT*,'training period II (recommend 2000 cycles)?'
      READ(5,*) MC1

      PRINT*,'record every mth cycle:'
      READ(5,*) mth
	
      PRINT*,'iseed1, iseed2:'
      READ(5,*) iseed1,iseed2

      NCMC=0; MCT=0
      F_MIN=1.D10
      ACC_P=0.
        
      print*,'setall'
      CALL SETALL(ISEED1,ISEED2)

      CALL BTRAFO(NDIM,X0,SU0)
      CALL FUNC(SU0,NDIM,F00); NCT=NCT+1
      F0=F00

      CALL GNUPLOT_INI(NPAR)

      N_ACC=0

c --- compute variances

	DEL=.001d0; MEANV=0.

 	CALL HESSIAN(NPAR,NDIM,IVN,X0,Y,DF,DDF,W,Z,NCT)

c        PRINT*,(DF(I),I=1,NPAR)
c        PRINT*,(DDF(I,I),I=1,NPAR)
c        pause 'check 1'
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

c --- initialize covariance matrix

          SM=0.; C_SUM=0.d0

C --- start 1. training M chain -------------------

        DO N=1,MCT0

c --- USE MULTIVARIATE NORMAL NOW !compute width of random uniform proposal kernel

        ACC_P=0.                !by default reject
        CALL GENMN(PARM,X,WORK)

c --- compute proposal x' (=y). If unacceptable, goto 4  

        X00=X0; Y=X0
	DO I=1,NPAR
	Y(IVN(I))=X0(IVN(I))+1.*X(I)
	IF(ABS(Y(IVN(I))).GT.XINF) then
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

        ACC_P=MIN(1.D0,DEXP((-FP+F0)*beta))
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

c	if(mod(n,mth).eq.0) then
c	WRITE(8,1001) F0,(SU0(I),I=1,NPAR)  !,sum1/float(n),
c     &  sum2/float(n),sum3/float(n),dexp(-f0*(1.d0-beta))
c	endif

c ---------------------------------------

        if(mod(n,100).eq.0) then
        WRITE(6,1007) N,N_ACC,F0,(SU0(IVN(I)),I=1,NPAR)
	N_ACC=0  !resetting the acceptance counter to zero
        endif
	
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
         
	 PRINT 999; PRINT*,'covariance matrix:'

          SM=0.D0; COVM=C_SUM/FLOAT(MCT0) 

         DO I=1,NPAR
	  PRINT 9050,(COVM(I,J)/SQRT(COVM(I,I))/
     &   SQRT(COVM(J,J)),J=1,NPAR)
         ENDDO
          
C ---	get eigenvalues and eigenvectors of covariance matrix

	NRC=NPAR; IE=0
        CALL RS(NMAX,NRC,C_SUM,W,0,S,fwork1,fwork2,ie)
	PRINT*,'RSM ierr: ',ierr

	PRINT 999; PRINT 9040,(W(I),I=1,NPAR); PRINT 999
        DO I=1,NPAR
	PRINT 9050,(S(I,J),J=1,NPAR)
        ENDDO
	PRINT 999

c --- initialize multivariate normal sampler (on original scale)

        CALL SETGMN(MEANV,COVM,NPAR,PARM)

c ---  procede with Q sampling from multivariate normal with COVM

        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
	WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
        close(8)

	MC_Smax=2; MC_S=0; N_ACC=0

	print 999; print*,'second pilot chain:'; print 999

c ---  reinitialize COV stuff

        print*,'scale for m.v. sampler (recommend 2.):'; read*,scl

        CALL GNUPLOT()

 66     SM=0.D0
	C_SUM=0.D0
	
C --- start revised chain 

        DO N=1,MC1

        ACC_P=0.                !by default reject
 	CALL GENMN(PARM,X,WORK)

c --- compute proposal x' (=y). If unacceptable, got 14  

        SU00=SU0; Y=SU0
	DO I=1,NPAR
	Y(IVN(I))=SU0(IVN(I))+scl*X(I)
	IF(Y(IVN(I)).LT.S1(IVN(I))) THEN
           PRINT*,'down move rejected!',ivn(i)
           GOTO 14
         ELSEIF(Y(IVN(I)).GT.S2(IVN(I))) THEN
           PRINT*,'up move rejected!',ivn(i)
           GOTO 14
        ENDIF
	ENDDO

c ---   get neg. log likelihood or energy value

	IERR=0
	CALL FUNC(Y,NDIM,FP); NCT=NCT+1

c --- boundary checks from within subroutine func()

        IF(IERR.EQ.1) GOTO 14 !checks on positivity problems

C ---	R ratio 

        ACC_P=MIN(1.D0,DEXP((-FP+F0)*beta))

 14	IF(ACC_P.EQ.1.D0) THEN
	  SU0=Y; F0=FP; n_acc=n_acc+1
	ELSE
	 UD01=GENUNF(0.,1.)  !RAN2(ISEED)
	 IF(UD01.LE.ACC_P) THEN
	  SU0=Y; F0=FP; n_acc=n_acc+1
	 ENDIF
	ENDIF

        ESU0=ESU0+SU0
        NCMC=NCMC+1

c --- update minimum in neg. log likelihood (approx. mle's)

        IF(ACC_P.EQ.1.) THEN
         IF(FP.LT.F_MIN) THEN
           SU0_MC=SU0; F_MIN=FP
         ENDIF
        ENDIF

c ---------------------------------------------

        if(mod(n,100).eq.0) then

c	DO I=1,NPAR
c	if(a(i).lt.25. and. n.le.500) w(i)=.25*w(i) !balancing acceptance rates
c	if(a(i).gt.75. and. n.le.500) w(i)= 4.*w(i) !balancing acceptance rates
c	ENDDO
        WRITE(6,1007) N,N_ACC,F0,(SU0(IVN(I)),I=1,NPAR)
	N_ACC=0  !resetting the acceptance counter to zero

        endif
	
c --- write

	IF(MOD(N,MTH).EQ.0) THEN
        OPEN(8,FILE='mcmc.out',STATUS='unknown',POSITION='append')
	WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
        CLOSE(8)
	ENDIF

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
          C_SUM(I,J)=C_SUM(I,J)+((SU0(IVI)-SU00(IVI))-THETA(I))*
     &                          ((SU0(IVJ)-SU00(IVJ))-THETA(J))
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
	
        CALL FTRAFO(NDIM,X0,SU0_MC)

        PRINT 1003
	DO I=1,NDIM
        SU(I)=SU0_MC(I)
	PRINT 1002,LABELS(I),ESU0(I)/FLOAT(NCMC),SU(I),RNF(I)
	ENDDO
        
        F0=F_MIN
        
        CLOSE(8)

        IDUM=SYSTEM('killall -9 gnuplot')
        IDUM=SYSTEM('killall -9 gnuplot_x11')

	RETURN
 999    FORMAT(/)
 1000	FORMAT(I5,F12.4,30E12.4)
 1001	FORMAT(30E16.8)
 1002   FORMAT(A10,2E12.5,5x,A1)
 1003  FORMAT('pram',t15,' mc mean',t24,' appr. mode')
 1007  FORMAT(2I5,F12.4,30E12.4)
 9040  FORMAT(/,' eigenvalues:  ',/,30E12.4,/)
 9050  FORMAT(30E14.6)
	END


      EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_INI(NPAR)
      INTEGER NPAR
c --- initialize gnuplot graphics

      INCLUDE '/usr/pgi/linux86/include/lib3f.h'
      IDUM=SYSTEM('if [ -f mcmc.out ]; then mv -f mcmc.out mcmc.out.old; fi')
      IDUM=SYSTEM('if [ -f mcmc.het.out ]; then mv -f mcmc.out mcmc.het.out.old; fi')

        OPEN(20,FILE='mcmc.gnu',STATUS='UNKNOWN')
        WRITE(20,1004) 
        DO J=3,NPAR+2
        WRITE(20,1005) J-2,J
        ENDDO
        WRITE(20,1006)
        CLOSE(20)

 1004   FORMAT(' set title "MC Markov chain (-log L)"',/,
     1  'plot "mcmc.out" using 2 w l',/,
     2  'pause 10')
 1005   FORMAT(
     3  ' set title "MC Markov chain (parameter ',i2,')"',/,
     4  'plot "mcmc.out" using ',i2,' w l',/,
     5  'pause 10')
 1006   FORMAT('reread')
        END SUBROUTINE GNUPLOT_INI


      EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT()
      INCLUDE '/usr/pgi/linux86/include/lib3f.h'
      IDUM=SYSTEM('gnuplot mcmc.gnu &')
      END SUBROUTINE GNUPLOT
