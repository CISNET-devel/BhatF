c *** MARKOV-CHAIN MONTE-CARLO METHOD SEARCH ALGORITHM
c     
c     INCLUDES: mcmc.block
      SUBROUTINE GLOBAL(NDIM,NPAR,NCT,X0,F0,COVM)
      
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(NMAX=100,EPS=.5d-1)

      CHARACTER*10 LABELS(NMAX),LABELS_ORG(NMAX),ANSWER,BLOCK
      CHARACTER RNF(NMAX)

      DIMENSION Y(NMAX),X0(NMAX),IVN(NMAX),SU0_MC(NMAX)
      DIMENSION SQ(NMAX),ISQ(NMAX),JSQ(NMAX),SU0(NMAX),SU00(NMAX)
      DIMENSION DI(NMAX),N_ACCEPTED(NMAX),N_TRIED(NMAX),SM(NMAX)
      DIMENSION C_SUM(NMAX,NMAX),THETA(NMAX),ACC(NPAR),IR(1)
      REAL MEANV(NPAR),COVM(NPAR,NPAR),X(NPAR),WORK(NPAR)
      REAL PARM(NPAR * (NPAR+3)/2+1), genunf
      EXTERNAL genunf

      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /DELTA/DEL(NMAX)
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)

      DATA BETA/1.d0/ N_ACCEPTED/NMAX*0/

! -----------------------------------------------------------------
!!      INCLUDE '/usr/pgi/linux86/include/lib3f.h'

!! 	INTERFACE
!! 	 EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_GLOBAL_INI(NPAR_)
!!          INTEGER, INTENT(IN) :: NPAR_
!! 	 END SUBROUTINE GNUPLOT_GLOBAL_INI
!! 	END INTERFACE

!! 	INTERFACE
!! 	 EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_GLOBAL()
!! 	 END SUBROUTINE GNUPLOT_GLOBAL
!! 	END INTERFACE
! -----------------------------------------------------------------

!      PRINT*,'number of cycles?'
!      READ(5,*) MC
      MC=NPAR*400

c  --- scale dsc: used in defining proposal

      DSC=.01

!      PRINT*,'adjustment  period?'
!      READ(5,*) MCT0
      MCT0=200

      PRINT*,'number of cycles for search   ',MC
      PRINT*,'adjusting acceptance rate for ',MCT0,' cycles'

      IF(MCT0.GT.MC) MCT0=MC

!      PRINT*,'inverse temperature for heating likelihood'
!      READ(5,*) BETA
      BETA=1.

!      PRINT*,'iseed1, iseed2:'
!      READ(5,*) iseed1,iseed2

      iseed1=66666
      iseed2=77777

      MCT=0
      F_MIN=1.D10
        
      CALL SETALL(ISEED1,ISEED2)

c      PRINT*,'block updates:'
c      READ(5,9) BLOCK
      BLOCK='n'
        
        CALL BTRAFO(NDIM,X0,SU0)
      CALL FUNC(SU0,NDIM,F00); NCT=NCT+1
        F0=F00

      CALL GNUPLOT_GLOBAL_INI(NPAR)

        OPEN(8,FILE='global.out',STATUS='unknown',POSITION='append')
      WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
      WRITE(8,1001) 1.00000001*F0,(SU0(IVN(I)),I=1,NPAR)
        close(8)

        ac2=1.-(.25) !**(1./float(npar)) !1-target acceptance rate

c --- create suitable random uniform proposal kernel (width)

      DO J=1,NPAR
       IVJ=IVN(J)
       DI(J)=DSC*MIN(DABS(SU0(IVJ)-S1(IVJ)),DABS(SU0(IVJ)-S2(IVJ)))
      ENDDO

c     start plotting

        CALL GNUPLOT_GLOBAL()

C --- start M chain ------------------------------------------------------------

      SUM1=0.D0; SUM2=0.D0; SUM3=0.D0; N_ACCEPT2=0; N_TRIED=0; Y=SU0
        SM=0.D0; C_SUM=0.D0

        ACC=1.
        DO N=1,MCT0 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        SU00=SU0
        DO I=1,NPAR

      DH=DI(I)/2.d0; ACC_P1=0.; ACC_P2=0.

        IVI=IVN(I)
 3      SU(IVI)=SU0(IVI)+GENUNF(real(-DH),real(DH)) !proposal

        IF(SU(IVI).LE.S1(IVI)) GOTO 4 !these are rejected
        IF(SU(IVI).GE.S2(IVI)) GOTO 4

      CALL FUNC(SU,NDIM,FP); NCT=NCT+1
        IF(FP.EQ.nan .or. FP.EQ.-nan) GOTO 4
        IF(FP.EQ.inf .or. FP.EQ.-inf) GOTO 4
        
        ACC(I)=DEXP((-FP+F0)*BETA)
        ACC_P1=MAX(EPS,MIN(1.D0,ACC(I))) !max(eps,...)  avoids freezing at initialization 
        
        N_TRIED(I)=N_TRIED(I)+1

 4      IF(ACC_P1.EQ.1.D0) THEN
       SU0(IVI)=SU(IVI); F0=FP
       n_accept2=n_accept2+1
         N_ACCEPTED(I)=N_ACCEPTED(I)+1
         SU0_MC=SU0 !update the best likelihood estimate (ble)
         F_MIN=FP
         GOTO 3
      ELSE
       UD01=genunf(0.,1.) !RAN2(ISEED)
       IF(UD01.LE.ACC_P1) THEN
        SU0(IVI)=SU(IVI); F0=FP
        n_accept2=n_accept2+1
         N_ACCEPTED(I)=N_ACCEPTED(I)+1
       ENDIF
      ENDIF

        ENDDO

c ---   outputs

        IF(MOD(N,100).EQ.0) THEN

c ---   update proposal width array DI every 100 updates for acceptance near 50 %

        WRITE(6,'("Parameter",t12,"no. tried",t24,'//
     $          '"no. accepted",t38,"ratio")')
        DO I=1,NPAR
        FREQ=FLOAT(N_ACCEPTED(I))/FLOAT(N_TRIED(I))
        WRITE(6,'(I5,T12,I5,T24,I5,T36,F6.2)')
     $          I,N_TRIED(I),N_ACCEPTED(I),FREQ

      DI(I)=MAX(1.d-10,DI(I)*(FREQ+.5)**2)
        ENDDO
      N_ACCEPTED=0; N_TRIED=0  !resetting the acceptance counter to zero

        WRITE(6,1000) N,N_ACCEPT2,FP,(SU0(IVN(I)),I=1,NPAR)
        N_ACCEPT2=0

        ENDIF
      
        OPEN(8,FILE='global.out',STATUS='unknown',POSITION='append')
      WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
        CLOSE(8)

c     --- monitor sample covariance 

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


      ENDDO !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        MEANV(1:NPAR)=0.; COVM=C_SUM/FLOAT(MCT0) 
        SCL=2. !scale increment
        CALL SETGMN(MEANV,COVM,NPAR,PARM)

        
        DO N=1,MC !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c ---   now use empirical cov of increments

        ACC_P=0.                !by default reject
        CALL GENMN(PARM,X,WORK)

        SU00=SU0; ACC_P1=0.; ACC_P2=0.
        SU(IVN(1:NPAR))=SU0(IVN(1:NPAR))+SCL*X !proposal

        DO I=1,NPAR
        IVI=IVN(I)

        IF(SU(IVI).LE.S1(IVI)) THEN
          PRINT*,'down move rejected!',IVI
          GOTO 5 
        ELSEIF (SU(IVI).GE.S2(IVI)) THEN
          PRINT*,'UP move rejected!',IVI
          GOTO 5 
        ENDIF
        ENDDO

      CALL FUNC(SU,NDIM,FP); NCT=NCT+1
        IF(FP.EQ.nan .or. FP.EQ.-nan) GOTO 5
        IF(FP.EQ.inf .or. FP.EQ.-inf) GOTO 5
        
        ACC(I)=DEXP((-FP+F0)*BETA)
        ACC_P1=MAX(EPS,MIN(1.D0,ACC(I))) !max(eps,...)  avoids freezing at initialization 
        
        N_TRIED(1)=N_TRIED(1)+1

 5      IF(ACC_P1.EQ.1.D0) THEN
       SU0=SU; F0=FP
       n_accept2=n_accept2+1
         N_ACCEPTED(1)=N_ACCEPTED(1)+1
         SU0_MC=SU0 !update the best likelihood estimate (ble)
         F_MIN=FP
      ELSE
       UD01=genunf(0.,1.) !RAN2(ISEED)
       IF(UD01.LE.ACC_P1) THEN
        SU0=SU; F0=FP
        n_accept2=n_accept2+1
         N_ACCEPTED(1)=N_ACCEPTED(1)+1
       ENDIF
      ENDIF

c ---   outputs

        IF(MOD(N,100).EQ.0) THEN

c ---   update COVM scale SCL every 100 cycles for acceptance near 50 %

        WRITE(6,'("no. tried",t12,"no. accepted",t26,"ratio")')
        FREQ=FLOAT(N_ACCEPTED(1))/FLOAT(N_TRIED(1))
        WRITE(6,'(I5,T12,I5,T24,F6.2)') N_TRIED(1),N_ACCEPTED(1),FREQ
      SCL=SCL*(FREQ+.75)**2
      N_ACCEPTED=0; N_TRIED=0  !resetting the acceptance counter to zero

        WRITE(6,1000) N,N_ACCEPT2,FP,(SU0(IVN(I)),I=1,NPAR)
        N_ACCEPT2=0

        ENDIF
      
        OPEN(8,FILE='global.out',STATUS='unknown',POSITION='append')
      WRITE(8,1001) F0,(SU0(IVN(I)),I=1,NPAR)
        CLOSE(8)

        ENDDO !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        F0=F_MIN; CALL FTRAFO(NDIM,X0,SU0_MC)

      RETURN
 1000	FORMAT(2I5,F12.4,30E12.4)
 1001	FORMAT(30E16.8)
 1002   FORMAT(I5,2E12.5,2F12.6)
 1003  FORMAT(' pram',t10,' mc mean',t19,' appr. mode',t32,'acc. rates',
     &  t46,'l. width')
      END


!!      EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_GLOBAL_INI(NPAR)
      SUBROUTINE GNUPLOT_GLOBAL_INI(NPAR)
      INTEGER NPAR
c --- initialize gnuplot graphics

!!      INCLUDE '/usr/pgi/linux86/include/lib3f.h'
      IDUM=SYSTEM('if [ -f global.out ]; '//
     $            'then mv -f global.out global.out.old; fi')

        OPEN(20,FILE='global.gnu',STATUS='UNKNOWN')
        WRITE(20,1004) 
        WRITE(20,1006)
        CLOSE(20)

 1004   FORMAT(' set title "GLOBAL (-log L)"',/,
     1  'plot "global.out" using 1 w l',/,
     2  'pause 10')
 1006   FORMAT('reread')
        END SUBROUTINE GNUPLOT_GLOBAL_INI


!!      EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_GLOBAL()
      SUBROUTINE GNUPLOT_GLOBAL()
!!      INCLUDE '/usr/pgi/linux86/include/lib3f.h'
      IDUM=SYSTEM('gnuplot global.gnu &')
      END SUBROUTINE GNUPLOT_GLOBAL
      
