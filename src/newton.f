      SUBROUTINE NEWTON(NCT,X0,F0,EPS,I_HESSE,I_RETRACT,INEW,NWT_MAX)
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NMAX=100,NMAX2=NMAX*NMAX,SENS=0.005,XINF=32.)
      INTEGER  NRC, IE
      CHARACTER*28 STAMP,STATUS
      CHARACTER*10 LABELS(NMAX)
      CHARACTER    RNF(NMAX)
      CHARACTER*10 STRING1,ANSWER

      DIMENSION Y(NMAX),X0(NMAX),SE(NMAX),IVN(NMAX),DSTEPS(NMAX)
      DIMENSION XN(NMAX,5000),F(5000),DF(NMAX),DDF(NMAX,NMAX)
     &,FDD(NMAX,NMAX),DISC(NMAX),WAREA(NMAX),ENS(NMAX,NMAX)
     &,GGF(NMAX,NMAX),BIDY(NMAX,NMAX),DISC0(NMAX),X00(NMAX)
      DIMENSION W(NMAX),Z(NMAX,NMAX),FWORK1(NMAX),FWORK2(NMAX)
      DIMENSION SC(NMAX),ETY(nmax,nmax)
      DIMENSION SUL(NMAX),SUU(NMAX),X0L(NMAX),X0U(NMAX)

      COMMON /BOUNDS/SU(NMAX),S_E(NMAX),S1(NMAX),S2(NMAX)
      COMMON /NPRODCTS/NDIM,NPAR,NP,MP,NP2,NDIM2,NPAR2,NLM  
      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /IFLAGS/INDEX,IMCMC,IBOOT,ISEED,IERR  
      COMMON /DELTA/DEL(NMAX)

      DATA DF/NMAX*0.D0/ DDF/NMAX2*1.D3/ 
      DATA GGF/NMAX2*0.D0/ ENS/NMAX2*0.D0/
      DATA SMALL/-1.D-3/

C ***	INITIALIZE ITERATION LOOP

!!      INCLUDE '/usr/pgi/linux86/include/lib3f.h'

!!	INTERFACE
!! 	 EXTRINSIC(HPF_SERIAL) SUBROUTINE BTIME(STAMP_)
!!          CHARACTER*28, INTENT(OUT) :: STAMP_
!! 	 END SUBROUTINE BTIME
!! 	END INTERFACE

        CALL BTIME(STAMP)

        CALL BTRAFO(NDIM,X0,SU)
        CALL FUNC(SU,NDIM,F0); NCT=NCT+1

        DO 10 N=1,NWT_MAX
		
	PRINT*,' ' 
	PRINT*,'FUNCTION VALUE: ',F0,'     CALL: ',NCT
	PRINT*,' ' 

C -------------------------------------------------------------------

	IF(I_HESSE.EQ.1)	THEN	!compute Hessian to high accuracy:

	CALL HESSIAN(NPAR,NDIM,IVN,X0,Y,DF,DDF,W,Z,NCT)
   
        BIDY=0.D0; GGF=DDF
    
        IF(INEW.EQ.0) PRINT 9001,(((GGF(I,J)),I=1,NPAR),J=1,NPAR)	
	
	DO I=1,NPAR
	BIDY(I,I)=1.D0
	ENDDO

	CALL AXEB(GGF,NPAR,NMAX,BIDY,NPAR,NMAX,0,WAREA,IE)

	ELSE

C *** INITIALIZE NEWTON BY CALCULATION OF A
C     NLM POINT SIMPLEX GEOMETRY IN N DIM PARAMETER SPACE

        DO 20 IN=1,NPAR
        DO 20 JN=1,NLM
20      XN(IN,JN)=X0(IVN(IN))
C
C ***   ON AXES

        IF(N.EQ.1) CALL DQSTEP(NCT,NDIM,NPAR,X0,SENS,XINF,DSTEPS)

        DO 30 I=1,NPAR

C ***	compute delta distance in each direction

	IF(N.EQ.1) DEL(I)=DSTEPS(I)       ! .0001D0+.001d0*dabs(X0(I))

                J=2*I
                K=J-1
                XN(I,J)=XN(I,J)-DEL(I)
30              XN(I,K)=XN(I,K)+DEL(I)

C ***   OFF AXIS
        IF(NPAR.GE.2) THEN
        MC=NP2+1
        DO 40 I=2,NPAR
                DO 40 J=1,I-1
                XN(I,MC)=XN(I,MC)+DEL(I)
                XN(J,MC)=XN(J,MC)+DEL(J)
40              MC=MC+1
        ENDIF

        Y=X0
        DO I=1,NLM
        Y(IVN(1:NPAR))=XN(1:NPAR,I)
        CALL BTRAFO(NDIM,Y,SU)
        CALL FUNC(SU,NDIM,F(I)); NCT=NCT+1
        ENDDO

C ***   FIRST AND DIAGONAL SECOND DERIVATIVES

        DO I=1,NPAR
        I2=I*2
        I1=I2-1
        DF(I)=(F(I2)-F(I1))/2.D0/DEL(I)
        DDF(I,I)=(F(I2)+F(I1)-2.D0*F0)/DEL(I)/DEL(I)
        ENDDO

C ***   SECOND DERIVATIVES
        IF(NPAR.GE.2) THEN
        MC=NP2+1
        DO 80 I=2,NPAR
                DO 81 J=1,I-1
                IP=2*I-1
                JP=2*J-1
                DDF(I,J)=(F0+F(MC)-F(JP)-F(IP))/DEL(I)/DEL(J)
                DDF(J,I)=DDF(I,J)
                MC=MC+1
 81             CONTINUE
 80             CONTINUE
       ENDIF

C ***   MATRIX INVERSION

        BIDY=0.D0; FDD=DDF
        DO 200 I=1,NPAR
 200    BIDY(I,I)=1.D0

C ***	CHECK EIGENVALUES

	NRC=NPAR
        CALL RS(NMAX,NRC,FDD,W,0,Z,fwork1,fwork2,ie)

        FDD=DDF         !restore Hessian

	CALL AXEB(FDD,NPAR,NMAX,BIDY,NPAR,NMAX,0,WAREA,IE)

C ***   ACCURACY TEST
        ety=0.
        do i=1,npar
           do j=1,npar
              do k=1,npar
              ety(i,j)=ety(i,j)+ddf(k,i)*bidy(k,j)
              enddo
           enddo
        enddo
c        print*,'test 1'
c        print*,(ety(i,i),i=1,npar)
c        print*,ety(2,3),ety(3,2)

	IF(IBOOT.EQ.0) THEN
	PRINT 1000
	PRINT 9040,(W(I),I=1,NPAR)
C	PRINT 9050,((Z(J,I),J=1,NPAR),I=1,NPAR)
	ENDIF

	ENDIF

C ---	variable Newton stepping

        RNEPS=1.; F1=F0; NEPSC=0; DISC=0.d0; X00=X0

        DO J=1,NPAR
          DISC(J)=DISC(J)+SUM(BIDY(J,1:NPAR)*DF(1:NPAR))
        ENDDO
        DISC0=DISC

c --- check for convergence

        NZ=0
        DO J=1,NPAR
           IF(DABS(DF(J)).LE.EPS) NZ=NZ+1
        ENDDO

c --- update current position

 100    DISC=RNEPS*DISC0

        DO J=1,NPAR
           X0(IVN(J))=X00(IVN(J))+DISC(J)
        ENDDO

        CALL BTRAFO(NDIM,X0,SU)
        CALL FUNC(SU,NDIM,F0); NCT=NCT+1
        
        IF(I_RETRACT.EQ.1) THEN
        IF(F0.GT.F1 .AND. NEPSC.LT.8) THEN
           RNEPS=.75*RNEPS
           PRINT*,'reducing step size:',RNEPS
           NEPSC=NEPSC+1
           GOTO 100
        ENDIF
        IF(NEPSC.EQ.8) PRINT*,'problem finding lower F0, will continue.'
        ENDIF

C ---   STOP-CRITERION FOR ITERATION --------------------------------------

        IF(NZ.EQ.NPAR.AND.W(1).GT.SMALL) THEN
           
        PRINT*,' '
        PRINT*,'CONVERGENCE: ------------------------------------------'
        STATUS='converged'

        IF(IBOOT.NE.2) THEN

C ---	compute Hessian symmetrically:

	CALL HESSIAN(NPAR,NDIM,IVN,X0,Y,DF,DDF,W,Z,NCT)

        BIDY=0.D0; GGF=DDF
	DO I=1,NPAR
	BIDY(I,I)=1.D0
	ENDDO

	CALL AXEB(GGF,NPAR,NMAX,BIDY,NPAR,NMAX,0,WAREA,IE)

	DO I=1,NPAR
	DO J=1,NPAR
	DO L=1,NPAR
	ENS(I,J)=ENS(I,J)+DDF(I,L)*BIDY(L,J)
	ENDDO
	ENDDO
	ENDDO
	
C        PRINT 9002,(((ENS(I,J)),I=1,NPAR),J=1,NPAR)
C        PRINT 9000,(((BIDY(I,J)),I=1,NPAR),J=1,NPAR)	

        X0L=X0; X0U=X0
 	DO I=1,NPAR
           IF(BIDY(I,I).LE.0.D0) THEN
              SE(I)=0.D0
              PRINT*,'INFORMATION MATRIX NOT POS. DEF. OR SOME OTHER PROBLEM'
           ELSE
              SE(I)=DSQRT(BIDY(I,I))
           ENDIF

C	GGF(IVN(I),IVN(I))=se(i)*1.96
           X0L(IVN(I))=X0(IVN(I))-se(i)*1.96
           X0U(IVN(I))=X0(IVN(I))+se(i)*1.96
	ENDDO	

        CALL BTRAFO(NDIM,X0,SU)        
        CALL BTRAFO(NDIM,X0L,SUL)
        CALL BTRAFO(NDIM,X0U,SUU)
	
        WRITE(6,9009) STAMP,STATUS
        WRITE(6,9005)
        J=1
        IF(IVN(J).EQ.1) THEN
          WRITE(6,9022) N,F0,LABELS(1),SU(1),DF(1),SUL(1),SUU(1),NCT
          J=J+1
        ELSE
          WRITE(6,9023) N,F0,LABELS(1),SU(1),NCT
        ENDIF

	DO I=2,NDIM
        IF(IVN(J).EQ.I) THEN
          WRITE(6,9032) LABELS(I),SU(I),DF(I),SUL(I),SUU(I) 
          J=J+1
        ELSE
          WRITE(6,9033) LABELS(I),SU(I)
        ENDIF
	ENDDO

	ELSE	!IBOOT=2

	WRITE(9,6000) F0,(SU(I),I=1,NPAR),W(1)
	RETURN	!GO BACK TO BOOTSTRAPLOOP

	ENDIF

c --- graphical diagnostic
!!        include 'newton.graph.incl'

        RETURN
        END IF

C ***	DURING NON-CONVERGEING CYCLES

	IF(INEW.EQ.0) THEN

        I_HESSE=1
        INEW  =1
	ENDIF
	
c        IF(INEW.EQ.0) PRINT 9000,(((BIDY(I,J)),I=1,NDIM),J=1,NDIM)	

        X0L=X0; X0U=X0
 	DO I=1,NPAR
           IF(BIDY(I,I).LE.0.D0) THEN
              SE(I)=0.D0
              PRINT*,'INFORMATION MATRIX NOT POS. DEF. OR SOME OTHER PROBLEM'
           ELSE
              SE(i)=DSQRT(BIDY(I,I))
           ENDIF

           SE(I)=SE(I)*1.96
           X0L(IVN(I))=X0(IVN(I))-se(i)
           X0U(IVN(I))=X0(IVN(I))+se(i)
	ENDDO	
	
        CALL BTRAFO(NDIM,X0,SU)        
        CALL BTRAFO(NDIM,X0L,SUL)
        CALL BTRAFO(NDIM,X0U,SUU)
	
        STATUS='not converged'
        WRITE(6,9009) STAMP,STATUS
        WRITE(6,9005)
        J=1
        IF(IVN(J).EQ.1) THEN
          WRITE(6,9022) N,F0,LABELS(1),SU(1),DF(1),SUL(1),SUU(1),NCT
          J=J+1
        ELSE
          WRITE(6,9023) N,F0,LABELS(1),SU(1),NCT
        ENDIF

	DO I=2,NDIM
        IF(IVN(J).EQ.I) THEN
          WRITE(6,9032) LABELS(I),SU(I),DF(I),SUL(I),SUU(I) 
          J=J+1
        ELSE
          WRITE(6,9033) LABELS(I),SU(I)
        ENDIF
	ENDDO

C -----------------------------------------------------------------    GET OUT	
	IF(INEW.EQ.0.AND.IBOOT.LE.1) THEN
        PRINT 1000
        PRINT*,'DO YOU WANT TO LEAVE NEWTON PROCEDURE?'
        READ(5,7010) ANSWER
        IF(ANSWER.EQ.'Y')  RETURN
	ELSEIF(INEW.EQ.0.AND.IBOOT.EQ.2) THEN
	PRINT*,'UNSUCCESSFULL NEWTON SEARCH - TRY GRADIENT SEARCH AGAIN'
	INEW=3
	RETURN
	ENDIF
C ----------------------------------------------------------------- 
	IF(MOD(N,10).EQ.0) INEW=0	! resetting inew after 10 iterations

10      CONTINUE
        PRINT*,'NEWTON SEARCH NOT CONVERGED!'
	
	RETURN
6000  FORMAT(F15.5,30(E18.8),E16.6)
7010	FORMAT(A10)
7020	FORMAT('ENTER FURTHER OPTIONS (ERROR, CONTOUR , SAVE ETC.): ')
8000  FORMAT('set parametric',/,'set dgrid3d 20,20',/,
     1       'set isosamples 20',/,
     1       'set data style line',/,
     1       'set nokey',/,
     2       'set hidden3d',/,
     3       'set title "log-L for',a4,'and ',a4,'"',/,
     4       'splot "tmp_"',/,
     4       'set cntrparam levels auto 10',/,
     5       'set contour',/,
     6       'replot',/,
     7       'pause -1 "Hit return to continue"',/,
     8       'clear')
8010  FORMAT(/,'GRAPHING  LABELS:   ',2A10)
8020  FORMAT(i5,' ',a5)
9000  FORMAT(/,' VARIANCE MATRIX:',//,30(30(E11.4),/))
9001  FORMAT(/,' INFORM.  MATRIX:',//,30(30(E11.4),/))
9002  FORMAT(/,' TEST     MATRIX:',//,30(30(E11.4),/))
9005  FORMAT(/,' IT',T11,'VALUE',T18,'LBL',T26,'ESTIMATES',T39,
     & 'DERIVATIVES',T55,'95% LOWER',T67,
     & '95% UPPER',T78,'CALLS',/)
9009  FORMAT(/,'BHAT RUN: ',A28,T42,'STATUS: ',A15)
9010  FORMAT(/,' IT',T11,'VALUE',T18,'LBL',T26,'ESTIMATES',T39,
     & 'LOG-ESTIMATES',T61,'DERIVATIVES',T77,'CURVE',T87,'CALLS',/)
9020  FORMAT(I3,T5,F11.5,T18,A5,E12.5,T40,E12.5,T60,F12.2,E10.2,I10)
9030  FORMAT(T18,A5,E12.5,T40,E12.5,T60,F12.2,E10.2)
 9022 FORMAT(I3,T5,F11.5,T18,A5,E12.5,T38,F12.2,T52,E12.4,E12.4,I7)
 9023 FORMAT(I3,T5,F11.5,T18,A5,E12.5,T52,'  ...... PARAMETER FIXED',I7)
 9032 FORMAT(T18,A5,E12.5,T38,F12.2,T52,E12.4,E12.4)
 9033 FORMAT(T18,A5,E12.5,T52,'  ...... PARAMETER FIXED',I7)
9040	FORMAT('APPROX. INTERNAL EIGENVALUES:',/,6E12.4,/)
9050	FORMAT(30(E12.4))
9051	FORMAT(4F12.8)
c9052	FORMAT(F12.4)
9053	FORMAT(A10)
9054	FORMAT('xlabel "parameter ',I2,'"',/,'ylabel "parameter ',I2,'"',
     &  /,'xticks 10',/,'yticks 10',/,'equalscale = off',/,
     &  'toplabel "Xmin CONTOUR PLOT"',/)
 
1000    FORMAT(/)
	END

