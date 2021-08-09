C ------------------------------------------------------------------------
C                                           
C                                           
C
C       PROGRAM:     BHAT.f
C
C       FLOW CONTROL FOR MODULES GLOBAL,MCMC,SIMPLEX,NEWTON,MIGRAD,CI        
C                    
C
C ------------------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NMAX=100)
      PARAMETER (PI=3.1415926535898D0)
      CHARACTER*10 LABELS(NMAX),FLBL
      CHARACTER RNF(NMAX)
      CHARACTER*10 STRING1,ANSWER
      CHARACTER*6 Bhat
      REAL MEANV(NMAX),COVM(NMAX,NMAX),XP(NMAX)
      REAL*8 LAMBDA,LKH

      DIMENSION X(NMAX),YS(NMAX+1),XS0(NMAX)
     &,X0(NMAX),PS(NMAX+1,NMAX),IV(NMAX),IVN(NMAX)

      DIMENSION S_WORK(300),I_S_WORK(300), scale(300)

      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)
      COMMON /NPRODCTS/NDIM,NPAR,NP,MP,NP2,NDIM2,NPAR2,NLM  
      COMMON /LABELS/LABELS,RNF,IVN
      COMMON /IFLAGS/IFLAG,IMCMC,IBOOT,ISEED,IERR  
      COMMON /CHANGE/VAR(NMAX)
      COMMON /DELTA/DEL(NMAX)

!!      INCLUDE '/usr/pgi/linux86/include/lib3f.h'

!!      INTERFACE 
!!         INTEGER FUNCTION SYSTEM (COMMANDA) 
!!         CHARACTER(LEN=*) COMMANDA 
!!         END FUNCTION SYSTEM 
!!      END INTERFACE

      PRINT*,'---------------------------------------------------------'
      PRINT*,' '
      PRINT*,' '
      PRINT*,'Bhat Version 1.0  dated  Tue Oct 26 12:27:44 PDT 1999'
      PRINT*,' '
      PRINT*,' '
      PRINT*,'---------------------------------------------------------'
        PRINT*,' '
        
        Bhat='Bhat> '

        OPEN(7,FILE="bhat.inp",status="unknown")
      PRINT*,'INPUT FILE: bhat.inp'
  
        IVN=0; I=1; J=1
 1     READ (7,*,END=10,ERR=10) LABELS(I),SU(I),RNF(I),SE(I),S1(I),S2(I)
        IF(LABELS(I).EQ."END") GOTO 10
        RL=I
        PRINT 1002,RL,LABELS(I),SU(I),RNF(I),SE(I),S1(I),S2(I)
        IF(RNF(I).NE."S") THEN
           IVN(J)=I; J=J+1
        ENDIF
        I=I+1
        GOTO 1

 10     NDIM=I-1; NPAR=J-1
        PRINT*,'NDIM:',NDIM,' NPAR:',NPAR

C ---	PRODUCTS OF NDIM AND NPAR

      NP=NPAR
      MP=NP+1
      NP2=2*NPAR
      NDIM2=NDIM*NDIM
      NPAR2=NPAR*NPAR
      NLM=(NPAR*NPAR+3*NPAR)/2

      CLOSE(7)

      PRINT*,'---------------------------------------------------------'
      PRINT*,' '

C --- PARAMETER TRANSFORMATION

      CALL FTRAFO(NDIM,X0,SU)

        XS0 = X0
      NCT=0			! ... counts no. of FUNC calls

C *** FIRST FUNCTION CALL 

        IFLAG=1
        IMCMC=0
      IBOOT=0
        COVM=0.

        CALL FUNC(SU,NDIM,F_START); NCT=NCT+1
      PRINT 1000
      PRINT*,'first call with function value: ',F_START
        IFLAG=2   

      PRINT 1000
 380    WRITE(6, FMT='(A6)', ADVANCE='NO') Bhat
        READ(5, FMT='(A10)', ADVANCE='YES')  STRING1

C ---------------------------------------------------------------------- IF HELP

 	IF(STRING1.EQ.'HELP'.OR.STRING1.EQ.'help') THEN

      PRINT*,'main menu of commands:'
      PRINT*,'---------------------------------------------------------'
      PRINT*,'global :	Runs a Markov-Chain MC search'
      PRINT*,'simplex:	Runs a Simplex algorithm'
      PRINT*,'grad   :	Runs a Fletcher-Powell-Davidson algorithm'
      PRINT*,'newton :	Runs a Newton-Raphson algorithm'
      PRINT*,' '
      PRINT*,'(recommended sequence: simplex, grad, newton)'
      PRINT*,' '
      PRINT*,'mcmc   :	Runs a Markov-Chain MC algorithm'
      PRINT*,'raneff :	Runs a Markov-Chain MC with random effects'
      PRINT*,'boot   :	Run bootstrap
     $ (assumes that sampling is implemented)'
      PRINT*,'ci     :        Compute profile likelihood
     $ CIs of par <lbl>'
      PRINT*,' '
      PRINT*,'reset  :	Read input file bhat.inp'
      PRINT*,'perturb:	Change parameter with label <lbl>'
      PRINT*,'fix    :	Fix parameter with label <lbl>'
      PRINT*,'rel    :	Release parameter with label <lbl>'
        PRINT*,'log    :        invoke log entry'
      PRINT*,' '
      PRINT*,'help   :	Gives this menu'
      PRINT*,'exit   :	stop Bhat and save parameters in bhat.out'
      PRINT*,'---------------------------------------------------------'
      PRINT*,'Yes or NO answer with Y or N'

      PRINT 1000
      GOTO 380


C ---------------------------------------------------------------------- IF RESET

 	ELSEIF(STRING1.EQ.'RESET'.OR.STRING1.EQ.'reset') THEN

        OPEN(7,FILE="bhat.inp",status="unknown")

      PRINT 1000
      PRINT*,'resetting parameters to values in bhat.inp:'
      PRINT*,'-------------------------------------------'

        IVN=0; I=1; J=1
 2     READ (7,*,END=11,ERR=11) LABELS(I),SU(I),RNF(I),SE(I),S1(I),S2(I)
        IF(LABELS(I).EQ."END") GOTO 11
        RL=I
        PRINT 1002,RL,LABELS(I),SU(I),RNF(I),SE(I),S1(I),S2(I)
        IF(RNF(I).NE."S") THEN
           IVN(J)=I; J=J+1
        ENDIF
        I=I+1
        GOTO 2

 11     NDIM=I-1; NPAR=J-1
        PRINT*,'NDIM:',NDIM,' NPAR:',NPAR

        NCT=0
      NP=NPAR
      MP=NP+1
      NP2=2*NPAR
      NDIM2=NDIM*NDIM
      NPAR2=NPAR*NPAR
      NLM=(NPAR*NPAR+3*NPAR)/2

      CLOSE(7)

        CALL FTRAFO(NDIM,X0,SU)  
        XS0 = X0

        CALL FUNC(SU,NDIM,F_START); NCT=NCT+1
        PRINT*,' '
      PRINT*,'NEW VALUE:',F_START

      PRINT 1000
C --------------------------------------------------------------------- END RESET
      GOTO 380
C ---------------------------------------------------------------------- IF PERTURB

 	ELSEIF(STRING1.EQ.'PERTURB'.OR.STRING1.EQ.'perturb') THEN
           
 22        WRITE(6, FMT='(A7)', ADVANCE='NO') 'Label: '
           READ(5, FMT='(A10)', ADVANCE='YES', ERR=23) FLBL

           IF(FLBL.EQ.'') GOTO 23
           DO I=1,NDIM
              IF(FLBL.EQ.LABELS(I)) THEN
                 WRITE(6,FMT='(A20,2X,A10)') 'perturbing parameter',FLBL
                 IF(RNF(I).EQ.'S')
     $                PRINT*,'note, this is a STOPPED parameter'
                 IF(RNF(I).EQ.'R')
     $                PRINT*,'note, this is a RUNNING parameter'
                 WRITE(6, FMT='(A12)', ADVANCE='NO') 'what value: '
                 READ(5, FMT='(F10.0)', ADVANCE='YES', ERR=25) SU(I)
                 IF(SU(I).GE.S2(I) .OR. SU(I).LE.S1(I)) THEN
                    PRINT*,'value is out of selected bounds'//
     $                   ' (use reset to fix)'
                    GOTO 23
                 ENDIF
                 CALL FTRAFO(NDIM,X0,SU)  
                 XS0 = X0

                 CALL FUNC(SU,NDIM,F_START); NCT=NCT+1
                 PRINT*,' '
               PRINT*,'NEW VALUE:',F_START

                 GOTO 22
              ENDIF
           ENDDO
           PRINT*,'PARAMETER NOT FOUND!'
 23        CONTINUE
           
C ---------------------------------------------------------------------- END PERTURB
        GOTO 380
C ---------------------------------------------------------------------- IF FIX

 	ELSEIF(STRING1.EQ.'FIX'.OR.STRING1.EQ.'fix') THEN
           
 24        WRITE(6, FMT='(A7)', ADVANCE='NO') 'Label: '
           READ(5, FMT='(A10)', ADVANCE='YES', ERR=25) FLBL

           IF(FLBL.EQ.'') GOTO 25
           DO I=1,NPAR
              IF(FLBL.EQ.LABELS(IVN(I))) THEN
                 WRITE(6,FMT='(A16,2X,A10)') 'fixing parameter',FLBL
                 RNF(IVN(I))='S'
                 DO J=I,NPAR-1  !slide back seq. ids
                    IVN(J)=IVN(J+1)
                 ENDDO
                    IVN(NPAR)=0
                    
                 NPAR=NPAR-1
                 NP=NPAR
                 MP=NP+1
                 NP2=2*NPAR
                 NPAR2=NPAR*NPAR
                 NLM=(NPAR*NPAR+3*NPAR)/2
                 GOTO 24
              ENDIF
           ENDDO
           PRINT*,'PARAMETER NOT FOUND!'
 25        CONTINUE
           
C ---------------------------------------------------------------------- END FIX
        GOTO 380
C ---------------------------------------------------------------------- IF REL

 	ELSEIF(STRING1.EQ.'REL'.OR.STRING1.EQ.'rel') THEN
           
 26        WRITE(6, FMT='(A7)', ADVANCE='NO') 'Label: '
           READ(5, FMT='(A10)', ADVANCE='YES', ERR=25) FLBL

           IF(FLBL.EQ.'') GOTO 27
           DO I=1,NDIM
              IF(FLBL.EQ.LABELS(I)) THEN
                 WRITE(6,FMT='(A19,2X,A10)') 'releasing parameter ',FLBL
                 RNF(I)='R'
                 DO J=1,NPAR
                    IF(IVN(J).GT.I) THEN
                       IVN((J+1):(NPAR+1))=IVN(J:NPAR)
                       RNF(IVN((J+1):(NPAR+1)))=RNF(IVN(J:NPAR))
                       IVN(J)=I
                       EXIT
                    ENDIF
                 ENDDO
                    
                 NPAR=NPAR+1
                 NP=NPAR
                 MP=NP+1
                 NP2=2*NPAR
                 NPAR2=NPAR*NPAR
                 NLM=(NPAR*NPAR+3*NPAR)/2
                 GOTO 26
              ENDIF
           ENDDO
           PRINT*,'PARAMETER NOT FOUND!'
 27        CONTINUE
           
C ---------------------------------------------------------------------- END FIX
        GOTO 380
C ---------------------------------------------------------------------- IF VARY

 	ELSEIF(STRING1.EQ.'VARY'.OR.STRING1.EQ.'vary') THEN

C ---	CHANGE only variables that are 'stopped' and work on with 'running' variables
      
      PRINT 1000
        CALL BTRAFO(NDIM,X0,SU)
        CALL FUNC(SU,NDIM,F_START); NCT=NCT+1
      PRINT*,'NEW FUNCTION VALUE: ',F_START

C ----------------------------------------------------------------------   END VARY
      GOTO 380
C ----------------------------------------------------------------------   IF SIMPLEX

 	ELSEIF(STRING1.EQ.'SIMPLEX'.OR.STRING1.EQ.'simplex') THEN

C        PRINT*,'SIMPLEX DISPLACEMENT (.05 IS GOOD CHOICE)? '
C        READ*,XLA
C        PRINT*,'INITIAL FRACTIONAL CONVERGENCE FTOL? (1.D-4 TO START)'
C        READ*,FTOL

      XLA=.5
      FTOL=1.D-5

C ---   INITIALIZE SIMPLEX
c ---	SELECTION OF PARAMETERS --> IV()

      IV(1:NP)=IVN(1:NP); X(1:NDIM)=X0(1:NDIM)

        DO I=1,MP
        PS(I,1:NP)=X0(IV(1:NP))
        IF(I.NE.1) PS(I,I-1)=PS(I,I-1)+XLA !can be improved using vector xla
        ENDDO

        CALL BTRAFO(NDIM,X,SU)
        CALL FUNC(SU,NDIM,YS(1)); NCT=NCT+1
        
        DO I=2,MP
           DO J=1,NP
            X(IV(J))=PS(I,J)
           ENDDO
        CALL BTRAFO(NDIM,X,SU)
        CALL FUNC(SU,NDIM,YS(I)); NCT=NCT+1
        ENDDO

        CALL SIMPLEX(NCT,X0,IV,PS,YS,MP,NP,NPAR,NDIM,FTOL,ITER)
        
        WRITE(6,8010)
      WRITE(6,8020) NCT,YS(1),LABELS(1),SU(1),RNF(1),FTOL
      DO I=2,NDIM
      WRITE(6,8021) LABELS(I),SU(I),RNF(I)
      ENDDO
        
        CALL FTRAFO(NDIM,X0,SU)
        F_START=YS(1)

c ------------------------------------------------------------------	END SIMPLEX 
        GOTO 380

C ----------------------------------------------------------------------   IF SUBPLX

 	ELSEIF(STRING1.EQ.'SUBPLX'.OR.STRING1.EQ.'subplx') THEN

!!!!!!!!!!!!!!!!!!!!! needs to be worked over !!!!!!!!!!!!!!!!!!!!!!!!!!!
        PRINT*,'removed SUBPLX option'
        
        GOTO 380

c ---	SELECTION OF PARAMETERS --> IV()

        NP=NPAR
      IV(1:NP)=IVN(1:NP)
        X(1:NP)=X0(IV(1:NP))

C ---   INITIALIZE SUBPLX

        ic=0

        CALL BTRAFO(NPAR,X0,SU)
        CALL FUNC(SU,NPAR,FX)
      NCT=NCT+1

c ----------------------------------------------------------
c --- Set subplx's operating mode.

c --- Variables that define the sequence of tolerances.
c --- If tol = 0 is used, subplx will optimize to the limits
c --- of machine precision.
c --- See subplx comments for description of tol.

      tol1=.1
      tol2=1.e-3
      tolfac=.1

c --- Variables that define the sequence of maximum number of
c --- function evaluations.
c --- See subplx comments for description of maxnfe.

      nf1=100
      nf2=1000
      nfinc=100

c --- Set initial stepsizes for optimization.
c --- See subplx comments for description of scale.

      tol = tol1
      maxnfe = nf1

c ----------------------------------------------------------

       scl=.1
       scale(1) = -abs(scl)

c ---  First call to subplx so continuation mode is off.
        mdcont = 0

c ---  Using default options so user options mode is off.
        mduser = 0

c ---  Using optimization so single-step mode is off.
        mdsing = 0

        write (6,1010) 
 1010   format ('  max calls',t15,' tolerance',t30
     *  ,'        value',t45,' calls',t60,' flag',/)

   20   continue

        mode = 4*mdsing + 2*mduser + mdcont

cccc        call subplx (np,tol,maxnfe,mode,scale,x,ys(1),nfe,s_work,i_s_work,ifl)

c Print intermediate results.

        write (6,1020) maxnfe,tol,ys(1),nfe,ifl
 1020   format (t5,i6,t15,e10.2,t30,e13.5,t45,
     *          i6,t60,i5)

c Check ifl to see if done or which termination
c test needs to be reset before resuming optimization.

        if (ifl .eq. -1) then
          if (maxnfe .ge. nf2) go to 30
          maxnfe = maxnfe+nfinc
        else if (ifl .eq. 0) then
          if (tol .le. tol2) go to 30
          tol = tol*tolfac
        else
          go to 30
        end if

c Resume optimization in continuation mode.
c
        mdcont = 1
        ic=ic+1
        go to 20

c Print optimization results.
   30   continue

        PRINT*,' CYCLES:',IC
        WRITE(6,8010)
      WRITE(6,8020) NFE,YS(1),LABELS(1),(SU(1)),TOL
      DO I=2,NDIM
      WRITE(6,8021) LABELS(I),SU(I)
      ENDDO
        
c --- update estimates

        CALL FTRAFO(NDIM,X0,SU); NCT=NCT+NFE
        X0=X(IV); F_START=YS(1)

c ------------------------------------------------------------------	END SUBPLX 
      
      GOTO 380

c ------------------------------------------------------------------	GLOBAL

      ELSEIF(STRING1.EQ.'GLOBAL'.OR.STRING1.EQ.'global') THEN
      CALL GLOBAL(NDIM,NPAR,NCT,X0,F_START,COVM)
      
        PAUSE 'pausing before removing plot'
        CALL SYSTEM('killall gnuplot')
c ------------------------------------------------------------------	END GLOBAL

        GOTO 380

c ------------------------------------------------------------------	BERKSON     

      ELSEIF(STRING1.EQ.'BERK'.OR.STRING1.EQ.'berk') THEN
        CALL BERKSON(NDIM,NPAR,NCT,X0,F_START,MEANV,COVM,XP)
      IMCMC=0

c     ------------------------------------------------------------------	END BERKSON    

        GOTO 380

c ------------------------------------------------------------------	MCMC       
c
c *** to be incorporated into mcmc_raneff!
c
c	ELSEIF(STRING1.EQ.'MCMC'.OR.STRING1.EQ.'mcmc') THEN
c        CALL MCMC(NDIM,NPAR,NCT,X0,F_START,MEANV,COVM,XP)
c	IMCMC=0
c ------------------------------------------------------------------	END MCMC     
c
c        GOTO 380

c ------------------------------------------------------------------	RANEFF       

      ELSEIF(STRING1.EQ.'RANEFF'.OR.STRING1.EQ.'raneff') THEN

        PRINT*,'note: raneff algorithm uses a mv normal propsoal'

        PRINT*,'gamma (G) or log-normal (L) distributions:'
      READ (5, FMT='(A1)', ADVANCE='YES') ANSWER
        
        L=0
        DO I=1,NDIM
           if(rnf(i).eq.'H') l=l+1
        enddo

        IF(L.EQ.0) THEN
         PRINT*,'unable to perform hetero-mcmc, procede with plain mcmc'
         CALL MCMC_HET_GAM(NDIM,NPAR,L,NCT,X0,F_START,MEANV,COVM,XP)
        ELSEIF(ANSWER.EQ.'G') THEN
         CALL MCMC_HET_GAM(NDIM,NPAR,L,NCT,X0,F_START,MEANV,COVM,XP)
        ELSE
         CALL MCMC_HET_LNORM(NDIM,NPAR,L,NCT,X0,F_START,MEANV,COVM,XP)
        ENDIF
       IMCMC=0
      
c ------------------------------------------------------------------	END RANEFF     

        GOTO 380

c ------------------------------------------------------------------	GRADIENT       

      ELSEIF(STRING1.EQ.'GRAD'.OR.STRING1.EQ.'grad') THEN
      CALL MIGRAD(NDIM,NPAR,NCT,X0,F_START)
      
c ------------------------------------------------------------------	END GRADIENT     

      GOTO 380

c ------------------------------------------------------------------	PROFILE LKH CI

      ELSEIF(STRING1.EQ.'CI'.OR.STRING1.EQ.'ci') THEN

      CALL PLKHCI(X0,F_START,NPAR,NDIM,NCT)
      
c ------------------------------------------------------------------	END PROFILE LKH CI     

      GOTO 380

c ------------------------------------------------------------------	BEGIN NEWTON    
      ELSEIF(STRING1.EQ.'NEWTON'.OR.STRING1.EQ.'newton') THEN

      PRINT 1000
        WRITE(6,8003,ADVANCE='NO')
        READ(5,*) EPS, NWT_MAX,I_RETRACT
        IF(EPS.EQ.0.) THEN
           EPS=0.01; NWT_MAX=10;I_HESSE=1
        ENDIF
      I_HESSE=0
      INEW=0

      CALL NEWTON(NCT,X0,F_START,EPS,I_HESSE,I_RETRACT,INEW,NWT_MAX)	
      GOTO 380
      
c ------------------------------------------------------------------	END NEWTON   


      ELSEIF(STRING1.EQ.'BOOT'.OR.STRING1.EQ.'boot') THEN

c ------------------------------------------------------------------	BOOTSTRAP
      OPEN(9,FILE='boot_',STATUS='unknown')

        PRINT 1000
c	PRINT*,'GO THROUGH NEWTON ONCE BEFORE BOOTSTRAPPING (TEST)?'

        PRINT 8004
      READ (5, FMT='(A10)', ADVANCE='YES') ANSWER
      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y') THEN
        IBOOT=0
      IFLAG=2
        EPS=.001D0               !TOLERANCE ON DERIVATIVES
      PRINT 1000
      WRITE(6,8003,ADVANCE='NO')
        read*,EPS, NWT_MAX, I_RETRACT
      I_HESSE=0
      INEW=0

      CALL NEWTON(NCT,X0,F_START,EPS,I_HESSE,I_RETRACT,INEW,NWT_MAX)	
        ENDIF

        PRINT 1000
      PRINT*,'BEGIN BOOTSTRAP:'

        PRINT 1000
        PRINT 8002
        READ*,IBSIZE

      DO IB=1,IBSIZE

      PRINT 1000
      PRINT*,'--- BOOTSTRAP SAMPLE NO ',IB

        XS0 = X0

      NCT=0
      IFLAG=2
      IBOOT=1 			!CREATE BOOTSTRAP SAMPLE
        CALL BTRAFO(NDIM,XS0,SU)
        CALL FUNC(SU,NDIM,F_START)

      PRINT*,'STARTING AT:',F_START
      NCT=NCT+1
      IBOOT=2				!SAVE SAMPLE, WRITE OUTPUT ON UNIT 9

c33	CALL MIGRAD(NDIM,NPAR,NCT,XS0,F_START)
      CALL MIGRAD(NDIM,NPAR,NCT,XS0,F_START)
      
        IFLAG=3
        CALL BTRAFO(NDIM,XS0,SU)
        CALL FUNC(SU,NDIM,F_START)

c	IHESSE=1
c	INEW=1
c	CALL NEWTON(NCT,XS0,F_START,EPS,IHESSE,INEW,NWT_MAX)
c	IF(INEW.EQ.3) GOTO 33
      
      PRINT 1000
      PRINT*,'--- END OF BOOTSTRAP SAMPLE COMPUTATION NO ',IB
      ENDDO
      CLOSE(9)
      GOTO 380
       
c ------------------------------------------------------------------	FIT       

 	ELSEIF(STRING1.EQ.'FIT'.OR.STRING1.EQ.'fit') THEN
        PRINT 1000
        IFLAG=3
        CALL BTRAFO(NDIM,X0,SU)
      CALL FUNC(SU,NDIM,F0); NCT=NCT+1
      IFLAG=2
      GOTO 380

c ------------------------------------------------------------------	END FIT       

c ------------------------------------------------------------------	LOG      
 	ELSEIF(STRING1.EQ.'log'.OR.STRING1.EQ.'LOG') THEN
      

        CALL BTRAFO(NDIM,X0,SU)   
        CALL FUNC(SU,NDIM,F_START)
      NCT=NCT+1

        WRITE(6,FMT='(A28)') 'adding to log file: bhat.log'

        OPEN(7,FILE='bhat.log',status='unknown',position='append')

        DO I=1,NDIM 
         WRITE(7, 1003)  LABELS(I),SU(I),RNF(I),SE(I),S1(I),S2(I)
        END DO 
         WRITE(7, FMT='(A4,2x,F16.6)') 'END!',F_START

      CLOSE(7)

      GOTO 380
c ------------------------------------------------------------------	END LOG      

c ------------------------------------------------------------------	EXIT      
 	ELSEIF(STRING1.EQ.'EXIT'.OR.STRING1.EQ.'exit'.OR.STRING1.EQ.
     1  'quit'.OR. STRING1.EQ.'QUIT') THEN
      

        CALL BTRAFO(NDIM,X0,SU)   
        CALL FUNC(SU,NDIM,F_START)
      NCT=NCT+1

        PRINT*,'exiting and creating file: bhat.out'

        OPEN(7,FILE="bhat.out",status="unknown")

        DO I=1,NDIM 
         WRITE(7, 1003)  LABELS(I),SU(I),RNF(I),SE(I),S1(I),S2(I)
        END DO 
         WRITE(7, FMT='(A4,2x,F16.6)') 'END!',F_START

      CLOSE(7)

      STOP
c ------------------------------------------------------------------	END EXIT      

      ELSE

      PRINT*,'PLEASE ENTER A COMMMAND (or type help for a menu):' 
      GOTO 380
  
      ENDIF	!SIMPLEX, GRAD, NEWTON

        STOP

 1000   FORMAT(/)
 1001   FORMAT (F10.0,A5,4F10.0)
 1002   FORMAT (F3.0,2x,A5,E12.5,A5,3E12.5)
 1003   FORMAT (2x,A5,E12.5,A5,3E12.5)

8002	FORMAT('BOOTSTRAP SAMPLE SIZE (B)? ')
8003    FORMAT('DERIV. LIMIT (0.01), ITERATIONS (10),
     $ RETRACTION(0/1): ')
8004    FORMAT('TEST CONVERGENCE WITH NEWTON/RAPHSON (Y/N): ')
8005    FORMAT('WHICH PARAMETER: ')
8010  FORMAT(/,'CALLS',T11,'VALUE',T18,'LBL',T26,'ESTIMATES',T40,
     & 'STATE',T48,'  TOLERANCE')
8020  FORMAT(/,I4,T6,E12.5,T20,A5,E12.5,T44,A1,T47,E12.4)
8021  FORMAT(T20,A5,E12.5,T44,A1)
8030  FORMAT(/,T5,E12.4,T20,30E12.6) 
8072  FORMAT(/,' TRANSF. POINTS:',T20,30F8.4)
        END

!!      EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT_INI(NPAR)
      SUBROUTINE GNUPLOT_INI(NPAR)
      INTEGER NPAR
c --- initialize gnuplot graphics

!!      INCLUDE '/usr/pgi/linux86/include/lib3f.h'
      CALL SYSTEM('if [ -f mcmc.out ]; '//
     $     'then mv -f mcmc.out mcmc.out.old; fi')
      CALL SYSTEM('if [ -f mcmc.het.out ]; '//
     $     'then mv -f mcmc.out mcmc.het.out.old; fi')

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


!!      EXTRINSIC(HPF_SERIAL) SUBROUTINE GNUPLOT()
      SUBROUTINE GNUPLOT()
!!      INCLUDE '/usr/pgi/linux86/include/lib3f.h'
      CALL SYSTEM('gnuplot mcmc.gnu &')
      END SUBROUTINE GNUPLOT
