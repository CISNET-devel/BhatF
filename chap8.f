      PROGRAM GOLDSEAR
C     NUMERICAL METHODS: FORTRAN Programs, (c) John H. Mathews 1994
C     To accompany the text:
C     NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
C     Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
C     This free software is complements of the author.
C
C     Algorithm 8.1 (Golden Search for a Minimum).
C     Section 8.1, Minimization of a Function, Page 413
C
      PARAMETER(Delta=1E-5,Epsilon=1E-5)
      INTEGER K
      REAL A,A0,B,B0,P
      CHARACTER ANS*1
      EXTERNAL F
10    CALL INPUTS(A0,B0,A,B)
      CALL SEARCH(F,A,B,P,K,Delta,Epsilon)
      CALL RESULT(F,A0,B0,A,B,P,Delta,Epsilon,K)
      WRITE(9,*)' '
      WRITE(9,*)'WANT TO FIND THE MINIMUM ON ANOTHER INTERVAL ? '
      READ(9,'(A)') ANS
      IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') GOTO 10
      STOP
      END

      REAL FUNCTION F(X)
      REAL X
      F=X*X-SIN(X)
      RETURN
      END

      SUBROUTINE PRINTFUN
      WRITE(9,*)'F(X)  =  X*X - SIN(X)'
      RETURN
      END

      SUBROUTINE SEARCH(F,A,B,P,K,Delta,Epsilon)
      INTEGER K
      REAL A,B,C,D,Delta,Epsilon,H,P,Rone,Rtwo,YP,YA,YB,YC,YD
      EXTERNAL F
      WRITE(9,*) '               A                C                B'
      Rone=(SQRT(5.0)-1)/2
      Rtwo=Rone*Rone
      H=B-A
      YA=F(A)
      YB=F(B)
      C=A+Rtwo*H
      D=A+Rone*H
      YC=F(C)
      YD=F(D)
      K=1
      WHILE (ABS(YB-YA).GT.Epsilon .OR. H.GT.Delta)
        IF (YC.LT.YD) THEN
          B=D
          YB=YD
          D=C
          YD=YC
          H=B-A
          C=A+Rtwo*H
          YC=F(C)
        ELSE
          A=C
          YA=YC
          C=D
          YC=YD
          H=B-A
          D=A+Rone*H
          YD=F(D)
        ENDIF
        K=K+1
        WRITE(9,1000) A,C,B
      REPEAT
      P=A
      YP=YA
      IF (YB.LT.YA) THEN
        P=B
        YP=YB
      ENDIF
      PAUSE
      RETURN
1000  FORMAT(5X,3F17.7)
      END

      SUBROUTINE XSEARCH(F,A,B,P,K,Delta,Epsilon)
C     This subroutine uses labeled DO loop(s).
      INTEGER K
      REAL A,B,C,D,Delta,Epsilon,H,P,Rone,Rtwo,YP,YA,YB,YC,YD
      EXTERNAL F
      WRITE(9,*) '               A                C                B'
      Rone=(SQRT(5.0)-1)/2
      Rtwo=Rone*Rone
      H=B-A
      YA=F(A)
      YB=F(B)
      C=A+Rtwo*H
      D=A+Rone*H
      YC=F(C)
      YD=F(D)
      K=1
10    IF (ABS(YB-YA).GT.Epsilon .OR. H.GT.Delta) THEN
        IF (YC.LT.YD) THEN
          B=D
          YB=YD
          D=C
          YD=YC
          H=B-A
          C=A+Rtwo*H
          YC=F(C)
        ELSE
          A=C
          YA=YC
          C=D
          YC=YD
          H=B-A
          D=A+Rone*H
          YD=F(D)
        ENDIF
        K=K+1
        WRITE(9,1000) A,C,B
        GOTO 10
      ENDIF
      P=A
      YP=YA
      IF (YB.LT.YA) THEN
        P=B
        YP=YB
      ENDIF
      PAUSE
      RETURN
1000  FORMAT(5X,3F17.7)
      END

      SUBROUTINE INPUTS(A0,B0,A,B)
      INTEGER I
      REAL A0,B0,A,B
      DO 10 I=1,16
        WRITE(9,*)' '
10    CONTINUE
      WRITE(9,*)'    GOLDEN SEARCH FOR FINDING THE MINIMUM OF THE'
      WRITE(9,*)' '
      WRITE(9,*)'UNIMODAL FUNCTION  F(X)  OVER THE INTERVAL  [A,B].'
      WRITE(9,*)' '
      WRITE(9,*)'SEQUENCES {A } AND {B } CONSTRUCTED, WHICH BRACKET
     + THE MINIMUM,'
      WRITE(9,*)'            K        K'
      WRITE(9,*)' '
      WRITE(9,*)'I.E.  A  < X   < B .  TERMINATION OCCURS WHEN  B - A  <
     + DELTA.'
      WRITE(9,*)'        K   MIN   K                             N   N'
      WRITE(9,*)' '
      WRITE(9,*)'THE GIVEN FUNCTION F(X) IS:'
      WRITE(9,*)' '
      CALL PRINTFUN
      WRITE(9,*)' '
      WRITE(9,*)'ENTER LEFT  ENDPOINT  A = '
      READ(9,*) A
      A0=A
      WRITE(9,*)' '
      WRITE(9,*)'ENTER RIGHT ENDPOINT  B = '
      READ(9,*) B
      B0=B
      RETURN
      END

      SUBROUTINE RESULT(F,A0,B0,A,B,P,Delta,Epsilon,K)
      INTEGER I,K
      REAL A0,B0,A,B,H,Delta,Epsilon,P
      EXTERNAL F
      WRITE(9,*)' '
      WRITE(9,*)'     THE GOLDEN SEARCH WAS USED TO FIND THE MINIMUM'
      WRITE(9,*)' '
      WRITE(9,*)'OF THE FUNCTION  F(X)  OVER THE INTERVAL  [A,B].'
      WRITE(9,*)' '
      WRITE(9,*)' '
      WRITE(9,*)'     IT TOOK ',K,' ITERATIONS TO MINIMIZE'
      WRITE(9,*)' '
      CALL PRINTFUN
      WRITE(9,*)' '
      WRITE(9,*)'OVER THE INTERVAL  [',A0,',',B0,'  ].'
      WRITE(9,*)' '
      WRITE(9,*)'THE  MINIMUM  VALUE F(P)'
      WRITE(9,*)' '
      WRITE(9,*)'          OCCURS  AT  P  = ',P
      WRITE(9,*)' '
      H=B-A
      WRITE(9,*)'    THE  ACCURACY  IS    +-',H
      WRITE(9,*)' '
      WRITE(9,*)'THE FUNCTION VALUE IS:'
      WRITE(9,*)' '
      WRITE(9,*)'F(',P,'  )  = ',F(P)
      RETURN
      END














      PROGRAM NELDMEAD
C     NUMERICAL METHODS: FORTRAN Programs, (c) John H. Mathews 1994
C     To accompany the text:
C     NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
C     Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
C     This free software is complements of the author.
C
C     Algorithm 8.2 (Nelder-Mead's Minimization Method).
C     Section 8.1, Minimization of a Function, Page 414
C
      PARAMETER(MaxN=6,Epsilon=1E-5)
      INTEGER Count,Hi,J,K,Lo,N
      REAL Norm,V,Y
      CHARACTER ANS*1
      DIMENSION V(0:MaxN,1:MaxN),Y(1:MaxN)
      EXTERNAL F
10    CALL INPUT(F,V,Y,N)
      CALL NELDER(F,V,Y,N,Epsilon,Norm,Lo,Hi,Count)
      CALL OUTPUT(V,Y,N,Norm,Lo,Hi,Count)
      WRITE(9,*)' '
      WRITE(9,*)'WANT TO TRY NEW STARTING VERTICES ? '
      READ(9,'(A)') ANS
      IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') GOTO 10
      END

      REAL FUNCTION F(P)
      PARAMETER(MaxN=6)
      REAL P,X,Y,Z,U,V,W
      DIMENSION P(1:MaxN)
        X=P(1)
        Y=P(2)
        Z=P(3)
        U=P(4)
        V=P(5)
        W=P(6)
        F=X*X - 4*X + Y*Y - Y - X*Y
      RETURN
      END

      SUBROUTINE PRINTFUN(N)
      INTEGER K,N
      CHARACTER C*2,FUN*80,L*6,R*4
      DIMENSION C(1:6)
        N = 2
      FUN ='X*X - 4*X + Y*Y - Y - X*Y'
      C(1)='(X'
      C(2)=',Y'
      C(3)=',Z'
      C(4)=',U'
      C(5)=',V'
      C(6)=',W,'
        L ='     F'
        R =') = '
      WRITE(9,999) L,C(1),(C(K),K=2,N),R,FUN
999   FORMAT(9A)
      RETURN
      END

      SUBROUTINE NELDER(F,V,Y,N,Epsilon,Norm,Lo,Hi,Count)
      PARAMETER(MaxN=6,Min=10,Max=200)
      INTEGER Count,Hi,Ho,J,K,Li,Lo,N
      REAL C,E,Epsilon,M,Norm,R,S,V,Y,YM,YC,YE,YR,Z
      DIMENSION C(1:MaxN),E(1:MaxN),M(1:MaxN),R(1:MaxN)
      DIMENSION V(0:MaxN,1:MaxN),Y(1:MaxN),Z(1:MaxN)
      EXTERNAL F
      Count=0
      Lo=0                       !Order the vertices:
      Hi=0
      DO J=1,N
        IF (Y(J).LT.Y(Lo)) Lo=J
        IF (Y(J).GT.Y(Hi)) Hi=J
      ENDDO
      Li=Hi
      Ho=Lo
      DO J=0,N
        IF (J.NE.Lo .AND. Y(J).LT.Y(Li)) Li=J
        IF (J.NE.Hi .AND. Y(J).GT.Y(Ho)) Ho=J
      ENDDO
200   WHILE (Y(Hi).GT.(Y(Lo)+Epsilon) .AND.
     +      (Count.LT.Max .OR. Count.LT.Min))
        DO K=1,N                 !The main loop starts here:
          S=0                    !Form the new points:
          DO J=0,N
            S=S+V(J,K)
          ENDDO
          M(K)=(S-V(Hi,K))/N
        ENDDO
        DO K=1,N
          R(K)=2*M(K)-V(Hi,K)
        ENDDO
        YR=F(R)
        IF (YR.LT.Y(Ho)) THEN    !Improve the simplex:
          IF (Y(Li).LT.YR) THEN
            DO K=1,N             !Replace a vertex:
              V(Hi,K)=R(K)
            ENDDO
            Y(Hi)=YR
          ELSE
            DO K=1,N
              E(K)=2*R(K)-M(K)
            ENDDO
            YE=F(E)
              IF (YE.LT.Y(Li)) THEN
                DO K=1,N
                  V(Hi,K)=E(K)
                ENDDO
                Y(Hi)=YE
              ELSE
                DO K=1,N         !Replace a vertex:
                  V(Hi,K)=R(K)
                ENDDO
                Y(Hi)=YR
              ENDIF
          ENDIF
        ELSE
          IF (YR.LT.Y(Hi)) THEN
            DO K=1,N             !Replace a vertex:
              V(Hi,K)=R(K)
            ENDDO
            Y(Hi)=YR
          ENDIF
          DO K=1,N
            C(K)=(V(Hi,K)+M(K))/2
          ENDDO
          YC=F(C)
          IF (YC.LT.Y(Hi)) THEN
            DO K=1,N
              V(Hi,K)=C(K)
            ENDDO
            Y(Hi)=YC
          ELSE
            DO J=0,N             !Shrink the simplex:
              IF (J.NE.Lo) THEN
                DO K=1,N
                  V(J,K)=(V(J,K)+V(Lo,K))/2
                  Z(K)=V(J,K)
                ENDDO
                Y(J)=F(Z)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
        Count=Count+1
        Lo=0                     !Order the vertices:
        Hi=0
        DO J=1,N
          IF (Y(J).LT.Y(Lo)) Lo=J
          IF (Y(J).GT.Y(Hi)) Hi=J
        ENDDO
        Li=Hi
        Ho=Lo
        DO J=0,N
          IF (J.NE.Lo .AND. Y(J).LT.Y(Li)) Li=J
          IF (J.NE.Hi .AND. Y(J).GT.Y(Ho)) Ho=J
        ENDDO
      REPEAT                     !End of the main loop.
      Norm=0                     !Determine the size of the simplex:
      DO J=0,N
        S=0
        DO K=1,N
          S=S+(V(Lo,K)-V(J,K))*(V(Lo,K)-V(J,K))
        ENDDO
        IF (S.GT.Norm) Norm=S
      ENDDO
      Norm=SQRT(Norm)
      RETURN
      END

      SUBROUTINE XNELDER(F,V,Y,N,Epsilon,Norm,Lo,Hi,Count)
      PARAMETER(MaxN=6,Min=10,Max=200)
C     This subroutine uses labeled DO loop(s).
      INTEGER Count,Hi,Ho,J,K,Li,Lo,N
      REAL C,E,Epsilon,M,Norm,R,S,V,Y,YM,YC,YE,YR,Z
      DIMENSION C(1:MaxN),E(1:MaxN),M(1:MaxN),R(1:MaxN)
      DIMENSION V(0:MaxN,1:MaxN),Y(1:MaxN),Z(1:MaxN)
      EXTERNAL F
      Count=0
      Lo=0                       !Order the vertices:
      Hi=0
      DO 100 J=1,N
        IF (Y(J).LT.Y(Lo)) Lo=J
        IF (Y(J).GT.Y(Hi)) Hi=J
100   CONTINUE
      Li=Hi
      Ho=Lo
      DO 110 J=0,N
        IF (J.NE.Lo .AND. Y(J).LT.Y(Li)) Li=J
        IF (J.NE.Hi .AND. Y(J).GT.Y(Ho)) Ho=J
110   CONTINUE
200   IF (Y(Hi).GT.(Y(Lo)+Epsilon) .AND.
     +   (Count.LT.Max .OR. Count.LT.Min)) THEN
        DO 310 K=1,N             !The main loop starts here:
          S=0                    !Form the new points:
          DO 300 J=0,N
            S=S+V(J,K)
300       CONTINUE
          M(K)=(S-V(Hi,K))/N
310     CONTINUE
        DO 320 K=1,N
          R(K)=2*M(K)-V(Hi,K)
320     CONTINUE
        YR=F(R)
        IF (YR.LT.Y(Ho)) THEN    !Improve the simplex:
          IF (Y(Li).LT.YR) THEN
            DO 400 K=1,N         !Replace a vertex:
              V(Hi,K)=R(K)
400         CONTINUE
            Y(Hi)=YR
          ELSE
            DO 410 K=1,N
              E(K)=2*R(K)-M(K)
410         CONTINUE
            YE=F(E)
              IF (YE.LT.Y(Li)) THEN
                DO 420 K=1,N
                  V(Hi,K)=E(K)
420             CONTINUE
                Y(Hi)=YE
              ELSE
                DO 430 K=1,N     !Replace a vertex:
                  V(Hi,K)=R(K)
430             CONTINUE
                Y(Hi)=YR
              ENDIF
          ENDIF
        ELSE
          IF (YR.LT.Y(Hi)) THEN
            DO 440 K=1,N         !Replace a vertex:
              V(Hi,K)=R(K)
440         CONTINUE
            Y(Hi)=YR
          ENDIF
          DO 450 K=1,N
            C(K)=(V(Hi,K)+M(K))/2
450       CONTINUE
          YC=F(C)
          IF (YC.LT.Y(Hi)) THEN
            DO 460 K=1,N
              V(Hi,K)=C(K)
460         CONTINUE
            Y(Hi)=YC
          ELSE
            DO 480 J=0,N         !Shrink the simplex:
              IF (J.NE.Lo) THEN
                DO 470 K=1,N
                  V(J,K)=(V(J,K)+V(Lo,K))/2
                  Z(K)=V(J,K)
470             CONTINUE
                Y(J)=F(Z)
              ENDIF
480         CONTINUE
          ENDIF
        ENDIF
        Count=Count+1
        Lo=0                     !Order the vertices:
        Hi=0
        DO 500 J=1,N
          IF (Y(J).LT.Y(Lo)) Lo=J
          IF (Y(J).GT.Y(Hi)) Hi=J
500     CONTINUE
        Li=Hi
        Ho=Lo
        DO 510 J=0,N
          IF (J.NE.Lo .AND. Y(J).LT.Y(Li)) Li=J
          IF (J.NE.Hi .AND. Y(J).GT.Y(Ho)) Ho=J
510     CONTINUE
        GOTO 200
      ENDIF                      !End of the main loop.
      Norm=0                     !Determine the size of the simplex:
      DO 610 J=0,N
        S=0
        DO 600 K=1,N
          S=S+(V(Lo,K)-V(J,K))*(V(Lo,K)-V(J,K))
600     CONTINUE
        IF (S.GT.Norm) Norm=S
610   CONTINUE
      Norm=SQRT(Norm)
      RETURN
      END

      SUBROUTINE INPUT(F,V,Y,N)
      PARAMETER(MaxN=6)
      INTEGER I,J,K,N
      REAL V,Y,Z
      DIMENSION V(0:MaxN,1:MaxN),Y(1:MaxN),Z(1:MaxN)
      EXTERNAL F
      DO 10 I=1,16
        WRITE(9,*)' '
10    CONTINUE
      WRITE(9,*)'   THE NELDER-MEAD SIMPLEX METHOD OR "POLYTOPE METHOD"
     + IS'
      WRITE(9,*)' '
      WRITE(9,*)'USED FOR FINDING THE MINIMUM OF THE FUNCTION  F(V)'
      WRITE(9,*)' '
      WRITE(9,*)'FOR CONVENIENCE, THE FUNCTION F(V) CAN BE EXPRESSED USI
     +NG'
      WRITE(9,*)' '
      WRITE(9,*)'THE VARIABLES X = v , Y = v , Z = v , U = v , V = v , W
     + = v .'
      WRITE(9,*)'                   1       2       3       4       5
     +    6 '
      WRITE(9,*)' '
      CALL PRINTFUN(N)
      WRITE(9,*)' '
      WRITE(9,*)'YOU MUST SUPPLY',N+1,' LINEARLY INDEPENDENT'
      WRITE(9,*)' '
      IF (N.EQ.2) THEN
        WRITE(9,*)'STARTING POINTS   V = (v   ,v   )   FOR  J=0,1,3'
        WRITE(9,*)'                   J    J,1  J,2'
      ELSEIF (N.EQ.3) THEN
        WRITE(9,*)'STARTING POINTS   V = (v   ,v   ,v   )   FOR  J=0,1,3
     +,4'
        WRITE(9,*)'                   J    J,1  J,2  J,3'
      ELSE
        WRITE(9,*)'STARTING POINTS   V = (v   ,v   ,...,v   )   FOR  J=0
     +,1,...,N'
        WRITE(9,*)'                   J    J,1  J,2      J,N'
        WRITE(9,*)'WHERE N =',N
      ENDIF
      DO 30 J=0,N
        WRITE(9,*)' '
        WRITE(9,*)'GIVE COORDINATES OF POINT  V'
        WRITE(9,*)'                  ',J
      DO 20 K=1,N
        WRITE(9,*)'V(',J,',',K,') = '
        READ(9,*) V(J,K)
        Z(K)=V(J,K)
20    CONTINUE
      Y(J)=F(Z)
30    CONTINUE
      RETURN
      END

      SUBROUTINE OUTPUT(V,Y,N,Norm,Lo,Hi,Count)
      PARAMETER(MaxN=6,Epsilon=1E-5)
      INTEGER Count,Hi,I,K,Lo,N
      REAL Norm,V,Y
      DIMENSION V(0:MaxN,1:MaxN),Y(1:MaxN)
      DO 10 I=1,16
        WRITE(9,*)' '
10    CONTINUE
      WRITE(9,*)' '
      WRITE(9,*)'   THE NELDER-MEAD METHOD WAS USED TO FIND THE MINIMUM
     + OF THE FUNCTION'
      WRITE(9,*)' '
      CALL PRINTFUN(N)
      WRITE(9,*)' '
      WRITE(9,*)'IT TOOK ',Count,' ITERATIONS TO FIND AN APPROXIMATION F
     +OR'
      WRITE(9,*)' '
      IF (N.EQ.2) THEN
        WRITE(9,*)'THE COORDINATES OF THE LOCAL MINIMUM  P = (p ,p )'
        WRITE(9,*)'                                            1  2'
      ELSEIF (N.EQ.3) THEN
        WRITE(9,*)'THE COORDINATES OF THE LOCAL MINIMUM  P = (p ,p ,p )'
        WRITE(9,*)'                                            1  2  3'
      ELSE
        WRITE(9,*)'THE COORDINATES OF THE LOCAL MINIMUM  P = (p ,p ,...,
     +p )'
        WRITE(9,*)'                                            1  2
     + N'
        WRITE(9,*)'WHERE N =',N
      ENDIF
      DO 20 K=1,N
        WRITE(9,*)'P(',K,')  =',V(Lo,K)
20    CONTINUE
      WRITE(9,*)' '
      WRITE(9,*)'THE MAXIMUM DISTANCE TO THE OTHER VERTICES OF THE SIMP
     +LEX IS'
      WRITE(9,*)' '
      WRITE(9,*)'            DP  =',Norm
      WRITE(9,*)' '
      WRITE(9,*)'THE FUNCTION VALUE AT THE MINIMUM POINT IS'
      WRITE(9,*)' '
      WRITE(9,*)'F(P)  =',Y(Lo)
      WRITE(9,*)' '
      WRITE(9,*)'  DF  =',Y(Hi)-Y(Lo),'  IS AN ESTIMATE FOR THE ACCURACY
     +.'
      RETURN
      END














      PROGRAM QMINUM
C     NUMERICAL METHODS: FORTRAN Programs, (c) John H. Mathews 1994
C     To accompany the text:
C     NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
C     Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
C     This free software is complements of the author.
C
C     Algorithm 8.3 (Local Minimum Search: Quadratic Interpolation).
C     Section 8.1, Minimization of a Function, Page 416
C
      PARAMETER(Delta=1E-5,Epsilon=1E-7,Jmax=50,Kmax= 50)
      INTEGER Cond,Count
      REAL Error,H,P0,Pmin,Pold,Y0,Ymin
      CHARACTER ANS*1
      EXTERNAL F,F1
10    CALL INPUTS(P0,Pold)
      CALL QMINIM(F,F1,P0,Delta,Epsilon,Jmax,Kmax,Pmin,H,Ymin,
     +Error,Cond,Count)
      CALL RESULT(Pold,Pmin,H,Ymin,Error,Epsilon,Cond,Count)
      WRITE(9,*)' '
      WRITE(9,*)'WANT TO TRY ANOTHER STARTING POINT ? '
      READ(9,'(A)') ANS
      IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') GOTO 10
      END

      REAL FUNCTION F(X)
      REAL X
      F=X*X*X-X
      RETURN
      END

      REAL FUNCTION F1(X)
      REAL X
      F1=3*X*X-1
      RETURN
      END

      SUBROUTINE PRINTFUN
      WRITE(9,*)'F(X)  =  X*X*X-X'
      RETURN
      END

      SUBROUTINE QMINIM(F,F1,P0,Delta,Epsilon,Jmax,Kmax,Pmin,H,Ymin,
     +Error,Cond,Count)
      PARAMETER(Big=1E9)
      INTEGER Cond,Count,J,Jmax,K,Kmax
      REAL Delta,Epsilon,Error,H,P0,Pmin,Ymin
      REAL D,E0,E1,E2,H0,H1,H2,Hmin,P1,P2,Y0,Y1,Y2
      EXTERNAL F,F1
      K=0
      Error=1
      Cond=0
      Count=0
      H=1
      IF (ABS(P0).GT.1E4) H=ABS(P0)/1E4
      WHILE (K.LT.Kmax .AND. Error.GT.Epsilon .AND. Cond.NE.5)
        IF (F1(P0).GT.0) H=-ABS(H)
        P1=P0+H
        P2=P0+2*H
        Pmin= P0
        Y0=F(P0)
        Y1=F(P1)
        Y2=F(P2)
        Ymin=Y0
        Cond=0
        J=0
        WHILE (J.LT.Jmax .AND. ABS(H).GT.Delta .AND. Cond.EQ.0)
          IF (Y0.LE.Y1) THEN
            P2=P1
            Y2=Y1
            H=H/2
            P1=P0+H
            Y1=F(P1)
          ELSE
            IF (Y2.LT.Y1) THEN
              P1=P2
              Y1=Y2
              H=2*H
              P2=P0+2*H
              Y2=F(P2)
            ELSE
              Cond=-1
            ENDIF
          ENDIF
          J=J+1
          Count=Count+1
          IF (ABS(H).GT.Big .OR. ABS(P0).GT.Big) Cond=5
          WRITE(9,1000) Count,P1,Y1
        REPEAT
        IF (Cond.EQ.5) THEN
          Pmin=P1
          Ymin=F(P1)
        ELSE
          D=4*Y1-2*Y0-2*Y2                      !Start of a long block:
          IF (D.LT.0) THEN
            Hmin=H*(4*Y1-3*Y0-Y2)/D
          ELSE
            Hmin=H/3
            Cond=4
          ENDIF
          Pmin=P0+Hmin
          Ymin=F(Pmin)
          H=ABS(H)
          H0=ABS(Hmin)
          H1=ABS(Hmin-H)
          H2=ABS(Hmin-2*H)
          IF (H0.LT.H) H=H0
          IF (H1.LT.H) H=H1
          IF (H2.LT.H) H=H2
          IF (H.EQ.0)  H=Hmin
          IF (H.LT.Delta) Cond=1
          IF (ABS(H).GT.Big .OR. ABS(Pmin).GT.Big) Cond=5
          E0=ABS(Y0-Ymin)
          E1=ABS(Y1-Ymin)
          E2=ABS(Y2-Ymin)
          IF (E0.NE.0 .AND. E0.LT.Error) Error=E0
          IF (E1.NE.0 .AND. E1.LT.Error) Error=E1
          IF (E2.NE.0 .AND. E2.LT.Error) Error=E2
          IF (E0.EQ.0 .AND. E1.EQ.0 .AND. E2.EQ.0) Error=0
          IF (Error.LT.Epsilon) Cond=2
          P0=Pmin
          K=K+1
          Count=Count+1
        ENDIF                                   !End of the long block.
        WRITE(9,1000) Count,Pmin,Ymin
        IF (Cond.EQ.2 .AND. H.LT.Delta) Cond=3
      REPEAT
      RETURN
1000  FORMAT(5X,I4,5X,F18.7,5X,F18.7)
      END

      SUBROUTINE XMINIM(F,F1,P0,Delta,Epsilon,Jmax,Kmax,Pmin,H,Ymin,
     +Error,Cond,Count)
C     This subroutine uses labeled DO loop(s).
      PARAMETER(Big=1E9)
      INTEGER Cond,Count,J,Jmax,K,Kmax
      REAL Delta,Epsilon,Error,H,P0,Pmin,Ymin
      REAL D,E0,E1,E2,H0,H1,H2,Hmin,P1,P2,Y0,Y1,Y2
      EXTERNAL F,F1
      K=0
      Error=1
      Cond=0
      Count=0
      H=1
      IF (ABS(P0).GT.1E4) H=ABS(P0)/1E4
10    IF (K.LT.Kmax .AND. Error.GT.Epsilon .AND. Cond.NE.5) THEN
        IF (F1(P0).GT.0) H=-ABS(H)
        P1=P0+H
        P2=P0+2*H
        Pmin= P0
        Y0=F(P0)
        Y1=F(P1)
        Y2=F(P2)
        Ymin=Y0
        Cond=0
        J=0
20      IF (J.LT.Jmax .AND. ABS(H).GT.Delta .AND. Cond.EQ.0) THEN
          IF (Y0.LE.Y1) THEN
            P2=P1
            Y2=Y1
            H=H/2
            P1=P0+H
            Y1=F(P1)
          ELSE
            IF (Y2.LT.Y1) THEN
              P1=P2
              Y1=Y2
              H=2*H
              P2=P0+2*H
              Y2=F(P2)
            ELSE
              Cond=-1
            ENDIF
          ENDIF
          J=J+1
          Count=Count+1
          IF (ABS(H).GT.Big .OR. ABS(P0).GT.Big) Cond=5
          WRITE(9,1000) Count,P1,Y1
          GOTO 20
        ENDIF
        IF (Cond.EQ.5) THEN
          Pmin=P1
          Ymin=F(P1)
        ELSE
          D=4*Y1-2*Y0-2*Y2                      !Start of a long block:
          IF (D.LT.0) THEN
            Hmin=H*(4*Y1-3*Y0-Y2)/D
          ELSE
            Hmin=H/3
            Cond=4
          ENDIF
          Pmin=P0+Hmin
          Ymin=F(Pmin)
          H=ABS(H)
          H0=ABS(Hmin)
          H1=ABS(Hmin-H)
          H2=ABS(Hmin-2*H)
          IF (H0.LT.H) H=H0
          IF (H1.LT.H) H=H1
          IF (H2.LT.H) H=H2
          IF (H.EQ.0)  H=Hmin
          IF (H.LT.Delta) Cond=1
          IF (ABS(H).GT.Big .OR. ABS(Pmin).GT.Big) Cond=5
          E0=ABS(Y0-Ymin)
          E1=ABS(Y1-Ymin)
          E2=ABS(Y2-Ymin)
          IF (E0.NE.0 .AND. E0.LT.Error) Error=E0
          IF (E1.NE.0 .AND. E1.LT.Error) Error=E1
          IF (E2.NE.0 .AND. E2.LT.Error) Error=E2
          IF (E0.EQ.0 .AND. E1.EQ.0 .AND. E2.EQ.0) Error=0
          IF (Error.LT.Epsilon) Cond=2
          P0=Pmin
          K=K+1
          Count=Count+1
        ENDIF                                   !End of the long block.
        WRITE(9,1000) Count,Pmin,Ymin
        GOTO 10
      ENDIF
      IF (Cond.EQ.2 .AND. H.LT.Delta) Cond=3
      PAUSE
      RETURN
1000  FORMAT(5X,I4,5X,F18.7,5X,F18.7)
      END

      SUBROUTINE INPUTS(P0,Pold)
      INTEGER I
      REAL P0
      WRITE(9,*)' '
      WRITE(9,*)'     THE QUADRATIC SEARCH METHOD IS USED'
      WRITE(9,*)' '
      WRITE(9,*)'TO FIND A LOCAL MINIMUM OF THE FUNCTION: '
      WRITE(9,*)' '
      CALL PRINTFUN
      WRITE(9,*)' '
      WRITE(9,*)'ONE INITIAL APPROXIMATION P0 IS NEEDED.'
      WRITE(9,*)' '
      WRITE(9,*)'ENTER P0 = '
      READ(9,*) P0
      Pold=P0
      WRITE(9,*)' '
      RETURN
      END

      SUBROUTINE RESULT(Pold,Pmin,H,Ymin,Error,Epsilon,Cond,Count)
      INTEGER Cond,Count,I
      REAL Error,Epsilon,H,Pold,Pmin,Ymin
      DO 10 I=1,18
        WRITE(9,*)' '
10    CONTINUE
      WRITE(9,*)'     THE QUADRATIC SEARCH METHOD WAS USED '
      WRITE(9,*)' '
      WRITE(9,*)'TO FIND A LOCAL MINIMUM OF THE FUNCTION: '
      WRITE(9,*)' '
      CALL PRINTFUN
      WRITE(9,*)' '
      WRITE(9,*)'STARTING WITH THE APPROXIMATION  P0 =',Pold
      WRITE(9,*)' '
      WRITE(9,*)'AFTER ',Count,' ITERATIONS AN APPROXIMATION FOR THE MIN
     +IMUM IS:'
      WRITE(9,*)' '
      WRITE(9,*)'     P  =',Pmin
      WRITE(9,*)' '
      WRITE(9,*)'    DP  =',ABS(H),'  IS THE ESTIMATED ACCURACY FOR P.'
      WRITE(9,*)' '
      WRITE(9,*)'       F(',Pmin,'  )  =',Ymin
      WRITE(9,*)' '
      IF (Cond.EQ.0) THEN
        WRITE(9,*)'CONVERGENCE IS DOUBTFUL BECAUSE THE'
        WRITE(9,*)' '
        WRITE(9,*)'MAXIMUM NUMBER OF ITERATIONS WAS EXCEEDED.'
      ELSEIF (Cond.EQ.1) THEN
        WRITE(9,*)'CONVERGENCE HAS BEEN ACHIEVED BECAUSE '
        WRITE(9,*)' '
        WRITE(9,*)'CONSECUTIVE ABSCISSAS ARE CLOSER THAN ',Delta
      ELSEIF (Cond.EQ.2) THEN
        WRITE(9,*)'CONVERGENCE HAS BEEN ACHIEVED BECAUSE '
        WRITE(9,*)' '
        WRITE(9,*)'CONSECUTIVE ORDINATES ARE CLOSER THAN ',Epsilon
      ELSEIF (Cond.EQ.3) THEN
        WRITE(9,*)'CONVERGENCE HAS BEEN ACHIEVED BECAUSE'
        WRITE(9,*)' '
        WRITE(9,*)'CONSECUTIVE ABSCISSAS ARE CLOSER THAN ',H
        WRITE(9,*)' '
        WRITE(9,*)'CONSECUTIVE ORDINATES ARE CLOSER THAN ',Error
      ELSEIF (Cond.EQ.4) THEN
        WRITE(9,*)'CONVERGENCE IS DOUBTFUL BECAUSE DIVISION BY ZERO WAS
     +ENCOUNTERED.'
        IF (Error.LT.Epsilon/100) THEN
          WRITE(9,*)'HOWEVER, CONSECUTIVE ORDINATES ARE CLOSE.'
        ENDIF
      ELSEIF (Cond.EQ.5) THEN
        WRITE(9,*)'CONVERGENCE IS DOUBTFUL, H IS TOO LARGE, H  =',H
        WRITE(9,*)' '
        WRITE(9,*)'PERHAPS A DIFFERENT STARTING VALUE SHOULD BE USED.'
        WRITE(9,*)' '
        WRITE(9,*)'IT IS POSSIBLE THAT THERE IS NO LOCAL MINIMUM.'
      ENDIF
      RETURN
      END














      PROGRAM GRADMETH
C     NUMERICAL METHODS: FORTRAN Programs, (c) John H. Mathews 1994
C     To accompany the text:
C     NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
C     Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
C     This free software is complements of the author.
C
C     Algorithm 8.4 (Steepest Descent  or Gradient Method).
C     Section 8.1, Minimization of a Function, Page 418
C
      PARAMETER(MaxN=6,Jmax=200,Max=50,Delta=1E-5,Epsilon=1E-7)
      INTEGER Cond,Count,N
      REAL Error,H,P0,P1,P2,Pmin,S,Y0,Ymin
      CHARACTER ANS*1
      DIMENSION P0(1:MaxN),P1(1:MaxN),P2(1:MaxN),Pmin(1:MaxN),S(1:MaxN)
      EXTERNAL F
10    CALL INPUTS(P0,Y0,N)
      CALL GRADSR(F,Y0,P0,P1,P2,Pmin,S,N,Jmax,Max,
     +Delta,Epsilon,H,Ymin,Error,Cond,Count)
      CALL RESULT(P0,Pmin,N,H,Ymin,Error,Epsilon,Cond,Count)
      WRITE(9,*)' '
      WRITE(9,*)'WANT TO TRY ANOTHER STARTING POINT ? '
      READ(9,'(A)') ANS
      IF (ANS.EQ.'Y' .OR. ANS.EQ.'y') GOTO 10
      END

      REAL FUNCTION F(P,N)
      PARAMETER(MaxN=6)
      INTEGER N
      REAL P,X,Y,Z,U,V,W
      DIMENSION P(1:MaxN)
      X=P(1)
      Y=P(2)
      Z=P(3)
      U=P(4)
      V=P(5)
      W=P(6)
      F=X*X-4*X+Y*Y-Y-X*Y
      RETURN
      END

      SUBROUTINE PRINTFUN(N)
      INTEGER K,N
      CHARACTER C*2,FUN*80,L*6,R*4
      DIMENSION C(1:6)
        N = 2
      FUN ='X*X-4*X+Y*Y-Y-X*Y'
      C(1)='(X'
      C(2)=',Y'
      C(3)=',Z'
      C(4)=',U'
      C(5)=',V'
      C(6)=',W,'
        L ='     F'
        R =') = '
      WRITE(9,999) L,C(1),(C(K),K=2,N),R,FUN
999   FORMAT(9A)
      RETURN
      END

      SUBROUTINE GRADVECT(P,S,N)
      PARAMETER(MaxN=6)
      INTEGER K,N
      REAL G,Length,P,S,Sum,X,Y,Z,U,V,W
      DIMENSION G(1:MaxN),P(1:MaxN),S(1:MaxN)
        X=P(1)
        Y=P(2)
        Z=P(3)
        U=P(4)
        V=P(5)
        W=P(6)
        G(1)=2*X-4-Y
        G(2)=2*Y-1-X
        Sum=0
        DO K=1,N
          Sum=Sum+G(K)*G(K)
        ENDDO
        Length=SQRT(Sum)
        IF (Length.EQ.0) THEN
          Length=1
        ENDIF
        DO K=1,N
          G(K)=-G(K)
          S(K)=G(K)/Length
        ENDDO
        RETURN
      END

      SUBROUTINE QMIN(F,P0,Y0,N,Jmax,Delta,Epsilon,
     &                P1,P2,Pmin,S,H,Ymin,Error,Cond,Count)
      PARAMETER(MaxN=6)
      INTEGER Cond,Count,J,Jmax,K,N
      REAL Delta,Epsilon,Error,H,P0,P1,P2,Pmin,S,Y0,Ymin
      REAL D,E0,E1,E2,H0,H1,H2,Hmin,Length,Y1,Y2
      DIMENSION P0(1:MaxN),P1(1:MaxN),P2(1:MaxN),Pmin(1:MaxN),S(1:MaxN)
      EXTERNAL F
      PARAMETER (Big=1E8)
      DO K=1,N
        P1(K)=P0(K)+H*S(K)
        P2(K)=P0(K)+2*H*S(K)
      ENDDO
      Y1=F(P1,N)
      Y2=F(P2,N)
      Cond=0
      J=0
      WHILE (J.LT.Jmax .AND. Cond.EQ.0)
        Sum=0
        DO K=1,N
          Sum=Sum+P0(K)*P0(K)
        ENDDO
        Length=SQRT(Sum)
        IF (Y0.LE.Y1) THEN
          DO K=1,N
            P2(K)=P1(K)
          ENDDO
          Y2=Y1
          H=H/2
          DO K=1,N
            P1(K)=P0(K)+H*S(K)
          ENDDO
          Y1=F(P1,N)
        ELSE
          IF (Y2.LE.Y1) THEN
            DO K=1,N
              P1(K)=P2(K)
            ENDDO
            Y1=Y2
            H=2*H
            DO K=1,N
              P2(K)=P0(K)+2*H*S(K)
            ENDDO
            Y2=F(P2,N)
          ELSE
              Cond=-1
          ENDIF
        ENDIF
        J=J+1
        IF (H.LT.Delta) Cond=1
        IF (ABS(H).GT.Big .OR. Length.GT.Big) Cond=5
      REPEAT
      Count=Count+J
      IF (Cond.EQ.5) THEN
        DO K=1,N
          Pmin(K)=P1(K)
        ENDDO
        Ymin=Y1
      ELSE
        D=4*Y1-2*Y0-2*Y2             !Start of a long block:
        IF (D.LT.0) THEN
          Hmin=H*(4*Y1-3*Y0-Y2)/D
        ELSE
          Cond=4
          Hmin=H/3
        ENDIF
        DO K=1,N
          Pmin(K)=P0(K)+Hmin*S(K)
        ENDDO
        Ymin=F(Pmin,N)
        H0=ABS(Hmin)
        H1=ABS(Hmin-H)
        H2=ABS(Hmin-2*H)
        IF (H0.LT.H) H=H0
        IF (H1.LT.H) H=H1
        IF (H2.LT.H) H=H2
        IF (H.EQ.0)  H=Hmin
        IF (H.LT.Delta) Cond=1
        E0=ABS(Y0-Ymin)
        E1=ABS(Y1-Ymin)
        E2=ABS(Y2-Ymin)
        IF (E0.NE.0 .AND. E0.LT.Error) Error=E0
        IF (E1.NE.0 .AND. E1.LT.Error) Error=E1
        IF (E2.NE.0 .AND. E2.LT.Error) Error=E2
        IF (E0.EQ.0 .AND. E1.EQ.0 .AND. E2.EQ.0) Error=0
        IF (Error.LT.Epsilon) Cond=2
        IF (Cond.EQ.2 .AND. H.LT.Delta) Cond=3
      ENDIF      !End of the long block.
      RETURN
      END

      SUBROUTINE GRADSR(F,Y0,P0,P1,P2,Pmin,S,N,Jmax,Max,
     &                  Delta,Epsilon,H,Ymin,Error,Cond,Count)
      PARAMETER(MaxN=6)
      INTEGER Cond,Count,Jmax,K,Max,N
      REAL Delta,Epsilon,Error,H,Length,P0,P1,P2,Pmin,S,Y0,Ymin
      DIMENSION P0(1:MaxN),P1(1:MaxN),P2(1:MaxN),Pmin(1:MaxN),S(1:MaxN)
      EXTERNAL F
      H=1
      Sum=0
      DO K=1,N
        Sum=Sum+P0(K)*P0(K)
      ENDDO
      Length=SQRT(Sum)
      IF (Length.GT.1E4) THEN
        H=Length/1E4
      ENDIF
      Error=1
      Count=0
      Cond=0
      WHILE ( Count.LT.Max .AND. Cond.NE.5 .AND.
     &   (H.GT.Delta .OR. Error.GT.Epsilon) )
        CALL GRADVECT(P0,S,N)
        CALL QMIN(F,P0,Y0,N,Jmax,Delta,Epsilon,
     &            P1,P2,Pmin,S,H,Ymin,Error,Cond,Count)
        DO K=1,N
          P0(K)=Pmin(K)
        ENDDO
        Y0=Ymin
        Count=Count+1
      REPEAT
      RETURN
      END

      SUBROUTINE XGRADVECT(P,S,N)
C     This subroutine uses labeled DO loop(s).
      PARAMETER(MaxN=6)
      INTEGER K,N
      REAL G,Length,P,S,Sum,X,Y,Z,U,V,W
      DIMENSION G(1:MaxN),P(1:MaxN),S(1:MaxN)
        X=P(1)
        Y=P(2)
        Z=P(3)
        U=P(4)
        V=P(5)
        W=P(6)
        G(1)=2*X-4-Y
        G(2)=2*Y-1-X
        Sum=0
        DO 10 K=1,N
          Sum=Sum+G(K)*G(K)
10      CONTINUE
        Length=SQRT(Sum)
        IF (Length.EQ.0) THEN
          Length=1
        ENDIF
        DO 20 K=1,N
          G(K)=-G(K)
          S(K)=G(K)/Length
20      CONTINUE
        RETURN
      END

      SUBROUTINE XQMIN(F,P0,Y0,N,Jmax,Delta,Epsilon,
     +P1,P2,Pmin,S,H,Ymin,Error,Cond,Count)
C     This subroutine uses labeled DO loop(s).
      PARAMETER(MaxN=6)
      INTEGER Cond,Count,J,Jmax,K,N
      REAL Delta,Epsilon,Error,H,P0,P1,P2,Pmin,S,Y0,Ymin
      REAL D,E0,E1,E2,H0,H1,H2,Hmin,Length,Y1,Y2
      DIMENSION P0(1:MaxN),P1(1:MaxN),P2(1:MaxN),Pmin(1:MaxN),S(1:MaxN)
      EXTERNAL F
      PARAMETER (Big=1E8)
      DO 10 K=1,N
        P1(K)=P0(K)+H*S(K)
        P2(K)=P0(K)+2*H*S(K)
10    CONTINUE
      Y1=F(P1,N)
      Y2=F(P2,N)
      Cond=0
      J=0
20    IF (J.LT.Jmax .AND. Cond.EQ.0) THEN
        Sum=0
        DO 30 K=1,N
          Sum=Sum+P0(K)*P0(K)
30      CONTINUE
        Length=SQRT(Sum)
        IF (Y0.LE.Y1) THEN
          DO 40 K=1,N
            P2(K)=P1(K)
40        CONTINUE
          Y2=Y1
          H=H/2
          DO 50 K=1,N
            P1(K)=P0(K)+H*S(K)
50        CONTINUE
          Y1=F(P1,N)
        ELSE
          IF (Y2.LE.Y1) THEN
            DO 60 K=1,N
              P1(K)=P2(K)
60          CONTINUE
            Y1=Y2
            H=2*H
            DO 70 K=1,N
              P2(K)=P0(K)+2*H*S(K)
70          CONTINUE
            Y2=F(P2,N)
          ELSE
              Cond=-1
          ENDIF
        ENDIF
        J=J+1
        IF (H.LT.Delta) Cond=1
        IF (ABS(H).GT.Big .OR. Length.GT.Big) Cond=5
        GOTO 20
      ENDIF
      Count=Count+J
      IF (Cond.EQ.5) THEN
        DO 80 K=1,N
          Pmin(K)=P1(K)
80      CONTINUE
        Ymin=Y1
      ELSE
        D=4*Y1-2*Y0-2*Y2             !Start of a long block:
        IF (D.LT.0) THEN
          Hmin=H*(4*Y1-3*Y0-Y2)/D
        ELSE
          Cond=4
          Hmin=H/3
        ENDIF
        DO 90 K=1,N
          Pmin(K)=P0(K)+Hmin*S(K)
90      CONTINUE
        Ymin=F(Pmin,N)
        H0=ABS(Hmin)
        H1=ABS(Hmin-H)
        H2=ABS(Hmin-2*H)
        IF (H0.LT.H) H=H0
        IF (H1.LT.H) H=H1
        IF (H2.LT.H) H=H2
        IF (H.EQ.0)  H=Hmin
        IF (H.LT.Delta) Cond=1
        E0=ABS(Y0-Ymin)
        E1=ABS(Y1-Ymin)
        E2=ABS(Y2-Ymin)
        IF (E0.NE.0 .AND. E0.LT.Error) Error=E0
        IF (E1.NE.0 .AND. E1.LT.Error) Error=E1
        IF (E2.NE.0 .AND. E2.LT.Error) Error=E2
        IF (E0.EQ.0 .AND. E1.EQ.0 .AND. E2.EQ.0) Error=0
        IF (Error.LT.Epsilon) Cond=2
        IF (Cond.EQ.2 .AND. H.LT.Delta) Cond=3
      ENDIF      !End of the long block.
      RETURN
      END

      SUBROUTINE XGRADSR(F,Y0,P0,P1,P2,Pmin,S,N,Jmax,Max,
     +Delta,Epsilon,H,Ymin,Error,Cond,Count)
C     This subroutine uses labeled DO loop(s).
      PARAMETER(MaxN=6)
      INTEGER Cond,Count,Jmax,K,Max,N
      REAL Delta,Epsilon,Error,H,Length,P0,P1,P2,Pmin,S,Y0,Ymin
      DIMENSION P0(1:MaxN),P1(1:MaxN),P2(1:MaxN),Pmin(1:MaxN),S(1:MaxN)
      EXTERNAL F
      H=1
      Sum=0
      DO 10 K=1,N
        Sum=Sum+P0(K)*P0(K)
10    CONTINUE
      Length=SQRT(Sum)
      IF (Length.GT.1E4) THEN
        H=Length/1E4
      ENDIF
      Error=1
      Count=0
      Cond=0
20    IF ( Count.LT.Max .AND. Cond.NE.5 .AND.
     &   (H.GT.Delta .OR. Error.GT.Epsilon) ) THEN
        CALL GRADVECT(P0,S,N)
        CALL QMIN(F,P0,Y0,N,Jmax,Delta,Epsilon,
     &            P1,P2,Pmin,S,H,Ymin,Error,Cond,Count)
        DO 30 K=1,N
          P0(K)=Pmin(K)
30      CONTINUE
        Y0=Ymin
        Count=Count+1
        GOTO 20
      ENDIF
      RETURN
      END

      SUBROUTINE INPUTS(P0,Y0,N)
      PARAMETER(MaxN=6)
      INTEGER I,J,K,N
      REAL P0,Y0,Z
      CHARACTER RESP*40
      DIMENSION P0(1:MaxN),Z(1:MaxN)
      EXTERNAL F
      DO 10 I=1,18
        WRITE(9,*)' '
10    CONTINUE
      WRITE(9,*)'   THE GRADIENT METHOD OR "METHOD OF STEEPEST DESCENT" I
     +S USED FOR FINDING'
      WRITE(9,*)' '
      WRITE(9,*)'THE MINIMUM OF THE FUNCTION  F(P)  WHERE  P = (p ,p ,..
     +.,p )  AND  N<=6.'
      WRITE(9,*)'                                                1  2
     +   N'
      WRITE(9,*)' '
      WRITE(9,*)'FOR CONVENIENCE, THE FUNCTION F(P) CAN BE EXPRESSED USI
     +NG THE VARIABLES'
      WRITE(9,*)' '
      WRITE(9,*)'X = p , Y = p , Z = p , U = p , V = p , W = p .'
      WRITE(9,*)'     1       2       3       4       5       6 '
      WRITE(9,*)' '
      WRITE(9,*)'    THE FUNCTION TO BE MINIMIZED IS:'
      WRITE(9,*)' '
      CALL PRINTFUN(N)
      WRITE(9,*)' '
      IF (N.EQ.2) THEN
        WRITE(9,*)'GIVE THE STARTING POINT P = (p ,p )'
        WRITE(9,*)'                              1  2'
      ELSEIF (N.EQ.3) THEN
        WRITE(9,*)'GIVE THE STARTING POINT P = (p ,p ,p )'
        WRITE(9,*)'                              1  2  3'
      ELSE
        WRITE(9,*)'GIVE THE STARTING POINT P = (p ,p ,...,p )'
        WRITE(9,*)'                              1  2      N'
        WRITE(9,*)'WHERE N=',N
      ENDIF
      WRITE(9,*)' '
      DO 20 K=1,N
        WRITE(9,*)'P(',K,') = '
        READ(9,*) P0(K)
20    CONTINUE
      Y0=F(P0,N)
      RETURN
      END

      SUBROUTINE RESULT(P0,Pmin,N,H,Ymin,Error,Epsilon,Cond,Count)
      PARAMETER(MaxN=6)
      INTEGER Cond,Count,I,K,N
      REAL Epsilon,Error,H,P0,Pmin,Ymin
      DIMENSION P0(1:MaxN),Pmin(1:MaxN)
      EXTERNAL F
      DO 10 I=1,16
        WRITE(9,*)' '
10    CONTINUE
      WRITE(9,*)'   THE GRADIENT SEARCH METHOD WAS USED TO FIND THE MINI
     +MUM OF'
      WRITE(9,*)' '
      CALL PRINTFUN(N)
      WRITE(9,*) Count,' ITERATIONS WERE REQUIRED'
      WRITE(9,*)' '
      DO 20 K=1,N
        WRITE(9,*)'P(',K,')  = ',Pmin(K)
20    CONTINUE
      WRITE(9,*)' '
      WRITE(9,*)'  DP  = ',H   ,'  IS AN ESTIMATE FOR THE ACCURACY.'
      WRITE(9,*)' '
      WRITE(9,*)'THE MINIMUM VALUE OF THE FUNCTION IS'
      WRITE(9,*)' '
      WRITE(9,*)'F(P)  = ',Ymin
      WRITE(9,*)' '
      WRITE(9,*)'  DF  = ',Error,'  IS AN ESTIMATE FOR THE ACCURACY.'
      WRITE(9,*)' '
      IF (Cond.EQ.0) THEN
        WRITE(9,*)'CONVERGENCE IS DOUBTFUL BECAUSE THE'
        WRITE(9,*)'MAXIMUM NUMBER OF ITERATIONS WAS EXCEEDED.'
      ELSEIF (Cond.EQ.1) THEN
        WRITE(9,*)'CONVERGENCE HAS BEEN ACHIEVED BECAUSE '
        WRITE(9,*)'CONSECUTIVE POINTS  ARE  CLOSER  THAN ',H
      ELSEIF (Cond.EQ.2) THEN
        WRITE(9,*)'CONVERGENCE HAS BEEN ACHIEVED BECAUSE '
        WRITE(9,*)'CONSECUTIVE FUNCTION VALUES DIFFER BY ',Error
      ELSEIF (Cond.EQ.3) THEN
        WRITE(9,*)'CONVERGENCE HAS BEEN ACHIEVED BECAUSE '
        WRITE(9,*)'CONSECUTIVE POINTS  ARE  CLOSER  THAN ',H
        WRITE(9,*)'CONSECUTIVE FUNCTION VALUES DIFFER BY ',Error
      ELSEIF (Cond.EQ.4) THEN
        WRITE(9,*)'CONVERGENCE IS DOUBTFUL BECAUSE DIVISION BY ZERO WAS
     +ENCOUNTERED.'
        IF (Error.LT.Epsilon/100.0) THEN
          WRITE(9,*)'HOWEVER, CONSECUTIVE POINTS ARE CLOSE.'
        ENDIF
      ELSEIF (Cond.EQ.5) THEN
        WRITE(9,*)'CONVERGENCE IS DOUBTFUL, H IS TOO LARGE, H = ',H
        WRITE(9,*)'IT IS POSSIBLE THAT THERE IS NO LOCAL MINIMUM.'
      ENDIF
      RETURN
      END
  
