C-------------------------------------------------------------------------
C   PURPOSE             - MATRIX DECOMPOSITION, LINEAR EQUATION
C                           SOLUTION - SPACE ECONOMIZER SOLUTION -
C                           COMPLEX MATRICES
C
C   USAGE               - CALL AXEB (A,N,IA,B,M,IB,IJOB,WA,IER)
C
C   ARGUMENTS    A      - INPUT REAL N BY N MATRIX CONTAINING THE
C                           REAL COEFFICIENTS OF THE EQUATION AX = B
C                         ON OUTPUT, A CONTAINS THE L-U DECOMPOSITION O
C                           A ROWWISE PERMUTATION OF THE INPUT MATRIX A
C                N      - ORDER OF MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT REAL N BY M MATRIX CONTAINING THE
C                           M REAL VALUED RIGHT HAND SIDES OF THE
C                           EQUATION AX = B.
C                         ON OUTPUT, THE SOLUTION MATRIX X REPLACES B.
C                           IF IJOB=1, B IS NOT USED.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS SPECIFIE
C                           IN THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                IJOB   - INPUT OPTION PARAMETER.  IJOB=I IMPLIES WHEN
C                           I=0, FACTOR THE MATRIX AND SOLVE THE
C                             EQUATION AX=B.
C                           I=1, FACTOR THE MATRIX A.
C                           I=2, SOLVE THE EQUATION AX=B.  THIS
C                             OPTION IMPLIES THAT LEQT1C HAS ALREADY
C                             BEEN CALLED USING IJOB=0 OR 1 SO THAT
C                             THE MATRIX HAS ALREADY BEEN FACTORED.  IN
C                             THIS CASE OUTPUT MATRIX A MUST HAVE BEEN
C                             SAVED FOR REUSE IN THE CALL TO LEQT1C.
C                WA     - WORK AREA OF LENGTH N CONTAINING THE PIVOT
C                           INDICES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR.  (SEE THE
C                             CHAPTER L PRELUDE.)
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  WHEN IJOB=1, ARGUMENTS B, M AND IB ARE NOT USED BY
C                LEQT1C.
C            2.  INPUT MATRIX A IS DESTROYED WHEN IJOB=0 OR 1. WHEN
C                IJOB=0 OR 2, B IS REPLACED WITH THE SOLUTION X.
C            3.  LEQT1C CAN BE USED TO COMPUTE THE INVERSE OF A REAL
C                MATRIX. THIS IS DONE BY CALLING LEQT1C WITH M=N,
C                B=THE N BY N IDENTITY MATRIX AND IJOB=0. WHEN N IS
C                LARGE, IT MAY BE MORE PRACTICAL TO COMPUTE THE INVERSE
C                A COLUMN AT A TIME. TO DO THIS, FIRST CALL LEQT1C WITH
C                IJOB=1 TO FACTOR A. MAKE SUCCEEDING CALLS WITH M=1, B
C                =A COLUMN OF THE IDENTITY MATRIX AND IJOB=2. B WILL BE
C                REPLACED BY THE CORRESPONDING COLUMN OF A INVERSE.
C            4.  THE DETERMINANT OF A CAN BE COMPUTED AFTER LEQT1C HAS
C                BEEN CALLED AS FOLLOWS
C
C                  DET = 1.D0
C                  DO 5 I = 1,N
C                     IPVT = WA(I)
C                     IF (IPVT .NE. I) DET = -DET
C                     DET = DET*A(I,I)
C                5 CONTINUE
C
C
C
C----------------------------------------------------------------------
C
cccccccccccccccccccccccccccccccccccccccccc
c     TO BE REPLACED BY LINPACK ROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE AXEB (A,N,IA,B,M,IB,IJOB,WA,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION   A(IA,N),B(IB,M),WA(N),T(2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      EQUIVALENCE        (SUM1,T(1))
      DATA               ZERO/0.0D0/,ONE/1.D0/
C                                  INITIALIZATION
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IJOB .EQ. 2) GO TO 75
      RN = N
C                                  FIND EQUILIBRATION FACTORS
      DO 10 I=1,N
         BIG = ZERO
         DO 5 J=1,N
            TEMP = A(I,J)
            P = DABS(TEMP)
            IF (P .GT. BIG) BIG = P
    5    CONTINUE
         IF (BIG .EQ. ZERO) GO TO 105
         WA(I) = ONE/BIG
   10 CONTINUE
C                                  L-U DECOMPOSITION
      DO 70 J = 1,N
         JM1 = J-1
         IF (JM1 .LT. 1) GO TO 25
C                                  COMPUTE U(I,J), I=1,...,J-1
         DO 20 I=1,JM1
            SUM1 = A(I,J)
            IM1 = I-1
            IF (IM1 .LT. 1) GO TO 20
            DO 15 K=1,IM1
               SUM1 = SUM1-A(I,K)*A(K,J)
   15       CONTINUE
            A(I,J) = SUM1
   20    CONTINUE
   25    P = ZERO
C                                  COMPUTE U(J,J) AND L(I,J), I=J+1,...
         DO 45 I=J,N
            SUM1 = A(I,J)
            IF (JM1 .LT. 1) GO TO 40
            DO 35 K=1,JM1
               SUM1 = SUM1-A(I,K)*A(K,J)
   35       CONTINUE
            A(I,J) = SUM1
   40       Q = WA(I)*DABS(SUM1)
            IF (P .GE. Q) GO TO 45
            P = Q
            IMAX = I
   45    CONTINUE
C                                  TEST FOR ALGORITHMIC SINGULARITY
         Q = RN+P
         IF (Q .EQ. RN) GO TO 105
         IF (J .EQ. IMAX) GO TO 60
C                                  INTERCHANGE ROWS J AND IMAX
         DO 50 K=1,N
            TEMP = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = TEMP
   50    CONTINUE
         WA(IMAX) = WA(J)
   60    WA(J) = IMAX
         JP1 = J+1
         IF (JP1 .GT. N) GO TO 70
C                                  DIVIDE BY PIVOT ELEMENT U(J,J)
         TEMP = A(J,J)
         DO 65 I = JP1,N
            A(I,J) = A(I,J)/TEMP
   65    CONTINUE
   70 CONTINUE
   75 IF (IJOB .EQ. 1) GO TO 9005
      DO 103 K = 1,M
C                                  SOLVE UX = Y FOR X
         IW = 0
         DO 90 I = 1,N
            IMAX = WA(I)
            SUM1 = B(IMAX,K)
            B(IMAX,K) = B(I,K)
            IF (IW .EQ. 0) GO TO 85
            IM1 = I-1
            DO 80 J = IW,IM1
               SUM1 = SUM1-A(I,J)*B(J,K)
   80       CONTINUE
            GO TO 88
   85       IF (T(1) .NE. ZERO .OR. T(2) .NE. ZERO) IW = I
   88       B(I,K) = SUM1
   90    CONTINUE
C                                  SOLVE LY = B FOR Y
         N1 = N+1
         DO 100 IW = 1,N
            I = N1-IW
            JP1 = I+1
            SUM1 = B(I,K)
            IF (JP1 .GT. N) GO TO 98
            DO 95 J = JP1,N
               SUM1 = SUM1-A(I,J)*B(J,K)
   95       CONTINUE
   98       B(I,K) = SUM1/A(I,I)
  100    CONTINUE
  103 CONTINUE
      GO TO 9005
C                                  ALGORITHMIC SINGULARITY
  105 IER = 129
 9000 CONTINUE
C                                  PRINT ERROR
      PRINT*,'SOME ERROR HAS OCCURED'
      STOP
C      CALL UERTST(IER,6HLEQT1C)
 9005 RETURN
      END
