      SUBROUTINE zmfINIT(U, UPRIME, X, JJ,NATJ,NEQ)
      use var
      IMPLICIT NONE  
C
C*****precision > double
c     IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

C*****END precision > double
      DOUBLE PRECISION ::  U(NATJ*JJ),UPRIME(NATJ*JJ),X(JJ),XJ
      INTEGER :: K,J,JJ,NEQ,NATJ,I,M

C       JJ    - Number of Nodes ( In CVCM reactor JJ =1 )
C       NATJ  - Number of Dependent Variables
C       U     - Dependent Variable Matrix
C             dimension(NATJ*JJ)       
C       UPRIME(*) -- Set this array to the initial values of the NATJ*JJ first
C               derivatives of the solution components at the initial
C               point.  You must dimension YPRIME at least NATJ*JJ in your
C               calling program.   
C       U(NZ) - Mixture Fraction 
      M=JJ-2
        DO 10 J = 0,M+1
          XJ = J*(X(2)-X(1))
          I = J + 1
          U(I) = 16.0D0*XJ*(1.0D0-XJ)
 10       CONTINUE

      DO I=1,NEQ
      UPRIME(I)=0.0
      END DO
      END SUBROUTINE zmfinit
