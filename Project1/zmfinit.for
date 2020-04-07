      SUBROUTINE zmfINIT(U, UPRIME, X, JJ,NATJ,NEQ)
      use var
      IMPLICIT NONE  
C
C*****precision > double
c     IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

C*****END precision > double
      DOUBLE PRECISION ::  U(NATJ*JJ),UPRIME(NATJ*JJ),X(JJ)
      INTEGER :: K,J,JJ,NEQ,NATJ,I

C       JJ    - Number of Nodes ( In CVCM reactor JJ =1 )
C       NATJ  - Number of Dependent Variables
C       U     - Dependent Variable Matrix
C             dimension(NATJ*JJ)       
C       UPRIME(*) -- Set this array to the initial values of the NATJ*JJ first
C               derivatives of the solution components at the initial
C               point.  You must dimension YPRIME at least NATJ*JJ in your
C               calling program.   
C       U(NZ) - Mixture Fraction 


      DO I=1,NEQ
      UPRIME(I)=0.0
      END DO
      END SUBROUTINE zmfinit
