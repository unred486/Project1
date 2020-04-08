      SUBROUTINE CVCMINIT(ICKWRK,RCKWRK,U, UPRIME, WDOT,RPAR,KSYM,h,
     1               cp,WT,KK,JJ,NATJ,LOUT)
      use var
      IMPLICIT NONE  
C
C*****precision > double
c     IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

C*****END precision > double
      DOUBLE PRECISION ::  U(NATJ*JJ),WT(*),RCKWRK(*),
     1 UPRIME(NATJ*JJ),RPAR(NATJ),WDOT(KK),h(KK),cp(KK),Xi(KK)

      INTEGER :: ICKWRK(*)
      INTEGER :: K,J,JJ,LOUT,KK,NATJ
      DOUBLE PRECISION :: RU,sumwdot,sumhwdot,sumdenom,sumX,dummy,
     1                      WTM,rho,P,cvm
      character(len=32) :: my_fmt,my_fmt2
      CHARACTER*(*) :: KSYM(*)
C       JJ    - Number of Nodes ( In CVCM reactor JJ =1 )
C       NATJ  - Number of Dependent Variables
C       U     - Dependent Variable Matrix
C             dimension(NATJ*JJ)       
C       UPRIME(*) -- Set this array to the initial values of the NATJ*JJ first
C               derivatives of the solution components at the initial
C               point.  You must dimension YPRIME at least NATJ*JJ in your
C               calling program.   
C       U(CM) - Pressure 
C               cgs units - dynes/cm**2
C       U(NT) - Chamber Temperature
C               unit - K
C       U(NY) - Mass fraction
C       RU     - Universal gas constant.
C                   cgs units - 8.314E7 ergs/(mole*K)
C                   Data type - real scalar
C       CP   - Specific heats at constant pressure in mass units
C              for the species.
C                   cgs units - ergs/(mole*K)
C                   Data type - real array
C                   Dimension CP(*) at least KK, the total number of
C                   species.
C        h    - Enthalpies in molar units for the species.
C                   cgs units - ergs/mole
C                   Data type - real array
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
c   ickwrk - integer chemkin work space.
c              dimensioning - see chemkin documentation.
c   rckwrk  - floating point chemkin work space.
c              dimensioning - see chemkin documentation.

C  SUBROUTINE CKRP 
C     Returns universal gas constants and the pressure of one standard
C     atmosphere            
      call CKRP   (ICKWRK, RCKWRK, RU, dummy, dummy)
      
      call CKMMWY (U(NY), ICKWRK, RCKWRK, WTM)
      
      P=U(NM)*RU/WTM*U(NT)
C  SUBROUTINE CKYTCR
C     Returns the molar concentrations given the mass density,
C     temperature and mass fractions
      CALL CKYTCR (P,U(NT), U(NY), ICKWRK, RCKWRK, Xi)  
C  SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
C     Returns the molecular weights of the species
      CALL CKWT   (ICKWRK, RCKWRK, WT)   


C  SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C     Returns the molar production rates of the species given the
C     pressure, temperature and mass fractions 
      call ckwyp (P, U(NT), U(NY), ickwrk, rckwrk, wdot)
      


C  SUBROUTINE CKHML
C     Returns the enthalpies in molar units
      call CKHML  (U(NT), ICKWRK, RCKWRK, h)
      
C  SUBROUTINE CKCPML (T, ICKWRK, RCKWRK, CPML)
C     Returns the specific heats at constant pressure in molar units
      call CKCVBS (U(NT),U(NY), ICKWRK, RCKWRK, cvm)
      
C       see pg 188 of an introduction to combustion by Turns 
C       for reference of UPRIME(NT) and U(NM) expression 

C       Initializing UPRIME of kth species
      DO K=1, KK
      wdot(K)=wdot(K)
      UPRIME(NYS+K)=wdot(K)*WT(K)/U(NM)

      END DO
      
      sumwdot=0.0
      sumhwdot=0.0
      sumdenom=0.0
      sumX=0.0

      do K=1,KK
      sumwdot= sumwdot + wdot(K)
      sumhwdot= sumhwdot + h(K)*wdot(K)
c      sumdenom= sumdenom + Xi(K)*(cp(K)-RU)
      sumX= sumX + Xi(K)
      end do
C       Initializing temperature UPRIME
      UPRIME(NT)= (RU*U(NT)*sumwdot-sumhwdot)/(u(nm)*cvm)

C       Initializing pressure UPRIME
      UPRIME(NM)=0! RU*U(NT)*sumwdot+RU*UPRIME(NT)*sumX

C////////   checking if C matrix is being imported correctly //////////     
      WRITE(LOUT,*) 'Initialization at t = 0s'
      WRITE(my_fmt , '(I0)')  KK+3
      my_fmt = '('//trim(my_fmt)//'(1PE12.4))' 
      WRITE(my_fmt2 , '(I0)')  KK+2
      my_fmt2 = '('//trim(my_fmt2)//'(1PE12.4))' 
      
      WRITE(LOUT,73) (trim(KSYM(K)),K=1,KK)
      WRITE(LOUT,my_fmt) (U(NYS+K),K=1,KK), U(NT),U(NM),P
      WRITE(LOUT,74) (trim(KSYM(K)),K=1,KK)
      WRITE(LOUT,my_fmt2) (UPRIME(NYS+K),K=1,KK),UPRIME(NT),UPRIME(NM)
      write(LOUT,*)
      
73    FORMAT(<KK> (' X'A''), ' Temperature Pressure') 
     
74    FORMAT( <KK> (' dX'A'/dt'), ' dT/dt dP/dt') 
c//////////////////////////////////////////////////////////////////////
      
c      WRITE(LOUT,*) UPRIME(NT),(UPRIME(NYS+K),K=1,KK),UPRIME(NM)
      END SUBROUTINE CVCMinit
