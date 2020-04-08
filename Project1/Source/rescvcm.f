C////////////////// RESCVCM                    ////////////////////////
C//////////////////////////////////////////////////////////////////////
C       Sets Up Boundary Condition and interior residuals
C       CVCM reactor only requires initial conditions
C       Therefore, simply calls FCVCM which calculates interior residuals
      SUBROUTINE RESCVCM(T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR)
      use var
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C     need for UPRIME calculation      
      DIMENSION :: U(*),UPRIME(*),DELTA(*),RPAR(*),IPAR(*)
      COMMON /SP/ KK,NATJ
      
      Call FCVCM(U,DELTA, RPAR)
      DO J=1,NATJ
      DELTA(J)=UPRIME(J)-DELTA(J)
      END DO
      END SUBROUTINE  RESCVCM 
      

C////////////////// FCVCM                      ////////////////////////
C//////////////////////////////////////////////////////////////////////
C       Residual consists of algebraic source term( ex : reaction term)
C       and differential term ( ex : tranport term)
C       DDASPK can use preconditions which differentiates between 
C       differential and algebraic terms to reduce the memory usage
C       FCVCM calculates differential residuals and calls WCVCM to 
C       calculate the reaction source terms in RPAR array
      SUBROUTINE FCVCM(U,CRATE,RPAR)     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION :: U(*),CRATE(*),RPAR(*)
      COMMON /SP/ KK,NATJ
      
      CALL WCVCM(T,U,RPAR)
      DO J=1,NATJ
      CRATE(J) = RPAR(J)
      END DO
      
      END SUBROUTINE FCVCM
      
      
C//////////////////////////////////////////////////////////////////////
C////////////////// WCVCM                      ////////////////////////
C//////////////////////////////////////////////////////////////////////      
C       WCVCM
C       calculate the reaction source terms in CRATE array
C       CRATE = object array ( like actual crate box)

      SUBROUTINE WCVCM(T,U,CRATE)
      use var 
      use varray
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DOUBLE PRECISION :: WTM
      DIMENSION :: U(*), CRATE(*),Xi(KK)
      COMMON /SP/ KK,NATJ
      
C       NEQ - the number of equations in the system.
c       NHBW - Half Band Width
C       JJ    - Number of Nodes ( In CVCM reactor JJ =1 )
C       KK    - Total Number of Species
C       NATJ  - Number of Dependent Variables
C       U     - Dependent Variable Matrix
C             dimension(NATJ*JJ)    
C       U(NM) - mixture density 
C               cgs units - gm/cm**3
C       U(NT) - Chamber Temperature
C               unit - K
C       U(NY) - Mass fraction
C       UPRIME(*) -- Set this array to the initial values of the NATJ*JJ first
C               derivatives of the solution components at the initial
C               point.  You must dimension YPRIME at least NATJ*JJ in your
C               calling program.   
C       P     - Pressure 
C               unit - dynes/cm**2
C       WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C        h    - Enthalpies in molar units for the species.
C                   cgs units - ergs/mole
C                   Data type - real array
C       CVM   - Specific heats at constant pressure in mass units
C              for the micture.
C                   cgs units - ergs/(gm*K)
C                   Data type - real array
C       WT     - Molecular weights of the species.
C                   cgs units - gm/mole  
C     WTM    - Mean molecular weight of the species mixture.
C                   cgs units - gm/mole
C       RU     - Universal gas constant.
C                   cgs units - 8.314E7 ergs/(mole*K)
C                   Data type - real scalar
C       rho - density 
C               cgs units - gm/cm^3
C       INFO   - DDASPK control array  
C       RW   - REAL WORK SPACE for DDASPk
C       IW   - INTEGER WORK SPACE for DDASPK
C       RCKWRK - REAL WORK SPACE for Chemkin Library
C       ICKWRK - INTEGER WORK SPACE for Chmekin Library
C       ATOL   - ABSOLUTE CONVERGENCE CRITERIA
C       RTOL   - RELATIVE CONVERGENCE CRITERIA
C       IDID -- reports what the code did, Refer to DDASPK.f for more info
C       T      - Start Time for current simulation  unit : s
C       TOUT   - END TIME for current simulation    unit : s    
C       DT1    - time step  
CC/////////////////////////////////////////////////////////////////////

C  SUBROUTINE CKRP 
C     Returns universal gas constants and the pressure of one standard
C     atmosphere            
      call CKRP   (ICKWRKALL, RCKWRKALL, RU, dummy, dummy)
      
      CALL CKMMWY (U(NY), ICKWRKALL, RCKWRKALL, WTM)      
C  SUBROUTINE CKRHOY 
C     Returns the mass density of the gas mixture given the pressure,
C     temperature and mass fractions
      P=U(NM)*RU/WTM*U(NT)
C      CALL CKRHOY (P, U(NT), U(NY), ICKWRKALL, RCKWRKALL, RHO)
      
C  SUBROUTINE CKYTCR
C     Returns the molar concentrations given the mass density,
C     temperature and mass fractions
c      CALL CKYTCR (U(NM),U(NT), U(NY), ICKWRKall, RCKWRKall, Xi)  
      CALL CKWT   (ICKWRKall, RCKWRKall, WTALL)  
      

C  SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C     Returns the molar production rates of the species given the
C     pressure, temperature and mass fractions 
      CALL ckwyp (P, U(NT), U(NY), ICKWRKALL, RCKWRKALL, WDOTALL)


      
C  SUBROUTINE CKHML
C     Returns the enthalpies in molar units
      CALL CKHML  (U(NT), ICKWRKALL, RCKWRKALL, hall)
      
C
C  SUBROUTINE CKCPML (T, ICKWRK, RCKWRK, CPML)
C     Returns the specific heats at constant pressure in molar units
      call CKCVBS (U(NT), U(NY), ICKWRKALL, RCKWRKALL, cvm)

C       see pg 188 of an introduction to combustion by Turns 
C       for reference of UPRIME(NT) and U(NM) expression 

C       Initializing UPRIME of kth species
      DO K=1, KK
      wdotALL(K)=wdotALL(K)
      CRATE(NYS+K)=wdotALL(K)*WTALL(K)/U(NM)
      END DO
      
c      write(6,10)
c      write(6,11) wdotall(5),u(6),crate(6),NYS
c10    format('dXO2/dt-wdotside XO2 dXO2/dt-crateside NYS')
c11    format(3(4X,E10.3), I2)  
c      write(6,120)(wdotALL(k),K=1,KK)
c120   format(57(E10.3))

      sumwdot=0.0
      sumhwdot=0.0
      sumdenom=0.0
      sumX=0.0
      
      do K=1,KK
      sumwdot= sumwdot + wdotALL(K)
      sumhwdot= sumhwdot + hALL(K)*wdotall(K)
c      sumdenom= sumdenom + Xi(K)*(cpALL(K)-RU)
c      sumX= sumX + Xi(K)
      end do
      
      
C      Temperature Residual Calculation 
      CRATE(NT)=(RU*U(NT)*sumwdot-sumhwdot)/(cvm*u(nm))
C      Pressure Residual Calculation
      CRATE(NM)=0!RU*U(NT)*sumwdot+RU*CRATE(NT)*sumX    
      
        
c           write(6,100) (cpall(K),K=1,KK) 
c100     format(57(4X,E10.3))
c      write(6,140) RU*U(NT)*sumwdot,sumwdot,sumhwdot,sumdenom,
c     1 CRATE(NT)
c140   format(5(4X,E10.3))




      END SUBROUTINE WCVCM
