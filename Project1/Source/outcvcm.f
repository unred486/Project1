C//////////////////////////////////////////////////////////////////////
C////////////////// OUTCVCM ////////////////////////
C//////////////////////////////////////////////////////////////////////  
C
      SUBROUTINE OUTCVCM(T,U,UPRIME,ICKWRK,RCKWRK,LOUT,IDID,IW,RW,
     1 NHBW,KSYM,KK) 
C       OUTCVCM prints out solution variables such as U and UPRIME
     
C       NEQ    - the number of equations in the system.
c       NHBW   - Half Band Width
C       JJ     - Number of Nodes ( In CVCM reactor JJ =1 )
C       KK     - Total Number of Species
C       NATJ   - Number of Dependent Variables
C       U      - Dependent Variable Matrix
C                dimension(NATJ*JJ)   
C       U(NM)  - mixture density 
C                cgs units - gm/cm**3
C       U(NT)  - Chamber Temperature
C                unit - K
C       U(NY)  - Mass fraction 
C    UPRIME(*) -- Set this array to the initial values of the NATJ*JJ first
C               derivatives of the solution components at the initial
C               point.  You must dimension YPRIME at least NATJ*JJ in your
C               calling program.   
C       P      - Pressure 
C               unit - dynes/cm**2
C       WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C        h     - Enthalpies in molar units for the species.
C                   cgs units - ergs/mole
C                   Data type - real array
C       CP     - Specific heats at constant pressure in mass units
C              for the species.
C                   cgs units - ergs/(mole*K)
C                   Data type - real array
C                   Dimension CP(*) at least KK, the total number of
C                   species.
C       WT     - Molecular weights of the species.
C                   cgs units - gm/mole  
C       WTM    - Mean molecular weight of the species mixture.
C                   cgs units - gm/mole
C       RU     - Universal gas constant.
C                   cgs units - 8.314E7 ergs/(mole*K)
C                   Data type - real scalar
C       rho    - density 
C               cgs units - gm/cm^3
C       INFO   - DDASPK control array  
C       RW     - REAL WORK SPACE for DDASPk
C       IW     - INTEGER WORK SPACE for DDASPK
C       RCKWRK - REAL WORK SPACE for Chemkin Library
C       ICKWRK - INTEGER WORK SPACE for Chmekin Library
C       ATOL   - ABSOLUTE CONVERGENCE CRITERIA
C       RTOL   - RELATIVE CONVERGENCE CRITERIA
C       IDID   - reports what the code did, Refer to DDASPK.f for more info
C       T      - Start Time for current simulation  unit : s
C       TOUT   - END TIME for current simulation    unit : s    
C       DT1    - time step  
      use var
      IMPLICIT NONE
      DOUBLE PRECISION :: U(*),UPRIME(*),RCKWRK(*),RW(*)
      INTEGER ::ICKWRK(*),IW(*)
      DOUBLE PRECISION :: T,dummy,WTM,P,RU,AVDIM
      INTEGER :: NHBW,IDID,LOUT, KK,J,JJ,K,MBAND,NST,NPE,PRE,
     1          NNI,NLI,NPS,NCFN,NCFL
      character(len=32) :: my_fmt,my_fmt2
      CHARACTER*(*) :: KSYM(*)
      
      WRITE(LOUT, 80) T
      WRITE(LOUT,*)

      MBAND = NHBW + NHBW + 1      
      NST = IW(11)
      NPE = IW(13)
      NRE = IW(12) + NPE*MBAND
      LIW = IW(17)
      LRW = IW(18)
      NNI = IW(19)
      NLI = IW(20)
      NPS = IW(21)
      IF (NNI .NE. 0) AVDIM = REAL(NLI)/REAL(NNI)
      NCFN = IW(15)
      NCFL = IW(16)
      WRITE (LOUT,90) LRW,LIW,NST,NRE,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL

      WRITE(my_fmt , '(I0)')  KK+3
      my_fmt = '('//trim(my_fmt)//'(1PE12.4))'     
      WRITE(my_fmt2 , '(I0)')  KK+2
      my_fmt2 = '('//trim(my_fmt2)//'(1PE12.4))'   
C       Calculates RU      
      call CKRP   (ICKWRK, RCKWRK, RU, dummy, dummy)     
C       Calculates WTM
      call CKMMWY (U(NY), ICKWRK, RCKWRK, WTM)      
      P=U(NM)*RU/WTM*U(NT)
      
      WRITE(LOUT,73) (trim(KSYM(K)),K=1,KK)
      WRITE(LOUT,my_fmt) (U(NYS+K),K=1,KK), U(NT),U(NM),P
      WRITE(LOUT,74) (trim(KSYM(K)),K=1,KK)
      WRITE(LOUT,my_fmt2)  (UPRIME(NYS+K),K=1,KK), UPRIME(NT),UPRIME(NM)
c      CALL CKRHOY (U(NM), U(NT), U(NY), ICKWRK, RCKWRK, RHO)
      write(LOUT,*) IDID
      write(LOUT,75)
73    FORMAT(<KK> (' X'A''), '  Temperature RHO Pressure')      
74    FORMAT(<KK> (' dX'A'/dt'), '  dT/dt dRHO/dt') 
75    FORMAT(/,'///////////////////////////////////////////////// ',/,
     2        'DDASPK TIME STEP INTEGRATION COMPLETE ',/,
     1          '///////////////////////////////////////////////// ',//)
  

      


80    FORMAT('At time t = ',ES16.5)
90    FORMAT(/' Final statistics for this run..'/
     1   '   RWORK size =',I5,'   IWORK size =',I4/
     1   '   Number of time steps ................ =',I5/
     1   '   Number of residual evaluations ...... =',I5/
     1   '   Number of preconditioner evaluations  =',I5/
     1   '   Number of preconditioner solves ..... =',I5/
     1   '   Number of nonlinear iterations ...... =',I5/
     1   '   Number of linear iterations ......... =',I5/
     1   '   Average Krylov subspace dimension =',F8.4/
     1   I5,' nonlinear conv. failures,',I5,' linear conv. failures'/)
     
      END SUBROUTINE OUTCVCM  
