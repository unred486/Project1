C   Written by Han Ju Lee
C   University of Maryland
C   Contact : hanjulee@umd.edu
C       
C   see DDASPK.f and examples in ddaspk for more details on how this 
C   subroutine uses dasppk solver. 

      SUBROUTINE cvcm(NATJ,JJ,LOUT,KK,ICKWRK,RCKWRK,RW,IW,WDOT,C,
     1              INFO,cp,h,KSYM,P,WT,DT1,RTOL,ATOL)
C       CVCM.f starts the constant volume and constant mass calculation
C       Calls CVCMINIT, SETPAR, DDASPK, PRINTCVCM, OUTCVCM subroutines

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
C       CP   - Specific heats at constant pressure in mass units
C              for the species.
C                   cgs units - ergs/(mole*K)
C                   Data type - real array
C                   Dimension CP(*) at least KK, the total number of
C                   species.
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
C*****precision > double
      use var
      use varray
      IMPLICIT none
      EXTERNAL RESCVCM,JACCVCM, PSOLCVCM
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER :: ICKWRK(*),IW(*),INFO(20),IPAR(2)
      DOUBLE PRECISION :: C(NATJ*JJ),WDOT(KK),h(*),cp(*),
     1  WT(*),RCKWRK(*),RW(*),CCPRIME(NATJ*JJ),RPAR(NATJ),
     2  RTOL(1),ATOL(1)
      CHARACTER*(*) :: KSYM(*)
      DOUBLE PRECISION :: P,DT1,TOUT,T,DT2
      INTEGER :: LOUT, NEQ, NATJ, NHBW,ML,MU,KK,JJ,K,J,NOUT,IOUT
     1           ,IDID,IMOD3 
      NEQ=NATJ
      NHBW=1   
C       Refer to DDASPK.F subroutine for information about INFO array
C       which controls how ddaspk solver behaves
      WRITE(LOUT,71)
      WRITE(LOUT,70) RTOL(1),ATOL(1),INFO(11),INFO(12)

C          **** Is this the first call for this problem ...
C                yes - set INFO(1) = 0     
      INFO(1) = 0
      IF (INFO(12) .EQ. 0) THEN
         INFO(6) = 1
         IW(1) = NHBW
         IW(2) = NHBW
         WRITE(LOUT,72)NHBW
      ENDIF
      
C       INFO(12) 
C       direct methods - set INFO(12) = 0.
C       Krylov method  - set INFO(12) = 1,
      IF (INFO(12) .EQ. 1) THEN
c        JBG = 0
c        IPAR(2) = JBG
C First set the preconditioner choice JPRE.
C  JPRE = 1 means reaction-only (block-diagonal) factor A_R
C  JPRE = 2 means spatial factor (Gauss-Seidel) A_S
C  JPRE = 3 means A_S * A_R
C  JPRE = 4 means A_R * A_S
C Use IPAR to communicate JPRE to the preconditioner solve routine.
c        JPRE = 1
c        IPAR(1) = JPRE
c        WRITE(LOUT,100)JPRE
c 100    FORMAT(' Preconditioner flag is JPRE =',I3/
c     1      '  (1 = reaction factor A_R, 2 = spatial factor A_S,',
c     2      ' 3 = A_S*A_R, 4 = A_R*A_S )'/)
C Here call DMSET2 if JBG = 0 to set the
C mesh parameters and block-grouping data, and the IWORK segment ID
C indicating the differential and algebraic components.
c
c          CALL DMSET2 (MX, MY, NS, NP, 40, IWORK)
c          WRITE(LOUT,110)
c 110      FORMAT(' No block-grouping in reaction factor')
c        IF (JBG .EQ. 1) THEN
c          NXG = 5
c          NYG = 5
c          NG = NXG*NYG
c          CALL DGSET2 (MX, MY, NS, NP, NXG, NYG, 40, IWORK)
c          WRITE(LOUT,120)NG,NXG,NYG
c 120      FORMAT(' Block-grouping in reaction factor'/
c     1           ' Number of groups =',I5,
c     2           '   (NGX by NGY, NGX =',I3,',  NGY =',I3,')')
c         ENDIF
      ML = NHBW
      MU = NHBW
      IPAR(1) = ML
      IPAR(2) = MU
C        INFO(15) - used when INFO(12) = 1 (Krylov methods).
C               When using preconditioning in the Krylov method,
C               you must supply a subroutine, PSOL, which solves the
C               associated linear systems using P.
C               The usage of DDASPK is simpler if PSOL can carry out
C               the solution without any prior calculation of data.
C               However, if some partial derivative data is to be
C               calculated in advance and used repeatedly in PSOL,
C               then you must supply a JAC routine to do this,
C               and set INFO(15) to indicate that JAC is to be called
C               for this purpose.  For example, P might be an
C               approximation to a part of the matrix A which can be
C               calculated and LU-factored for repeated solutions of
C               the preconditioner system.  The arrays WP and IWP
C               (described under JAC and PSOL) can be used to
C               communicate data between JAC and PSOL.
      INFO(15) = 1
          CALL DMSET2 (NEQ, 40, IW)
          WRITE(LOUT,110)
 110      FORMAT(' No block-grouping in reaction factor')
       ENDIF
        
C Set the initial T and TOUT, and call CINIT to set initial values.  
      T = 0.0D0
      TOUT = DT1
      CALL CVCMINIT (ICKWRK,RCKWRK,C, CCPRIME, WDOT,RPAR,KSYM,h,cp,WT,
     1  KK,JJ,NATJ,LOUT)  
c      WRITE (LOUT,40)
c      WRITE (*,40)      

      IOUT = 0
      DO WHILE(1)
C       SETPAR sets pointers to arrays needed for calculation
C       in varray module for use by RESCVCM subroutine which calculates 
C       residual 
        CALL SETPAR(ICKWRK,RCKWRK,h,cp,wdot,WT,KK,NATJ,P,LOUT)
C       Calling main solver, DDASPK subroutine
100     CONTINUE
        CALL DDASPK (RESCVCM, NEQ, T, C, CCPRIME, TOUT, INFO, RTOL, 
     1            ATOL, IDID, RW,LRW, IW,LIW, RPAR, IPAR, JACCVCM,
     2            PSOLCVCM)
c        IF (IDID .EQ. -1) GO TO 100
C       Printing out the solver status 
c        CALL PRINTCVCM(IW,RW,LOUT,T)   
C       Printing out the solution variables
        IMOD3 = IOUT - 10000*(IOUT/10000)
        IF (IMOD3 .EQ. 0) THEN 
            CALL OUTCVCM (TOUT, C,CCPRIME,ICKWRK,RCKWRK,LOUT,IDID,IW,RW
     1                  ,NHBW,KSYM,KK) 
c            WRITE (LOUT,40)
c            WRITE (*,40)
        ENDIF
        TOUT=TOUT+DT1
        IOUT=IOUT+1
C       EXIT CONDITION
C       IF TOUT REACHES MAX TIME, program exits
C    or IF ERROR CONDITION BY DDASPK is reached, program exits
        IF (IDID .LT. -1) THEN            
            WRITE(LOUT,*)'error condition has been reached'
            CALL OUTCVCM (TOUT, C,CCPRIME,ICKWRK,RCKWRK,LOUT,IDID,IW,RW
     1                  ,NHBW,KSYM,KK) 
            WRITE(LOUT,73) IDID
            GO TO 210            
        END IF
        
        IF (TIMMX .LT. TOUT) THEN
C       Last Time step of DDASPK before reporting
            CALL SETPAR(ICKWRK,RCKWRK,h,cp,wdot,WT,KK,NATJ,P,LOUT)
            CALL DDASPK (RESCVCM, NEQ, T, C, CCPRIME, TOUT, INFO, RTOL,
     1             ATOL,IDID, RW,LRW, IW,LIW, RPAR, IPAR, JACCVCM,
     2            PSOLCVCM)
            CALL OUTCVCM (TOUT, C,CCPRIME,ICKWRK,RCKWRK,LOUT,IDID,IW,RW,
     1                     NHBW,KSYM,KK) 
            WRITE(LOUT,160)TOUT 
            GO TO 210
        ENDIF     
      END DO
        

40    FORMAT(5X,'t',12X,'NQ',8X,'H',8X,'STEPS',
     1       5X,'NNI',5X,'NLI')
70    FORMAT(' Tolerance parameters.. RTOL =',E10.2,'  ATOL =',E10.2,/
     1   ' Internal I.C. calculation flag INFO(11) =',I2,
     2   '   (0 = off, 1 = on)'/
     3   ' Linear Solver Method Flag = ', I2,
     4   '   ( Direct Method - 0  Krylov Method - 1)') 
71    FORMAT('DASPK OUTPUT')
72    FORMAT(' Difference-quotient banded Jacobian,'
     1         ' half-bandwidths =',I4)
73    FORMAT('ERROR CONDITION (IDID) = ', I4)
80    FORMAT(/,'DDASPK SIMULATION STARTED',/)  
160   FORMAT(//' Final time reached =',E12.4//) 
210   CONTINUE
        
      RETURN
      contains
C////////////////////////////////////////////////////////////////
C       SUBROUTINES UNDER CONTAINS///////////////////////////////
C////////////////////////////////////////////////////////////////   
      
      SUBROUTINE PRINTCVCM(IW,RW,LOUT,T)  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION IW(*),RW(*)      
      HU = RW(7)
      NQU = IW(8)
      NST = IW(11)
      NNI = IW(19)
      NLI = IW(20)
c      WRITE (LOUT,60) T,NQU,HU,NST,NNI,NLI
      WRITE (*,60) T,NQU,HU,NST,NNI,NLI
60      FORMAT(E15.5,I5,E14.3,I7,I9,I8)
      END SUBROUTINE PRINTCVCM  
      END SUBROUTINE  cvcm
   
      

