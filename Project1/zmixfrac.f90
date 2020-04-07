!C   Written by Han Ju Lee
!C   University of Maryland
!C   Contact : hanjulee@umd.edu
!C       
!C   see DDASPK.f and examples in ddaspk for more details on how this 
!C   subroutine uses dasppk solver. 

      SUBROUTINE zmf(NATJ,JJ,LOUT,NMAX,X,RW,IW,C,INFO,DT1,RTOL,ATOL,DMIX, &
          XST,ZST,XSTR,XEND)
!C       zmixfrac.f90 starts the constant volume and constant mass calculation
!C       Calls zmfINIT, SETPARzmf, DDASPK,  OUTzmf subroutines

!C       NEQ - the number of equations in the system.
!c       NHBW - Half Band Width
!C       JJ    - Number of Nodes ( In CVCM reactor JJ =1 )
!C       KK    - Total Number of Species
!C       NATJ  - Number of Dependent Variables
!C       U     - Dependent Variable Matrix
!C             dimension(NATJ*JJ)    
!C       U(NM) - mixture density 
!C               cgs units - gm/cm**3
!C       U(NT) - Chamber Temperature
!C               unit - K
!C       U(NY) - Mass fraction
!C       UPRIME(*) -- Set this array to the initial values of the NATJ*JJ first
!C               derivatives of the solution components at the initial
!C               point.  You must dimension YPRIME at least NATJ*JJ in your
!C               calling program.   
!C       P     - Pressure 
!C               unit - dynes/cm**2
!C       WDOT   - Chemical molar production rates of the species.
!C                   cgs units - moles/(cm**3*sec)
!C                   Data type - real array
!C                   Dimension WDOT(*) at least KK, the total number of
!C                   species.
!C        h    - Enthalpies in molar units for the species.
!C                   cgs units - ergs/mole
!C                   Data type - real array
!C       CP   - Specific heats at constant pressure in mass units
!C              for the species.
!C                   cgs units - ergs/(mole*K)
!C                   Data type - real array
!C                   Dimension CP(*) at least KK, the total number of
!C                   species.
!C       WT     - Molecular weights of the species.
!C                   cgs units - gm/mole  
!C     WTM    - Mean molecular weight of the species mixture.
!C                   cgs units - gm/mole
!C       RU     - Universal gas constant.
!C                   cgs units - 8.314E7 ergs/(mole*K)
!C                   Data type - real scalar
!C       rho - density 
!C               cgs units - gm/cm^3
!C       INFO   - DDASPK control array  
!C       RW   - REAL WORK SPACE for DDASPk
!C       IW   - INTEGER WORK SPACE for DDASPK
!C       RCKWRK - REAL WORK SPACE for Chemkin Library
!C       ICKWRK - INTEGER WORK SPACE for Chmekin Library
!C       ATOL   - ABSOLUTE CONVERGENCE CRITERIA
!C       RTOL   - RELATIVE CONVERGENCE CRITERIA
!C       IDID -- reports what the code did, Refer to DDASPK.f for more info
!C       T      - Start Time for current simulation  unit : s
!C       TOUT   - END TIME for current simulation    unit : s    
!C       DT1    - time step  
!C*****precision > double
      use var
      use varray
      use f90_module
      IMPLICIT none
      EXTERNAL RESZMF,DBANJA, DBANPS
!C*****END precision > double
!C
!C*****precision > single
!C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
!C*****END precision > single
!C
      INTEGER :: IW(*),INFO(20),IPAR(4)
      DOUBLE PRECISION :: C(NATJ*JJ),RW(*),CCPRIME(NATJ*JJ),RPAR(NATJ),  &
       RTOL(1),ATOL(1),X(*),DMIX(*)
      DOUBLE PRECISION :: DT1,TOUT,T,DT2,XST,ZST,XSTR,XEND
      INTEGER :: LOUT, NEQ, NATJ, NHBW,ML,MU,JJ,J,NOUT,IOUT &
                ,IDID,IMOD3,NMAX 
      NEQ=NATJ*JJ
      NHBW=1   
!C       Refer to DDASPK.F subroutine for information about INFO array
!C       which controls how ddaspk solver behaves
      WRITE(LOUT,71)
      WRITE(LOUT,70) RTOL(1),ATOL(1),INFO(11),INFO(12)

!C          **** Is this the first call for this problem ...
!C                yes - set INFO(1) = 0     
      INFO(1) = 0
      IF (INFO(12) .EQ. 0) THEN
         INFO(6) = 1
         IW(1) = NHBW
         IW(2) = NHBW
         WRITE(LOUT,72)NHBW
      ENDIF
      
!C       INFO(12) 
!C       direct methods - set INFO(12) = 0.
!C       Krylov method  - set INFO(12) = 1,
      IF (INFO(12) .EQ. 1) THEN
!c        JBG = 0
!c        IPAR(2) = JBG
!C First set the preconditioner choice JPRE.
!C  JPRE = 1 means reaction-only (block-diagonal) factor A_R
!C  JPRE = 2 means spatial factor (Gauss-Seidel) A_S
!C  JPRE = 3 means A_S * A_R
!C  JPRE = 4 means A_R * A_S
!C Use IPAR to communicate JPRE to the preconditioner solve routine.
!c        JPRE = 1
!c        IPAR(1) = JPRE
!c        WRITE(LOUT,100)JPRE
!c 100    FORMAT(' Preconditioner flag is JPRE =',I3/
!c     1      '  (1 = reaction factor A_R, 2 = spatial factor A_S,',
!c     2      ' 3 = A_S*A_R, 4 = A_R*A_S )'/)
!C Here call DMSET2 if JBG = 0 to set the
!C mesh parameters and block-grouping data, and the IWORK segment ID
!C indicating the differential and algebraic components.
!c
!c          CALL DMSET2 (MX, MY, NS, NP, 40, IWORK)
!c          WRITE(LOUT,110)
!c 110      FORMAT(' No block-grouping in reaction factor')
!c        IF (JBG .EQ. 1) THEN
!c          NXG = 5
!c          NYG = 5
!c          NG = NXG*NYG
!c          CALL DGSET2 (MX, MY, NS, NP, NXG, NYG, 40, IWORK)
!c          WRITE(LOUT,120)NG,NXG,NYG
!c 120      FORMAT(' Block-grouping in reaction factor'/
!c     1           ' Number of groups =',I5,
!c     2           '   (NGX by NGY, NGX =',I3,',  NGY =',I3,')')
!c         ENDIF
      ML = NHBW
      MU = NHBW
      IPAR(1) = ML
      IPAR(2) = MU
      IW(27) = LENWP
      IW(28) = NATJ*NMAX
!C        INFO(15) - used when INFO(12) = 1 (Krylov methods).
!C               When using preconditioning in the Krylov method,
!C               you must supply a subroutine, PSOL, which solves the
!C               associated linear systems using P.
!C               The usage of DDASPK is simpler if PSOL can carry out
!C               the solution without any prior calculation of data.
!C               However, if some partial derivative data is to be
!C               calculated in advance and used repeatedly in PSOL,
!C               then you must supply a JAC routine to do this,
!C               and set INFO(15) to indicate that JAC is to be called
!C               for this purpose.  For example, P might be an
!C               approximation to a part of the matrix A which can be
!C               calculated and LU-factored for repeated solutions of
!C               the preconditioner system.  The arrays WP and IWP
!C               (described under JAC and PSOL) can be used to
!C               communicate data between JAC and PSOL.
      INFO(15) = 1
!          CALL DMSET2 (NEQ, 40, IW)
!          WRITE(LOUT,110)
! 110      FORMAT(' No block-grouping in reaction factor')
       ENDIF
        
!C Set the initial T and TOUT, and call CINIT to set initial values.  
      T = 0.0D0
      TOUT = DT1
      CALL zmfINIT (C, CCPRIME, X, JJ,NATJ,NEQ)  
!c      WRITE (LOUT,40)
!c      WRITE (*,40)      
        IPAR(3)=NEQ
        IPAR(4)=JJ
      IOUT = 0
      DO WHILE(1)
!C       SETPAR sets pointers to arrays needed for calculation
!C       in varray module for use by RESCVCM subroutine which calculates 
!C       residual 
!        CALL updategridzmf(NATJ,JJ,NMAX,C,X,DMIX,XST,ZST,XSTR,XEND,IPAR)    
        CALL SETPARzmf(X,JJ,DMIX)
!C       Calling main solver, DDASPK subroutine
100     CONTINUE
        CALL DDASPK (RESZMF, NEQ, T, C, CCPRIME, TOUT, INFO, RTOL, &
                 ATOL, IDID, RW,LRW, IW,LIW, RPAR, IPAR,           & 
                 DBANJA, DBANPS)

!c        IF (IDID .EQ. -1) GO TO 100
!C       Printing out the solver status 
!c        CALL PRINTCVCM(IW,RW,LOUT,T)   
!C       Printing out the solution variables
        
        IMOD3 = IOUT - 1000*(IOUT/1000)
        IF (IMOD3 .EQ. 0) THEN 
            CALL OUTZMF (T,C,CCPRIME,LOUT,IDID,IW,RW,JJ,NHBW,X,ZST) 
        ENDIF
        TOUT=TOUT+DT1
        IOUT=IOUT+1
!C       EXIT CONDITION
!C       IF TOUT REACHES MAX TIME, program exits
!C    or IF ERROR CONDITION BY DDASPK is reached, program exits
        IF (IDID .LT. -1) THEN            
            WRITE(LOUT,*)'error condition has been reached'
            CALL OUTZMF (T,C,CCPRIME,LOUT,IDID,IW,RW,JJ,NHBW,X,ZST) 
            WRITE(LOUT,73) IDID
            GO TO 210            
        END IF
        
        IF (TIMMX .LT. TOUT) THEN
!C       Last Time step of DDASPK before reporting
!            CALL updategridzmf(NATJ,JJ,NMAX,C,X,DMIX,XST,ZST,XSTR,XEND, &
!                IPAR)

            CALL SETPARzmf(X,JJ,DMIX)
            CALL DDASPK (RESZMF, NEQ, T, C, CCPRIME, TOUT, INFO, RTOL,  &
                  ATOL,IDID, RW,LRW, IW,LIW, RPAR, IPAR,      &
                 DBANJA, DBANPS)
            CALL OUTZMF (T,C,CCPRIME,LOUT,IDID,IW,RW,JJ,NHBW,X,ZST) 
            WRITE(LOUT,160)TOUT 
            GO TO 210
        ENDIF     
      END DO
        

40    FORMAT(5X,'t',12X,'NQ',8X,'H',8X,'STEPS', &
            5X,'NNI',5X,'NLI')
70    FORMAT(' Tolerance parameters.. RTOL =',E10.2,'  ATOL =',E10.2,/ &
        ' Internal I.C. calculation flag INFO(11) =',I2,        &
        '   (0 = off, 1 = on)'/                                 &
        ' Linear Solver Method Flag = ', I2,                    &
        '   ( Direct Method - 0  Krylov Method - 1)') 
71    FORMAT('DASPK OUTPUT')
72    FORMAT(' Difference-quotient banded Jacobian,'    &
              ' half-bandwidths =',I4)
73    FORMAT('ERROR CONDITION (IDID) = ', I4)
80    FORMAT(/,'DDASPK SIMULATION STARTED',/)  
160   FORMAT(//' Final time reached =',E12.4//) 
210   CONTINUE
        
      RETURN 
      END SUBROUTINE  zmf
   
      

