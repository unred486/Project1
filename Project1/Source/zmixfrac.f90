!C   Written by Han Ju Lee
!C   University of Maryland
!C   Contact : hanjulee@umd.edu
!C       
!C   see DDASPK.f and examples in ddaspk for more details on how this 
!C   subroutine uses dasppk solver. 

      SUBROUTINE zmf(NATJ,JJ,LOUT,NMAX,X,RW,IW,C,INFO,DT1,RTOL,ATOL,DMIX, &
          XST,ZST,XSTR,XEND,ZBUR,ZAMB,LRCRVR,TOUT,FRZTIM,DMAX)
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
      DOUBLE PRECISION :: C(NMAX),RW(5*NMAX + 2*((NMAX/3) + 1)+91+19*NMAX) &
          ,CCPRIME(NATJ*NMAX),RPAR(4),  RTOL(1),ATOL(1),X(NMAX),DMIX(NMAX)
      DOUBLE PRECISION :: DT1,TOUT,T,DT2,XST,ZST,XSTR,XEND,ZBUR,ZAMB, &
          XSTOLD,FRZTIM,DMAX
      INTEGER :: LOUT, NEQ, NATJ, NHBW,ML,MU,JJ,J,NOUT,IOUT &
                ,IDID,IMOD3,NMAX,JOLD,LRCRVR,IREFLAG
      XSTOLD=XST
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
      ML = NHBW
      MU = NHBW
      IPAR(1) = ML
      IPAR(2) = MU
      IW(27) = LENWP-NEQ
      IW(28) = NEQ
      INFO(15) = 1
    ENDIF
        
!C Set the initial T and TOUT, and call CINIT to set initial values.  
IF (LRSTRT) then
    T=TOUT
    TOUT=T+DT1
else
      T = 0.0D0
      TOUT = DT1
endif

      CALL zmfINIT (C, CCPRIME, X, JJ,NATJ,NEQ,RPAR,ZBUR,ZAMB)  
!c      WRITE (LOUT,40)
!c      WRITE (*,40)      
        IPAR(3)=NEQ
        IPAR(4)=JJ
      IOUT = 0
!      CALL updategridzmf(NATJ,JJ,NMAX,C,X,DMIX,XST,ZST,XSTR,XEND,&
!          IPAR, IW,NEQ)  
      CALL SETPARzmf(X,JJ,DMIX)
      DO WHILE(1)
!C       SETPAR sets pointers to arrays needed for calculation
!C       in varray module for use by RESCVCM subroutine which calculates 
!C       residual 
!C       Calling main solver, DDASPK subroutine
100 CONTINUE
        JOLD=JJ
        CALL DDASPK (RESZMF, NEQ, T, C, CCPRIME, TOUT, INFO, RTOL, &
                 ATOL, IDID, RW,LRW, IW,LIW, RPAR, IPAR,           & 
                 DBANJA, DBANPS)
        TOUT=TOUT+DT1
        CALL updategridzmf(NATJ,JJ,NMAX,C,X,DMIX,XST,ZST,XSTR,XEND,&
          IPAR, IW,NEQ,TOUT,XSTOLD,IREFLAG,FRZTIM,DMAX)    
        IF ((JOLD .NE. JJ) .OR. (IREFLAG .eq. 1)) THEN 
            CALL OUTZMF (T,C,CCPRIME,LOUT,IDID,IW,RW,JJ,NHBW,X,ZST,DMIX) 
                WRITE (LRCRVR) NATJ, JJ,TOUT,DMAX
                WRITE (LRCRVR) (X(J), J=1,JJ)
                WRITE (LRCRVR) (C(J), J=1,JJ)
                write(LOUT,*) 'save file stored'
            go to 210
        end if
        CALL SETPARzmf(X,JJ,DMIX)
        
        IMOD3 = IOUT - 100*(IOUT/100)
        IF (IMOD3 .EQ. 0) THEN 
            CALL OUTZMF (T,C,CCPRIME,LOUT,IDID,IW,RW,JJ,NHBW,X,ZST,DMIX) 
        ENDIF

        
        IOUT=IOUT+1
        IF (IDID .LE. -1) THEN            
            WRITE(LOUT,*)'error condition has been reached'
            CALL OUTZMF (T,C,CCPRIME,LOUT,IDID,IW,RW,JJ,NHBW,X,ZST,DMIX) 
            WRITE(LOUT,73) IDID
            GO TO 210            
        END IF
        
        IF (TIMMX .LT. TOUT) THEN
            CALL DDASPK (RESZMF, NEQ, T, C, CCPRIME, TOUT, INFO, RTOL,  &
                  ATOL,IDID, RW,LRW, IW,LIW, RPAR, IPAR,      &
                 DBANJA, DBANPS)
            CALL OUTZMF (T,C,CCPRIME,LOUT,IDID,IW,RW,JJ,NHBW,X,ZST,DMIX) 
                WRITE (LRCRVR) NATJ, JJ,TOUT,DMAX
                WRITE (LRCRVR) (X(J), J=1,JJ)
                WRITE (LRCRVR) (C(J), J=1,JJ)
                write(LOUT,*) 'save file stored'            
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