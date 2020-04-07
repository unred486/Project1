C//////////////////////////////////////////////////////////////////////
C////////////////// OUTZMF ////////////////////////
C//////////////////////////////////////////////////////////////////////  
C
      SUBROUTINE OUTZMF(T,U,UPRIME,LOUT,IDID,IW,RW,JJ,NHBW,X,ZST) 
      use f90_module
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
      DOUBLE PRECISION :: U(*),UPRIME(*),RW(*),X(*)
      INTEGER :: IW(*)
      DOUBLE PRECISION :: T,dummy,AVDIM,XST,ZST
      INTEGER :: NHBW,IDID,LOUT,J,JJ,K,MBAND,NST,NPE,PRE,
     1          NNI,NLI,NPS,NCFN,NCFL,LFTEMPTIME,LSPACE
      LOGICAL :: exist_ftemp_time

      LFTEMPTIME = 2617   ! Flame Temperature vs time
      LSPACE =2618        ! Spatial mixture fraction
      OPEN(UNIT=LSPACE,FORM='FORMATTED',FILE='SpatialProfile.tec',
     1 POSITION='append')
      write(LSPACE,8026)
      write(LSPACE,80277)
      write(LSPACE,8029) T
      inquire(file='Temptimex.tec', exist=exist_ftemp_time)
      
      if (exist_ftemp_time) then
        OPEN(UNIT=LFTEMPTIME, file='Temptimex.tec', Status='OLD', 
     1 POSITION='append', FORM='FORMATTED', ACTION='READWRITE')
      else
        OPEN(UNIT=LFTEMPTIME, file='Temptimex.tec', Status='NEW', 
     1  FORM='FORMATTED',ACTION='READWRITE')
      !Write names and tecplot related variables only when the file is created

      WRITE (LFTEMPTIME, 8026)
      WRITE (LFTEMPTIME, 80276)
      WRITE (LFTEMPTIME, 8029) T
      endif 
      
      CALL XFIND(X,U,ZST,XST)
      WRITE (LFTEMPTIME, 80298) T,XST
      
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
      WRITE (LOUT,*)
      WRITE (LOUT,50)
      
      DO J=1,JJ
          WRITE(LOUT, 51) X(J),U(J)
          WRITE(LSPACE, 51) X(J),U(J)
      END DO
      
      write(LOUT,*) IDID
      write(LOUT,75)
      
      close(LFTEMPTIME)
      close(LSPACE)
75    FORMAT(/,'///////////////////////////////////////////////// ',/,
     2        'DDASPK TIME STEP INTEGRATION COMPLETE ',/,
     1          '///////////////////////////////////////////////// ',//)
  
50    FORMAT('X Zmix')
51    FORMAT(2ES14.3)      
8026  FORMAT('VARIABLES=')      
80298 FORMAT (2 (1PE12.4),2X)
80276 FORMAT ('Time(s)' 2X 'XST' )
80277 FORMAT ('X' 2X 'Z')
8029  FORMAT('ZONE T="',1PE11.3,' s"')
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
     
      END SUBROUTINE OUTZMF
