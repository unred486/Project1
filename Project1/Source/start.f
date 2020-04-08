C///////////////////////////////////////////////////////////////////
C        ONE DIMENSIONAL LAMINAR SPHERICAL DIFFUSION FLAME
C           (MODIFIED FROM PREMIX VERSION 2.55)
C        WRITTEN BY:
C            ZHEN SUN
C            DEPARTMENT OF MECHANICAL ENGINEERING
C            WASHINGTON UNIVERSITY IN ST. LOUIS
C            ST. LOUIS, MO 63130
C            ZS3@CEC.WUSTL.EDU
C            Version 2.x modified by
C            Beei-Huan Chao and Karl Santa
C            Department of Mechanical Engineering
C            University of Hawaii at Manoa
C            Honolulu, HI 96822
C            
C            NEW SPHDIFF PROGRAM
C            UPDATED BY HANJU LEE
C            UNIVERSITY OF MARYLAND
C            hanjulee@umd.edu
C///////////////////////////////////////////////////////////////////
C        
C
C SUBROUTINE START
C CALLED BY FLDRIV 
C INITIALIZED THE VECTOR SOLUTION BASED ON THE INPUT FILE COMMANDS.
C Updates to initialize differently to include homogeneous reactor
C
C CALL SUBROUTINE LINWMX
C CALL SUBROUTINE TEMP
C CALL SUBROUTINES FROM CHEMKIN LIBRARY


      SUBROUTINE START (KK, NMAX, NATJ, JJ, LOUT,
     1                  NREAC, NINTM, NPROD, REAC, XINTM, PROD, KR, KI,
     2                  KP, FLRT, NTEMP, XGIVEN, TGIVEN, XSTR, XCEN,
     3                  XEND, WMIX, ICKWRK, RCKWRK, Y, SI, EPS, JFIXT,
     4                  TFIXT, X, S, EPA, TBUR,TAMB,KSYM,
     5                  HC,C_P,RK,phi,P,ZGIVEN,NZMIX,NDMIX,XXZ,XXD,DM,
     6                  DMIX,ZST )
       USE f90_module
       USE var
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      

c      INCLUDE 'LOCS.cmn'
      INCLUDE 'TEST.cmn' 
      INCLUDE 'LOGICAL.cmn'
C
      DIMENSION REAC(*), XINTM(*), PROD(*), KR(*), KI(*), KP(*),
     1          ICKWRK(*), RCKWRK(*), Y(*), SI(*), EPS(*),
     2          XGIVEN(*), TGIVEN(*), X(*), S(NATJ, *), EPA(*),
     3          ZGIVEN(*),XXZ(NMAX),XXD(*),DM(*),DMIX(*),XNEW(NMAX),
     4          SS(NATJ,NMAX)   
     
      CHARACTER*(*) KSYM(*)
      character(len=32) :: my_fmt
C
      LOGICAL LBURNR, LTIME, LTIME2, LMOLE, LENRGY, LTDIF,
     1        LRSTRT, LCNTUE, LASEN, LRSEN, LESEN, LMSEN, LPSEN, LMULTI,
     2        LVARMC, RSTCNT, ERROR, FUNCTN, JACOBN, REENTR, SOLVE,
     3        STORE, SUCCES, ADAPT, SHOW, UPDATE, 
     4        ENERGY, BURNER, LUSTGV, IERR, LVCOR, LHSEN,
     5        LADAPT, LSTEADY, ADIAMF, SOOTP, freeMolOnly,LPOST,LREF,
     6        LGLOB
c      LOGICAL LMOLE, LUMESH, LBURNR, SOOTP, LGLOB,LHOM,LTEST,LSOOTINIT

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT-
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   NMAX   - MAXIMUM NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT, AND
C            THE EXACT FIRST DIMENSION OF S(NATJ,*), SN(NATJ,*), AND
C            F(NATJ,*).  NATJ=KK+2.
C   JJ     - NUMBER OF MESH POINTS IN THE STARTING MESH.
C   LMOLE  - IF LMOLE=.TRUE. THEN THE INPUT AND OUTPUT IS IN TERMS OF
C            MOLE FRACTIONS.
C   LUMESH - IF LUMESH=.TRUE. THEN A UNIFORM STARTING MESH IS USED.
C   LBURNR - IF LBURNR=.TRUE. THIS IS A BURNER STABILIZED FLAME PROBLEM.
C            IF LBURNR=.FALSE. THIS IS A FREELY PROPAGATING ADIABATIC
C            FLAME.
C   LHOM   - IF LHOM= .TRUE. THEN HOMOGENEOUS REACTOR
C   NREAC  - NUMBER OF REACTANT SPECIES SPECIFIED.
C   NINTM  - NUMBER OF INTERMEDIATE SPECIES SPECIFIED.
C   NPROD  - NUMBER OF PRODUCT SPECIES SPECIFIED.
C   REAC   - ARRAY OF REACTANT SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION REAC(*) AT LEAST KK.
C   XINTM  - ARRAY OF INTERMEDIATE SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION XINTM(*) AT LEAST KK.
C   PROD   - ARRAY OF PRODUCT SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION PROD(*) AT LEAST KK.
C   KR     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE REACTANTS AS
C            SPECIFIED IN THE REAC ARRAY.
C            DIMENSION KP(*) AT LEAST KK.
C   KI     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE INTERMEDIATES AS
C            SPECIFIED IN THE XINTM ARRAY.
C            DIMENSION KI(*) AT LEAST KK.
C   KP     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE PRODUCTS AS
C            SPECIFIED IN THE PROD ARRAY.
C            DIMENSION KP(*) AT LEAST KK.
C   FLRT   - THE MASS FLOW RATE IN A BURNER STABILIZED PROBLEM.
C              CGS UNITS - GM/(CM**2-SEC)
C   NTEMP  - THE NUMBER OF TEMPERATURE-POSITION PAIRS GIVEN FOR BURNER
C            STABILIZED PROBLEMS WHERE THE ENERGY EQUATION IS NOT
C            SOLVED, AND THE TEMPERATURES ARE INTERPOLATED FROM FIXED
C            PROFILES.
C   XGIVEN - THE DISTANCES FROM THE BURNER AT WHICH TEMPERATURES ARE
C            SPECIFIED.
C              DIMENSION XGIVEN(*) AT LEAST NTEMP.
C   TGIVEN - THE TEMPERATURES AT THE XGIVEN LOCATIONS.
C              DIMENSION TGIVEN(*) AT LEAST NTEMP.
C   XSTR   - BEGINNING X POSITION FOR THE MESH.
C              CGS UNITS - CM
C   XCEN   - X POSITION WHERE THE INITIAL STARTING ESTIMATES ARE
C            CENTERED.
C              CGS UNITS - CM
C   XEND   - ENDING POSITION FOR THE MESH.
C              CGS UNITS - CM
C   WMIX   - WIDTH OF THE MIXING REGION OVER WHICH THE STARTING
C            ESTIMATES ARE FIT.
C              CGS UNITS - CM
C   TFIXT  - TEMPERATURE THAT IS SPECIFIED TO BE HELD FIXED FOR FREE
C            FLAMES.
C              CGS UNITS - K
C   KSYM   - ARRAY OF CHARACTER*(*) CHEMKIN SPECIES NAMES.
C              DIMENSION KSYM LEAST KK.
C
C WORK AND SCRATCH SPACE-
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C   Y      - ARRAY OF MASS FRACTIONS.
C              DIMENSION Y(*) AT LEAST KK.
C   SI     - ARRAY OF THE SUM OF THE INTERMEDIATES.  USED FOR
C            DECREMENTING THE REACTANTS AND PRODUCTS TO KEEP
C            THE MOLE FRACTIONS SUMMED TO 1.
C              DIMENSION SI(*) AT LEAST JJ.
C
C OUTPUT-
C   EPS    - INLET MASS FLUX FRACTONS.
C              DIMENSION EPS(*) AT LEAST KK.
C   JFIXT  - MESH POINT THAT HAS A FIXED TEMPERATURE IN AN ADIABATIC
C            FLAME PROBLEM.
C   X      - STARTING MESH LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ.
C   S      - STARTING SOLUTION ESTIMATES.  S(N,J) IS THE STARTING
C            ESTIMATE FOR COMPONENT N AT NODE J.
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION
C              AND AT LEAST JJ FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          write(*,*) LZMIX
      IF (LHOM) THEN
C       Initializing mass fraction 
        DO K =1, KK
            S(NYS+K,1) = 0.0
        END DO
                
        DO N=1,NREAC
            S(NYS+KR(N),1)=REAC(KR(N))
        END DO

C       Initializing temperature
            S(NT,1)    = TGIVEN(1)



C       Converting to mass fraction if given mole fraction 
        IF (LMOLE) THEN
         CALL CKXTY (S(NY,1), ICKWRK, RCKWRK, Y)
        ENDIF

        DO K=1,KK
            S(NYS+K,1)=Y(K)
        END DO
        
C       Initializing Pressure 
C       S(NM,1) is a pressure variable for CVCM reactor case
      CALL CKRHOY (P, S(NT,1), S(NY,1), ICKWRK, RCKWRK, RHO)
            S(NM,1)    = RHO
C       Limiting homogeneous closed reactor node to just one node
C       Total number of grid is 1
        JJ=1
c////////////////////////////////////////////////////////////////        
C       Testing initialization
c////////////////////////////////////////////////////////////////
c        WRITE(my_fmt , '(I0)')  KK+2
c       my_fmt = '('//trim(my_fmt)//'(1PE12.4))'       
c        Write(LOUT,99) 
c 99     format('Initial solution')
c        WRITE(LOUT,98) (trim(KSYM(K)),K=1,KK)
c 98     format('TEMP(K)',2X,'FLOWRATE',<KK> (' X'A'') )  
c        WRITE(LOUT,my_fmt) S(NT,1), S(NM,1), (S(NYS+K,1),K=1,KK)
c////////////////////////////////////////////////////////////////

      ELSEIF (LZMIX) THEN 
C
C           SET UNIFORM X MESH COORDINATES
C
         DX = (XEND-XSTR) / (JJ-1)
         DO J = 1,JJ
            X(J) = XSTR + DX*(J-1)
         ENDDO
      
      DO J = 1, JJ
         CALL TEMP (NZMIX, X(J), XXZ, ZGIVEN, S(NZ,J))
         CALL TEMP (NDMIX, X(J), XXD, DM,DMIX(J))
      ENDDO    
      
!      IRMAX = 400
!      CALL XFIND(XXZ,ZGIVEN,ZST,XST)
!      CALL GRIDGEN(XSTR,XEND,XST,XNEW,IRMAX,IREFINE)
!      CALL COMPRESSZMF(X,S,JJ,XNEW,SS,IRMAX,NATJ,XST)
!!      do J=1,JJ
!!      write(*,*) S(NATJ,J)
!!      end do
!      JJ=IRMAX   
!      CALL DCOPY (NATJ*JJ, SS, 1, S, 1)
!      CALL DCOPY (JJ,XNEW,1,X,1) 
!      DO J = 1, JJ
!         CALL TEMP (NDMIX, X(J), XXD, DM,DMIX(J))
!      ENDDO 
      ELSE
cC        Burner or Premix flame initialization
cC        These flames have transport terms in their governing equations      
cC
cC
cC        INITIALIZE MASS FLUX FRACTIONS AND MASS FRACTIONS TO 0.
cC
c      DO K = 1, KK
c         EPS(K) = 0.0
c	 EPA(K) = 0.0
c      ENDDO 

c      DO K = 1, KK
c         DO J = 1, NMAX
c            S(NYS+K, J) = 0.0
c         ENDDO
c      ENDDO
cC
c      DO J = 1, JJ
c     	 DO L = 1, N_moments
c      	  	 S(NSM+L,J) = 0.0
c      	 ENDDO 
c      ENDDO

c      IF (LUMESH) THEN
cC
cC           SET UNIFORM X MESH COORDINATES
cC
c         DX = (XEND-XSTR) / (JJ-1)
c         DO J = 1,JJ
c            X(J) = XSTR + DX*(J-1)
c         ENDDO
c      ENDIF
cC
cC          FOR FREE FLAMES, ADD THE FIXED TEMPERATURE POINT TO THE MESH
cC
c      IF (.NOT. LBURNR) THEN
cC
c         DO 300 N = 1, NTEMP
c            NFIXT = N
c            IF (TGIVEN(N) .GE. TFIXT) GO TO 350
c300      CONTINUE
c         WRITE (LOUT,*) ' ERROR...NO USABLE MESH LOCATION FOR TFIXT'
c         STOP

c350      CONTINUE
cC
c         IF (TGIVEN(NFIXT) .EQ. TFIXT) THEN
c            XFIXT = XGIVEN(NFIXT)
c         ELSE
c            XFIXT = XGIVEN(NFIXT-1) + (TFIXT-TGIVEN(NFIXT-1)) *
c     1              (XGIVEN(NFIXT) - XGIVEN(NFIXT-1)) /
c     2              (TGIVEN(NFIXT) - TGIVEN(NFIXT-1))
c         ENDIF
cC
c         DO 400 J = 1, JJ
c            IF (XFIXT .EQ. X(J)) THEN
c               JFIXT = J
c               GO TO 700
c            ENDIF
c400      CONTINUE
cC
c         DO J = 1, JJ-1
c            IF (XFIXT.GT.X(J) .AND. XFIXT.LT.X(J+1)) JFIXT = J+1
c         ENDDO
cC
c         JJ = JJ + 1
c         DO 600 J = JJ, JFIXT+1, -1
c            X(J) = X(J-1)
c600      CONTINUE
cC
c         X(JFIXT) = XFIXT
cC
c700      CONTINUE
cC
c      ENDIF
cC
cC         SET INTERMEDIATE GAUSSIANS PROFILE
cC
   
c      DO N = 1,NINTM

c         GBAS = 0.0
c         GM   = XINTM(KI(N))
c         GMIX = 0.15*GM + GBAS
c         W    = -LOG((GMIX-GBAS)/GM)/(WMIX/2.)**2

c         DO J = 1, JJ
c            	S(NYS+KI(N), J) = GM*EXP(-W*(X(J)-XCEN)**2) + GBAS
c         ENDDO
c      ENDDO
      
cC
c      IF (LSOOTINIT.OR.LTEST) CALL InitSootPrfl (NATJ,JJ,X,S,ICKWRK,
c     1                           RCKWRK,WT,REAC,NREAC,KR,KK,KSYM)

cC            SUM THE INTERMEDIATES AT EACH J
cC
c      DO J = 1, JJ
c         SI(J) = 0.

c         DO  N = 1, NINTM
c            SI(J) = SI(J) + S(NYS+KI(N), J)
c         ENDDO

c      ENDDO
cC
cC            SET STARTING SPECIES PROFILES
cC
c      DO 1400 J = 1, JJ
c         CALL LINWMX (WMIX, XCEN, X(J), XRE, XPD)
c         FAC = 1.0 - SI(J)
c         DO 1100 N = 1, NREAC
c            S(NYS+KR(N), J) = (XPD*PROD(KR(N)) + XRE*REAC(KR(N))) * FAC
c1100     CONTINUE
c         DO 1300 N = 1, NPROD
c            DO 1200 L = 1, NREAC
c               IF (KP(N) .EQ. KR(L)) GO TO 1300
c1200        CONTINUE
c            S(NYS+KP(N), J) = (XPD*PROD(KP(N)) + XRE*REAC(KP(N))) * FAC
c1300     CONTINUE
c1400  CONTINUE
cC
cC             SET THE MASS FLUX FRACTION BOUNDARY CONDITIONS
cC
c      IF (LMOLE) THEN
c         CALL CKXTY (REAC, ICKWRK, RCKWRK, EPS)
c      ELSE
c         DO  N = 1, NREAC
c            EPS(N) = REAC(N)
c         END DO
c      ENDIF


c      IF (LMOLE) THEN
c         CALL CKXTY (PROD, ICKWRK, RCKWRK, EPA)
c      ELSE
c         DO N = 1, NPROD
c            EPA(N) = PROD(N)
c         END DO
c      ENDIF
cC
cC             CONVERT STARTING ESTIMATES TO MASS FRACTION, IF NEEDED
cC
c      IF (LMOLE) THEN
c         DO J = 1,JJ
c            CALL CKXTY (S(NY, J), ICKWRK, RCKWRK, Y)
c            DO K = 1, KK
c               S(NYS+K, J) = Y(K)
c            ENDDO
c         ENDDO
c      ENDIF
cC     Initialization of species character location variables     
c      iO2=0
c      iFUEL=0
c      IF (LGLOB) THEN
cC       Storing fuel mass fraction at the burner surface and the 
cC       oxidizer mass fraction at the ambient side      
c        CALL CKCOMP ('O2'  , KSYM, KK, iO2  )
c        CALL CKCOMP ('CH4'  , KSYM, KK, iFUEL  )
c        RFUEL  =   S(NYS+iFUEL, 1)
c        ROXI   =   S(NYS+iO2,  JJ)
c!       nn           number of data points for initialization   
c        nn=100
c        CALL INIDROP(TBUR,TAMB,XCEN,XSTR,XEND,RFUEL,ROXI,HC,C_p,
c     1             RK,PHI,RR,T,nn,Lout)
c!   TEMP2 initializes the temperature profile created from INIDROP subroutine
c!   in XGIVEN and TGIVEN matrix which are the matrices that this program
c!   uses for initialization
c        Call TEMP2(RR,T,nn,XGIVEN,TGIVEN)
c        NTEMP=nn
c      END IF
        
cC
cC        SET THE TEMPERATURE AND FLOW RATE PROFILES
cC
c      DO J = 1, JJ
c         CALL TEMP (NTEMP, X(J), XGIVEN, TGIVEN, S(NT,J))
c         S(NM,J) = FLRT
c      ENDDO
c      IF (.NOT. LBURNR) S(NT,JFIXT)=TFIXT
cC
c      RETURN
      END IF
      END
C
