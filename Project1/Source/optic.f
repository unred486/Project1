c       /////////////////////////////////
c	//this is the function of inside integral for radiation
c	//at burner surface
c	//given: array of XX,TJ (radius and its temperature)
c	//       JJ,KK: number of array elements grids and species
C	//       YJ: array of molar fraction of specis at J
c	//       RB,X: burner outer surface radius and any radius
c	////////////////////////////////////
	FUNCTION FUNC(CCKWRK,ICKWRK,RCKWRK,LOUT,
     1                XX,TJ,YJ,KKGAS,JJ,X,RB,TB)

        CHARACTER*(*) CCKWRK(*)
	INTEGER ICKWRK(*)
	DOUBLE PRECISION RCKWRK(*)
	DOUBLE PRECISION X,T,RB,TB,FUNC,CON1,CON2,
     1                   XX(*), TJ(*), YJ(*), SKAPA
	INTEGER LOUT,KKGAS,JJ
	CON1=(1.-SQRT(1.-(RB/X)**2))*X*X
	CALL FSGEINT(XX,TJ,JJ,X,T)
	CON2=T**4-TB**4
	CALL FSRCOEF(CCKWRK,ICKWRK,RCKWRK,LOUT,
     1               YJ,KKGAS,T,SKAPA)
	FUNC=SKAPA*CON1*CON2
	RETURN
	END


c	//////////////////////////////////
c	//this is the function of Plank Mean Absorption 
c	//coefficient of CO2
c	//given: T, temperature
c	//return: FUNCO2KP, the coefficient
c	///////////////////////////////

	FUNCTION FUNCO2KP(T)
	DOUBLE PRECISION T, TR

	IF(T .LE. 554) THEN
	  FUNCO2KP = 3.e-4*T+0.2203
	ELSEIF(T .GT. 758) THEN
	  FUNCO2KP = 1.3495 * EXP(-0.0016 * T)
        ELSE
          FUNCO2KP = -0.00000185 * T * T + 0.0025 * T - 0.4307
        END IF

	RETURN
	END


c	//////////////////////////////////
c	//this is the function of Plank Mean Absorption 
c	//coefficient of H2O
c	//given: T, temperature
c	//return: FUNH2OKP, the coefficient
c	///////////////////////////////

	FUNCTION FUNH2OKP(T)
	DOUBLE PRECISION T

        FUNH2OKP =  439539 * T **(-2.3046) 
	IF(FUNH2OKP .GT. 0.172784) FUNH2OKP = 0.172784

	RETURN
	END


c	////////////////////////////////
c	//this subroutine gives linear intepolation
c	//given: XX,TJ, array of radius and its temperature
c	//       JJ, number of array elements
c	//       X, any radius
c	//return: T, temperature at X
c	/////////////////////////////////////
	SUBROUTINE FSGEINT(XX,TJ,JJ,X,T)
	DOUBLE PRECISION XX(*),TJ(*), X,T
	INTEGER J,JJ,NL,NU, JL, JU, NFD
	NL=0
	NU=0
	NFD=0
	DO 20 J=1,JJ,1
  	IF(X .GT. XX(J)) THEN
 	  NL=1
	  JL=J
	ELSE IF(X .LT. XX(J)) THEN
	  NU=1
	  JU=J
	ENDIF
	IF(NL. EQ. 1 .AND. NU. EQ. 1) THEN
	  T=TJ(JU)+(X-XX(JU))*(TJ(JL)-TJ(JU))/(XX(JL)-XX(JU))
	  NL=0
	  NU=0
	  NFD=1
	ENDIF
20	CONTINUE
        IF(NFD .EQ. 0) THEN
          WRITE(*,*)"ERROR: ",X," IS NOT WITHIN THE REGION."
	  STOP
        ENDIF

	RETURN
	END

c	/////////////////////////////////////////
c	//this is GAUSS QUADRATURAL intergation
c	//given: NPAIR of intergation PAIR, FUNC, array of 
c	//       radius and temp XX,TJ, number of array elements,JJ
c	//       and burner surface radius RB and its TB
c	//Return: SS, the integral of the function FUNC between
c	//        all NPAIR regions, by ten-point Gauss Legendre
c	//        integration: the function is evaluated exactly ten 
c	//        times as interior points in the range.
c	/////////////////////////////////////////////

	SUBROUTINE FSQGAUS(CCKWRK,ICKWRK,RCKWRK,LOUT,YJ,
     1                     FUNC,XX,TJ,KK,JJ,RB,TB,PAIR,NPAIR,SS)

      common /OPTHIN/ SIGMA, SLEPS, IFXTM, TFXTM 
	DOUBLE PRECISION PAIR(100,2),SS,SSTEMP,FUNC,RB,TB, SIGMA, SLEPS
	EXTERNAL FUNC
	INTEGER I,J,NPAIR,KK,JJ,LOUT,ICKWRK(*), IRADI
        CHARACTER*(*) CCKWRK(*)
c	the abscissas and weights
	DOUBLE PRECISION DX, XM,XR,X1,X2,W(5),X(5)
	DOUBLE PRECISION XX(*),TJ(*),RCKWRK(*), YJ(*)
	SAVE W, X
	DATA W/.2955242247, .2692667193, .2190863625,
     *         .1494513491, .0666713443/
	DATA X/.1488743389, .4333953941, .6794095682, 
     *         .8650633666, .9739065285/
 	SS=0
	DO 30 I=1,NPAIR,1

	XM=0.5*(PAIR(I,2)+PAIR(I,1))
	XR=0.5*(PAIR(I,2)-PAIR(I,1))

c	SS will be twice the average value of the function 
c	since the ten weights (five numbers above each used twice)
c	sum to 2
	SSTEMP=0
	DO 11 J=1, 5
	  DX=XR*X(J)
	  X1 = XM + DX
	  X2 = XM - DX
	  F1=FUNC(CCKWRK,ICKWRK,RCKWRK,LOUT,
     1            XX,TJ,YJ,KK,JJ,X1,RB,TB)
	  F2=FUNC(CCKWRK,ICKWRK,RCKWRK,LOUT,
     1            XX,TJ,YJ,KK,JJ,X2,RB,TB)
	  SSTEMP=SSTEMP+W(J)*(F1+F2)
11	CONTINUE
c	scale the answer to the range of integration
	SSTEMP=XR*SSTEMP
	SS=SS+SSTEMP
30	CONTINUE
        SS=SS * 2. * SIGMA /RB /RB
	RETURN
	END


c       ////////////////////////////////////////////////////
	SUBROUTINE FSRE(RB,JB,TB,XJ,JJ,TJ,NPAIR,PAIR)
c       /////////////////////////////////////////////////////
c	//given  XJ, TJ, array of location and temperature
c	//       RB, JB, TB, value of burner radius, grid number, 
c	//	 and burner temperature
c	//return NPAIR, number of pairs regions with temp higher than TB
c	//	 PAIR, array of pairs with starting and ending position 
c	//	 within which temperature is higher than TB
c       ///////////////////////////////////////////////////////
        DIMENSION PAIR(100, 2),	XJ(*),TJ(*)
	DOUBLE PRECISION RB,TB,PAIR, REL, REU, XJ, TJ
	INTEGER JB,JJ,NPAIR,J, NL, NU, INC,DEC
	INC=0
	DEC=0
	NPAIR=0
	REL=RB
	NL=1
	NU=0
	DO 10 J=JB+1, JJ, 1
	  IF(TJ(J). GT. TB) INC=1
	  IF(TJ(J). LT. TB) DEC=1
	  IF(TJ(J) .LT. TB .AND. INC.EQ.1) THEN
	    REU=XJ(J)+(TB-TJ(J))*(XJ(J-1)-XJ(J))/(TJ(J-1)-TJ(J))
	    NU=1
	    INC=0	
	    IF(NL.EQ.1 .AND. NU.EQ.1) THEN
	      NPAIR=NPAIR+1
	      PAIR(NPAIR, 1)=REL
	      PAIR(NPAIR, 2)=REU
	      NL=0
	      NU=0
	    ENDIF
	  ELSE IF(TJ(J) .GT. TB .AND. DEC.EQ.1) THEN
	    REL=XJ(J)+(TB-TJ(J))*(XJ(J-1)-XJ(J))/(TJ(J-1)-TJ(J))
	    NL=1
	    DEC=0
	    IF(NL.EQ.1 .AND. NU.EQ.1) THEN
	      NPAIR=NPAIR+1
	      PAIR(NPAIR, 1)=REL
	      PAIR(NPAIR, 2)=REU
	      NL=0
	      NU=0
	    ENDIF
	  ENDIF
10	CONTINUE

	RETURN
	END

c	//////////////////////////////
c	//calculate sum of radiation coefficient
c	////////////////////////////////
	SUBROUTINE FSRCOEF(CCKWRK,ICKWRK,RCKWRK,LOUT,
     1                     YJ,KKGAS,T,SKAPA)
	DOUBLE PRECISION SKAPA, T, YJ(200), MJ(200)
	INTEGER K, KKGAS, JCO2, JH2O, LOUT
	CHARACTER*16 KNAME(200)
        CHARACTER*(*) CCKWRK(*)
	INTEGER ICKWRK(*)
	DOUBLE PRECISION RCKWRK(*)
	LOGICAL KERR
	CALL CKSYMS(CCKWRK, LOUT, KNAME, KERR)
	CALL CKYTX(YJ, ICKWRK, RCKWRK, MJ)
	DO K = 1, KKGAS
	  IF(KNAME(K) .EQ. 'CO2') JCO2 = K
	  IF(KNAME(K) .EQ. 'H2O') JH2O = K
	ENDDO
	SKAPA = MJ(JCO2) * FUNCO2KP(T) + MJ(JH2O) * FUNH2OKP(T)
c	skapa = 0.0
	IF(SKAPA . LT. 0.) THEN
	WRITE(4,*) "Kapa is a negative number"
	write(4,222) MJ(JH2O),MJ(JCO2), T
222     format("*************",3(f13.6, 3x))
c	STOP
	ENDIF
	RETURN
	END

c	/////////////////////////////////
c	//find grid number JB at RB(burner surface)
c	//given: XX,JJ array of radius, number of grids
c	//       RB, burner surface radius
c	//return: JB, the grid number of RB
c	//////////////////////////////////
	SUBROUTINE FSJB(XX, RB, JJ, JB)
	DOUBLE PRECISION RB, XX(*)
	INTEGER JB,JJ,J

	DO J=1,JJ,1
	   IF(XX(J) .EQ. RB) JB=J
	ENDDO

	IF(JB .EQ. 0) THEN
	   WRITE(*,*)"ERROR:  NO GRID ON BURNER SURFACE,
     *                QUIT AND ADD GRID AT",RB
           STOP
 	ENDIF

	RETURN
	END

C
C         OPTICAL THICK RADIATION MODEL
C            DISCRETE ORDINATE METHOD
C            STATISTIC NARROW BAND MODEL
C
C
C----------------------------------------------------------------------
C
      SUBROUTINE TABLECALC (RK, IWN, TEMP1, TABLEWN, TABLETEMP, TABLE, 
     1           NROW, NCOL)
C
C             CALCULATE THE VALUE GIVEN THE WN AND TEMPERATURE 
C             BY INTERPOLATING THE TABLE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C       RK ----- THE VARIABLE TO BE CALCULATED
C       IWN ---- INDEX OF BAND
C       TEMP --- TEMPERATURE AT WHICH THE RK TO BE VALUED
C       TABLEWN -- COLUMN IN THE TABLE STORE THE WAVE NUMBERS
C       TABLETEMP -- ROWS IN THE TABLE STORE THE TEMPERATURE
C       TABLE   ---- TABLE VALUES
C       NROW   ---- NUMBER OF ROWS
C       NCOL  ----- NUMBER OF COLUMNS
C
      DIMENSION TABLEWN(*), TABLETEMP(*), TABLE(NROW,*)

      IF(TEMP1 .LT. TABLETEMP(1)) THEN
	   TEMP = TABLETEMP(1)
	ELSE IF (TEMP1 .GE. TABLETEMP(NCOL)) THEN
	   TEMP = TABLETEMP(NCOL)
	ELSE
	   TEMP = TEMP1
	END IF

      I = 2

100   CONTINUE

	IF (TEMP .EQ. TABLETEMP(I)) THEN
	    RK = TABLE(IWN,I)
	    GOTO 200
	ELSEIF (TEMP .LT. TABLETEMP(I)) THEN
	   RK = (TEMP - TABLETEMP(I-1))/(TABLETEMP(I)-TABLETEMP(I-1))
     1         *(TABLE(IWN,I)-TABLE(IWN,I-1))+TABLE(IWN,I-1)
	   GOTO 200
	ELSE 
	   
! LECOUSTRE MODIFICATION
	   IF (I.GT.NCOL) THEN 
		RK = TABLE(IWN,NCOL)
		GOTO 200
	   ENDIF
! END OF MODIF
           I = I+1	
	   GOTO 100
	END IF

200   CONTINUE
      

	RETURN
	END 
C
C----------------------------------------------------------------------
C
      SUBROUTINE READTABLE (LTABLE, TABLETEMP, TABLEWN, TABLE, 
     1           NROW, NCOL, LOUT)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION TABLEWN(NROW), TABLETEMP(NCOL), TABLE(NROW,NCOL)

      DIMENSION VALUE(100)
	CHARACTER LINE*1000
	LOGICAL IERR

C        CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
C        KERR = KERR.OR.IERR
C        IFXTM = INT(VALUE(1))
C        TFXTM = VALUE(2)

      IERR = .FALSE.
      LINE = ' '
      READ(LTABLE, 1000) LINE
C      READ(LTABLE,1000) LINE
C	CALL CKXNUM(LINE,NCOL,LOUT,NVAL,VALUE,IERR)
      READ(LTABLE, 2000) (TABLETEMP(I), I = 1, NCOL)

C	DO I = 1,NCOL
C	   TABLETEMP(I) = VALUE(I)
C	END DO
      
10    CONTINUE
      
	DO I = 1, NROW
C	   LINE = ''
C	   READ(LTABLE, 1000) LINE
C	   CALL CKXNUM(LINE,NCOL+1,LOUT,NVAL,VALUE,IERR)
C	   TABLEWN(I) = VALUE(1)
C	   DO J = 1, NCOL	
C	      TABLE(I,J) = VALUE(J+1)
C	   END DO
         READ(LTABLE, 3000) ITABLEWN, (TABLE(I,J), J = 1, NCOL)
	   TABLEWN(I) = ITABLEWN
	END DO

      CLOSE(LTABLE)

1000  FORMAT(A)
2000  format(7x, 1p20e10.3)
3000  format(i5, 2x, 1p20e10.3)

	RETURN
	END        

C
C----------------------------------------------------------------------
C
      SUBROUTINE TRANSMIT_K (RK, TABLEWN, TABLETEMP, NROW,     
     1           NCOL, TABLE, TEMP, XCON,IBAND, P, NX) 
         
C      
C
C      CALCULATE K AT A GIVEN WAVE NUMBER FOR A GIVEN SPECIES
C
C         RK ---   ABSOBPTION COEFFICIENT
C         TABLEWN --
C         TABLETEMP --
C         TABLEK ---
C         NROW ---
C         NCOL ---
C         IBAND
C         TEMP ---
C         XCON ---
C         P ---
C         NX ---
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION TABLEWN(NROW), TABLETEMP(NCOL), TABLE(NROW,NCOL),
     1          RK(NX), TEMP(NX), XCON(NX)
C
C     SUBROUTINE TABLECALC (RK, IWN, TEMP, TABLEWN, TABLETEMP, TABLE, 
C     1           NROW, NCOL)
C
      AV = 6.0221367E23 !AVOGARDO'S Number
      TS = 296.0
      PS = 1.01325E6

      DO I = 1, NX-1
	   SZTEMP = 0.5*(TEMP(I)+TEMP(I+1))
	   CALL TABLECALC(RK(I),IBAND, SZTEMP,TABLEWN, TABLETEMP, 
     1                  TABLE, NROW,NCOL)
C
C            
         SZXCON = 0.5*(XCON(I)+XCON(I+1)) 
         RK(I) = AV*SZXCON*RK(I)
      END DO         

      RETURN
      END
     
C
C----------------------------------------------------------------------
C
      SUBROUTINE TRANSMSOOT(RKSOOT, TEMP, SOOTM1, NX, WN)

C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
!	WN : Wave number = 1/Lambda

	INCLUDE 'soot.h'

	DOUBLE PRECISION n, k, Kext
	DIMENSION TEMP(NX), SOOTM1(NX), RKSOOT(NX)
!       EXTERNAL sootVolumeFraction

!	DO I = 1, NX-1

!	   SZTEMP = 0.5*(TEMP(I)+TEMP(I+1))
!	   SM1TEMP = 0.5*(SOOTM1(I)+(SOOTM1(I+1)))

!	   FV = sootVolumeFraction(SM1TEMP)

!         RKSOOT(I) = 1600.*FV*SZTEMP/(9300.0)
!	   WRITE(*,*) RKSOOT(I), FV
!	   IF (RKSOOT(I).GT.1.0) RKSOOT(I) = 1.0

!       END DO         

! COMPUTATION OF THE EXTINCTION COEFFICIENT USING RAYLEIGH THEORY:
! ASSUMPTION: SPHERICAL PARTICLES, WITH DIAMETER SMALLER THAN CONSIDERED WAVELENGHT.
! Kext = -(PI/2)*(2*PI/LAMBDA)*(CD1)**3*M1*Im{(m**2-1)/(m**2+2)}
! WITH m = n-i*k (COMPLEX REFRACTIVE INDEX OF SOOT PARTICLE), WITH i**2=-1
! LAMBDA = WAVELENGHT OF INCIDENT RAYLIGHT.
! COMPLEXE REFRACTIVE INDEX (m = n-i*k):
!                           m = 1.57-0.56*i , from (W.H. Dalzell, A.F. Sarofim: J.Heat Transfer 91, 100 (1969))
!			    m = 1.98 - 0.93*i, from (D.W. Mackowski: Combust. Sci. and Tech. 53, 399 (1987))
!                           m = 2.10(+-0.12) - 0.48*i(+-0.06), from (P. Van-Hulle, M. Talbaut, M. Weill and A. Coppall: Meas. Sci. Techno 13 (2002))
!                           USE WAVELENGHT REFRACTIVE INDEX FROM: Chang and Charalampopoulos, TT Determination of the wavelenght dependence of the refractive indices of flame soot, Proc. of the Royal Society (London) Ser. A 1990 430:577-591
!                           n = 1.8110+0.1263*ln(Lambda)+0.027*ln^2(Lambda)+0.0417*ln^3(Lambda)
!                           k = 0.5821+0.1213*ln(Lambda)+0.2309*ln^2(Lambda)-0.01*ln^3(Lambda) (Lambda in micrometer)
! NOTE : LAMBDA (micrometer) = 10^4/k, with k the inverse of lambda, in cm^-1

! NOTE : -Im{(m**2-1)/(m**2+2)} = (6*n*k)/((n**2+k**2)**2+4*(n**2-k**2)+4)

! NOTE USE UNIT USED BY THE CODE. REMOVE FROM COMMENTS
! DEFINE REFRACTIVE INDEX
        n = 1.57
        k = 0.56
! 
!        n = 1.8110-0.1263*LOG(WN*1.D-4)+0.027*(LOG(WN*1.D-4))**2
!     1  -0.0417*(LOG(WN*1.D-4))**3
!        k = 0.5821-0.1213*(LOG(WN*1.D-4))+0.2309*(LOG(WN*1.D-4))**2
!     1  +0.01*(LOG(WN*1.D-4))**3 

        Kext = (PI/2.)*(2.*PI*WN)*(CD1)**3
        Kext = Kext*(6.*n*k)/((n**2+k**2)**2+4.*(n**2-k**2)+4.)
! 
        DO I = 1, NX-1
! 
 	SM1TEMP = 0.5*(SOOTM1(I)+(SOOTM1(I+1)))
!  
         RKSOOT(I) = Kext*SM1TEMP
! 
        END DO         




       RETURN
      END


C
C----------------------------------------------------------------------
C

      SUBROUTINE GEOMETRY (ALFA, X, NX, NORD, RMIU, WT)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C      X -----  SPATIAL GRIDS
C      NX ----  TOTAL NUMBER OF GRIDS
C      NORD --- NUMBER OF DIRECTION COSIN
C      ALFA ----   DIRECTIONAL TRANSFER COEFFICIENT
C
C
	DIMENSION ALFA(NX,NORD+1), X(NX), RMIU(NORD+1), WT(NORD+1)

      DO J = 1, NX-1
	   ALFA(J,1) = 0
	   XMID = X(J) + X(J+1)
	   CALL AREA(X(J), AREAJ)
	   CALL AREA(X(J+1), AREAJP1)
         DO M = 2, NORD+1
            ALFA(J,M) = ALFA(J,M-1) - WT(M)*RMIU(M)*(AREAJP1-AREAJ)
	   END DO
	END DO

	RETURN
	END
C
C----------------------------------------------------------------------
C
      FUNCTION PLANK(WN, T)
C
C          PLANK FUNCTION
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single	
C
C       WN --- WAVE NUMBER
C       T ---- TEMPERATURE
C	   
      C1 = 3.740E-5
      C2 = 1.4387
      PLANK = C1*WN**3/(EXP(C2*WN/T)-1)

	RETURN 
	END
C
C----------------------------------------------------------------------
C
C
      SUBROUTINE DISORD (NORD, Q, TEMP, P, XFRAC, TABLEWN, 
     1           TABLETEMP, NROW, NCOL, XH2O, XCO2, 
     2           XCO, NX, X, ALFA, RKH2O, RKCO2, RKCO, TABLEH2O, 
     3           TABLECO2, TABLECO, LTHICK,SOOTM1)
C
C          CALCULATE THE RADIATION HEAT FLUX Q FOR A GIVEN RADIATING
C          GAS SPECIES, USUALLY H20, CO2 AND CO
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C      NBAND       TOTAL NUMBER OF NARROW BANDS
C      NORD        TOTAL NUMBER OF ORDINATE
C
C
      EXTERNAL PLANK

      DIMENSION Q(NX), TEMP(NX), RKH2O(NX), ALFA(NX,NORD+1), X(NX),  
     1          XFRAC(NX), TABLEWN(NROW), TABLETEMP(NCOL), 
     2          TABLEH2O(NROW,NCOL), XH2O(NX),  XCO2(NX), XCO(NX), 
     3          RKCO2(NX), RKCO(NX), TABLECO2(NROW,NCOL),
     4          TABLECO(NROW, NCOL), RKSOOT(NX),SOOTM1(NX)

	DIMENSION RI(NORD+1,NX), RMIU(NORD+1), WT(NORD+1), RIBAND(NX,400),
     1          RIMIUM(NORD+1+1,NX), RIBJP(NX), RIPLUS(NORD+1,NX)

	LOGICAL   LTHICK
   

      IF (LTHICK) THEN 
	    ITHICK = 1
	ELSE 
	    ITHICK = 0
	END IF

      PI = 3.1415926 
      CALL GETORD(NORD,RMIU,WT)
C
C      SUBROUTINE GEOMETRY (ALFA, X, NX, NORD, RMIU, WT)
C
      CALL GEOMETRY(ALFA, X, NX, NORD, RMIU, WT)
C
      DWN = TABLEWN(2) - TABLEWN(1)
      DO L = 1, NROW 
	 WN = TABLEWN(L)
C
C      SUBROUTINE TRANSMIT_K (RK, TABLEWN, TABLETEMP, NROW,     
C     1           NCOL, TABLEK, TEMP, XCON,IBAND, P, NX) 
C         

        CALL TRANSMIT_K(RKH2O, TABLEWN, TABLETEMP, NROW,     
     1           NCOL, TABLEH2O, TEMP, XH2O, L, P, NX)          

        CALL TRANSMIT_K(RKCO2, TABLEWN, TABLETEMP, NROW,     
     1           NCOL, TABLECO2, TEMP, XCO2, L, P, NX) 

        CALL TRANSMIT_K(RKCO, TABLEWN, TABLETEMP, NROW,     
     1           NCOL, TABLECO, TEMP, XCO, L, P, NX)

	 CALL TRANSMSOOT(RKSOOT, TEMP, SOOTM1, NX, WN)

	                   
C
C             STARTING FROM RMIU= -1.0, NEGATIVE DIRECTION
C
        DO J = 1, NX-1
	      RIBJP(J) = PLANK(WN,(TEMP(J)+TEMP(J+1))/2)/PI
	 END DO

        M = 1
	 RI(M,NX) = 0.6*PLANK(WN,TEMP(NX))/PI
C       RI(M,NX) = 0
	 DO J = NX-1, 1, -1
		XMID = (X(J) + X(J+1))/2

		CALL AREA(X(J+1), AERAP)
		CALL AREA(X(J), AERAJ)

		AMIUP = ABS(RMIU(M))* (AERAP + AERAJ)
C	       VOLJ = (X(J+1) - X(J)) * AERAJ

		VOLJ = 4./3.*PI*(X(J+1)**3-X(J)**3)
		VOLKJ = VOLJ * (RKH2O(J)+RKCO2(J)+RKCO(J)+RKSOOT(J))

		RIPLUS(M,J) = (AMIUP * RI(M,J+1) + VOLKJ*RIBJP(J))
     1        /(AMIUP + VOLKJ*ITHICK)

		RI(M,J) = 2*RIPLUS(M,J) - RI(M,J+1)
		RIMIUM(M+1,J) = RIPLUS(M,J)

	 END DO
C
C                 RMIU < 0 NEGATIVE DIRECTION
C
	 DO M = 2, (NORD)/2+1
		RI(M,NX) = 0.6*PLANK(WN,TEMP(NX))/PI
	      	DO  J = NX-1,1,-1
      	     		XMID = (X(J) + X(J+1))/2
               	CALL AREA(X(J+1), AERAP)
	         	CALL AREA(X(J), AERAJ)
	         	AMIUP = ABS(RMIU(M))* (AERAP+AERAJ)
C	         VOLJ = (X(J+1) - X(J)) * AERAJ
               	VOLJ = 4./3.*PI*(X(J+1)**3-X(J)**3)
	         	VOLKJ = VOLJ * (RKH2O(J)+RKCO2(J)+RKCO(J)+RKSOOT(J))
	         	ALFAWT = (ALFA(J,M-1)+ALFA(J,M))/WT(M)
	         	RIPLUS(M,J) = (AMIUP * RI(M,J+1)+ALFAWT*RIMIUM(M,J) 
     1          	+ VOLKJ*RIBJP(J))/(AMIUP + ALFAWT + VOLKJ*ITHICK)
               	RIMIUM(M+1,J) = 2 * RIPLUS(M,J) - RIMIUM(M,J)
	         	RI(M,J) = 2 * RIPLUS(M,J) - RI(M,J+1)
	      	END DO
	 END DO
C
C                RMIU > 0 POSITIVE DIRECTION
C
	 DO M = (NORD)/2+2, NORD+1
		RI(M,1) = 0.6 * PLANK(WN,TEMP(1))/PI

	      	DO  J = 1, NX-1
      	    	 XMID = (X(J) + X(J+1))/2
               	 CALL AREA(X(J), AERAJ)
	         CALL AREA(X(J+1), AERAP)
	         AMIUJ = ABS(RMIU(M))* (AERAJ+AERAP)
C	         VOLJ = (X(J+1) - X(J)) * AERAP
               	 VOLJ = 4./3.*PI*(X(J+1)**3-X(J)**3)
	         VOLKJ = VOLJ * (RKH2O(J)+RKCO2(J)+RKCO(J)+RKSOOT(J))
               	 ALFAWT = (ALFA(J,M-1)+ALFA(J,M))/WT(M)
	         RIPLUS(M,J) = (AMIUJ * RI(M,J) + ALFAWT*RIMIUM(M,J) 
     1          	+ VOLKJ*RIBJP(J))/(AMIUJ + ALFAWT + VOLKJ*ITHICK)
               	 RIMIUM(M+1,J) = 2 * RIPLUS(M,J) - RIMIUM(M,J)
	         RI(M,J+1) = 2 * RIPLUS(M,J) - RI(M,J)

	      	END DO
	 END DO
C
C          INTEGRATE ON ALL ORDINATES
C       
         DO I = 1, NX-1
            RIBAND(I,L) = 0
            DO N = 1, NORD+1         
	         RIBAND(I,L) = RIBAND(I,L) + RMIU(N)*WT(N)*RIPLUS(N,I)
	    END DO
         END DO
C
C         END OF THE LOOP OF BANDS
C      
      END DO		   

C
C         SUMMATION OF ALL THE WAVE LENGTH
C
      DO I = 1, NX-1			   
        Q(I) = 0
        DO L = 1, NROW
		Q(I) = Q(I) + RIBAND(I,L)*DWN
	END DO
      END DO  
	
      RETURN 
      END   

C
C----------------------------------------------------------------------
C
C
      SUBROUTINE GETORD (NORD, RMIU, WT)
C
C         RETURN THE MIU AND WEITHT GIVEN THE ORDER OF ORDINATE 
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C      NORD ----- NUMBER OF DISCREATE BEAMS
C      RMIU ----- MIU
C      WT ------- WEIGHT
C
C

	DIMENSION RMIU(NORD+1), WT(NORD+1)

	PI = 3.1415926

        IF (NORD .EQ.4) THEN

        RMIU(1) = -1.0
	RMIU(2) = -0.8819171
	RMIU(3) = -0.3333333
	RMIU(4) = 0.3333333
	RMIU(5) = 0.8819171


	WT(1) = 0 
	WT(2) = 0.1666667*4*PI 
        WT(3) = 0.3333333*4*PI 
        WT(4) = 0.3333333*4*PI 
	WT(5) = 0.1666667*4*PI 	

	ELSE IF (NORD .EQ. 20) THEN

        RMIU(1) = -1.0

        RMIU(2)  = -0.99313
        RMIU(3)  = -0.96397
        RMIU(4)  = -0.91223        
        RMIU(5)  = -0.83912
        RMIU(6)  = -0.74633
        RMIU(7)  = -0.63605 
        RMIU(8)  = -0.51087       
        RMIU(9)  = -0.37371
        RMIU(10) = -0.22779 
        RMIU(11) = -0.07653

        RMIU(21) = 0.99313
        RMIU(20) = 0.96397
        RMIU(19) = 0.91223        
        RMIU(18) = 0.83912
        RMIU(17) = 0.74633
        RMIU(16) = 0.63605 
        RMIU(15) = 0.51087       
        RMIU(14) = 0.37371
        RMIU(13) = 0.22779 
        RMIU(12) = 0.07653

        WT(1)  = 0 
	WT(2)  = 0.15275*2*PI 
	WT(3)  = 0.14917*2*PI 
        WT(4)  = 0.14210*2*PI 
        WT(5)  = 0.13169*2*PI 
	WT(6)  = 0.11819*2*PI 
	WT(7)  = 0.10193*2*PI 
	WT(8)  = 0.08328*2*PI 
        WT(9)  = 0.06267*2*PI 
        WT(10) = 0.04060*2*PI 
	WT(11) = 0.01761*2*PI 

	WT(21) = 0.15275*2*PI 
	WT(20) = 0.14917*2*PI 
        WT(19) = 0.14210*2*PI 
        WT(18) = 0.13169*2*PI
	WT(17) = 0.11819*2*PI 
	WT(16) = 0.10193*2*PI 
	WT(15) = 0.08328*2*PI 
        WT(14) = 0.06267*2*PI 
        WT(13) = 0.04060*2*PI 
	WT(12) = 0.01761*2*PI 

      END IF
		
      RETURN
      END  





