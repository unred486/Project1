C SUBROUTINE POINT 
C DEFINE THE POINTERS USED FOR THE RESIDUAL VECTOR, STORAGE OF SOLUTION, INDEX OF CHEMISTRY AND TRANSPORT PROPERTIES
C FIRST FUNCTION CALLED BY SPHDIFB
C
! determines the size of arrays used by other subroutines
! Points to the location of first element in an array
      SUBROUTINE POINT (LINKCK, LINKMC, LENIWK, LENRWK, LENCWK, NMAX,
     1                  LOUT, LSAVE, KK, II, NATJ, LENTWP, LTOT, ITOT,
     2                  NTOT, ICTOT, I, R, C, NMXROW, NMXCOL, NORD,LIN)
      USE ini
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
      PARAMETER (LENPFAC = 10)
      PARAMETER (LENPLUFAC = 2)
      PARAMETER (IPREMETH = 2)  ! =1 means ILUT preconditioner used
                                ! =2 means ILUTP preconditioner used
      PARAMETER (LFILILUT = 10)
      PARAMETER (IREORDER = 1)
      PARAMETER (ISRNORM = 1)
      PARAMETER (NORMTYPE = 2)
      PARAMETER (JACOUT = 0)
      PARAMETER (JSCALCOL = 1)
      PARAMETER (TOLILUT = 0.001)
      PARAMETER (PERMTOL = 0.01)
      
      DIMENSION I(*), R(*)
c      INTEGER, DIMENSION(:), POINTER :: ICKWW
c      DOUBLE PRECISION, DIMENSION(:), POINTER :: NCKWW
      CHARACTER C(*)*(*)
      INTEGER :: NMAX
C
c      INCLUDE 'FLFLFL.cmn'
C
c      INCLUDE 'LOCS.cmn'
C
c      INCLUDE 'OPTHICK.cmn'

C          SET THE POINTERS INTO THE SOLUTION VECTOR
C

      N_moments  = 0
C      NM  = KK+2
C
      CALL CKLEN (LINKCK, LOUT, LENICK, LENRCK, LENCCK)
      CALL MCLEN (LINKMC, LOUT, LENIMC, LENRMC)

C
C     real chemkin work space
      NCKW = 1
C     real transport work space
      NMCW = NCKW + LENRCK
      NTOT = NMCW + LENRMC
C     integer chemkin work space
      ICKW = 1
C     integer transport work space
      IMCW = ICKW + LENICK 
      ITOT = IMCW + LENIMC
C     character chemkin work space
      ICC = 1
      ICTOT = ICC + LENCCK
C
      IF (ITOT.LT.LENIWK .AND. NTOT.LT.LENRWK .AND. ICTOT.LT.LENCWK)
     1   THEN
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, I, R, C)
         CLOSE (LINKCK)
         CALL CKINDX (I, R, MM, KK, II, NFIT)
C
         CALL CKMXTP (I, MAXTP)
         NTR = MAXTP - 1
C
         CALL MCINIT (LINKMC, LOUT, LENIMC, LENRMC, I(IMCW), R(NMCW))
         CLOSE (LINKMC)
         REWIND LSAVE
         CALL PRSAVE (I, R, C, I(IMCW), R(NMCW), LOUT, LSAVE)
      ENDIF

C          Read INIKEY subroutine for allocating memory space for arrays
C          used in other solver subroutines  
C          NTOTG - total allowable grid01 
C          NTOTG not declared in this subroutine explicited but declared from
C          ini module
        
      CALL INIKEY(LIN,LOUT) 
      
      REWIND(LIN)
      
      IF (LHOM) THEN
      NMAX=NTOTG
C
      NSM = 0
      NT  = N_moments + NSM + 1 	
      NYS = NT
      NY  = NYS + 1
      NM  = NY + KK 
      NATJ = NM
!       NM  =  KK+2
!       NSM  = NM     ! NSM POINTER SOOT MOMENT. THE FIRST ELEMENT IS LOCATED IN NSM+1
!       NATJ = NSM+N_moments
      
      LENTWP = (7*NATJ+2)*NMAX
      
C          APPORTION THE BALANCE OF THE FLOATING POINT SPACE
C
      NEPS = NTOT
      NEPA = NEPS + KK
      NWT  = NEPA + KK
      NRE  = NWT  + KK
      NPD  = NRE  + KK
      NSCH = NPD  + KK
      NX   = NSCH + KK*6
      NCON = NX   + NMAX
      NVISM= NCON + NMAX
! REPLACE NTGV = NCON + NMAX BY NTGV = NVISM + NMAX
      NTGV = NVISM+ NMAX
      NXGV = NTGV + NMAX
      ND   = NXGV + NMAX
      NDKJ = ND   + KK*NMAX
      NTDR = NDKJ + KK*KK*NMAX
      NYV  = NTDR + KK*NMAX
      NSCL = NYV  + KK*NMAX
      NABV = NSCL + NATJ
      NBLW = NABV + NATJ*NMAX
      NBUF = NBLW + NATJ*NMAX
      NTWP = NBUF + NATJ*NMAX
      NS   = NTWP + LENTWP
      NSN  = NS   + NATJ*NMAX
      NF   = NSN  + NATJ*NMAX
      NFN  = NF   + NATJ*NMAX
      NDS  = NFN  + NATJ*NMAX
      NSSV = NDS  + NMAX
C     thermodynamic coefficients a6
      NKA6 = NSSV + NMAX
      NA   = NKA6 + NTR
      NTOT = NA   + (6*NATJ-2) * (NATJ*NMAX) - 1
C
C         APPORTION THE SPACE FOR RADIATION

      NTBTEMP =  NTOT
      NTBWN   =  NTBTEMP + NMXCOL
      NZM     =  NTBWN   + NMXROW
      NTABLE  =  NZM   + NMAX
      NALFA   =  NTABLE  + 6*NMXROW*NMXCOL
      NSPACE  =  NALFA   + NMAX*(NORD+1)
      NTOT    =  NSPACE  + 10*NMAX

! TRANSIENT REFINE GRID ARRAY POINTERS
      NXNEW = NTOT
      NXINT = NXNEW+NMAX
      NSS  = NXINT+NMAX
      NDMIXP = NSS+NMAX*NMAX
      NTOT = NDMIX+NMAX
      
C    When INFO(12) = 0, with INFO(5) = 0, INFO(6) = 1:
C    The length required for RWORK is
C    50 + (2*ML+MU+11)*NEQ + 2*(NEQ/(ML+MU+1) + 1) .
C        If INFO(12) = 0 (standard direct method), the base value is
C               base = 50 + max(MAXORD+4,7)*NEQ.
C               The default value is MAXORD = 5 (see INFO(9)).  With the
C               default MAXORD, base = 50 + 9*NEQ.
C               Additional storage must be added to the base value for
C               any or all of the following options:
C                 if INFO(6) = 0 (dense matrix), add NEQ**2
C                 if INFO(6) = 1 (banded matrix), then
C                    if INFO(5) = 0, add (2*ML+MU+1)*NEQ + 2*(NEQ/(ML+MU+1)+1),
C                    if INFO(5) = 1, add (2*ML+MU+1)*NEQ,
C                 if INFO(16) = 1, add NEQ.
C
C        When INFO(12) = 1(Krylow method), with INFO(5) = 0, INFO(6) = 1:
C        The length required for RWORK is 
C              91 + 19*NEQ + LENWP,
C     where LENWP is given by
C       LENWP =  length of RWORK segment WP = 
C              2*LENPFAC*NEQ + LENPLUFAC*NEQ + ISRNORM*NEQ + 2*(NEQ+1)
        
C         NEQ=NATJ*NMAX but NMAX for cvcm is 1
      NEQ=NATJ*NMAX

      IF (INFOT .EQ. 0) THEN
      
      LRW=NTOT
      NTOT=LRW+50+9*NEQ+NEQ+2*(NEQ+1)
      
      ELSE IF (INFOT .EQ. 1) THEN
      
      LRW=NTOT
      LENWP=2*LENPFAC*NEQ+LENPLUFAC*NEQ+ISRONORM*NEQ+2*(NEQ+1)
      NTOT=LRW+91+19*NEQ+LENWP
      
      END IF
      
c homogeneous reactor pointers for cvcm subroutine which uses daspk
      
C    The length required fo r IWORK is
C              40 + NEQ + LENIWP,
C    where LENIWP is given by
C       LENIWP = length of IWORK segment IWP =
C              4*(NEQ+1) + 3*LENPFAC*NEQ + 2*LENPLUFAC*NEQ
C                 + 2*IREORDER*NEQ + 2*(IPREMETH-1)*NEQ
      
C
C           APPORTION THE BALANCE OF THE INTEGER SPACE
C
      IKR  = ITOT
      IKI  = IKR  + KK
      IKP  = IKI  + KK
      IIP  = IKP  + KK
      ITOT = IIP  + NATJ*NMAX - 1
      LIW  = ITOT
      
      
      IF (INFOT .EQ. 0) THEN
      ITOT = LIW + 40 + 2*NEQ
      ELSE IF (INFOT .EQ. 1) THEN
      LENIWP = 4*(NEQ+1) + 3*LENPFAC*NEQ + 2*LENPLUFAC*NEQ
     1            + 2*IREORDER*NEQ + 2*(IPREMETH-1)*NEQ
      ITOT = LIW + 40+ NEQ+ LENIWP
      END IF
      IINFO=ITOT
      ITOT=IINFO+20
C
C           APPORTION THE LOGICAL SPACE
C
      LAC  = 1
      LMK  = LAC  + NATJ
      LTOT = LMK  + NMAX - 1

C           APPORTION THE BALANCE OF THE CHARACTER SPACE
C
      IKS = ICTOT
      ICTOT = IKS + KK - 1
      
      ELSE IF (LZMIX) THEN

      NMAX=NTOTG
C
      NSM = 0
      NZ  = 1
      NATJ = NZ
      
      LENTWP = (7*NATJ+2)*NMAX
      
C          APPORTION THE BALANCE OF THE FLOATING POINT SPACE
C
      NEPS = NTOT
      NEPA = NEPS + KK
      NWT  = NEPA + KK
      NRE  = NWT  + KK
      NPD  = NRE  + KK
      NSCH = NPD  + KK
      NX   = NSCH + KK*6
      NCON = NX   + NMAX
      NVISM= NCON + NMAX
! REPLACE NTGV = NCON + NMAX BY NTGV = NVISM + NMAX
      NTGV = NVISM+ NMAX
      NXGV = NTGV + NMAX
      ND   = NXGV + NMAX
      NDKJ = ND   + KK*NMAX
      NTDR = NDKJ + KK*KK*NMAX
      NYV  = NTDR + KK*NMAX
      NSCL = NYV  + KK*NMAX
      NABV = NSCL + NATJ
      NBLW = NABV + NATJ*NMAX
      NBUF = NBLW + NATJ*NMAX
      NTWP = NBUF + NATJ*NMAX
      NS   = NTWP + LENTWP
      NSN  = NS   + NATJ*NMAX
      NF   = NSN  + NATJ*NMAX
      NFN  = NF   + NATJ*NMAX
      NDS  = NFN  + NATJ*NMAX
      NSSV = NDS  + NMAX
C     thermodynamic coefficients a6
      NKA6 = NSSV + NMAX
      NA   = NKA6 + NTR
      NTOT = NA   + (6*NATJ-2) * (NATJ*NMAX) - 1
C
C         APPORTION THE SPACE FOR RADIATION

      NTBTEMP =  NTOT
      NTBWN   =  NTBTEMP + NMXCOL
      NTABLE=  NTBWN   + NMXROW
      NALFA   =  NTABLE  + 6*NMXROW*NMXCOL
      NSPACE  =  NALFA   + NMAX*(NORD+1)
      NTOT    =  NSPACE  + 10*NMAX

! TRANSIENT REFINE GRID ARRAY POINTERS
      NXNEW = NTOT
      NXINT = NXNEW+NMAX
      NSS  = NXINT+NMAX
      NDMIXP = NSS+NMAX*NMAX
      NDGIVEN= NDMIXP+NMAX
      NXXT= NDGIVEN+NMAX
      NZM = NXXT+NMAX
      NTOT = NZM+NMAX
      
C    When INFO(12) = 0, with INFO(5) = 0, INFO(6) = 1:
C    The length required for RWORK is
C    50 + (2*ML+MU+11)*NEQ + 2*(NEQ/(ML+MU+1) + 1) .
C        If INFO(12) = 0 (standard direct method), the base value is
C               base = 50 + max(MAXORD+4,7)*NEQ.
C               The default value is MAXORD = 5 (see INFO(9)).  With the
C               default MAXORD, base = 50 + 9*NEQ.
C               Additional storage must be added to the base value for
C               any or all of the following options:
C                 if INFO(6) = 0 (dense matrix), add NEQ**2
C                 if INFO(6) = 1 (banded matrix), then
C                    if INFO(5) = 0, add (2*ML+MU+1)*NEQ + 2*(NEQ/(ML+MU+1)+1),
C                    if INFO(5) = 1, add (2*ML+MU+1)*NEQ,
C                 if INFO(16) = 1, add NEQ.
C
C        When INFO(12) = 1(Krylow method), with INFO(5) = 0, INFO(6) = 1:
C        The length required for RWORK is 
C              91 + 19*NEQ + LENWP,
C     where LENWP is given by
C       LENWP =  length of RWORK segment WP = 
C              2*LENPFAC*NEQ + LENPLUFAC*NEQ + ISRNORM*NEQ + 2*(NEQ+1)
        
C         NEQ=NATJ*NMAX but NMAX for cvcm is 1
      NEQ=NATJ*NMAX

      IF (INFOT .EQ. 0) THEN
      
      LRW=NTOT
      NTOT=LRW+50+9*NEQ+NEQ+2*(NEQ+1)
      
      ELSE IF (INFOT .EQ. 1) THEN
      
      LRW=NTOT
      LENWP=5*NEQ + 2*((NEQ/3) + 1)
      NTOT=LRW+91+19*NEQ+LENWP
      
      END IF
      
c homogeneous reactor pointers for cvcm subroutine which uses daspk
      
C    The length required fo r IWORK is
C              40 + NEQ + LENIWP,
C    where LENIWP is given by
C       LENIWP = length of IWORK segment IWP =
C              4*(NEQ+1) + 3*LENPFAC*NEQ + 2*LENPLUFAC*NEQ
C                 + 2*IREORDER*NEQ + 2*(IPREMETH-1)*NEQ
      
C
C           APPORTION THE BALANCE OF THE INTEGER SPACE
C
      IKR  = ITOT
      IKI  = IKR  + KK
      IKP  = IKI  + KK
      IIP  = IKP  + KK
      ITOT = IIP  + NATJ*NMAX - 1
      LIW  = ITOT
      
      
      IF (INFOT .EQ. 0) THEN
      ITOT = LIW + 40 + 2*NEQ
      ELSE IF (INFOT .EQ. 1) THEN
      LENIWP = 4*(NEQ+1) + 3*LENPFAC*NEQ + 2*LENPLUFAC*NEQ
     1            + 2*IREORDER*NEQ + 2*(IPREMETH-1)*NEQ
      ITOT = LIW + 40+ NEQ  !+ LENIWP
      END IF
      IINFO=ITOT
      ITOT=IINFO+20
C
C           APPORTION THE LOGICAL SPACE
C
      LAC  = 1
      LMK  = LAC  + NATJ
      LTOT = LMK  + NMAX - 1

C           APPORTION THE BALANCE OF THE CHARACTER SPACE
C
      IKS = ICTOT
      ICTOT = IKS + KK - 1
      
      ENDIF
      RETURN
      END
C
C----------------------------------------------------------------------
