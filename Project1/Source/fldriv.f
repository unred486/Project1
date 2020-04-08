C SUBROUTINE FLDRIV 
C CORE FUNCTION. INITIALIZED DATA, CALLS DIFFERENT METHODS
C SECOND AND LAST FUNCTION CALLED BY SPHDIFB
C
      SUBROUTINE FLDRIV (LENTWP, LIN, LOUT, LREST, LSAVE, LRCRVR, 
     1                   KK, II, NATJ, NMAX, RCKWRK, 
     2                   RMCWRK, EPS, WT, REAC, SCRCHK, X, COND, XXD,
     3                   TGIVEN, XGIVEN, D, DKJ, TDR, YV, SCAL, ABOVE, 
     4                   BELOW, BUFFER, TWPWK, S, SN, F, FN, DS, SSAVE, 
     5                   A6, A, ICKWRK, IMCWRK, KSYM, CCKWRK,
     6                   KR, KI, KP, IP, ACTIVE, MARK, PROD, EPA, LRATE,
     7                   ALFA, RSPACE, TABLEWN, TABLETEMP, TABLE,
     8                   NMXROW, NMXCOL, NORD, LH2O,LCO2, LCO,XNEW,XINT,
     9                   SS,RW,ZGIVEN,DMIX,DGIVEN,XXT,IW,INFO,LTESTT)

c      USE f90_module
      USE var

C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
      INCLUDE 'TEST.cmn'
      INCLUDE 'LOGICAL.cmn'

C
      DIMENSION ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*),RW(*),IW(*)
      DIMENSION EPS(*),  REAC(*), X(*),D(KK,*),COND(*),
     1          TDR(KK,*), YV(KK,*), SCAL(*), S(NATJ,*), SN(NATJ,*),
     2          F(NATJ,*), XGIVEN(*), TGIVEN(*), SCRCHK(KK,*), KR(*),
     3          KI(*), KP(*), DKJ(KK,KK,*), A6(*), PROD(*), EPA(*),
     4          RHO(NMAX), XNEW(*), XINT(*),SS(NATJ,*),TJ(NMAX),INFO(*),
     5          ATOL(1),RTOL(1),ZGIVEN(*),DMIX(*),XXT(*),XXD(*),
     6          DGIVEN(*),WT(*)

C
      CHARACTER*(*) CCKWRK(*), KSYM(*)
      CHARACTER*16 ICHR, ISOLUT, ISENSI, ICKLNK, IMCLNK, IHSENS
      CHARACTER*16 KNAME(KK)
      CHARACTER*12 PAH
      character(len=32) :: my_fmt
C
C         DIMENSION NEWTON SPACE
C
      DIMENSION ABOVE(NATJ,*), BELOW(NATJ,*), BUFFER(NATJ,*), TWPWK(*),
     1          A(6 * NATJ - 2, *), IP(NATJ,*), FN(NATJ,*), DS(*),
     2          SSAVE(*), LEVEL(2)
C
      DIMENSION ALFA(NMAX,NORD+1), RSPACE(NMAX,10), 
     1          TABLEWN(NMXROW), TABLETEMP(NMXCOL), 
     2          TABLE(NMXROW,NMXCOL,3)

      LOGICAL LBURNR, LTIME, LTIME2, LMOLE, LENRGY, LTDIF,
     1        LRSTRT, LCNTUE, LASEN, LRSEN, LESEN, LMSEN, LPSEN, LMULTI,
     2        LVARMC, RSTCNT, ERROR, FUNCTN, JACOBN, REENTR, SOLVE,
     3        STORE, SUCCES, ADAPT, SHOW, SAVE, UPDATE, ACTIVE(*),
     4        MARK(*), ENERGY, BURNER, LUSTGV, IERR, LVCOR, LHSEN,
     5        LADAPT, LSTEADY, ADIAMF, SOOTP, freeMolOnly,LPOST,LREF,
     6        LGLOB
      
C
      INTEGER CALL, CALLS
C
      DATA ADAFLR /1.0E-08/, LCNTUE/.FALSE./
      DATA ISOLUT/'SOLUTION        '/, ISENSI/'SENSITIVITY     '/,
     1     ICKLNK/'CKLINK          '/, IMCLNK/'MCLINK          '/,
     2     IHSENS/'HSENSITIVITY    '/

      DOUBLE PRECISION M0_cutoff
C
200   CONTINUE ! NEW PROBLEM

C------------------------------------------------------------------------
C         CALL THE KEYWORD INPUT
C------------------------------------------------------------------------

      IF (.NOT.LPOST) THEN
      CALL RDKEY      (KK, NMAX, LIN, LOUT, KSYM, PATM,
     1                  MFILE, IPRNT, UFAC, DTMAX, DTMIN,
     2                   JJ, NTOT, NADP, X, SCAL, NREAC, NINTM,
     3                  NPROD, REAC, PROD, KR, KI, KP, XSTR,
     4                  XCEN, XEND, WMIX, FLRT, GRAD, CURV, SFLR,
     5                  NTEMP, XGIVEN, TGIVEN, TFIXT, ATIM, RTIM,
     6                  NJAC, ITJAC, NINIT, NUMDT, IRETIR, DT1,
     7                  NUMDT2, DT2, WNDFAC, GFAC, SPOS, N1CALL, TBUR,
     8                  TAMB, IRADI,DVCF,RADCOEFF,P,
     9                  PAH,coeff_K, IREFN, RFUEL,HC,CP,RK,PHI,INFO,
     %                  ZGIVEN,XXT,XXD,DGIVEN,ATOL,RTOL,NZMIX,NDMIX,ZST)
      
      ENDIF
c//////////// for testing purposes delete when done updating///////////
c//////////////////////////////////////////////////////////////////////
c       WRITE(LOUT, 9999) LBURNR, LTIME, LTIME2, LMOLE, LENRGY, LTDIF, 
c     1        LRSTRT, LCNTUE, LUSTGV, LUMESH,LMULTI, LVCOR, 
c     3        LADAPT, ADIAMF, SOOTP, freeMolOnly,LPREMIX, LCARTE, 
c     4        LHYDRO, LTEST, LSOOTINIT, LNUCL, LSURF,LOXID, 
c     5         LCOAG,LPOST,LREF,LGLOB,LHOM
     
c9999  format('LBURNR = ', 1L5, 2X,'LTIME = ',1L5,2X,'LTIME2 = ',1L5,2X,
c     1       'LMOLE = ',1L5, 2X,'LENRGY = ', 1L5,2X,'LTDIF = ',1L5,2X, 
c     2       'LRSTRT = ',1L5, 2X,'LCNTUE = ',1L5, 2X,'LUSTGV = ',1L5,2X,
c     3       'LUMESH = ',1L5,2X,'LMULTI = ',1L5,2X,'LVCOR = ',1L5,2X,
c     4       'LADPT = ',1L5,2X,'LADIAMF = ',1L5,2X,
c     5       'LSOOTP = ',1L5,2X,'freemol = ',1L5,2X,
c     6       'LPREMIX = ',1L5,2X,'LCARTE = ',1L5,2X,'LHYDRO = ',1L5,2X,
c     7       'LTEST = ',1L5,2X,'LSOOTINIT = ',1L5,2X,'LNUCL = ',1L5,2X,
c     8       'LSURF = ',1L5,2X,'LOXID = ',1L5,2X,'LCOAG = ',1L5,2X,
c     9       'LPOST = ',1L5,2X,'LREF = ',1L5,2X,'LGLOB = ',1L5,2X,
c     &       'LHOM = ',1L5,2X)      
c/////////////////////////////////////////////////////////////////////

        CALL START (KK, NMAX, NATJ, JJ, LOUT,
     1               NREAC, NINTM, NPROD, REAC, SCRCHK(1,2),
     2               PROD, KR, KI, KP, FLRT, NTEMP, XGIVEN,
     3               TGIVEN, XSTR, XCEN, XEND, WMIX, ICKWRK, RCKWRK,
     4               SCRCHK(1,4), COND,EPS, JFIXT, TFIXT, 
     5		         X, S, EPA, TBUR,TAMB,KSYM,HC,CP,RK,PHI,
     6               P,ZGIVEN,NZMIX,NDMIX,XXT,XXD,DGIVEN,DMIX,ZST)


c////////////////////////////////for testing////////////      
c      WRITE(my_fmt , '(I0)')  KK+2
c      my_fmt = '('//trim(my_fmt)//'(1PE12.4))'        
c      WRITE(LOUT,9997) (trim(KSYM(K)),K=1,KK)
c      DO J=1,JJ
c      WRITE(LOUT,my_fmt) (S(NYS+K,J),K=1,KK), S(NT,J), S(NM,J)
c      END DO
c9997  FORMAT(<KK> (' X'A''), ' Temperature MassFlowRate') 
c      WRITE(LOUT,5000)
c5000  FORMAT('XX ZMIX')
c5010  FORMAT (2 (1PE12.4))
c      DO J=1,JJ
c          WRITE(LOUT, 5010) X(J), S(NZ,J)
c      END DO


c//////////////////////////////////////////////////////////////

c call steady state calculation
c  --> using TWOPNT
c    --> move all TWOPNT instruction variables to this subroutine
c    --> pick out the variables to move into this subroutine
c    --> 
c call subroutines for closed constant volue constant mass homogeneous reactor calculation
c  --> using DASPK
c  --> using 4th order RUNGE KUTTA  

c calling a subroutine that initializes/call daspk solver for cvcm_reactor
C direct solver

      if (LHOM) THEN
        call cvcm(NATJ,JJ,LOUT,KK,ICKWRK,RCKWRK,RW,IW,SCRCHK(1,2),
     1           S,INFO,SCRCHK(1,3),SCRCHK(1,4),KSYM,P,WT,DT1,RTOL,
     2           ATOL)
      END IF
      
      if (LZMIX) THEN 
          call zmf(NATJ,JJ,LOUT,NMAX,X,RW,IW,S,INFO,DT1,RTOL,ATOL,DMIX,
     1     XST,ZST,XSTR,XEND)       
      end if 
      
      END SUBROUTINE
