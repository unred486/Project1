C///////////////////////////////////////////////////////////////////
C        NEW ONE DIMENSIONAL LAMINAR SPHERICAL DIFFUSION FLAME(NSphD)
C           (MODIFIED FROM SphDiff)
C        WRITTEN BY:
C            Han Ju Lee
C            University of Maryland
C            Department of Aerospace Engineering
C            hanjulee@umd.edu
C
C///////////////////////////////////////////////////////////////////
C       Update Progress
C       1.1 Modifying Sphdiff code to test homogeneous reactor
C           
C///////////////////////////////////////////////////////////////////
C
C       To do list
C       1. organize common list for an easier use
C       2. 

      SUBROUTINE SPHDIFB (NMAX, LIN, LOUT, LINKCK, LINKMC, LREST, LSAVE,
     1                   LRCRVR, LENLWK, L, LENIWK, I, LENRWK,R,
     2                   LENCWK, C, LRATE, NMXROW, NMXCOL, NORD, LH2O,
     3                   LCO2, LCO,LTESTT)
       USE var     
c       use varray
        
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C


      DIMENSION I(*),R(*)
      LOGICAL   L(*)
      CHARACTER C(*)*(*)
      
      COMMON /SP/ KK,NATJ
C
c      INCLUDE 'FLFLFL.cmn'
c      INCLUDE 'LOCS.cmn'
c      INCLUDE 'INDEXSP.cmn'
c      INCLUDE 'OPTHICK.cmn'
C
C
C          WRITE VERSION NUMBER
C
      WRITE (LOUT,15)
   15 FORMAT (
     1/' NSphDiff:  One-dimensional diffusion laminar flame code',
     2/'          Version 1.0, Nov 2019',
C***** precision > double
     3/'          DOUBLE PRECISION')
C*****END precision > double
C
C 
  


C          SET UP INTERNAL WORK POINTERS
C
C          Initial point subroutine used from previous version
      

      CALL POINT (LINKCK, LINKMC, LENIWK, LENRWK, LENCWK, NMAX, LOUT,
     1            LSAVE, KK, II, NATJ, LENTWP, LTOT, ITOT, NTOT,
     2            ICTOT, I, R, C, NMXROW, NMXCOL, NORD,LIN)


C
C           CHECK FOR ENOUGH SPACE
C
      WRITE (LOUT, 7000) LENLWK, LTOT, LENIWK, ITOT, LENRWK, NTOT,
     1                   LENCWK, ICTOT
7000  FORMAT (/,'                WORKING SPACE REQUIREMENTS',
     1        /,'                 PROVIDED        REQUIRED ',
     2        /,' LOGICAL  ' , 2I15,
     3        /,' INTEGER  ' , 2I15,
     4        /,' REAL     ' , 2I15,
     5        /,' CHARACTER' , 2I15,/)

      
      
C
      IF (LTOT.GT.LENLWK .OR. ITOT.GT.LENIWK .OR. NTOT.GT.LENRWK
     1                   .OR. ICTOT.GT.LENCWK) THEN
         WRITE (LOUT, *) '  FATAL ERROR, NOT ENOUGH WORK SPACE PROVIDED'
         STOP
      ENDIF
         
         
         
      
C
C
      CALL FLDRIV (LENTWP, LIN, LOUT, LREST, LSAVE, LRCRVR, KK, II, 
     1           NATJ, NMAX, R(NCKW), R(NMCW), R(NEPS), R(NWT),
     2           R(NRE), R(NSCH), R(NX),R(NCON),R(NVISM), R(NTGV),
     3           R(NXGV),R(ND), R(NDKJ), R(NTDR), R(NYV), R(NSCL), 
     4           R(NABV),R(NBLW), R(NBUF), R(NTWP), R(NS), R(NSN), 
     5           R(NF), R(NFN), R(NDS), R(NSSV), R(NKA6), R(NA), 
     6           I(ICKW), I(IMCW), C(IKS), C(ICC), I(IKR), I(IKI),
     7           I(IKP), I(IIP), L(LAC),  L(LMK), R(NPD), R(NEPA),
     8           LRATE, R(NALFA), R(NSPACE), R(NTBWN), R(NTBTEMP),
     9           R(NTABLE), NMXROW, NMXCOL, NORD, LH2O, LCO2, LCO,
     &           R(NXNEW),R(NXINT),R(NSS),R(LRW),R(NZM), R(NDMIXP),
     $           R(NDGIVEN),R(NXXT),I(LIW),I(IINFO),LTESTT)
C
      RETURN
      END

