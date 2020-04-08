      SUBROUTINE PRSAVE (ICKWRK, RCKWRK, CCKWRK, 
     1                   IMCWRK, RMCWRK, LOUT, LSAVE)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*)
      CHARACTER*(*) CCKWRK(*)
      CHARACTER*16 ILINK
C
      ILINK = 'CKLINK          '
      WRITE (LSAVE) ILINK
      CALL CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
C
      ILINK = 'MCLINK          '
      WRITE (LSAVE) ILINK
      CALL MCSAVE (LOUT, LSAVE, IMCWRK, RMCWRK)
C
      RETURN
      END
C