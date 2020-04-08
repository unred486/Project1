      SUBROUTINE PSOLCVCM(NEQ, T, CC, CCPRIME, SAVR, WK, CJ, WT, WP, 
     1                   IWP, B, EPLIN, IER, RPAR, IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(*), CCPRIME(*), SAVR(*), WK(*), WP(*), IWP(*), B(*),
     1          RPAR(*), IPAR(*)
C       Check what IER is
C       Need to modify this code if JPRE and JBG value from cvcm.f changes
      IER = 0
      CALL DRBDPS (B, WP, IWP)
     
      END SUBROUTINE PSOLCVCM
