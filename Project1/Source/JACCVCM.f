      SUBROUTINE JACCVCM(RES, IRES, NEQ, T, CC, CCPRIME, REWT, SAVR, WK,
     1                  H, CJ, WP, IWP, IER, RPAR, IPAR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      EXTERNAL WCVCM
      DIMENSION CC(*), CCPRIME(*), REWT(*), SAVR(*), WK(*), WP(*),
     1          IWP(*), RPAR(*), IPAR(*)
C       Need to modify this code if JPRE and JBG value from cvcm.f changes
        CALL DRBDJA (T, CC, RPAR, WCVCM, WK, REWT, CJ, WP, IWP, IER)

      END SUBROUTINE JACCVCM
