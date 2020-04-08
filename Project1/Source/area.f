C
      SUBROUTINE AREA (X,A)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
      include 'LOGICAL.cmn'
C
      IF (LCARTE) THEN
	A = 1.0
      ELSE
	A = 4.* 3.1415926 * X*X
      ENDIF

      RETURN

      END
C
