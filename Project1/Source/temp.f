C
C----------------------------------------------------------------------
C
      SUBROUTINE TEMP (NPTS, X, XX, TT, T)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION XX(NPTS), TT(NPTS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     THIS SUBROUTINE USES BISECTION TO LINEARLY INTERPOLATE
C     AN ARRAY OF XX,TT PAIRS.  GIVEN AN XX,TT PAIR THIS ROUTINE
C     RETURNS THE INTERPOLATED VALUE OF THE T AT THE POINT X.
C
C INPUT-
C   NPTS   - NUMBER OF XX,TT PAIRS.
C   X      - LOCATION AT WHICH INTERPOLATED T IS DESIRED.
C   XX     - ARRAY OF X POINTS AT WHICH TT ARE GIVEN.
C   TT     - ARRAY OF T VALUES AT THE XX LOCATIONS.
C
C OUTPUT-
C   T     - INTERPOLATED T AT POINT X
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C             check for x outside (1,npts)
C
      IF (X .LE. XX(2)) THEN
         N = 2
         S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
      ELSEIF (X .GE. XX(NPTS-1)) THEN
         N = NPTS-1
         S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
      ELSE
         NLO = 1
         NHI = NPTS
         S   = 0.0
C
C        bisect interval
C
50       CONTINUE
         N = (NLO+NHI)/2
         IF (X .LT. XX(N)) THEN
            IF (X .LT. XX(N-1)) THEN
               NHI = N
               GO TO 50
            ELSEIF (X .EQ. XX(N-1)) THEN
               N = N-1
            ELSE
               S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
            ENDIF
         ELSEIF (X .GT. XX(N)) THEN
            IF (X .GT. XX(N+1)) THEN
               NLO = N
               GO TO 50
            ELSEIF (X .EQ. XX(N+1)) THEN
               N = N + 1
            ELSE
               S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
            ENDIF
         ENDIF
      ENDIF
C
  100 CONTINUE
      T      = TT(N) + S * (X - XX(N))
      RETURN
      END
C--------------------------------------------------------------------