      MODULE INI
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)      
      INTEGER :: NTOTG,INFOT
      contains 
      SUBROUTINE INIKEY(LIN,LOUT)
C*****precision > double
      use var
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)

      
      DIMENSION :: VALUE(5)
      CHARACTER :: KEYWRD*4, LINE*100
      LOGICAL :: IERR, KERR
      INTEGER :: INFO2
      LZMIX = .FALSE.
      LHOM  = .FALSE. 
C
C         READ NEXT INPUT LINE
C         KEYWRD gets read line by line and is initiated accordingly
C         by reading 4 letter keywords
C
90    CONTINUE
      KEYWRD = ' '
      LINE = ' '
      IERR = .FALSE.
      READ  (LIN,  7000) KEYWRD, LINE
c      WRITE (LOUT, 8000) KEYWRD, LINE
      CALL UPCASE (KEYWRD)
C
C               IS THIS A KEYWORD COMMENT?
C
      IF (KEYWRD(1:1) .EQ. '.' .OR. KEYWRD(1:1) .EQ. '/') GO TO 90
      IND = INDEX(LINE,'(')
      IF (IND .GT. 0) LINE(IND:) = ' '
C
C--------------PROBLEM TYPE KEYWORDS--------------------
C          
C     TOTAL NUMBER OF ALLOWABLE GRID 
      IF (KEYWRD .EQ. 'NTOT')     THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NTOTG = INT(VALUE(1))
C     Controls closed homogeneous reactor's solver method
C     INFO(12) which is a control matrix for ddaspk is set to be 0 if direct method 
C     ELSE INFO(12) is 1 if krylov method         
      ELSEIF (KEYWRD .EQ. 'DIRC')  THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IF (INT(VALUE(1)) .EQ. 0) THEN
         INFOT = 0
         ELSE
         INFOT = 1
         END IF
      ELSEIF (KEYWRD .EQ. 'ZMIX') THEN
          LZMIX = .TRUE.
      ELSEIF (KEYWRD .EQ. 'LHOM') THEN
          LHOM  = .TRUE.
      ELSEIF (KEYWRD .EQ. 'END ') THEN
         GO TO 6000
      ELSE

C
C--------------END OF KEYWORDS--------------------
C
C        TO GET HERE, AN INVALID KEYWORD WAS READ
C
      ENDIF
      GO TO 90
6000  CONTINUE
      RETURN
      
7000  FORMAT (A4, A)
8000  FORMAT (10X, A4, A76)      
      END SUBROUTINE
      END MODULE INI
