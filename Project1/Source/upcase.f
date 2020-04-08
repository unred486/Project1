C Subroutine to convert lower case to upper case
       SUBROUTINE UPCASE(ISTR)
      CHARACTER ISTR*(*), LCASE(26)*1, UCASE(26)*1
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      DO 10 J = 1, LEN(ISTR)
         DO 10 N = 1, 26
            IF (ISTR(J:J) .EQ. LCASE(N)) ISTR(J:J) = UCASE(N)
   10 CONTINUE
      RETURN
      END
