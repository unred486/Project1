subroutine restart(LRCRVR,X,S,JJ,TOUT,NMAX,NATJ,DMIX,ZST,XSTR,XST,DMAX)
    use var
    use f90_module
    implicit none
    integer :: LRCRVR,JJ,NMAX,J,NATJ
    double precision :: X(NMAX),S(NMAX),DMIX(NMAX)
    double precision :: TOUT,ZST,XST,DMAX,XSTR
    
    READ (LRCRVR) NATJ, JJ,TOUT,DMAX
    READ (LRCRVR) (X(J), J=1,JJ)
    READ (LRCRVR) (S(J), J=1,JJ)
    REWIND(LRCRVR)

    CALL XFIND(X,S,ZST,XST)
    do J=1,JJ
          DMIX(J)=0.0
    end do    
      CALL DMIXSETUP(X,S,JJ,DMAX,XST,DMIX,XSTR)
    end subroutine restart