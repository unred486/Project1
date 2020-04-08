        module varray
c        use var
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
        POINTER :: hALL(:),cpALL(:),wdotALL(:),ICKWRKALL(:),
     1              RCKWRKALL(:),WTall(:),XP(:),DMIXall(:)
        POINTER :: PALL
        DOUBLE PRECISION :: RU
        contains
        subroutine SETPAR(ICKWRK,RCKWRK,h,cp,wdot,WT,KK,NATJ,P,LOUT)
            use var

            IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)    
           DOUBLE PRECISION, TARGET :: h(*),cp(*),wdot(*),
     1                                  RCKWRK(*),P,WT(*)
           INTEGER,TARGET ::  ICKWRK(*)
           hALL => h(1:KK)
           cpALL=> cp(1:KK)
           wdotALL=> wdot(1:KK)
           ICKWRKALL=>ICKWRK(1:LENICK)
           RCKWRKALL=>RCKWRK(1:LENRCK)
           PALL=>P
           WTALL=>WT(1:KK)

            
        end subroutine SETPAR
        
        subroutine SETPARZMF(X,JJ,DMIX)
            use var

            IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)    
           DOUBLE PRECISION, TARGET :: X(*),DMIX(*)
           INTEGER :: JJ
           
           XP => X(1:JJ)
           DMIXall => DMIX(1:JJ) 
        end subroutine SETPARZMF
        end module varray
