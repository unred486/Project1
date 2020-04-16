MODULE f90_module
CONTAINS      
      SUBROUTINE GRIDGEN(XSTR,XEND,XFLAME,XNEW,JJ,IREFINE)
      ! Written By HAN JU LEE  1/7/2019
      ! A subroutine that calculates a new grid according to the 
      ! algorithm written in Vivien Lacoustre's PhD Dissertation
      ! For Region 3, instead of the algorithm, a geometric series is used
      ! XSTR = Position of Starting Element
      ! XEND = Position of Ending Element
      ! XFLAME = Center of Reaction
      ! XNEW = New Grid for Update
      ! JJ = NUmber of Mesh Points
      ! A, Alpha = Algorithm Generation Constants
      ! Dx = Increment Size in Reaction Region
      ! IRN1 = Number of Points in Pre-Reaction Region
      ! IRN2 = Number of Points in Reaction Region
      ! IRN3 = Number of Points in Post - Reaction Region

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
            
      PARAMETER  (A = 1.0, Alpha = 0.1)
      PARAMETER  (Dx = 0.001 , RATIO = 1.05)
      DOUBLE PRECISION, DIMENSION(JJ) ::  XINT
      REAL ::  XMID
      DOUBLE PRECISION, DIMENSION(JJ) :: XNEW
      !Calculating how many points are needed in each region
      IRN1 = JJ/6
      IRN2 = JJ/3
      IRN3 = IRN1+IRN2
      
      !WRITE(LTEST,*) XFLAME
      XMID = 0.9!(XSTR + XFLAME)/2.0
      

      !Calculating the intermediate mesh x_i for region 1 
       DO J = 1, IRN1 
         XINT(J) = (REAL(J)-1.0)/(REAL(IRN1)-1.0)*(XFLAME-(Dx*IRN2/2) &
         -XSTR)+XSTR
       END DO


       !Calculating Region Sum For SUM ( A EXP(-alpha(xi-xmid)^2)) term
       XSUM = XSTR
       DO J=2,IRN1
         XSUM = A*EXP(-Alpha*(XINT(J)-XMID)**2) + XSUM
       END DO
      
       DO J = 1,JJ 
         IF (J .LE. IRN1) THEN
         !Calculating ith Sum for SUM ( A EXP(-alpha(xi-xmid)^2)) term
             XSUM2=XSTR
             IF (J .GE. 2) THEN
                DO K=2,J
                   XSUM2 = A*EXP(-Alpha*(XINT(K)-XMID)**2) + XSUM2
                END DO
             END IF       
                 XNEW(J)= XSTR+(XSUM2-XSTR)/(XSUM-XSTR)*(XFLAME- &
                 (Dx*IRN2/2)-XSTR)
        
         !Creating Geometric Sequence for Region 3
         ELSE IF (J .GT. IRN3) THEN
                 XNEW(J)= RATIO*(XNEW(J-1)-XNEW(J-2))+XNEW(J-1)
                 IF (XNEW(J) .GE. XEND) THEN
                 XNEW(J)=XEND
                 JJ = J
                 EXIT
                 END IF
        !Creating a uniform mesh near flame region
         ELSE
             DO K =1, IRN2
                 XNEW(J) = XNEW(J-1)+Dx
             END DO
         END IF         
       END DO 
       DO J = 1 , JJ
       !WRITE(LTEST,*) XNEW(J)
       END DO
       
      IREFINE=0
      RETURN
      
      END SUBROUTINE
!/////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////
! INIDROP SUBROUTINE
!/////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////      
SUBROUTINE INIDROP(T_s,T_inf,r_f,r_s,r_inf,YF_s,YO_inf,Hc,cp,&
                    k,phi,rrr,TT,n,LOUT)
IMPLICIT none
!     INIDROP calculates the temperature profile obtained from a simple
!     droplet burning analysis using SZ variable. This temperature 
!     profile is to be used for initialization of the simulation
real*8 :: T_s, T_inf, YO_inf,a,YF,YF_s,T_f,T,rr
real*8 :: k, mdot,r_s,r_inf,r_f,phi,B,cp,Hc,Latent,r_s_new,r_inf_new
integer :: n,i,j,LOUT,error
real*8, dimension(*) :: TT,rrr
real*8, parameter :: YF_f=0, YO_f=0, YF_inf=0, YO_s=0
!       T_s         Surface Temperature         (K)
!       T_inf       Ambient Temperature         (K)
!       r_f         Flame Radius                (cm)
!       r_s         Burner Radius               (cm)
!       r_inf       ambient Radius              (cm)
!       Yf_s        Fuel mass fraction at Burner
!       YO_inf      Oxidizer mass fraction at Ambient
!       Hc          Heat of Combustion          (KJ/kg)
!       cp          thermal capacity            (KJ/K-kg)
!       k           thermal conductivity        (KW/m-K)
!       phi         fuel to oxi mass equivalence ratio 
!       T           temporary temperature calculation 
!       rr          temporary position calculation 
!
!       output:
!       rrr         position matrix             (cm)
!                   size is determined during calculation
!       TT           temperature matrix          (K)
!                   size is determined during calculation
!       n           number of data points for initialization       
!
! flame radius calculated from analytical method is known to produce 
! unrealistic flame diameter(often much bigger than the physical one)
! Therefore, flame radius is fed from user's input

! mass fraction at the surface is artificially set if it is set to be 
! too high as it causes the droplet model to function inappropriately
If (YF_s .GE. 0.90) then
YF_s = 0.873
end if
Latent=400
!variable latent is intended to be latent heat of vaporization in the 
!droplet model. while for our burner case, there is no latent heat of 
!vaporization, this term mathmatically smoothes the initial temperature
!profile so it is artificially added.

!converting units from cm to m
r_f=r_f/100
r_s_new=r_s/2000
!surface radius is manually set to r_s/2000 m because starting at the burner
!size of 0.0032 m produces unphysical temperature profile
r_inf_new=r_inf/100
B=-((YF_s)*Hc+cp*(T_s-T_inf))/((YF_s-1)*Hc+Latent)
! fuel flow rate is artificially calculated from the droplet evaporation 
! model
mdot=k/r_s_new/cp*log(1+B)
T_f=T_inf+(YO_inf)*((YF_s-1)*Hc+Latent)/(cp*phi*(YF_s-1))
a=mdot*r_s_new**2/(k/cp)

j=1
do i=1,n
    rr=r_s_new+(i-1)*(r_inf_new-r_s_new)/n
    if (rr.le.r_f) then
    YF=((1+B)*(exp(-a/rr)-1))*(YF_s-1)-YO_inf/phi
    else 
    YF=0
    end if   
T=T_inf+(((1+B)*(exp(-a/rr)-1))*(Latent+(YF_s-1)*Hc)-YF*Hc)/cp
rr=100*rr ! converting to cm from m 
!deleting the radius and temperature points before 0.32cm radius because
!the initial radius was artificially set to a very small number 

if (rr .gt. r_s) then
j=j+1
rrr(j)=rr
TT(j)=T
end if
end do
write(Lout,30)mdot,B, r_f,a,T_f
30 format('mdot= ',f9.2,3X,'B= ', f9.2,3X, 'r_f= ',f12.5,3X, &
            'a= ', f13.7, 'T_f= ',f12.5,3X)

rrr(1)=r_s
rrr(j+1)=r_inf
TT(1)=T_s
TT(j+1)=T_inf
n=j+1

write(LOUT,*)'Temperature Profile Generated from Modified Droplet model'
do i=1,n
write(LOUT,40) rrr(i),TT(i)
40 format('TEMP',3X,f10.4,3X,f10.4)
end do
RETURN
END SUBROUTINE

!/////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////
! TEMP2 SUBROUTINE
!/////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////    
    
SUBROUTINE TEMP2(rr,T,nn,XGIVEN,TGIVEN)
IMPLICIT none
DOUBLE Precision, dimension(nn) :: rr, T, XGIVEN, TGIVEN
integer :: nn,n
  
    DO n =1 , nn
        XGIVEN(n)=rr(n)
        TGIVEN(n)=T(n)
    END Do
    RETURN
END Subroutine


SUBROUTINE timekeeper(active)
implicit none
double precision :: active, active2 , dummy
integer :: error,IFTIMEKEEP,LFTEMPTIME,i
logical :: exist_ftemp_time, exist_timekeeper

LFTEMPTIME = 2617
IFTIMEKEEP = 2619

inquire(file='Temptimex.tec', exist=exist_ftemp_time)
inquire(file='Timekeeper.out', exist=exist_timekeeper)
error=0
if (exist_ftemp_time) then
    Open(UNIT=LFTEMPTIME, FILE='Temptimex.tec', ACTION='read', &
    STATUS='old', FORM='FORMATTED')
    OPEN(UNIT=IFTIMEKEEP, FILE='Timekeeper.out', ACTION='write',&
    STATUS='replace', FORM='FORMATTED')
    do i=1,3
    read(LFTEMPTIME,*)
    end do
    do while(error .ge. 0)
    READ(LFTEMPTIME,10,IOSTAT=error) active2,dummy, dummy, dummy
    end do
    WRITE(IFTIMEKEEP,*) active2
    close(LFTEMPTIME)
    close(IFTIMEKEEP)
else
    OPEN(UNIT=IFTIMEKEEP, FILE='Timekeeper.out', ACTION='write',&
    STATUS='replace', FORM='FORMATTED')
    WRITE(IFTIMEKEEP,*) 0.0
    close(IFTIMEKEEP)
end if
10 format (4 (1PE12.4),2X)
return
end subroutine

      
!//////////////////////////////////////////////////////////////////////
! DERIV calculates dT_flame/dt term and dX_flame/dt term
!//////////////////////////////////////////////////////////////////////
      subroutine DERIV(Tmax,Xflame,ACTIME3,Tderiv,Xderiv,LFTEMPTIME,&
      TIME)
                
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      ierror=0
      II=0
      write(*,*) LFTEMPTIME
      rewind(LFTEMPTIME)
      do i=1,3
      read(LFTEMPTIME,*)
      end do
      TIME=0.0
      do while(ierror .eq. 0)
      read(LFTEMPTIME,10,iostat=ierror)  TIME,Tmaxn,Xflamen,  &
      dummy,Tderiv,Xderiv,dummy
      II=II+1
      end do
      
      dt=ACTIME3-TIME
      Tderiv=(Tmax-Tmaxn)/dt
      Xderiv=(Xflame-Xflamen)/dt
      
      If (II .eq. 1) then
      Tderiv=0.0
      Xderiv=0.0
      END IF  
    
10 format(6 (1PE12.4),2X)
      return
      end subroutine
!//////////////////////////////////////////////////////////////////////
! DERIV2 calculates dX_Zst/dt term
!//////////////////////////////////////////////////////////////////////
      subroutine DERIV2(XZst,ACTIME3,XZstderi,LFTEMPTIME,&
      TIME)
                
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      ierror=0
      II=0
      write(*,*) LFTEMPTIME
      rewind(LFTEMPTIME)
      do i=1,3
      read(LFTEMPTIME,*)
      end do
      TIME=0.0
      do while(ierror .eq. 0)
      read(LFTEMPTIME,10,iostat=ierror)  TIME,dummy,dummy,  &
      Xstn,Tderiv,Xderiv
      II=II+1
      end do
      
      dt=ACTIME3-TIME
      XZstderi =(XZst-Xstn)/dt
      
      If (II .eq. 1) then
      XZstderi=0.0
      END IF  
    
10 format(6 (1PE12.4),2X)
      return
      end subroutine
      
!//////////////////////////////////////////////////////////////////////
! LYPREVIOUS calculates dT_flame/dt term and dX_flame/dt term
! calling LYPREVIOUS subroutine calculates transient terms in 
! continuity equation over spatial domain 
! it is trying to look at the effect of each term     
!//////////////////////////////////////////////////////////////////////
subroutine LYPREVIOUS(S,JJ,KK,KSYM,ACTIME3,X,RHO,NATJ,ICKWRK,RCKWRK, &
                      TIME)
IMPLICIT NONE
COMMON /LOCS/ NT, NM, NYS, NY, NTR, N_moments, NSM

integer, parameter :: LSOLPREV=100 , LdVARdt = 101
integer :: JJ, KK, K, J, NATJ, NT, NM, NYS, NY, NTR, N_moments, NSM, I
integer :: ICKWRK(*)
logical :: exist_LSOLPREV
Character :: KSYM(*)*(*)
character(len=32) :: my_fmt
real*8 :: ACTIME3,U,areaj,XN, SNT,RHON,SNM,UN,TIME,YSUM,wtm,dt,ydterm
real*8 :: tdtterm,YSUM2
real*8, dimension(JJ) :: X,RHO
real*8 ::  S(NATJ,JJ),SYN(KK), RCKWRK(*)
! JJ = total number of grid points
! KK = total number of species
! KSYM = ARRAY of CHEMKIN SPECIES


! if previous solution file, LSOLPREV.tec exists, then this is a 
! dvariable/dt term can be calculated and requires that dvar/dt terms 
! are stored in LdVardt.tec
!
! if previous solution file does not exists, then this subroutine
! simply creates LSOLPREV.tec file and saves the current solution data
! in it for the next time step data to use for the calculation of dvar/dt

! checking if previous solution file, LSOLPREV.tec exists
inquire(file ='LSOLPREV.tec', exist = exist_LSOLPREV)


WRITE(my_fmt , '(I0)')  KK+5
my_fmt = '('//trim(my_fmt)//'(1PE12.4))' 


if (exist_LSOLPREV) then
    open(unit=LSOLPREV, file = 'LSOLPREV.tec', action='readwrite', &
    status ='old', form = 'formatted')
    open(unit = LdVARdt, file = 'LdVardt.tec', action ='readwrite', &
    position ='append', form='formatted')
!  writing solution variable names to LdVardt.tec
    write(LdVardt,8026)
    WRITE(LdVardt, 80271)
    WRITE(LdVardt, 8028) ACTIME3,JJ  


!   reading the first 3 liness of LSOLPREV.tec that contains 
!   non-data related info
    do I =1 ,3
    read(LSOLPREV,*)
    end do
!   reading previous "timestep" data. It is the the real timestep that 
!   solution marched through. It is the timestep from when
!   the previous solution was saved to the current solution time.    
!   This loop also writes the dvariable/dt data in LdVardt.tec

    do J =1 ,JJ
    read(LSOLPREV,my_fmt) XN, SNT,RHON,SNM,UN, (SYN(K),K=1,KK)
!   calculating dydt sum of species
    call dysum(TIME,SYN,ICKWRK,RCKWRK,S,KK,YSUM,ACTIME3,J,JJ,NATJ,dt,&
              YSUM2)
    call ckmmwy (s(ny,j), ickwrk, rckwrk, wtm)
    ydterm=rho(j)*wtm*ysum
!   calculating dtdt term
    tdtterm=(rho(j)/s(nt,j)*(s(nt,j)-SNT)/dt)
    call area ( X(J) , areaj )
    U=S(NM,J)/RHO(J)/areaj
!   writing variables in LdVardt.tec file
    write(LdVardt,80299) X(J),tdtterm,ydterm,U,S(nm,j),rho(j),YSUM2 
    
    end do    
close(LSOLPREV)
close(LdVARdt)
!
!   overwriting the LSOLPREV.tec with the current solution for the next
!   timestep to use   
    open(unit=LSOLPREV, file = 'LSOLPREV.tec', action='readwrite', &
    status ='replace', form = 'formatted')   
    write(LSOLPREV,8026)
    WRITE(LSOLPREV, 80270) (trim(KSYM(K)), K=1,KK)
    WRITE(LSOLPREV, 8028) ACTIME3,JJ   
    
    DO J =1,JJ
    
    
    call area ( X(J) , areaj )
    U=S(NM,J)/RHO(J)/areaj
    write(LSOLPREV,my_fmt) X(J), S(NT,J),RHO(J), S(NM,J), U,   &
    (S(NYS+K,J),K=1,KK)
    END DO
    close(LSOLPREV)

else
    open(unit=LSOLPREV, file = 'LSOLPREV.tec', action='readwrite', &
    status ='new', form = 'formatted')
    open(unit = LdVARdt, file = 'LdVardt.tec', action ='readwrite', &
    status='new',form='formatted')
    close(LdVARdt)
    
    write(LSOLPREV,8026)
    WRITE(LSOLPREV, 80270) (trim(KSYM(K)), K=1,KK)
    WRITE(LSOLPREV, 8028) ACTIME3,JJ
    
    DO J =1,JJ
    call area ( X(J) , areaj )
    U=S(NM,J)/RHO(J)/areaj
    write(LSOLPREV,my_fmt) X(J), S(NT,J),RHO(J), S(NM,J), U,   &
    (S(NYS+K,J),K=1,KK)
    END DO
    close(LSOLPREV)
    
end if


8028  FORMAT('ZONE T="',1PE11.3,' s", I=',I4,' F=POINT')
80270 FORMAT('Distance_from_burner(cm) Temperature(K) &
                RHO M U(cm/s)', <KK> (' Y'A''))
80271 FORMAT('Distance_from_burner(cm) dT/dt dY/dt &
             U(cm/s) M rho YSUM')
8026  FORMAT('VARIABLES=')
80299 FORMAT (7(1PE12.4),2X)
end subroutine  

subroutine dysum(TIME,SYN,ICKWRK,RCKWRK,S,KK,YSUM,ACTIME3,J,JJ,NATJ,dt,&
                YSUM2)
implicit none
COMMON /LOCS/ NT, NM, NYS, NY, NTR, N_moments, NSM

integer :: KK,K,NT, NM, NYS, NY, NTR, N_moments, NSM,J,JJ, ICKWRK(*)
integer :: NATJ
real*8 ::  S(NATJ,JJ), RCKWRK(*), TIME,SYN(KK),WT(KK),YSUM,ACTIME3,dt, &
            YSUM2

CALL CKWT   (ICKWRK, RCKWRK, WT)
dt=ACTIME3-TIME
ysum=0.0
ysum2=0.0
DO k=1,kk
ysum=ysum+1/wt(k)*(s(nys+k,j)-syn(k))/dt
ysum2=ysum2+s(nys+k,j)
end do

end subroutine

Subroutine XFIND(XXZ,ZGIVEN,ZST,XST)
implicit none

DOUBLE PRECISION :: XXZ(*),ZGIVEN(1500)
DOUBLE PRECISION :: XST,ZST
INTEGER :: J
J=1
DO WHILE (ZGIVEN(J) .GE. ZST)
J=J+1
END DO
J=J-1
IF (ZGIVEN(J) .GT. ZST) THEN
XST=XXZ(J)+(XXZ(J+1)-XXZ(J))/(ZGIVEN(J+1)-ZGIVEN(J))*(ZST-ZGIVEN(J))
ELSE
XST=XXZ(J)
END IF

END Subroutine XFIND

subroutine updategridzmf(NATJ,JJ,NMAX,C,X,DMIX,XST,ZST,XSTR,XEND,IPAR,&
    IW,NEQ)
use var
implicit none

double precision :: XNEW(NMAX),C(NMAX),SS(NATJ,NMAX),X(NMAX),DMIX(NMAX),DM(NMAX) &
    ,XOLD(NMAX)
double precision :: XST,ZST,XSTR,XEND,DMAX
integer :: IPAR(*),IW(*)
integer :: IRMAX,IREFINE,JJ,J,NMAX,NATJ,JOLD,NEQ

        JOLD=JJ
!      DO J=1, JJ
!          DM(J)=DMIX(J)
!          XOLD(J)= X(J)
!      END DO
!      IRMAX = 400
      CALL XFIND(X,C,ZST,XST)
      IF (XST .GT. 0.35) THEN
      write(*,*)'hello'
      
      END IF
      !if (XST .GT. 0.35) then
      !CALL GRIDGEN(XSTR,XEND,XST,XNEW,IRMAX,IREFINE)
      !else 
      CALL GRIDGENTWO(XSTR,XEND,XST,XNEW,JJ,NMAX,IRMAX)
      !end if
      CALL COMPRESSZMF(X,C,JOLD,XNEW,SS,IRMAX,NATJ,XST)
!      do J=1,JJ
!      write(*,*) S(NATJ,J)
!      end do
      JJ=IRMAX   
      CALL DCOPY (NATJ*JJ, SS, 1, C, 1)
      !CALL DCOPY (JJ,XNEW,1,X,1)
      DO J=1,JJ
          X(J)=XNEW(J)
      END DO
      
      !DO J = 1, JJ
      !   CALL TEMP (JOLD, X(J), XOLD, DM,DMIX(J))
      !ENDDO 
      !temporary DMAX
      DMAX=3.06
      do J=1,JJ
          DMIX(J)=0.0
      end do    
      CALL DMIXSETUP(X,JJ,DMAX,XST,DMIX,XSTR)
      IPAR(3)=NATJ*JJ
      IPAR(4)=JJ
      NEQ=NATJ*JJ
!      IW(27) = LENWP-NEQ
!      IW(28) = NEQ
END SUBROUTINE updategridzmf             
    
SUBROUTINE GRIDGENTWO(XSTR,XEND,XST,XNEW,JJ,NMAX,IRMAX)
implicit none
double precision :: XSTR, XEND, XST,DX
double precision :: XNEW(NMAX)
integer :: IRMAX,JJ,ZONE_1,ZONE_2,I,J,NMAX,IREFINE

if (XST .GT. 0.35) then 
    JJ=1500
CALL GRIDGEN3(XSTR,XEND,XST,XNEW,JJ,IREFINE)
IRMAX=JJ
else
XNEW(1)=XSTR
ZONE_2=1000
DX = 0.0001
DO I=2,ZONE_2
XNEW(I)=XNEW(I-1)+DX
if (abs(xst-xnew(i)) .lt. dx/2.0) XNEW(I)=XST
END DO
I=ZONE_2
DO WHILE (XNEW(I) .LE. XEND)
    XNEW(I+1)=XNEW(I)+((XNEW(I)-XNEW(I-1))*1.03)
    I=I+1
END DO
XNEW(I)=XEND
IRMAX=I
ENDIF

END SUBROUTINE GRIDGENTWO

SUBROUTINE GRIDGEN3(XSTR,XEND,XFLAME,XNEW,JJ,IREFINE)
      ! Written By HAN JU LEE  1/7/2019
      ! A subroutine that calculates a new grid according to the 
      ! algorithm written in Vivien Lacoustre's PhD Dissertation
      ! For Region 3, instead of the algorithm, a geometric series is used
      ! XSTR = Position of Starting Element
      ! XEND = Position of Ending Element
      ! XFLAME = Center of Reaction
      ! XNEW = New Grid for Update
      ! JJ = NUmber of Mesh Points
      ! A, Alpha = Algorithm Generation Constants
      ! Dx = Increment Size in Reaction Region
      ! IRN1 = Number of Points in Pre-Reaction Region
      ! IRN2 = Number of Points in Reaction Region
      ! IRN3 = Number of Points in Post - Reaction Region

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
            
      PARAMETER  (A = 1.0, Alpha = 0.1)
      PARAMETER  (Dx = 0.0001 , RATIO = 1.005)
      DOUBLE PRECISION, DIMENSION(1500) ::  XINT
      double precision ::  XMID
      DOUBLE PRECISION, DIMENSION(1500) :: XNEW
      !Calculating how many points are needed in each region
      IRN1 = 200
      IRN2 = 100
      IRN3 = IRN1+IRN2
      
      !WRITE(LTEST,*) XFLAME
      XMID = (XSTR + XFLAME)/2.0
      

      !Calculating the intermediate mesh x_i for region 1 
       DO J = 1, IRN1 
         XINT(J) = (REAL(J)-1.0)/(REAL(IRN1)-1.0)*(XFLAME-(Dx*IRN2/2) &
         -XSTR)+XSTR
       END DO


       !Calculating Region Sum For SUM ( A EXP(-alpha(xi-xmid)^2)) term
       XSUM = XSTR
       DO J=2,IRN1
         XSUM = A*EXP(-Alpha*(XINT(J)-XMID)**2) + XSUM
       END DO
      
       DO J = 1,JJ 
         IF (J .LE. IRN1) THEN
         !Calculating ith Sum for SUM ( A EXP(-alpha(xi-xmid)^2)) term
             !XSUM2=XSTR
             !IF (J .GE. 2) THEN
             !   DO K=2,J
             !      XSUM2 = A*EXP(-Alpha*(XINT(K)-XMID)**2) + XSUM2
             !   END DO
             !END IF       
             !    XNEW(J)= XSTR+(XSUM2-XSTR)/(XSUM-XSTR)*(XFLAME- &
             !    (Dx*IRN2/2)-XSTR)
             XNEW(J)=(dble(J)-1.0)/(dble(IRN1)-1.0)*(XFLAME-(DX*dble(IRN2)/2.0) &
                 -XSTR)+XSTR
        
         !Creating Geometric Sequence for Region 3
         ELSE IF (J .GT. IRN3) THEN
                 XNEW(J)= RATIO*(XNEW(J-1)-XNEW(J-2))+XNEW(J-1)
                 IF (XNEW(J) .GE. XEND) THEN
                 XNEW(J)=XEND
                 JJ = J
                 EXIT
                 END IF
        !Creating a uniform mesh near flame region
         ELSE
             DO K =1, IRN2
                 XNEW(J) = XNEW(J-1)+Dx
             END DO
         END IF         
       END DO 
       DO J = 1 , JJ
       !WRITE(LTEST,*) XNEW(J)
       END DO
       
      IREFINE=0
      RETURN
end subroutine GRIDGEN3

SUBROUTINE DMIXSETUP(X,JJ,DMAX,XST,DMIX,XSTR)
implicit none

double precision :: x(1500),dmix(1500)
double precision :: dmax,xst,xstr,gap
integer :: jj,j

gap=0.35
! if (xst .gt. gap) then 
do j=1,jj
    if (x(j) .ge. xst) then
        dmix(j)=-(DMAX-0.178)/(2.0)*(x(j)-xst)+dmax
        if (dmix(j) .lt. 0.178) dmix(j)=0.178
    else
    dmix(j)=0.178
    endif
end do    
!else
!do j=1,jj
!    if (x(j) .gt. gap) then
!        dmix(j)=-(DMAX-0.178)/(2.0)*(x(j)-gap)+dmax
!        if (dmix(j) .lt. 0.178) dmix(j)=0.178
!    else
!    dmix(j)=0.178
!    endif
!end do    
!end if





end subroutine DMIXSETUP
END MODULE f90_module
