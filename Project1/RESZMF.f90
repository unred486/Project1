!C////////////////// RESZMF                    ////////////////////////
!C//////////////////////////////////////////////////////////////////////
!C       Sets Up Boundary Condition and interior residuals
!C       Therefore, simply calls FCVCM which calculates interior residuals
      SUBROUTINE RESZMF(T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR)
      use var
      use varray
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
!C     need for UPRIME calculation      
      DIMENSION :: U(1000),UPRIME(1000),DELTA(1000),RPAR(*),IPAR(*)
      COMMON /SP/ KK,NATJ
      
      NEQ=IPAR(3)
      JJ= IPAR(4)
      ZBUR=RPAR(1)
      ZAMB=RPAR(2)
      RHO = 1.187E-3
!     M = 10.0 mg/s
      DMFR=10.0e-3
      dxp  =     (xp(2) - xp(1)  )
      xmid =    0.5*(xp(2)+xp(1))
      call area ( xmid, Areap    )
      call area (xp(1),Ai)
      V = DMFR/(RHO*Ai)
      
      !!////////// Dirichlet BC 
      !DELTA(1)=U(1)-ZBUR
      !DELTA(JJ)=U(JJ)-ZAMB
      
     !/////////////// Diffusion into burner //////////////////
      !DELTA(1)  =   -1.0*DMIXall(1)*((U(2)-U(1))/dxp)+V*(U(1)-1.0)
      DELTA(1)  =   -(1.0-U(1))-DMIXall(1)*areap*((U(2)-U(1))/dxp)/V
      DELTA(JJ) =   U(JJ)-ZAMB   
      
      Call FZMF(U,DELTA,IPAR)
      DO I=2,NEQ-1
      DELTA(I)=UPRIME(I)-DELTA(I)
      END DO

      RETURN     
    
    END SUBROUTINE  RESZMF
      

!C///////////////////////// FZMF     /////////////////////////////////
!C//////////////////////////////////////////////////////////////////////
!C       Residual consists of algebraic source term( ex : reaction term)
!C       and differential term ( ex : tranport term)
!C       DDASPK can use preconditions which differentiates between 
!C       differential and algebraic terms to reduce the memory usage
!C       FZMF calculates differential residuals
      SUBROUTINE FZMF(U,CRATE,IPAR)     
      use varray
      use var
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION :: U(*),CRATE(*),IPAR(4),CONV(1000),DIFF(1000),DCHAN(1000)
      
      COMMON /SP/ KK,NATJ
        
      !     Hard coded constant density ==> rho = 1.187 e-3 gm/cm^3      
      RHO = 1.187E-3
!     M = 10.0 mg/s
      DMFR=10.0e-3
      D=1.0
      JJ= IPAR(4)
!    DO I=1,NATJ
      DO J=2,JJ-1 
          K=J
          DCHAN(K)=xp(K)-xp(K-1)
          xmid = 0.5 * ( xp(j)   + xp(j+1) )
          xmidm= 0.5 * ( xp(j-1) + xp(j)   )
          dxp  =     (xp(j+1) - xp(j)  )
          dxm  =     (xp(j)   - xp(j-1))
          dxav = 0.5*(xp(j+1) - xp(j-1))
          
          call area ( xmidm, aream )
          call area ( xmid , areap )
          call area ( Xp(J), Ai    )          
          V = DMFR/(RHO*Ai)  !velocity at ith grid          

          
         ! (cond(j)*areap*(s(nt,j+1)-s(nt,j))/dxp-cond(j-1)*aream*(s(nt,j)-s(nt,j-1))/dxm )
          CONV(K)=(V*(U(K)-U(K-NATJ))/dxm)
          FirstDIFF=areap*DmixALL(K)*(U(K+NATJ)-U(K))/dxp
          SecondDIFF=aream*DmixALL(K-NATJ)*(U(K)-U(K-NATJ))/dxm
          !FirstDIFF=D*(U(K+NATJ)-U(K))/dxp
          !SecondDIFF=D*(U(K)-U(K-NATJ))/dxm
          DIFF(K)=(1/Ai)*(FirstDIFF-SecondDIFF)/dxav
          !DIFF(K)=(FirstDIFF-SecondDIFF)/dxav
          
          !CRATE(K)=-conv(k)
          CRATE(K)=0.1*DIFF(K)-CONV(K)
          !CRATE(K)=DIFF(K)
      END DO
 !   END DO
      
      END SUBROUTINE FZMF
