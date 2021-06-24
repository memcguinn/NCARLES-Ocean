! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine estimates the stability functions for momentum, m,   !
!         and scalars, c, from input of the stability parameter zeta = z/L.    !
! ============================================================================ !
!
SUBROUTINE fzol(zeta,phim,phis,psim,psis)
!
    DATA c1/5./
    DATA a3,b3,a4,b4/1.258,8.382,-28.862,98.9545/
    DATA zetam,zetas/-0.2,-1.0/
!
    SAVE c1, a3, b3, a4, b4, zetam, zetas
!
! --------------------------------------------------------------------------- !
!
  psimu(Y)  = 1.571 + 2.0*(ALOG(0.5*(1.0 + Y)) - ATAN(Y)) + ALOG(0.5 + 0.5*Y**2)
  psisu(Y)  = 2.0*ALOG(0.5 + 0.5*Y)
  psicu(Y,G)= (1.0 - G)*ALOG(ABS(Y - 1.0))                  &
            + 0.5*(G + 2.0)*ALOG(ABS(Y**2 + Y + 1.0))       &
            - (2.0*G + 1.0) / SQRT(3.0) * ATAN((Y + 0.5)*2.0/SQRT(3.0))
  Xm(zol)   = (1.0 - 16.0*zol)**0.25
  Xs(zol)   = SQRT(1.0 - 16.0*zol)
  Xc(zol,f) = ABS(1.0 - f*zol)**(4.0/3.0)/(1.0 - f*zol)
!
! STABLE
  IF (zeta .GE. 0.0) THEN
    IF (zeta .LE. 1.0) THEN
      phim = 1.0 + c1 * zeta
      psim = - c1 * zeta
      phis = phim
      psis = psim
!
!   USE LIMITING FORM
    ELSE
      phim = c1 + zeta
      psim = (1.0 - c1)*(1.0 + ALOG(zeta) ) - zeta
      phis = phim
      psis = psim
    END IF
!
! UNSTABLE
  ELSE
!
    IF (zeta .GE. zetam) THEN
      phim = 1.0 / Xm(zeta)
      psim = psimu(Xm(zeta))
!
!   USE CONVECTIVE LIMIT FOR MOMENTUM
    ELSE
      X = (1.0 - b3/a3 * zeta)**(1.0/3.0)
      fm = a3**(-1.0/3.0)
      phim = fm / Xc(zeta,b3/a3)
      psim = psimu(Xm(zetam)) + psicu(Xc(zeta,b3/a3),fm) &
           - psicu(Xc(zetam,b3/a3),fm)
    END IF
!
    IF (zeta .GE. zetas) THEN ! UNSTABLE scalars
      phis = 1.0/Xs(zeta)
      psis = psisu(Xs(zeta))
!
!   USE CONVECTIVE LIMIT FOR SCALARS
    ELSE
      fs =   ABS(a4)**(-1.0/3.0)*ABS(a4)/a4
      phis = (a4 - b4*zeta)**(-1.0/3.0)
      psis = psisu(Xs(zetas)) + psicu(Xc(zeta,b4/a4),fs) &
           - psicu(Xc(zetas,b4/a4),fs)
    END IF
  END IF
!
RETURN
!
END SUBROUTINE fzol
