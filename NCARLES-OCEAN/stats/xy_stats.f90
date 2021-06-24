! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine solves for the statistics of the flow. We use mpi    !
!         to reduce over the z-directions. We are then able to get averages by !
!         setting stat(.,.) = 0 for all z per process.                         !
!                                                                              !
!         Note:                                                                !                                                       !
!            stat(.,1) = u*u = ups                                             !
!            stat(.,2) = v*v = vps                                             !
!            stat(.,3) = w*w = wps                                             !
!            stat(.,4) = w**3 = wcube                                          !
!            stat(.,5) = w**4 = wfour                                          !
!            stat(.,6) = resolved tke at zw = englez                           !
!            stat(.,7) = sgs e at zu = engsbz                                  !
!            stat(.,8) = sgs e at zw = eavg                                    !
!            stat(.,9) = resolved uw at zw = uwle                              !
!            stat(.,10) = resolved vw at zw = vwle                             !
!            stat(.,m1) = resolved scalar flux wt at zw = wtle                 !
!            stat(.,m2) = resolved scalar flux ut at zw = utle                 !
!            stat(.,m3) = resolved scalar flux vt at zw = vtle                 !
!            stat(.,m4) = scalar t*t at zu = tps                               !
!            stat(.,m5) = scalar t*t*t at zu = tcube                           !
! ============================================================================ !
!
SUBROUTINE xy_stats
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
!
    PARAMETER (js = 10,                         & ! Number of non-scalar stats
               ns = 5,                          & ! Number of scalar stats
               nstat = js + ns*nscl)
!
    REAL :: stat(1:nnz,nstat)
!
! ---------------------------------------------------------------------------- !
!
  DO i=1,nstat
    DO iz=1,nnz
      stat(iz,i) = 0.0
    END DO
  END DO
!
! DEFINE SCALAR INDICES
  m1 = js
  m2 = js + nscl
  m3 = js + 2*nscl
  m4 = js + 3*nscl
  m5 = js + 4*nscl
  sgn = 1.0
!
  IF (iupwnd .EQ. 1) sgn = -1.0
  DO iz=izs,ize
    izp2 = iz + 2
    izp1 = iz + 1
    izm1 = iz - 1
    DO iy=iys,iye
      DO ix=1,nnx
        stat(iz,1) = stat(iz,1) + (u(ix,iy,iz) - uxym(iz))**2
        stat(iz,2) = stat(iz,2) + (v(ix,iy,iz) - vxym(iz))**2
        stat(iz,3) = stat(iz,3) + (w(ix,iy,iz) - wxym(iz))**2
        stat(iz,4) = stat(iz,4) + (w(ix,iy,iz) - wxym(iz))**3
        stat(iz,5) = stat(iz,5) + (w(ix,iy,iz) - wxym(iz))**4
        stat(iz,6) = stat(iz,6) + ((w(ix,iy,iz)-wxym(iz))**2 +             &
              (0.5*(u(ix,iy,iz)-uxym(iz) + u(ix,iy,izp1)-uxym(izp1)))**2 + &
              (0.5*(v(ix,iy,iz)-vxym(iz) + v(ix,iy,izp1)-vxym(izp1)))**2)*0.5
        stat(iz,7) = stat(iz,7) + 0.5*(e(ix,iy,iz)+e(ix,iy,izm1))
        stat(iz,8) = stat(iz,8) + e(ix,iy,iz)
        stat(iz,9) = stat(iz,9) + (w(ix,iy,iz)-wxym(iz))*                  &
                   0.5*((u(ix,iy,iz)-uxym(iz))+(u(ix,iy,izp1)-uxym(izp1)))
        stat(iz,10) = stat(iz,10) + (w(ix,iy,iz)-wxym(iz))*                &
                   0.5*((v(ix,iy,iz)-vxym(iz))+(v(ix,iy,izp1)-vxym(izp1)))
      END DO
    END DO
!
!   GET SCALAR RESOLVED FLUXES AND VARIANCES
    DO l=1,nscl
      IF (iupwnd .NE. 1 .OR. iz .EQ. nnz) THEN
        DO iy=iys,iye
          DO ix=1,nnx
            stat(iz,m1+l)=stat(iz,m1+l) +                                      &
                      (w(ix,iy,iz)-wxym(iz))*0.5*(t(ix,iy,l,iz)-txym(iz,l) +   &
                      t(ix,iy,l,izp1)-txym(izp1,l))
          END DO
        END DO
      ELSE
!
!       MONOTONE FLUXES
        DO iy=iys,iye
          DO ix=1,nnx
          stat(iz,m1+l) = stat(iz,m1+l)                                      &
                      + amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,l,iz)             &
                      + rlim(t(ix,iy,l,izp1),t(ix,iy,l,iz),t(ix,iy,l,izm1))) &
                      +     amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,l,izp1)       &
                      + rlim(t(ix,iy,l,iz),t(ix,iy,l,izp1),t(ix,iy,l,izp2)))
          END DO
        END DO
      END IF
!
      stat(iz,m1+l)= sgn*stat(iz,m1+l)
!
!     GET HORIZONTAL SCALAR RESOLVED FLUXES
      DO iy=iys,iye
        DO ix=1,nnx
          stat(iz,m2+l) = stat(iz,m2+l)                                    &
                        + (u(ix,iy,iz)-uxym(iz))*(t(ix,iy,l,iz)-txym(iz,l))
               stat(iz,m3+l) = stat(iz,m3+l)+                              &
                             (v(ix,iy,iz)-vxym(iz))*(t(ix,iy,l,iz)-txym(iz,l))
        END DO
      END DO
!
!     FIND SCALAR VARIANCES AND HIGHER MOMENTS
      DO  iy=iys,iye
        DO ix=1,nnx
          stat(iz,m4+l) = stat(iz,m4+l) + (t(ix,iy,l,iz) - txym(iz,l))**2
          stat(iz,m5+l) = stat(iz,m5+l) + (t(ix,iy,l,iz) - txym(iz,l))**3
        END DO
      END DO
!
    END DO
  END DO
!
! FIND PARTIAL SUMS AND SEND
  CALL mpi_sum_z(stat(1,1),i_root,myid,nnz*nstat,1)
!
! FILL ARRAYS FOR PRINTOUT
  DO iz=1,nnz
    ups(iz)    = stat(iz,1)*fnxy
    vps(iz)    = stat(iz,2)*fnxy
    wps(iz)    = stat(iz,3)*fnxy
    wcube(iz)  = stat(iz,4)*fnxy
    wfour(iz)  = stat(iz,5)*fnxy
    englez(iz) = stat(iz,6)*fnxy
    engsbz(iz) = stat(iz,7)*fnxy
    eavg(iz)   = stat(iz,8)*fnxy
    uwle(iz)   = stat(iz,9)*fnxy
    vwle(iz)   = stat(iz,10)*fnxy
    uw_tot(iz) = uwle(iz) + uwsb(iz)
    vw_tot(iz) = vwle(iz) + vwsb(iz)
!
!   GET SCALAR RESOLVED FLUXES AND VARIANCES
    DO l=1,nscl
      wtle(iz,l)   = stat(iz,m1+l)*fnxy
      utle(iz,l)   = stat(iz,m2+l)*fnxy
      vtle(iz,l)   = stat(iz,m3+l)*fnxy
      tps(iz,l)    = stat(iz,m4+l)*fnxy
      tcube(iz,l)  = stat(iz,m5+l)*fnxy
      wt_tot(iz,l) = wtle(iz,l) + wtsb(iz,l)
    END DO
  END DO
!
RETURN
END SUBROUTINE xy_stats
