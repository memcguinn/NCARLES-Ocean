! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine sets the boundary condition on upper boundary        !
!         iz=nnz for special radiation.                                        !
!                                                                              !
!         Upper boundary conditions are denoted by index f(.,.,1)              !
! ============================================================================ !
!
SUBROUTINE upper
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    INTEGER istatus(mpi_status_size)
!
! --------------------------------------------------------------------------- !
!
  iz   = nnz
  izm1 = iz - 1
  izm2 = iz - 2
  izp1 = iz + 1
  izp2 = iz + 2
  IF (ibcu .EQ. 0) THEN
!
! FOR GRADIENT BOUNDARY CONDITIONS
    DO iy=iys,iye
      DO ix=1,nnx
        wbc(ix,iy,1) = 0.0
        ebc(ix,iy,1) = 0.0
        ubc(ix,iy,1) = u(ix,iy,iz)
        vbc(ix,iy,1) = v(ix,iy,iz)
        pbc(ix,iy,1) = 0.0
        pbc2(ix,iy,1)= 0.0
      END DO
    END DO
!
    DO iscl=1,nscl                              ! Get average scalar gradient
      dtdzf(iscl) = 0.0
      DO iy=iys,iye
        DO ix=1,nnx
          dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) - &
                        t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
        END DO
      END DO
      dtdzf(iscl) = dtdzf(iscl)*fnxy
    END DO
!
    CALL mpi_sum_xy(dtdzf,myid,iss,ise,nscl)
!
    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + dtdzf(iscl)*dzu(nnzp1)
        END DO
      END DO
    END DO
!
  ELSE IF (ibcu .EQ. 1) THEN
!
!   IRADUP BC TO ESTIMATE W FROM CONTINUITY AND LINEARIZED RELATION FOR PRESSURE
    xmeanp = 0.0
    grad_ug = ug(nnz) - ug((nnz-1))
!
    DO iy=iys,iye
      DO ix=1,nnx
        wbc(ix,iy,1) = w(ix,iy,izm1) - (ux(ix,iy,iz)+vy(ix,iy,iz))*dzw(iz)
        pbc(ix,iy,1) = .5*(w(ix,iy,izm1)+wbc(ix,iy,1))
        ebc(ix,iy,1) = 0.0
        ubc(ix,iy,1) = u(ix,iy,iz) + grad_ug
        vbc(ix,iy,1) = v(ix,iy,iz)
        pbc2(ix,iy,1)= 0.5*(u(ix,iy,iz)**2 + v(ix,iy,iz)**2) +    &
                       0.25*(w(ix,iy,izm1)**2 + wbc(ix,iy,1)**2)
        xmeanp = xmeanp + pbc2(ix,iy,1)
      END DO
    END DO
!
    CALL mpi_sum_xy(xmeanp,myid,iss,ise,1)
!
    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + dtdzf(iscl)*dzu(nnzp1)
        END DO
      END DO
    END DO
!
    xmeanp = xmeanp*fnxy
!
    DO iy=iys,iye
      DO ix=1,nnx
        pbc2(ix,iy,1) = pbc2(ix,iy,1) - xmeanp
      END DO
    END DO
!
  END IF
!
  DO iy=iys,iye
    DO ix=1,nnx
      w(ix,iy,iz)   = wbc(ix,iy,1)
      e(ix,iy,iz)   = ebc(ix,iy,1)
      r3(ix,iy,iz)  = 0.0
      r5(ix,iy,iz)  = 0.0
      u(ix,iy,izp1) = ubc(ix,iy,1)
      v(ix,iy,izp1) = vbc(ix,iy,1)
      r3(ix,iy,izp1)= 0.0
      r5(ix,iy,izp1)= 0.0
      wx(ix,iy,izp1) = 0.0
      wy(ix,iy,izp1) = 0.0
      ux(ix,iy,izp1) = 0.0
      uy(ix,iy,izp1) = 0.0
      vx(ix,iy,izp1) = 0.0
      vy(ix,iy,izp1) = 0.0
    END DO
  END DO
!
  DO iscl=1,nscl
    DO iy=iys,iye
      DO ix=1,nnx
        t(ix,iy,iscl,izp1) = tbc(ix,iy,iscl,1)
        t(ix,iy,iscl,izp2) = tbc(ix,iy,iscl,1)
      END DO
    END DO
  END DO
!
RETURN
END SUBROUTINE upper
