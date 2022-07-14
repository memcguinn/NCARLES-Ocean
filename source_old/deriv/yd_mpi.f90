! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine uses the fft routines, mpi and storage (a0, (a1,b1), !
!         (a2,b2),... to calculate y derivatives. It assumes that the first    !
!         wavenumber yk(1) = 0.0 and the wavenumbers are normalized by a       !
!         number of points, ny.                                                !
! ============================================================================ !
!
SUBROUTINE yd_mpi(ay,trigy,yk,nx,ny,ixs,ixe,ix_s,ix_e,    &
                                        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
!
!    USE fft
!
    REAL :: yk(ny), trigy(2*ny+15), ay(nx,iys:iye,iz1:iz2)
    REAL :: ayt(ny,ixs:ixe,iz1:iz2)
!
    INTEGER :: ix_s(0:np-1), ix_e(0:np-1), iy_s(0:np-1), iy_e(0:np-1)
!
! --------------------------------------------------------------------------- !
!
  CALL xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e, &
                                        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
!
  DO iz=iz1,iz2
    DO ix=ixs,ixe
      CALL rfftf(ny,ayt(1,ix,iz),trigy)
!
      ii = 1
      ayt(1,ix,iz)  = 0.0
      ayt(ny,ix,iz) = 0.0
!
      DO iy=2,ny-1,2
        ii              = ii + 1
        temp            = ayt(iy,ix,iz)
        ayt(iy,ix,iz)   = -yk(ii)*ayt(iy+1,ix,iz)
        ayt(iy+1,ix,iz) = yk(ii)*temp
      END DO
!
      CALL rfftb(ny,ayt(1,ix,iz),trigy)
!
    END DO
  END DO
!
  CALL ytox_trans(ayt,ay,nx,ny,ixs,ixe,ix_s,ix_e, &
                  iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
!
RETURN
END SUBROUTINE yd_mpi
