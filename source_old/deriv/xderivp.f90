! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine uses the fft routines and storage (a0, (a1,b1),      !
!         (a2,b2),... to calculate x derivatives. It assumes that the first    !
!         wavenumber xk(1) = 0.0 and the wavenumbers are normalized by a       !
!         number of points, nx.                                                !
! ============================================================================ !
!
SUBROUTINE xderivp(ax,trigx,xk,nnx,iys,iye)
!
!    USE fft
!
    REAL :: xk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)
!
! ---------------------------------------------------------------------------- !
!
  DO iy=iys,iye
    CALL rfftf(nnx,ax(1,iy),trigx)
    ii = 1
    ax(1,iy) = 0.0
    ax(nnx,iy) = 0.0
    DO ix=2,nnx-1,2
      ii          = ii + 1
      temp        = ax(ix,iy)
      ax(ix,iy)   = -xk(ii)*ax(ix+1,iy)
      ax(ix+1,iy) = xk(ii)*temp
    END DO
    CALL rfftb(nnx,ax(1,iy),trigx)
  END DO
!
RETURN
END SUBROUTINE xderivp
