! ============================================================================ !
! ABOUT:                                                                       !
!           This subroutine fills surface arrays on root processors.           !
! ============================================================================ !
!
SUBROUTINE f_suft2(rbuf,nnx,mxs,mxe,iys,iye,nscl,tau13m,tau23m,taut3m,t_grnd)
!
    REAL :: rbuf(2+2*nscl,mxs:mxe,iys:iye)
    REAL :: tau13m(nnx,iys:iye), tau23m(nnx,iys:iye)
    REAL :: taut3m(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl)
!
! --------------------------------------------------------------------------- !
!
  DO iy=iys,iye
    DO ix=mxs,mxe
      tau13m(ix,iy) = rbuf(1,ix,iy)
      tau23m(ix,iy) = rbuf(2,ix,iy)
    END DO
  END DO
  DO iscl=1,nscl
    DO iy=iys,iye
      DO ix=mxs,mxe
        taut3m(ix,iy,iscl) = rbuf(2+iscl,ix,iy)
        t_grnd(ix,iy,iscl) = rbuf(2+nscl+iscl,ix,iy)
      END DO
    END DO
  END DO

RETURN
END SUBROUTINE f_suft2