! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine sovles for pressure using a matrix transpose across  !
!         mpi tasks and tridiagonal solver. The transposed array is            !
!         dimensioned (0:nnz+1). Values (0, nnz+1) are not needed. They are    !
!         useful in the matrix transpose when we return. See send_ztox.        !
!         On exit, pressure is defined at all [izs-1:ize+1].                   !
! ============================================================================ !
!
SUBROUTINE pressure
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    INTEGER :: istatus(mpi_status_size)
!
    REAL :: pfft(nny,jxs:jxe,izs-1:ize+1), pt(0:nnz+1,jxs:jxe,jys:jye), &
            ptopfft(nny,jxs:jxe,1:2), psum(1:nnz)
!
! --------------------------------------------------------------------------- !
!
! FOURIER ANALYZE THE RHS AT ALL IZ = IZS,IZE - RESULTS IN PFFT
  CALL fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,      &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,        &
                 izs,ize,myid,ncpu_s,numprocs,-2)
!
! FOURIER ANALYZE THE RADIATION BC ARRAYS
  IF (ibcu .EQ. 1) THEN
    CALL fft2d_mpi(ptop(1,iys,1),ptopfft(1,jxs,1),trigx(1,1),trigc,  &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,1,2,    &
                 myid,ncpu_s,numprocs,-2)
  END IF
!
! TRANSPOSE FIRST AND LAST INDEX OF ARRAY THE ORDER OF PFFT IS Y,X,Z
  CALL xtoz_trans(pfft,pt,nny,nnz,jys,jye,jy_s,jy_e,jxs,jxe,izs,ize, &
                  iz_s,iz_e,myid,ncpu_s,numprocs)
  CALL solve_trid(pt, ptopfft)
!
! TRANSPOSE BACK
  CALL ztox_trans(pt,pfft,nny,nnz,jys,jye,jy_s,jy_e,jxs,jxe,izs,ize, &
                  iz_s,iz_e,myid,ncpu_s,numprocs)
!
  iz_ee = ize+1
!
  IF (ise .EQ. numprocs-1) THEN
    iz_ee = ize
  END IF
!
! INVERSE FFT AT ALL IZ = IZS, IZ_EE TO GET P, SEE Z INDICES
  CALL fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,nnx,nny, &
                 jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,izs,iz_ee,myid,    &
                 ncpu_s,numprocs,2)
!
! PARTIAL SUMS FOR PRESSURE
  DO iz=1,nnz
    psum(iz) = 0.0
  END DO
!
  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        psum(iz) = psum(iz) + p(ix,iy,iz)
      END DO
    END DO
!
    psum(iz) = psum(iz)*fnxy
  END DO
!
  CALL mpi_sum_z(psum,i_root,myid,nnz,1)
!
  DO iz=izs,iz_ee
    DO iy=iys,iye
      DO ix=1,nnx
        p(ix,iy,iz) = p(ix,iy,iz) - psum(iz)
      END DO
    END DO
  END DO
!
RETURN
!
END SUBROUTINE pressure
