SUBROUTINE dealias
! wave cutoff filter using 2d fft

    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats

    REAL wve(nny,jxs:jxe,izs:ize)
    REAL wves(nnxp2,iys:iye,izs:ize)

! --- sharp spectral cutoff, specific to current 2dfft
  ix_cut   = 2*INT(FLOAT(nnx)/3.) + 3
  iy_cut_l = INT(FLOAT(nny)/3.) + 2
  iy_cut_u = nnyp2 - iy_cut_l
! --- u-equation ---------------------------------------------------------------
  CALL fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc, &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,  &
                 izs,ize,myid,ncpu_s,numprocs,-2)
  CALL sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
  CALL fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc, &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,  &
                 izs,ize,myid,ncpu_s,numprocs,2)
! --- v-equation ---------------------------------------------------------------
  CALL fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc, &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,  &
                 izs,ize,myid,ncpu_s,numprocs,-2)
  CALL sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
  CALL fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc, &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,  &
                 izs,ize,myid,ncpu_s,numprocs,2)
! --- w-equation ---------------------------------------------------------------
  CALL fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc, &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,  &
                 izs,ize,myid,ncpu_s,numprocs,-2)
  CALL sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
  CALL fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc, &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,  &
                 izs,ize,myid,ncpu_s,numprocs,2)
! --- e-equation ---------------------------------------------------------------
  CALL fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc, &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,  &
                 izs,ize,myid,ncpu_s,numprocs,-2)
  CALL sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
  CALL fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc, &
                 nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,  &
                 izs,ize,myid,ncpu_s,numprocs,2)
! --- scalars, not stored in correct order
  DO iscl=1,nscl
    DO iz=izs,ize
      DO iy=iys,iye
        DO ix=1,nnx
          wves(ix,iy,iz) = t(ix,iy,iscl,iz)
        END DO
      END DO
    END DO
    CALL fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),         &
                   trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e, &
                   izs,ize,myid,ncpu_s,numprocs,-2)
    CALL sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
    CALL fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),         &
                   trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e, &
                   izs,ize,myid,ncpu_s,numprocs,2)
    DO iz=izs,ize
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,iscl,iz) = wves(ix,iy,iz)
        END DO
      END DO
    END DO
  END DO

RETURN
END SUBROUTINE dealias
