SUBROUTINE fft2d_mpi(ax,at,trigx,trigc,nx,ny,             &
                     jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e, &
                     iz1,iz2,myid,ncpu,np,isgn)
!
! -------- get 2d fft using fftpack routines and parallel mpi
!          use fftpack storage a0, (a1,b1), (a2,b2),...,
!
!         isgn = -1 do forward transform, get coefficients
!                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
!                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
!
!         isgn = -2 do forward transform, get coefficients
!                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
!                   outgoing array is at(ny,jxs:jxe,iz1:iz2)
!
!         isgn =  1 do inverse transform, move to physical space
!                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
!                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
!
!         isgn =  2 do inverse transform, move to physical space
!                   incoming array is at(ny,jxs:jxe,iz1:iz2)
!                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)

    INTEGER jx_s(0:np-1), jx_e(0:np-1), iy_s(0:np-1), iy_e(0:np-1)
    REAL ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2)
    REAL trigx(2*nx+15), trigc(4*ny+15), a2d(2,ny), a_wrk(nx)


  nxp2 = nx + 2
  IF (isgn .LT. 0) THEN
    fn   = 1.0/(FLOAT(nx)*FLOAT(ny))
! --- 1d fft in x over [iys,iye] for all z
    DO iz=iz1,iz2
      DO iy=iys,iye
        DO ix=1,nx
          a_wrk(ix) = ax(ix,iy,iz)*fn
        END DO
        CALL rfftf(nx,a_wrk(1),trigx(1))
        ax(1,iy,iz) = a_wrk(1)
        ax(2,iy,iz) = 0.0
        DO  ix=2,nx
          ax(ix+1,iy,iz) = a_wrk(ix)
        END DO
        ax(nx+2,iy,iz) = 0.0
      END DO
    END DO
    CALL xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e, &
                    iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

! --- 1d fft in y over [jxs,jxe] for all z
    DO iz=iz1,iz2
      DO ix=jxs,jxe,2
        DO iy=1,ny
          a2d(1,iy) = at(iy,ix,iz)
          a2d(2,iy) = at(iy,ix+1,iz)
        END DO
        CALL cfftf(ny,a2d(1,1),trigc(1))
        DO iy=1,ny
          at(iy,ix,iz)   = a2d(1,iy)
          at(iy,ix+1,iz) = a2d(2,iy)
        END DO
      END DO
    END DO
! --- decide whether to transpose back or leave as is
    IF (isgn .EQ. -1) THEN
      CALL ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e, &
                      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
    END IF
  ELSE
! --- decide whether to first transpose or leave as is
    IF (isgn .EQ. 1) THEN
      CALL xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e, &
                      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
    END IF
! --- 1d fft in y over [jxs,jxe] for all z
    DO iz=iz1,iz2
      DO ix=jxs,jxe,2
        DO iy=1,ny
          a2d(1,iy) = at(iy,ix,iz)
          a2d(2,iy) = at(iy,ix+1,iz)
        END DO
        CALL cfftb(ny,a2d(1,1),trigc(1))
        DO iy=1,ny
          at(iy,ix,iz)   = a2d(1,iy)
          at(iy,ix+1,iz) = a2d(2,iy)
        END DO
      END DO
    END DO
    CALL ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e, &
                    iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
! ---  1d fft in x over [iys,iye] for all z
    DO iz=iz1,iz2
      DO iy=iys,iye
        a_wrk(1) = ax(1,iy,iz)
        DO ix=2,nx
          a_wrk(ix) = ax(ix+1,iy,iz)
        END DO
        CALL rfftb(nx,a_wrk(1),trigx(1))
        DO ix=1,nx
          ax(ix,iy,iz) = a_wrk(ix)
        END DO
      END DO
    END DO
  END IF

RETURN
END SUBROUTINE fft2d_mpi
