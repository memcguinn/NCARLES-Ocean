! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine sets up the pressure solver.                         !
! ============================================================================ !
!
SUBROUTINE comp_p
!
    USE pars
    USE fftwk
    USE fields
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    REAL fnt1(nnx,iys:iye,izs:ize), fs(nnx,iys:iye,2), fr(nnx,iys:iye,2)
!
    INTEGER istatus(mpi_status_size)
!
! --------------------------------------------------------------------------- !
!
  gami = 1.0/dtgama
  nb = myid - ncpu_s
  nt = myid + ncpu_s
!
! SEND BOTH R3 AND UPDATED W (FROM COMP1) TO PROCESSOR ABOVE CURRENT MYID
  IF (iss .EQ. 0) THEN
    nb = mpi_proc_null
  END IF
!
  IF (ise .EQ. numprocs-1) THEN
    nt = mpi_proc_null
  END IF
!
  nsend = 2*nnx*(iye + 1 - iys)
  nrecv = nsend
!
  DO iy=iys,iye
    DO ix=1,nnx
      fs(ix,iy,1) = r3(ix,iy,ize)
      fs(ix,iy,2) = w(ix,iy,ize)
    END DO
  END DO
!
  CALL mpi_sendrecv(                         &
           fs(1,iys,1),nsend,mpi_real8,nt,2, &
           fr(1,iys,1),nrecv,mpi_real8,nb,2, &
           mpi_comm_world,istatus,ierr)
!
  IF (iss .NE. 0) THEN
    DO iy=iys,iye
      DO ix=1,nnx
        r3(ix,iy,izs-1) = fr(ix,iy,1)
        w(ix,iy,izs-1)  = fr(ix,iy,2)
      END DO
    END DO
  END IF
!
! SETUP GENERAL PRESSURE CALCULATION RELIES ON RHS FROM STEP N-1, INCLUDED
! IN VELOCITY ARRAYS
  DO iz=izs,ize
    izm1 = iz -1
!
    DO iy=iys,iye
      DO ix=1,nnx
        fnt1(ix,iy,iz) = u(ix,iy,iz)*gami + r1(ix,iy,iz)
      END DO
    END DO
!
    CALL xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
!
    IF (iz .EQ. 1) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          p(ix,iy,iz) = fnt1(ix,iy,iz) + ((w(ix,iy,iz) - wbc(ix,iy,2))*gami &
                      + r3(ix,iy,iz))*dzw_i(iz)
        END DO
      END DO
!
    ELSE IF (iz .EQ. nnz) THEN
!
      DO iy=iys,iye
        DO ix=1,nnx
          p(ix,iy,iz) = fnt1(ix,iy,iz) + ((wbc(ix,iy,1) - w(ix,iy,izm1))*gami &
                      - r3(ix,iy,izm1))*dzw_i(iz)
        END DO
      END DO
!
    ELSE
!
      DO iy=iys,iye
        DO ix=1,nnx
          p(ix,iy,iz) = fnt1(ix,iy,iz) + ((w(ix,iy,iz)  - w(ix,iy,izm1))*gami &
                      + r3(ix,iy,iz) - r3(ix,iy,izm1))*dzw_i(iz)
        END DO
      END DO
    END IF
  END DO
!
! CHECK FOR RADIATION BOUNDARY CONDITION, ALL PROCESSORS
  IF (ibcu .EQ. 1) THEN
    DO iy=iys,iye
      DO ix=1,nnx
        ptop(ix,iy,1) = pbc(ix,iy,1)
        ptop(ix,iy,2) = pbc2(ix,iy,1)
      END DO
    END DO
  END IF
!
! Y CONTRIBTION
  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        fnt1(ix,iy,iz) = v(ix,iy,iz)*gami + r2(ix,iy,iz)
      END DO
    END DO
  END DO
!
  CALL yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e,  &
              iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
!
  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        p(ix,iy,iz) = p(ix,iy,iz) + fnt1(ix,iy,iz)
      END DO
    END DO
  END DO
!
  CALL pressure
!
RETURN
END SUBROUTINE comp_p
