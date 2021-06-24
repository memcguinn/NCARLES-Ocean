! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine exchanges ghost points with mpi. Note that nb and nt !
!         are the destination and source nodes. Allows for 1z per CPU.         !
! ============================================================================ !
!
SUBROUTINE exchange
!
    USE pars
    USE fields
!
    INCLUDE 'mpif.h'
!
    INTEGER :: istatus(mpi_status_size)
!
    REAL :: fs(nnx,iys:iye,(4+nscl)),fr(nnx,iys:iye,(4+nscl))
!
! --------------------------------------------------------------------------- !
!
  nb = myid - ncpu_s
  nt = myid + ncpu_s
!
! ACCOUNT FOR ENDPOINTS
  IF (iss .EQ. 0) THEN
    nb = mpi_proc_null
  END IF
!
  IF (ise .EQ. numprocs-1) THEN
    nt = mpi_proc_null
  END IF
!
  nsend = nnx*(iye + 1 - iys)*(4+nscl)
  nrecv = nsend
!
! SEND TOP OF MYID, RECEIVE BOTTOM OF MYID - NCPU_S
  DO iy=iys,iye
    DO ix=1,nnx
      fs(ix,iy,1) = u(ix,iy,ize)
      fs(ix,iy,2) = v(ix,iy,ize)
      fs(ix,iy,3) = w(ix,iy,ize)
      fs(ix,iy,4) = e(ix,iy,ize)
    END DO
  END DO
!
  DO iscl=1,nscl
    jloc = 4 + iscl
    DO iy=iys,iye
      DO ix=1,nnx
        fs(ix,iy,jloc) = t(ix,iy,iscl,ize)
      END DO
    END DO
  END DO
!
  CALL mpi_sendrecv(                           &
           fs(1,iys,1),nsend,mpi_real8,nt,0,   &
           fr(1,iys,1),nrecv,mpi_real8,nb,0,   &
           mpi_comm_world,istatus,ierr)
!
  IF (iss .NE. 0) THEN
    izm1 = izs-1
    DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,izm1) = fr(ix,iy,1)
        v(ix,iy,izm1) = fr(ix,iy,2)
        w(ix,iy,izm1) = fr(ix,iy,3)
        e(ix,iy,izm1) = fr(ix,iy,4)
      END DO
    END DO
!
    DO iscl=1,nscl
      jloc = 4 + iscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,iscl,izm1) = fr(ix,iy,jloc)
        END DO
      END DO
!
    END DO
  END IF
!
! SEND BOTTOM OF MYID, RECIEVE BOTTOM FROM MYID + NCPU_S
  DO iy=iys,iye
    DO ix=1,nnx
      fs(ix,iy,1) = u(ix,iy,izs)
      fs(ix,iy,2) = v(ix,iy,izs)
      fs(ix,iy,3) = w(ix,iy,izs)
      fs(ix,iy,4) = e(ix,iy,izs)
    END DO
  END DO
!
  DO iscl=1,nscl
    jloc = 4 + iscl
    DO iy=iys,iye
      DO ix=1,nnx
        fs(ix,iy,jloc) = t(ix,iy,iscl,izs)
      END DO
    END DO
  END DO
!
  CALL mpi_sendrecv(                           &
           fs(1,iys,1),nsend,mpi_real8,nb,1,   &
           fr(1,iys,1),nrecv,mpi_real8,nt,1,   &
           mpi_comm_world,istatus,ierr)
!
  IF (ise .NE. numprocs-1) THEN
    izp1 = ize+1
    DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,izp1) = fr(ix,iy,1)
        v(ix,iy,izp1) = fr(ix,iy,2)
        w(ix,iy,izp1) = fr(ix,iy,3)
        e(ix,iy,izp1) = fr(ix,iy,4)
      END DO
    END DO
!
    DO iscl=1,nscl
      jloc = 4 + iscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,iscl,izp1) = fr(ix,iy,jloc)
        END DO
      END DO
!
    END DO
  END IF
!
! SEND EXTRA SCALAR POINTS
  nsend = nnx*(iye + 1 - iys)*nscl
  nrecv = nsend
!
! SEND TOP OF MYID, RECEIVE BOTTOM FROM MYID - NCPU_S
  izm1 = ize-1
  DO iscl=1,nscl
    DO iy=iys,iye
      DO ix=1,nnx
        fs(ix,iy,iscl) = t(ix,iy,iscl,izm1)
      END DO
    END DO
  END DO
!
  CALL mpi_sendrecv(                           &
           fs(1,iys,1),nsend,mpi_real8,nt,0,   &
           fr(1,iys,1),nrecv,mpi_real8,nb,0,   &
           mpi_comm_world,istatus,ierr)
!
  IF (iss .NE. 0) THEN
    izm2 = izs-2
    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,iscl,izm2) = fr(ix,iy,iscl)
        END DO
      END DO
    END DO
  END IF
!
! SEND BOTTOM OF MYID, RECEIVE BTTOM FROM MYID + NCPU_S
  izp1 = izs+1
  DO iscl=1,nscl
    DO iy=iys,iye
      DO ix=1,nnx
        fs(ix,iy,iscl) = t(ix,iy,iscl,izp1)
      END DO
    END DO
  END DO
!
  CALL mpi_sendrecv(                           &
           fs(1,iys,1),nsend,mpi_real8,nb,1,   &
           fr(1,iys,1),nrecv,mpi_real8,nt,1,   &
           mpi_comm_world,istatus,ierr)
!
  IF (ise .NE. numprocs-1) THEN
    izp2 = ize+2
    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,iscl,izp2) = fr(ix,iy,iscl)
        END DO
      END DO
    END DO
  END IF
!
RETURN
END SUBROUTINE exchange
