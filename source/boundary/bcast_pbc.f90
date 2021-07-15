! ============================================================================ !
! ABOUT:                                                                       !
!           This subroutine sends upper boundary conditions to other           !
!           processors for fft solutions of pressure.                          !
! ============================================================================ !
!
SUBROUTINE bcast_pbc
!
    USE pars
    USE fields
!
    INCLUDE 'mpif.h'
!
    INTEGER istatus(mpi_status_size)
!
! --------------------------------------------------------------------------- !
!
  IF (numprocs .EQ. 1) THEN
    CONTINUE
  ELSE
    irow_r = MOD(myid,ncpu_s)
    irow_t = is_s(numprocs-1) + irow_r
    num = nnx*(iye+1-iys)
!
! FIND MYID ROW
    IF (iss .NE. is_s(numprocs-1)) THEN
!
!   IF MYID NOT IN TOP ROW, RECIEVE FROM TOP ROW
      CALL mpi_recv(pbc(1,iys,1),num,mpi_real8,irow_t,1,  &
              mpi_comm_world,istatus,ierr)
    ELSE
!   IF MYID IN TOP ROW, SEND BELOW
      DO l=irow_r,irow_t-ncpu_s,ncpu_s
        CALL mpi_send(pbc(1,iys,1),num,mpi_real8,l,1,mpi_comm_world,ierr)
      END DO
    END IF
!
! REPEAT ABOVE FOR NEXT VARIABLE
    IF (iss .NE. is_s(numprocs-1)) THEN
      CALL mpi_recv(pbc2(1,iys,1),num,mpi_real8,irow_t,1, &
               mpi_comm_world,istatus,ierr)
    ELSE
      DO l=irow_r,irow_t-ncpu_s,ncpu_s
        CALL mpi_send(pbc2(1,iys,1),num,mpi_real8,l,1,mpi_comm_world,ierr)
      END DO
    END IF
  END IF

RETURN
END SUBROUTINE bcast_pbc
