! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine gets the sums on the root or all processors for all  !
!         z in vector f(i).                                                    !
! ============================================================================ !
!
SUBROUTINE mpi_sum_z(f,i_root,myid,nsend,iall)
!
    INCLUDE 'mpif.h'
!
    INTEGER :: istatus(mpi_status_size)
    REAL :: recv_b(nsend), f(nsend)
!
! --------------------------------------------------------------------------- !
!
! ROOT GETS RESULT
  IF (iall .NE. 1) THEN
!
    CALL mpi_reduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,i_root,  &
                    mpi_comm_world,ierr)
!
    IF (myid .EQ. i_root) THEN
      DO i=1,nsend
        f(i) = recv_b(i)
      END DO
    END IF
  ELSE
!
!   RESULT SENT OUT
    CALL mpi_allreduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,      &
                       mpi_comm_world,ierr)
!
    DO i=1,nsend
      f(i) = recv_b(i)
    END DO
  END IF
!
RETURN
END SUBROUTINE mpi_sum_z
