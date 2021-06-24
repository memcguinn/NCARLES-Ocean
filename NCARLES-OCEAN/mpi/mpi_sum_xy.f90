! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine gets the sum over a set of processors [iss:ise] for  !
!         vector f(i). f(i) is overwritten and skipped if working with a       !
!         single processor.                                                    !
! ============================================================================ !
!
SUBROUTINE mpi_sum_xy(f,myid,iss,ise,nsend)
!
    INCLUDE 'mpif.h'
!
    INTEGER :: istatus(mpi_status_size)
!
    REAL :: work(nsend,iss:ise), f(nsend)
!
! --------------------------------------------------------------------------- !
!
  IF (iss .EQ. ise) THEN
    CONTINUE
  ELSE
!
    DO j=1,nsend
      work(j,myid) = f(j)
      f(j)         = 0.0
    END DO
!
    DO i=iss,ise
      IF (i .NE. myid) THEN
        CALL mpi_sendrecv(work(1,myid),nsend,mpi_real8,i,1,   &
                    work(1,i),nsend,mpi_real8,i,1,mpi_comm_world,istatus,ierr)
      END IF
    END DO
!
    DO i=iss,ise
      DO j=1,nsend
        f(j) = f(j) + work(j,i)
      END DO
    END DO
  END IF
!
RETURN
END SUBROUTINE mpi_sum_xy
