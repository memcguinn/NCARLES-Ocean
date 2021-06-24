! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine sends the root results to above processors           !
! ============================================================================ !
!
SUBROUTINE mpi_send_root(fs,num,myid,np,ns)
!
    INCLUDE 'mpif.h'
!
    INTEGER :: istatus(mpi_status_size)
!
    REAL :: fs(num)
!
! --------------------------------------------------------------------------- !
!
  IF (np .EQ. 1) THEN
    CONTINUE
  ELSE
!
    irow_r = mod(myid,ns)
!
    IF (myid .GE. ns) THEN
      CALL mpi_recv(fs(1),num,mpi_real8,irow_r,1,mpi_comm_world,istatus,ierr)
    ELSE
!
      DO l=irow_r+ns,np-1,ns
        CALL mpi_send(fs(1),num,mpi_real8,l,1,mpi_comm_world,ierr)
      END DO
!
    END IF
  END IF
!
RETURN
END SUBROUTINE mpi_send_root
