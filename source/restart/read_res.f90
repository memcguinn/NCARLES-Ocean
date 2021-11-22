! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine reads restart file, including constant file changed  !
!         for iys:iye.                                                         !
! ============================================================================ !
!
SUBROUTINE read_res
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
    USE inputs
!
    INCLUDE 'mpif.h'
!
    INTEGER                       :: status(mpi_status_size), ierr
    INTEGER(kind=mpi_offset_kind) :: offset, disp
    INTEGER(kind=k8)              :: nsize, nsize2
!
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp
!
! --------------------------------------------------------------------------- !
!
  ALLOCATE (temp(nvar,nnx,iys:iye))
!
! OPEN FILE
  CALL mpi_file_open(mpi_comm_world, path_res,mpi_mode_create+mpi_mode_rdwr, &
                   mpi_info_null, nvel, ierr)
!
! SET FILE VIEW
  disp = 0
  CALL mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,'native',             &
                     mpi_info_null,ierr)
!
! READ 3D FIELDS
  nsize  = int(nvar,k8)*nnx*nny
  nsize2 = int(nvar,k8)*nnx*(iys-1)
  n_read = nvar*nnx*(iye+1-iys)
!
  DO k=izs,ize
    offset = int((k-1),k8)*nsize + nsize2
    CALL mpi_file_read_at_all(nvel,offset,temp,n_read,mpi_real8,status,ierr)
!
    IF (ierr .NE. 0) THEN
      WRITE (6,6100) nvel,iz
      CALL mpi_finalize(ierr)
      STOP
    END IF
!
    DO j=iys,iye
      DO i=1,nnx
        u(i,j,k) = temp(1,i,j)
        v(i,j,k) = temp(2,i,j)
        w(i,j,k) = temp(3,i,j)
        e(i,j,k) = temp(nvar,i,j)
      END DO
    END DO
!
    DO is = 1,nscl
      DO j = iys,iye
        DO i = 1,nnx
          t(i,j,is,k) = temp(3+is,i,j)
        END DO
      END DO
    END DO
!
  END DO
!
! CLOSE FILE
  CALL mpi_file_close(nvel, ierr)
  DEALLOCATE(temp)
!
! EVERY MPI READS CONSTANT FILE
  REWIND(nvelc)
  READ (nvelc,err=9993) c_c, c_hurr, c_s, case
  IF (iti .LE. 1) THEN
    time = 0.1
    dt = 0.1
  END IF
  CLOSE(nvelc)
!
  IF (l_root) WRITE (6,4001) case
!
! ---------------------------- RESTART CONDITIONS ---------------------------- !
! SET CASE NAME TO CASE INPUT
  case   = case_inp
  IF (l_root) WRITE (6,4002) case_inp, utau, utausv
!
! IF NEW VIS MODEL, SET OUTER GRID MATCH POINT
  nmatch = nnz
  utau = utausv
!
  IF (l_root) WRITE (6,4012) time
  IF (l_root) WRITE (6,4013) qstar(1) , nmatch, case
  CALL get_dz
!
RETURN
!
9993  CONTINUE
      WRITE (6,6200) nvelc
      CALL mpi_finalize(ierr)
      STOP
!
! ---------------------------------- FORMAT ---------------------------------- !
4001 FORMAT (' 4001, SR. RESTART: case from restart = ',a3)
4002 FORMAT (' 4002, SR. RESTART:',/,                                     &
             ' files will be saved with case name = ',a3,/,               &
             ' utau = ',e15.6,' utausv = ',e15.6)
4012 FORMAT (' SR. RESTART: restart completed at T=',e15.6)
4013 FORMAT ('    after restart qstar = ',e15.6,' nmatch = ',i5,' case = ',a3)
6100 FORMAT (' SR. READ_RES: error reading file on unit number = ',i2,/,  &
             '               at iz = ',i4)
6200 FORMAT (' SR. READ_RES:',/,                                          &
             '    error reading constant file on unit number = ',i2)
!------------------------------------------------------------------------------!
!
END SUBROUTINE read_res
