! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine saves the 3D fields.                                 !
! ============================================================================ !
!
SUBROUTINE save_v(it)
!
    USE pars
    USE fields
!
    INCLUDE 'mpif.h'
!
    INTEGER (kind=mpi_offset_kind) :: offset, disp
    INTEGER (kind=k8)              :: nsize, nsize2
    INTEGER                        :: status(mpi_status_size), ierr
!
    LOGICAL there
!
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp
!
! --------------------------------------------------------------------------- !
!
  ALLOCATE (temp(nvar,nnx,iys:iye))
!
! OPEN FILE
  CALL mpi_file_open(mpi_comm_world, path_sav_v,  &
          mpi_mode_create+mpi_mode_rdwr, mpi_info_null, nvel, ierr)
!
! SET FILE VIEW
  disp = 0
  CALL mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8, &
          'native',mpi_info_null,ierr)
!
! WRITE DATA
  nsize   = INT(nvar,k8)*nnx*nny
  nsize2  = INT(nvar,k8)*nnx*(iys-1)
  n_write = nvar*nnx*(iye+1-iys)
!
  DO k=izs,ize
    DO j = iys,iye
      DO i = 1,nnx
        temp(1,i,j)    = u(i,j,k)
        temp(2,i,j)    = v(i,j,k)
        temp(3,i,j)    = w(i,j,k)
        temp(nvar,i,j) = e(i,j,k)
      END DO
    END DO
!
    DO is = 1,nscl
      DO j = iys,iye
        DO i = 1,nnx
          temp(3+is,i,j) = t(i,j,is,k)
        END DO
      END DO
    END DO
!
    offset = INT((k-1),k8)*nsize + nsize2
    CALL mpi_file_write_at_all(nvel,offset,temp,n_write,mpi_real8,status,ierr)
!
    IF (ierr /= 0) THEN
      WRITE (6,6000) nvel, iz
      CALL mpi_finalize(ierr)
      STOP
    END IF
  END DO
!
! CLOSE FILE
  CALL mpi_file_close(nvel, ierr)
!
! CHECK FILE
  IF (l_root) THEN
    INQUIRE (file=path_sav_v,exist=there)
    IF (.NOT. there) THEN
      WRITE (6,8000) nvel,myid
      CALL mpi_finalize(ierr)
      STOP
    END IF
    WRITE (6,7000) it,path_sav_v
  END IF
!
  DEALLOCATE(temp)
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
6000 FORMAT (' SR. SAVE_V:',/,                                         &
             '    trouble cannot write restart file on unit = ',i2,/,  &
             '             at iz = ',i4)
7000 FORMAT (' **** DATA SET AT IT = ',I6,/,                           &
             '      VELOCITY DATA IS WRITTEN IN FILE  ',a80)
8000 FORMAT (' in SAVE_V: trouble writing file ',i5,'  myid = ',i5,    &
             ' at iz = ',i5)
!------------------------------------------------------------------------------!
!
END SUBROUTINE save_v
