! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine tells the root processes to write constant file      !
!         sequential Fortran binary.                                           !
! ============================================================================ !
!
SUBROUTINE save_c(it)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    CHARACTER :: options*8, passwd*1

    LOGICAL :: there
!
! --------------------------------------------------------------------------- !
!
! OPEN, WRITE, CLOSE FILE
  OPEN(nvelc,err=9992,file=path_sav_c,form='unformatted',status='unknown')
  WRITE (nvelc,err=9992) c_c, c_hurr, c_s, case
  CLOSE(nvelc)
!
  INQUIRE (file=path_sav_c,exist=there)
!
  IF (.NOT. there) THEN
    WRITE (6,8001) path_sav_c
    CALL mpi_finalize(ierr)
    STOP
  END IF
!
! CHECK: OUTPUT MESSAGE
  WRITE (6,7001) path_sav_c
!
RETURN
!
! ---------------------------------- ERRORS ---------------------------------- !
9992 CONTINUE
     WRITE (6,6100) nvelc
     CALL mpi_finalize(ierr)
     STOP
!
! ---------------------------------- FORMAT ---------------------------------- !
6100 FORMAT (' SR. SAVE_V:',/,  &
             '    trouble cannot open/write constant file on unit = ',i2)
7001 FORMAT ('      CONSTANT DATA IS WRITTEN IN FILE  ',a80)
8001 FORMAT (' SR. SAVE_C: Trouble constant file not in path =',a80)
!------------------------------------------------------------------------------!
!
END SUBROUTINE save_c
