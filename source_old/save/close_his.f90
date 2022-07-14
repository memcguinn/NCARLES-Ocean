! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine closes history files.                                !
! ============================================================================ !
!
SUBROUTINE close_his
!
    USE pars
!
    LOGICAL there
!
! --------------------------------------------------------------------------- !
!
! ROOT CLOSES AND CHECKS FILES 
  CLOSE(nhis1)
  CLOSE(nhisp)
!
  INQUIRE (file=path_sav_h,exist=there)
!
  IF (.NOT. there) THEN
    WRITE (6,8000) path_sav_h
    CALL mpi_finalize(ierr)
    STOP
  END IF
!
  INQUIRE (file=path_sav_hp,exist=there)
!
  IF (.NOT. there) THEN
    WRITE (6,8100) path_sav_hp
    CALL mpi_finalize(ierr)
    STOP
  END IF
!
  WRITE (6,7000) path_sav_h
  WRITE (6,7100) path_sav_hp
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
7000 FORMAT (' HISTORY DATA IS WRITTEN IN FILE  ',a80)
7100 FORMAT (' PROFILE HISTORY DATA IS WRITTEN IN FILE  ',a80)
8000 FORMAT (' SR. SAVE_HIS: Trouble history file not in path =',a80)
8100 FORMAT (' SR. SAVE_HIS: Trouble profile history file',' not in path =',a80)
!------------------------------------------------------------------------------!
!
END SUBROUTINE close_his
