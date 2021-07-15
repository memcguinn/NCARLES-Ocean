! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine sets the file path for restart, save directories,    !
!         history files.                                                       !
! ============================================================================ !
!
SUBROUTINE set_paths
!
    USE pars
!
! --------------------------------------------------------------------------- !
!
! PATH AND FILE NAME OF VELOCITY RESTART FILE
  path_res = './data/u.mp.30L00020'             ! Change file number to fit iti for restart
!
! SAVE PATH FOR 3D VOLUMES
  path_sav = './data'
!
! SAVE PATH FOR NEW HISTORY FILES
  path_his = './data'
!
RETURN
END SUBROUTINE set_paths
