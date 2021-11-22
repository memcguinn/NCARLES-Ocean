! ============================================================================ !
! ABOUT:                                                                       !
!         inputs contains the parameters most likely to be changed during the  !
!         initialization and runs of NCARLES-Ocean.                            !
! ============================================================================ !
!
MODULE inputs
!
CONTAINS
  SUBROUTINE input_params
    INTEGER ::    itn, iti
! --------------------------------- RESTART ---------------------------------- !
  iti = 0     ! Start iteration for restart (default 0)
  itn = 0     ! Number of velocity restart file
!
!  path_res = './data/u.mp.30L00000'   ! Change file number to fit itn for restart




! -------------------------------- CHEMISTRY --------------------------------- !
  RETURN
  END SUBROUTINE input_params
END MODULE inputs
