! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE get_units
!
    USE pars
!
! FILE UNIT NUMBERS
  nvel  = 20
  npre  = 30
  nhis1 = 40
  nvelc = 50
  nhisp = 60
!
! UNIT NUMBERS FOR STANDARD PRINT OUT FOR MPI TASKS
  nprt = 1
!
! OPEN UNIT FOR STNADARD PRINTOUT
  path_prt = case_inp(1:3)//'.mp.xxxxx.out'
  WRITE (path_prt(8:12),'(i5.5)') myid
  OPEN(nprt,file=path_prt,form='formatted')
!
RETURN
END SUBROUTINE get_units
