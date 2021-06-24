! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine opens history files by root isize determined in      !
!         con_stats.                                                           !
! ============================================================================ !
!
SUBROUTINE open_his(istep)
!
    USE pars
!
    INCLUDE 'mpif.h'
!
    CHARACTER cgrid*4, iblks*16
!
    LOGICAL there
!
! --------------------------------------------------------------------------- !
!
! BUILD CHARACTER STRINGS FOR ASCII HISTORY FILE NAMES
  cgrid = '.mp.'
!
  CALL blnk(iblks)
  WRITE (iblks(1:6),'(i6.6)') istep
!
  iblks(7:7) = '_'
!
  WRITE (iblks(8:13),'(i6.6)') (istep + itape)
!
  iblnk = index(path_his,' ')
!
  CALL blnk(path_sav_h)
!
  path_sav_h = path_his(1:iblnk-1)//'/his'//                 &
               cgrid(1:4)//case(1:3)//'.'//iblks(1:13)
!
! BUILD CHARACTER STRINGS FOR IEEE PROFILE HISTORY FILE SET RECORD COUNTER FOR
! DIRECT ACCESS FILE = 0
  krec = 0
  cgrid = '.mp.'
!
  CALL blnk(iblks)
  WRITE (iblks(1:6),'(i6.6)') istep
!
  iblks(7:7) = '_'
!
  WRITE (iblks(8:13),'(i6.6)') (istep + itape)
!
  iblnk = index(path_his,' ')
!
  CALL blnk(path_sav_hp)
!
  path_sav_hp = path_his(1:iblnk-1)//'/his'//                &
                cgrid(1:4)//case(1:3)//'.'//iblks(1:13)//'.ieee'
!
! SAVE DATA IN DIRECTORY
  IF (l_root) THEN
    CLOSE(nhis1)
    OPEN(nhis1,err=3000,file=path_sav_h,form='formatted')
    CLOSE(nhisp)
    OPEN(nhisp,err=4000,file=path_sav_hp,form='unformatted',   &
       access='direct',recl=isize*j_recl,status='unknown')
  END IF
!
RETURN
!
! ---------------------------------- ERRORS ---------------------------------- !
3000 CONTINUE
     WRITE (6,6301) nhis1, path_sav_h
     STOP
!
4000 CONTINUE
     WRITE (6,6302) nhisp, path_sav_hp
     STOP
!
! ---------------------------------- FORMAT ---------------------------------- !
6301 FORMAT (' 6301, SR. OPEN_HIS:',/,      &
       '    cannot open history1 file on unit = ',i2,/,'    path = ',a80)
6302 FORMAT (' 6302, SR. OPEN_HIS:',/,      &
       '    cannot open history profile file on unit = ',i2,/,'    path = ',a80)
!------------------------------------------------------------------------------!
!
END SUBROUTINE open_his
