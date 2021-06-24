! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine builds arrays for profiles.                          !
! ============================================================================ !
!
SUBROUTINE write_prof(nhisp,krec,num,f)
!
    INTEGER :: kode
!
    REAL   :: f(num)
    REAL*4 :: f32(num)
!
! --------------------------------------------------------------------------- !
!
! BUILD 32 BIT ARRAYS FOR PROFILES
  DO i=1,num
    f32(i) = f(i)
  END DO
!
  WRITE (nhisp,err=999,iostat=kode,rec=krec) (f32(i),i=1,num)
!
RETURN
!
! ---------------------------------- ERRORS ---------------------------------- !
999 CONTINUE
    WRITE (6,9000) num,krec,kode
    STOP
!
! ---------------------------------- FORMAT ---------------------------------- !
9000 FORMAT (' 9000, trouble in ','SR. save_prof cannot write profile data ', &
             /,' num = ',i8, ' krec = ',i6, ' kode = ',i6)
! ---------------------------------------------------------------------------- !
!
END SUBROUTINE write_prof
