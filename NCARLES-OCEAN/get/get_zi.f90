! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE get_zi(gradmax,gradout,len,itype)
!
    USE pars
!
    REAL :: gradMAX(*), gradout(*)
!
! --------------------------------------------------------------------------- !
!
  DO i=1,len,2
    IF (gradMAX(i) .GT. gradout(i)) THEN
      gradout(i)   = gradMAX(i)
      gradout(i+1) = gradMAX(i+1)
    END IF
  END DO
!
RETURN
END SUBROUTINE get_zi
