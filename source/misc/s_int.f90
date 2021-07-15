! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine solves for the integral using a mid-point rule,      !
!         from numerical recipes.                                              !
! ============================================================================ !
!
SUBROUTINE s_int(r_min,r_max,value)
!
! --------------------------------------------------------------------------- !
!
  iter = 10
  value = 0
!
  DO j=1,iter
    CALL midpnt(r_min,r_max,value,j)
  END DO
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
 1000 format(' j = ',i5,' value = ',e15.6)
!------------------------------------------------------------------------------!
!
END SUBROUTINE s_int
