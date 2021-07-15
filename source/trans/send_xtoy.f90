! ============================================================================ !
! ABOUT:                                                                       !
!          Grab and send section of array.                                     !
! ============================================================================ !
!
SUBROUTINE send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
!
    REAL f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)
!
! --------------------------------------------------------------------------- !
!
  DO k=izs,ize
    DO j=iys,iye
      DO i=ixs,ixe
        ft(i,j,k) = f(i,j,k)
      END DO
    END DO
  END DO

RETURN
END SUBROUTINE send_xtoy
