! ============================================================================ !
! ABOUT:                                                                       !
!          Grab and send section of array, skipping ghost points.              !
! ============================================================================ !
!
SUBROUTINE send_xtoz(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
!
    REAL f(nx,iys:iye,izs-1:ize+1), ft(ixs:ixe,iys:iye,izs:ize)
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
END SUBROUTINE send_xtoz
