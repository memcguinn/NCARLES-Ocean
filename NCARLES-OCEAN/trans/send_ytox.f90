! ============================================================================ !
! ABOUT:                                                                       !
!          Grab and send section of array, skipping ghost points.              !
! ============================================================================ !
!
SUBROUTINE send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
!
    REAL g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)
!
! --------------------------------------------------------------------------- !
!
  DO k=izs,ize
    DO i=ixs,ixe
      DO j=iys,iye
        gt(j,i,k) = g(j,i,k)
      END DO
    END DO
  END DO

RETURN
END SUBROUTINE send_ytox
