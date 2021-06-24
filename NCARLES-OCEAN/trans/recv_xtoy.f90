! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
!
    REAL g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)
!
! --------------------------------------------------------------------------- !
!
  DO k=izs,ize
    DO j=iys,iye
      DO i=ixs,ixe
        g(j,i,k) = gt(i,j,k)
      END DO
    END DO
  END DO

RETURN
END SUBROUTINE recv_xtoy
