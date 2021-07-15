! ============================================================================ !
! ABOUT:                                                                       !
!          Grab and send section of array, accounting for ghost points.        !
! ============================================================================ !
!
SUBROUTINE send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
!
    REAL g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)
!
! --------------------------------------------------------------------------- !
!
  DO j=iys,iye
    DO i=ixs,ixe
      DO k=izs-1,ize+1
         gt(k,j,i) = g(k,j,i)
      END DO
    END DO
  END DO

RETURN
END SUBROUTINE send_ztox
