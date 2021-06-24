! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine computes spacing for given vertical point            !
!         distribution.                                                        !
! ============================================================================ !
!
SUBROUTINE get_dz
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
! --------------------------------------------------------------------------- !
!
  DO iz=1,nnz+1
    dzw(iz) = z(iz) - z(iz-1)
  END DO
!
  dzw(0)     = dzw(1)
  dzw(nnz+2) = dzw(nnz+1)
!
  DO iz=0,nnz+2
    dzw_i(iz) = 1.0/dzw(iz)
  END DO
!
! BUILD Z GRID FOR U
  dzovr2 = dz*0.5
!
  DO iz=1,nnz+1
    zz(iz) = 0.5*(z(iz) + z(iz-1))
  END DO
!
  zz(0) = - zz(1)
!
  DO iz=1,nnz+1
    dzu(iz) = zz(iz) - zz(iz-1)
  END DO
!
  dzu(0)     = dzu(1)
  dzu(nnz+2) = dzu(nnz+1)
!
  DO iz=0,nnz+2
    dzu_i(iz) = 1.0/dzu(iz)
  END DO
!
RETURN
END SUBROUTINE get_dz
