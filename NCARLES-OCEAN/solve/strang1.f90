! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine splits the scalar reaction term:                     !
!           0.5 * react, advect, 0.5 * react for fast reactions (tau LEQ 1000).!
!         See RHS/RHS_SCL for slow reaction sources.                           !
!                                                                              !
!         The temp scalar arrays holds RHS from step n-1 and field variables   !
!         from step n.                                                         !
! ============================================================================ !
!
SUBROUTINE strang1(it)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
    USE reaction, ONLY: react_src
!
    INCLUDE 'mpif.h'
!
    REAL :: trhs(nnx,iys:iye,nscl,izs:ize)
    REAL,DIMENSION(nscl-1) :: tmp
!
! --------------------------------------------------------------------------- !
!
  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        tmp = react_src(ix,iy,1,iz)
        DO l=2,nscl
          trhs(ix,iy,l,iz) = tmp(l-1)
          IF (trhs(ix,iy,l,iz) .LE. 1.0e-20) THEN
            trhs(ix,iy,l,iz) = 1.0e-20
          END IF
        END DO
      END DO
    END DO
  END DO
!
  DO iz=izs,ize
    DO l=2,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
        END DO
      END DO
    END DO
  END DO
!
RETURN
END SUBROUTINE strang1
