! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine uses Schumann stability to correct length scales.    !
! ============================================================================ !
!
SUBROUTINE schu_vis(alk)
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
!
    REAL :: alk(nnx,iys:iye,izs-1:ize+1)
!
! ---------------------------------------------------------------------------- !
!
  DO iz=izs-1,ize+1
    izp1 = iz + 1
    dslk  = dsl_z(iz)
    IF (iz .GT. 0) dslk  = AMIN1(dsl_z(iz),vk*ABS(z(iz))/csmag)
    almin = almin_c*dsl_z(iz)
    IF (iz .EQ. 0 .OR. iz .EQ. nnz+1) THEN
      dfack = 1.0
    ELSE
      dfack = dfac(iz)
    END IF
!
!   NO STABILITY CORRECTED LENGTH SCALES
    DO j=iys,iye
      DO i=1,nnx
        alk(i,j,iz) = dslk
      END DO
    END DO
!
!   STABILITY CORRECTED FOR VERTICAL SCALAR FLUX
    DO j=iys,iye
      DO i=1,nnx
        vis_m(i,j,iz)  = ck*dslk*SQRT(e(i,j,iz))*dfack
        vis_s(i,j,iz)  = 3.0*vis_m(i,j,iz)
        stab = AMAX1(batag*(t(i,j,1,izp1) - t(i,j,1,iz))*dzu_i(izp1),0.0)
        vis_sv(i,j,iz) = vis_s(i,j,iz)*(e(i,j,iz)/ &
                                           (e(i,j,iz) + 0.3*stab*dslk**2))
      END DO
    END DO
!
!   FOR IZ = 1
    IF (iz .EQ. 1 .AND. ibcl .EQ. 0) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          vis_m(ix,iy,iz-1)  = vis_m(ix,iy,iz)
          vis_s(ix,iy,iz-1)  = vis_s(ix,iy,iz)
          vis_sv(ix,iy,iz-1) = vis_sv(ix,iy,iz)
        END DO
      END DO
    END IF
  END DO
!
RETURN
END SUBROUTINE schu_vis
