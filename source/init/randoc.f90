! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine defines random initial conditions.                   !
! ============================================================================ !
!
SUBROUTINE randoc
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    REAL :: psi(nnx,iys:iye), psix(nnx,iys:iye), psiy(nnx,iys:iye,izs:izs), &
            uxx(nnx,iys:iye), vyy(nnx,iys:iye,izs:izs)
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER :: n,s,sjk
!
! --------------------------------------------------------------------------- !
!
  zi  = z(izi)
  tmixed = 288.15                       ! Mixed layer temperature, [K]
!
  DO iz=izs,ize
    IF (iz.le.izi) THEN
      DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,iz)   = ugcont-ugal
        v(ix,iy,iz)   = vgcont
        w(ix,iy,iz)   = 0.0
        t(ix,iy,1,iz) = tmixed
        e(ix,iy,iz)   = 0.0
      END DO
      END DO
    END IF
!
    IF (iz.gt.izi) THEN
      DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,iz)   = ugcont-ugal
        v(ix,iy,iz)   = vgcont
        w(ix,iy,iz)   = 0.0
        t(ix,iy,1,iz) = tmixed + dtdzf(1)*(zz(iz)-zi)
        e(ix,iy,iz)   = 0.0
      END DO
      END DO
    END IF
!
    DO iy=iys,iye
    DO ix=1,nnx
      w(ix,iy,iz)    = 0.0
      r1(ix,iy,iz)   = 0.0
      r2(ix,iy,iz)   = 0.0
      r3(ix,iy,iz)   = 0.0
      r4(ix,iy,1,iz) = 0.0
      r5(ix,iy,iz)   = 0.0
    END DO
    END DO
  END DO
!
! INITIALIZE RANDOM NUMBER GENERATOR
  CALL random_seed(size=n)
  ALLOCATE(seed(n))
  CALL system_clock(s)
!
  sjk = 1010
  seed(:) = ABS(MOD((sjk*181)*((myid+1-83)*359), 104729) )
!
  CALL random_seed(put=seed)
!
! SET INITIAL FIELD TO BE DIVERGENCE FREE
  idum = -1
!
  DO iz=izs,ize
!
    IF (iz.le.8) THEN
      ampv = 0.01                     ! Max amplitude of random velocity field
      ampt = 0.01                     ! Max amplitude of random temperature field
!
!     RANDOM FIELD SCLAED BETWEEN 0,1
      DO iy=iys,iye
      DO ix=1,nnx
        CALL random_number(psi(ix,iy))
      END DO
      END DO
!
      DO iy=iys,iye
      DO ix=1,nnx
        psix(ix,iy) = psi(ix,iy)
        psiy(ix,iy,izs) = psi(ix,iy)
      END DO
      END DO
!
      CALL xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      CALL yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1), &
                  nnx,nny,ixs,ixe,ix_s,ix_e,        &
                  iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
!
      vmaxx = 0.0
      DO iy=iys,iye
      DO ix=1,nnx
        vmag = SQRT(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
        IF (vmag .GT. vmaxx) vmaxx = vmag
      END DO
      END DO
!
      facv = ampv/vmaxx
      WRITE (6,600) facv, MAXVAL(psiy)
!
      DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,iz) = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
        v(ix,iy,iz) = v(ix,iy,iz) + psix(ix,iy)*facv
        t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
        e(ix,iy,iz) = 0.0001
      END DO
      END DO
    END IF
!
!   CHECK: DIVERGENCE INITIAL FIELD
    DO iy=iys,iye
    DO ix=1,nnx
      uxx(ix,iy) = u(ix,iy,iz)
      vyy(ix,iy,izs) = v(ix,iy,iz)
    END DO
    END DO
!
    CALL xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
    CALL yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1), &
                nnx,nny,ixs,ixe,ix_s,ix_e,       &
                iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
!
    DO iy=iys,iye
    DO ix=1,nnx
      divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
    END DO
    END DO
!
    divz(iz) = divz(iz)*fnxy
!
  END DO
!
  CALL mpi_sum_z(divz(1),i_root,myid,nnz,1)
!
  WRITE (nprt,6000)
  WRITE (nprt,6100) (iz,divz(iz),iz=izs,ize)
!
  DO iz=izs,ize
    ug(iz)=ugcont
    vg(iz)=vgcont
  END DO
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
  600 FORMAT('facv = ',e15.8,' max psiy = ',e15.8)
 6000 FORMAT(' check of divergence for initial state',/,' iz ',5x,' divergence')
 6100 FORMAT(i5,e15.6)
!------------------------------------------------------------------------------!
!
END SUBROUTINE randoc
