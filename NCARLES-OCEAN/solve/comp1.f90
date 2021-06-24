! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine descrives the Or(3) Runge-Kutta time-step scheme     !
!         and monotone scalar fluxes in x,y,z. It is designed to use MPI in    !
!         the x and y directions.                                              !
! ============================================================================ !
!
SUBROUTINE comp1(istep,it)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    PARAMETER (js = 7, ns = 3, nstat = js + ns*nscl)
!
    INTEGER  :: istatus(mpi_status_size)
!
    REAL :: stat(1:nnz,nstat)
!
!   DEFINE TEMP ARRAYS TO HOLD RHS FROM STEP N-1 AND FIELD VARIABLES FROM STEP N
    REAL :: urhs(nnx,iys:iye,izs:ize), vrhs(nnx,iys:iye,izs:ize),  &
            wrhs(nnx,iys:iye,izs:ize), erhs(nnx,iys:iye,izs:ize),  &
            trhs(nnx,iys:iye,nscl,izs:ize)
!
! --------------------------------------------------------------------------- !
!
  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        urhs(ix,iy,iz) = u(ix,iy,iz) + dtzeta*r1(ix,iy,iz)
        vrhs(ix,iy,iz) = v(ix,iy,iz) + dtzeta*r2(ix,iy,iz)
        wrhs(ix,iy,iz) = w(ix,iy,iz) + dtzeta*r3(ix,iy,iz)
        erhs(ix,iy,iz) = e(ix,iy,iz) + dtzeta*r5(ix,iy,iz)
      END DO
    END DO
  END DO
!
  DO iz=izs,ize
    DO l=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          trhs(ix,iy,l,iz) = t(ix,iy,l,iz) + dtzeta*r4(ix,iy,l,iz)
        END DO
      END DO
    END DO
  END DO
!
! GET VISCOSITY AND RHS OF (E,U,V,W)
  CALL tke_vis(istep)
  CALL rhs_uvw(istep)
!
! EVALUATE RHS OF SCALAR EQUATIONS
  DO l=1,nscl
    CALL rhs_scl(istep,l,it)
  END DO
!
! GATHER STAT SUMS ON ROOT PROCESSOR, USE MPI_REDUCTION OVER ALL PROCESSORS
  IF (istep .EQ. 1) THEN
    DO j=1,nstat
      DO iz=1,nnz
        stat(iz,j) = 0.0
      END DO
    END DO
!
    DO iz=izs,ize
      stat(iz,1) = uwsb(iz)
      stat(iz,2) = vwsb(iz)
      stat(iz,3) = wwsb(iz)
      stat(iz,4) = tr_tau(iz)
      stat(iz,5) = triz(iz)
      stat(iz,6) = shrz(iz)
      stat(iz,7) = t_diss(iz)
    END DO
!
    m1 = js
    m2 = js + nscl
    m3 = js + 2*nscl
!
    DO l=1,nscl
      DO iz=izs,ize
        stat(iz,m1+l) = utsb(iz,l)
        stat(iz,m2+l) = vtsb(iz,l)
        stat(iz,m3+l) = wtsb(iz,l)
      END DO
    END DO
!
    CALL mpi_sum_z(stat(1,1),i_root,myid,nstat*nnz,1)
!
    DO iz=1,nnz
      uwsb(iz)   = stat(iz,1)
      vwsb(iz)   = stat(iz,2)
      wwsb(iz)   = stat(iz,3)
      tr_tau(iz) = stat(iz,4)
      triz(iz)   = stat(iz,5)
      shrz(iz)   = stat(iz,6)
      t_diss(iz) = stat(iz,7)
    END DO
!
    DO l=1,nscl
      DO iz=1,nnz
        utsb(iz,l) = stat(iz,m1+l)
        vtsb(iz,l) = stat(iz,m2+l)
        wtsb(iz,l) = stat(iz,m3+l)
      END DO
    END DO
!
    DO iz=1,nnz
      buyz(iz) = batag*wtsb(iz,1)
    END DO
!
  END IF
!
! SAVE OLD RHS IN FIELD VARIABLES FOR R-K ADVANCEMENT
  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,iz) = urhs(ix,iy,iz)
        v(ix,iy,iz) = vrhs(ix,iy,iz)
        w(ix,iy,iz) = wrhs(ix,iy,iz)
        e(ix,iy,iz) = erhs(ix,iy,iz)
      END DO
    END DO
  END DO
!
  DO iz=izs,ize
    DO l=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
        END DO
      END DO
    END DO
  END DO
!
RETURN
END SUBROUTINE comp1
