! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine sets the lower boundary condition for the entire     !
!         plane at iz = 1. You can use Businger or generalized formula         !
!         considering wind.                                                    !
!                                                                              !
!         Note that the index f(.,.,2) indicates the lower boundary.           !
! ============================================================================ !
!
SUBROUTINE lower(it)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    REAL :: sfc_flx(2+nscl)
    REAL :: u_level1(nnx,iys:iye,2+nscl), buf(2+2*nscl)
    REAL :: sbuf(2+2*nscl,mxs:mxe,iys:iye)
    REAL :: rbuf((2+2*nscl)*nnx*(iye+1-iys))
!
! ---------------------------------------------------------------------------- !
!
  IF (iss .EQ. 0 .AND. ifree .EQ. 0) THEN
    iz   = 1
    izm1 = iz - 1
    dz_i = dzu_i(1)
    DO iy=iys,iye
      DO ix=1,nnx
        ebc(ix,iy,2)  = AMAX1(e(ix,iy,iz),sml_eg)
        wbc(ix,iy,2)  = 0.0
        pbc(ix,iy,2)  = 0.0
        pbc2(ix,iy,2) = 0.0
      END DO
    END DO
!
    CALL sufto(it)
    DO iy=iys,iye
      DO ix=1,nnx
        tau13m(ix,iy) = -au13m
        tau23m(ix,iy) = -au23m
      END DO
    END DO
!
    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          taut3m(ix,iy,iscl) = aut3m(iscl)
        END DO
      END DO
    END DO
!
!   PARTIAL SUMS OF SURFACE FLUXES AND MEAN SCALAR
    sfc_flx(1) = 0.0
    sfc_flx(2) = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        sfc_flx(1) = sfc_flx(1) + tau13m(ix,iy)
        sfc_flx(2) = sfc_flx(2) + tau23m(ix,iy)
      END DO
    END DO
!
    DO iscl=1,nscl
      sfc_flx(2+iscl) = 0.0
      DO iy=iys,iye
        DO ix=1,nnx
          sfc_flx(2+iscl) = sfc_flx(2+iscl) + taut3m(ix,iy,iscl)
        END DO
      END DO
    END DO
!
    CALL mpi_sum_xy(sfc_flx,myid,iss,ise,(2+nscl))
    uwsfc = sfc_flx(1)*fnxy
    vwsfc = sfc_flx(2)*fnxy
    DO iscl=1,nscl
      wtsfc(iscl) = sfc_flx(2+iscl)*fnxy
    END DO
!
    DO iy=iys,iye
      DO ix=1,nnx
        dudz     = 2.*(u(ix,iy,iz) + ugal)*dz_i
        dvdz     = 2.*v(ix,iy,iz)*dz_i
        ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
        vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
      END DO
    END DO
!
    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
          tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
        END DO
      END DO
    END DO
!
!   INITIALIZE U,V,W,T AND DERIVATIVES AT IZM1
    DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,izm1)  = ubc(ix,iy,2)
        v(ix,iy,izm1)  = vbc(ix,iy,2)
        w(ix,iy,izm1)  = wbc(ix,iy,2)
        r3(ix,iy,izm1) = 0.0
        e(ix,iy,izm1)  = ebc(ix,iy,2)
        ux(ix,iy,izm1) = 0.0
        uy(ix,iy,izm1) = 0.0
        vx(ix,iy,izm1) = 0.0
        vy(ix,iy,izm1) = 0.0
        wx(ix,iy,izm1) = wbc(ix,iy,2)
        wy(ix,iy,izm1) = wbc(ix,iy,2)
      END DO
    END DO
!
!  NO NEED TO CALL DERIVATIVES SINCE WBC = 0, CHANGE TO GENERALIZE BC
    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
        END DO
      END DO
    END DO
!
!
ELSE IF (ifree .EQ. 1) THEN
!
!
!   SEND LEVEL 1 DATA OUT
    IF (iss .eq. 0) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          u_level1(ix,iy,1) = u(ix,iy,1)
          u_level1(ix,iy,2) = v(ix,iy,1)
        END DO
      END DO
      DO iscl=1,nscl
        DO iy=iys,iye
          DO ix=1,nnx
            u_level1(ix,iy,2+iscl) = t(ix,iy,iscl,1)
          END DO
        END DO
      END DO
    END IF
    num = nnx*(iye + 1 - iys)*(2+nscl)
!
!   SEND ROOT DATA TO OTHER PROCESSORS
    CALL mpi_send_root(u_level1(1,iys,1),num,myid,numprocs,ncpu_s)
!
!   TASKS ASSIGNED TO RESPECTIVE FLUXES AND SURFACE SCALARS
    CALL suft2(u_level1,it)
!
!   RETURN SURFACE SCALARS AND MOMENT FLUXES TO ROOT
    IF (numprocs .EQ. 1) THEN
      CONTINUE
    ELSE
      DO iy=iys,iye
        DO ix=mxs,mxe
          sbuf(1,ix,iy)  = tau13m(ix,iy)
          sbuf(2,ix,iy)  = tau23m(ix,iy)
        END DO
      END DO
      DO iscl=1,nscl
        DO iy=iys,iye
          DO ix=mxs,mxe
            sbuf(2+iscl,ix,iy)      = taut3m(ix,iy,iscl)
            sbuf(2+nscl+iscl,ix,iy) = t_grnd(ix,iy,iscl)
          END DO
        END DO
      END DO
      irow_r = MOD(myid,ncpu_s)
      IF (myid .GE. ncpu_s) THEN
        num = (2+2*nscl)*(mxe+1-mxs)*(iye+1-iys)
        CALL mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1, &
                      mpi_comm_world,ierr)
      ELSE
        DO l=irow_r+ncpu_s,numprocs-1,ncpu_s
          num = (2+2*nscl)*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
          CALL mpi_recv(rbuf(1),num,mpi_real8,l,1,mpi_comm_world,istatus,ierr)
          CALL f_suft2(rbuf,nnx,mx_s(l),mx_e(l),iys,iye,nscl,    &
                       tau13m,tau23m,taut3m,t_grnd)
        END DO
      END IF
    END IF
!
! IF ROOT ROW = 0, GET SUM OF SURFACE CONDITIONS AND SET SURFACE BOUNDARY CONDTITIONS
    IF (iss .EQ. 0) THEN
      buf(1) = 0.0
      buf(2) = 0.0
      DO iy=iys,iye
        DO ix=1,nnx
          buf(1) = buf(1) + tau13m(ix,iy)
          buf(2) = buf(2) + tau23m(ix,iy)
        END DO
      END DO
      DO iscl=1,nscl
        buf(2+iscl)      = 0.
        buf(2+nscl+iscl) = 0.
        DO iy=iys,iye
          DO ix=1,nnx
            buf(2+iscl)      = buf(2+iscl) + taut3m(ix,iy,iscl)
            buf(2+nscl+iscl) = buf(2+nscl+iscl) + t_grnd(ix,iy,iscl)
          END DO
        END DO
      END DO
      CALL mpi_sum_xy(buf,myid,iss,ise,2+2*nscl)
      uwsfc = buf(1)*fnxy
      vwsfc = buf(2)*fnxy
      DO iscl=1,nscl
        wtsfc(iscl) = buf(2+iscl)*fnxy
        tsfcc(iscl) = buf(2+nscl+iscl)*fnxy
      END DO
      iz   = 1
      izm1 = iz - 1
      dz_i = dzu_i(iz)
      DO iy=iys,iye
        DO ix=1,nnx
          ebc(ix,iy,2)=AMAX1(e(ix,iy,iz),sml_eg)
          wbc(ix,iy,2)= 0.0
          pbc(ix,iy,2) = 0.0
          pbc2(ix,iy,2) = 0.0
        END DO
      END DO
      DO iy=iys,iye
        DO ix=1,nnx
          dudz     = 2.*u(ix,iy,iz)*dz_i
          dvdz     = 2.*v(ix,iy,iz)*dz_i
          ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
          vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
        END DO
      END DO
      DO iscl=1,nscl
        DO iy=iys,iye
          DO ix=1,nnx
            dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
            tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
          END DO
        END DO
      END DO
!
!     INITIALIZE U,V,W,T AND DERIVATIVES AT IZM1
      DO iy=iys,iye
        DO ix=1,nnx
          u(ix,iy,izm1)  = ubc(ix,iy,2)
          v(ix,iy,izm1)  = vbc(ix,iy,2)
          w(ix,iy,izm1)  = wbc(ix,iy,2)
          r3(ix,iy,izm1) =  0.0
          e(ix,iy,izm1)  = ebc(ix,iy,2)
          ux(ix,iy,izm1) = 0.0
          uy(ix,iy,izm1) = 0.0
          vx(ix,iy,izm1) = 0.0
          vy(ix,iy,izm1) = 0.0
          wx(ix,iy,izm1) = wbc(ix,iy,2)
          wy(ix,iy,izm1) = wbc(ix,iy,2)
        END DO
      END DO
      DO iscl=1,nscl
        DO iy=iys,iye
          DO ix=1,nnx
            t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
          END DO
        END DO
      END DO
    ENDIF
!
  END IF
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
2345 FORMAT(' in lower 2345 uwsfc = ',e15.6,' vwsfc = ',e15.6,  &
             ' wtsfc = ',e15.6,' tsfcc = ',e15.6)
!------------------------------------------------------------------------------!
!
END SUBROUTINE lower
