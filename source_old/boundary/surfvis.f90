! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE surfvis(it)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    REAL :: xkvis(nnx,iys:iye), alwk(nnx,iys:iye), send(3), buf(3)
!
! --------------------------------------------------------------------------- !
!
  xksurf = 0.0
  viscon = 0.0
  vise   = 0.0
!
! ROOT PROCESS COMPUTE
  IF (iss .EQ. 0) THEN
    iz   = 1
    izm1 = iz - 1
    izp1 = iz + 1
    dz_i = dzu_i(izp1)
!
    CALL sufto(it)
    IF (qstar(1) .EQ. 0.) THEN
      zeta = 0.0
    ELSE
      zeta = ABS(z(1))/amonin
    END IF
!
    IF (ismlt .EQ. 1) THEN
      CALL busngr(zeta,phim,phis,psim,psis)
    ELSE
      CALL fzol(zeta,phim,phis,psim,psis)
    END IF
!
    viscon = vk*ABS(z(1))/(utau*phim)
    vise   = utau*vk*ABS(z(1))/phim
!
!   FIND VALUE AT Z1 TO MATCH SURFACE LAYER
    vws = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        uws = uws + 0.5*(u(ix,iy,iz)-u_mn(iz) +     &
              u(ix,iy,izp1) - u_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
        vws = vws + 0.5*(v(ix,iy,iz)-v_mn(iz) +     &
              v(ix,iy,izp1) - v_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
      END DO
    END DO
!
    uws = uws*fnxy
    vws = vws*fnxy
!
!   FIND AVERAGE FLUCTUATING EDDY VISCOSITY
    DO iy=iys,iye
      DO ix=1,nnx
        e(ix,iy,iz)=AMAX1(e(ix,iy,iz),sml_eg)
      END DO
    END DO
!
    dslk = AMIN1(dsl,vk*ABS(z(iz))/csmag)
    almin = almin_c*dsl_z(iz)
!
    DO iy=iys,iye
      DO ix=1,nnx
        alwk(ix,iy)=dslk
!
!       NO STABILITY CORRECTED LENGTH SCALES WHEN NEW EDDY VISCOSITY IS ON
        xkvis(ix,iy)=ck*alwk(ix,iy)*SQRT(e(ix,iy,iz))*dfac(1)
      END DO
    END DO
!
!   FIND AVERAGE VISCOSITY
    xkavg = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        xkavg = xkavg + xkvis(ix,iy)
      END DO
    END DO
!
    xkavg = xkavg*fnxy
    buf(1) = uws
    buf(2) = vws
    buf(3) = xkavg
!
    CALL mpi_sum_xy(buf,myid,iss,ise,3)
!
    uws   = buf(1)
    vws   = buf(2)
    xkavg = buf(3)
    xkz1 = vise - SQRT(uws**2 + vws**2)*viscon
    xksurf =  xkz1 - xkavg
    xksurf = AMAX1(xksurf,0.0)
    xksurf = AMIN1(xksurf,vise)
!
  END IF
!
! SEND VALUES TO OTHER PROCESSES
  send(1) = xksurf
  send(2) = viscon
  send(3) = vise
!
  CALL mpi_bcast(send,3,mpi_real8,i_root,mpi_comm_world,ierr)
!
  xksurf = send(1)
  viscon = send(2)
  vise   = send(3)
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
6000 FORMAT(' dfac = ',e12.4,' xkavg = ',e12.4,' xkz1 = ',e12.4,/,  &
             ' vise = ',e12.4,' xksurf = ',e12.4)
! ---------------------------------------------------------------------------- !
!
END SUBROUTINE surfvis
