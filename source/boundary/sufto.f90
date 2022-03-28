! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine uses a version of similarity theory which has been   !
!         adpated for ocean flows. The options are to use Businger or          !
!         Schumann, set in main/pars.f90.                                      !
! ============================================================================ !
!
SUBROUTINE sufto(it)
!
    USE inputs, ONLY: u10_windspeed
    USE pars
    USE fields
    USE con_data
    USE con_stats
    USE diurnal
!
    REAL buf(3+nscl)
!
! --------------------------------------------------------------------------- !
!
  iz    = 1
  izm1  = iz - 1
  izp1  = iz + 1
  z1_a  = ABS(z1)
  buf(1)  = 0.0
  buf(2)  = 0.0
  buf(3)  = 0.0
  tol     = 0.01
!
  DO iy=iys,iye
    DO ix=1,nnx
      buf(1) = buf(1) + u(ix,iy,iz)
      buf(2) = buf(2) + v(ix,iy,iz)
      wind(ix,iy) = SQRT((u(ix,iy,iz)+ugal)**2+v(ix,iy,iz)*v(ix,iy,iz))
      buf(3) = buf(3) + wind(ix,iy)
    END DO
  END DO
!
  DO iscl=1,nscl
    buf(3+iscl) = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
      END DO
    END DO
  END DO
!
! SUM OVER X-Y SLABS
  CALL mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
  u1xy  = buf(1)*fnxy + ugal
  v1xy  = buf(2)*fnxy
  windm = buf(3)*fnxy
!
  DO iscl=1,nscl
    t1xy(iscl) = buf(3+iscl)*fnxy
  END DO
!
  vsfc  = SQRT(u1xy*u1xy+v1xy*v1xy)
  windm = AMAX1(windm,ufree)
  vsfc  = AMAX1(vsfc,ufree)
  t10xy(1)=-qstar(1)/utau*zody*vk74in
!
! CHECK TEMPERATURE AT BOUNDARY CONDITIONS
  IF (isfc .EQ. 0 ) THEN
    tsfcc(1)=t1xy(1)-t10xy(1)
  END IF
!
! SURFACE WIND STRESS - tau = 0.0184n/m*m, rho = 1000kg/m^3)
  utau = SQRT(rho_a*(8.5e-4)*u10_windspeed*u10_windspeed/rho_w)
  utausv = utau
  utau2  = utau*utau
!
  IF (ibuoy .EQ. 0 .OR. qstar(1) .EQ. 0.) THEN
    amonin    = 1000.
    zeta      = 0.
    thstar(1) = 0.0
    t10xy(1)  = 0.0
  ELSE
    amonin = -utau2*utau/(batagk*qstar(1))
    zeta   = z1_a/amonin
  END IF
!
  IF (t10xy(1) .LT. 0. .AND. qstar(1) .LT. 0.) THEN
    WRITE(6,1234)u1xy,v1xy,t1xy(1),tsfcc(1),amonin,utau,it
    GO TO 9999
  END IF
!
! FIND DRIFT VELOCITY FOR STABLE/NEUTRAL/UNSTABLE PBL
  IF (ismlt .EQ. 1) THEN
    CALL busngr(zeta,phim,phis,psim,psis)
  ELSE
    CALL fzol(zeta,phim,phis,psim,psis)
  END IF
!
  udrift = windm + stokes(1) - stokess + utau*(zody-psim)*vkin
  vdrift = 0.0
  dnom   = (zody-psis)*vk74in
!
  IF (isfc .EQ. 1) THEN
    thstar(1) = (t1xy(1) - tsfcc(1))/dnom
    t10xy(1)  = thstar(1)*dnom
    qstar(1)  = - utau*thstar(1)
  ELSE
    thstar(1)  = -qstar(1)/utau
    tsfcc(1)   = t1xy(1)-thstar(1)*dnom
    t10xy(1)   = thstar(1)*dnom
  END IF
!
  zol = zeta
  hol = zol*zi/z1
  utau2 = utau*utau
  au13m = utau2
  au23m = 0.0
  aut3m(1)= wtsfc(1)
!
  RETURN
!
! ---------------------------------- ERROR ----------------------------------- !
! TROUBLE IN SL ROUTINE
  9999 CONTINUE
  WRITE (nprt,9000)
  CALL mpi_finalize(ierr)
!------------------------------------------------------------------------------!
!
! ---------------------------------- FORMAT ---------------------------------- !
1234 FORMAT(' ** check sfc u=',e12.3,' v=',e12.3,' t,ts=',2f10.3, &
             ' l=',e12.3,' u*=',e12.3,' at it=',i5)
9000 FORMAT(' Trouble in SR. sufto')
!------------------------------------------------------------------------------!
!
STOP
END SUBROUTINE sufto
