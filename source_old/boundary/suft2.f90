! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE suft2(u_level1,it)
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
!
    REAL :: u_level1(nnx,iys:iye,2+nscl)
!
  tol = 0.01
  ufree=0.07*(ABS(batag*qstar(1)*dzw(1)))**(1./3.)
  zeta_mn = -6.0
  zeta_mn_i = 1.0/zeta_mn
  iz   = 1
!
  DO iy=iys,iye
    DO ix=mxs,mxe
      utau = .001 ! first guess for utau
      t10xy(1) = -qstar(1)/utau*zody*vk74in
      tsfcc(1) = u_level1(ix,iy,3) - t10xy(1)
      vsfc2    = u_level1(ix,iy,1)**2 + u_level1(ix,iy,2)**2
      vsfc     = SQRT(vsfc2)
      windm    = ufree+vsfc
      utausv   = utau
      utau2    = utau*utau
      amonin   = -utau2*utau/(batagk*qstar(1))
      IF (amonin .EQ. 0.) THEN
        WRITE (6,5050) ix,iy,it,utau,amonin
        STOP
      END IF
!
!     FOR UNSTABLE, FREE CONVECTION PBL
      iter = 0
100   CONTINUE
!
!     LIMIT THE MINIMUM OF -L/Z TI ACCOMODATE STABLE FLOW
      zeta_i = AMIN1(amonin/z1,zeta_mn_i)
      zeta_a = 1.0/zeta_i
!
      IF (ismlt .EQ. 1) THEN
        CALL busngr(zeta_a,phim,phis,psim,psis)
      ELSE
        CALL fzol(zeta_a,phim,phis,psim,psis)
      END IF
!
      utau     = windm*vk/(zody-psim)
      thstar(1)= -qstar(1)/utau
      amonold  = amonin
      amonin   = utau*utau/(batagk*thstar(1))
      diff     = ABS(amonin - amonold)
      iter = iter+1
!
      IF (iter .GT. 10) GO TO 1000
      IF (diff .GT. ABS(tol*amonin)) GO TO 100
!
1000  CONTINUE
!
      IF (utau .GT. 10.) THEN
        WRITE(6,232)utau,windm
        STOP 9999
      END IF
!
      t10xy(1) = -qstar(1)/utau*vk74in*(zody-psis)
      t_grnd(ix,iy,1) = u_level1(ix,iy,3) - t10xy(1)
      zol = zeta_a
      hol = zol*zi/z1
      tep = u_level1(ix,iy,1)/windm
!
      IF (tep .GT. 1.) tep = 1.0
      IF (tep .LT. -1.) tep = -1.0
!
      thta  = ACOS(tep)
      utau2 = utau*utau
      tau13m(ix,iy)   = -utau2*COS(thta)
      tau23m(ix,iy)   = -utau2*SIN(thta)*SIGN(1.,u_level1(ix,iy,2))
      taut3m(ix,iy,1) = qstar(1)
!
    END DO
  END DO
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
 5050 FORMAT(' 5050, sr. suft2, trouble at ',/,' ix = ',i6,'iy = ', &
            i6,' it = ',i6,' utau = ',e15.6,' amonin = ',e15.6)
 232  FORMAT(' stop because utau=',e15.6,' windm=',e15.6)
!------------------------------------------------------------------------------!
!
END SUBROUTINE suft2
