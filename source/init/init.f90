! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE init
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
!
! --------------------------------------------------------------------------- !
!
  pi     = 4.0*ATAN(1.0)            ! Pi = 3.1417
  pi2    = 2.0*pi                   !
  d_to_r = pi2/360.0                ! Coversion factor - degree to radian
  grav   = 9.81                     ! Gravitational constant
  bfac   = 1.0                      !
!
  IF (ibuoy .EQ. 0) bfac = 0.
!
! CASE SPECIFIC INPUTS
  rho_a   = 1.0                     ! Density of air
  rho_w   = 1000.0                  ! Density of water
  t00     = 283                     !
  t00b    = 5000.0                  !
  cp_a    = 1.0057e03               ! Isobaric specific heat, air
  cp_w    = 4.20e03                 ! Isobaric specific heat, water - J/(kg K)
  gcp     = grav/cp_w               !
  batag   = bfac*grav/t00b          !
!
! STOKES DRIFT PARAMETERS
  cpou10  = 0.6                     !
  rlat    = 30                      ! Latitude, degrees North from equator
  fcor    = 2.0*pi2*SIN(rlat*d_to_r)/(24.0*3600.0)  ! Coriolis Parameter
  ugcont  = 0.0                     !
  vgcont  = 0.0                     !
  wtsfc(1) = 0.0                    ! Define surface cooling
  qstar(1) = wtsfc(1)               ! Heat flux (see wrsfc, surface cooling)
!
  dtdzf(1)= 0.010
  dtjump  = 0.
  divgls  = 0.
!
! PROBLEM DIMENSION INPUTS
  zo      = 0.0001                  ! Initial vertial point location
  zi      = -30                     ! Initial mixed layer depth
  xl      = 320.0                   ! Physical x-direction dimension
  yl      = 320.0                   ! Physical y-direction dimension
  zl      = -96.0                   ! Depth, z-direction dimension
  izi     = NINT((zi/zl)*nnz)       ! Vertical grid location, nearest bottom ML
!
! IF STRETCHED GRID, LOCATION OF INITIAL POINT
  ann      = 0.00615                ! Donelan parameters
  bnn      = 1.0                    !
  f2w      = 0.13                   !
  f_p      = f2w*grav/u_10          !
  npm      = 4                      !
  sigma_p  = pi2*f_p                !
  grav_w   = grav                   !
!
! RATIO OF K_1/K_P - PHILLIPS AND DONELAN
  r_kp     = grav*(f2w*pi2/u_10)**2
  r_k1     = r_fac*grav/(cd_10*u_10*u_10)
  rk_ratio = r_k1/r_kp
  time  = 0.0
!
! OUTERMOST COARSE GRID, INDICES ARE BOUNDS
  izlow = 1
  izup  = nnz
  dz    = zl/nnz
  dzg   = ABS(dz)
  IF(l_root) WRITE (6,4040) zl,nnz,dzg
!
! GENERATE Z GRIDS FOR MESH FROM IZ = 0,1,...,NNZ+1, ALLOWS INDEXING IZ(0)
  zwstrt = 0.0
!
! FOR UNIFORM VERTICAL SPACING
  IF (iz_space .EQ. 0) THEN
    DO iz=0,nnz+1
      z(iz) = dz*FLOAT(iz) + zwstrt   ! Build x grid for w points
    ENDDO
  ELSE
    CALL vgrid(zw1,zi,zl,nnz,z(0),l_root,l_debug)
  ENDIF
!
  CALL get_dz
!
  IF (l_root) THEN
    WRITE (6,8002) zwstrt
    WRITE (6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
  ENDIF
!
  nnzm1 = nnz-1
  dx    = xl/nnx                    ! Discretization size x-direction
  dy    = yl/nny                    ! Discretization size y-direction
  fnxy  = 1./FLOAT(nxy)
  dzdz  = dzw(1)*dzw(1)
  z1    = zz(1)
  c23   = 2.0/3.0
  dsl   = (dx*1.5*dy*1.5*ABS(dzw(1)))**(1./3.)
  dslg  = dsl
  cs    = 0.2
  vk    = 0.4
  batagk = batag*vk
  vkin   = 1./vk
  ttmean = 0.
  zody   = ALOG(ABS(z1/zo))
!
  WRITE (nprt, 9901) z1,zo,zody
!
  zodyin = 1./zody
  wstar  = ABS(batag*zi*wtsfc(1))**(1./3.)
!
  IF (ismlt .EQ. 1) THEN
!
!   SET CONSTANTS FOR BUSINGER SIMILARITY FUNCTIONS
    vk74   = vk*0.74
    vk74in = 0.74/vk
    zody74 = zody*0.74
!
  ELSE
!
!   SET CONSTANTS FOR LARGE SIMILARITY FUNCTIONS
    vk74    = vk
    vk74in  = 1.0/vk
    zody74  = zody
!
  ENDIF
!
  ugal   = ugcont*0.5
  cdbtm  = vk*vk/zody/zody
!
! SET SURFACE FRICTION VELOCITY, ALSO IN SURFTO
  utau  = SQRT(rho_a*(8.5e-4)*10*10/rho_w)
  utau2 = utau*utau
!
  IF (ibuoy .EQ. 0 .OR. qstar(1) .EQ. 0.) THEN
    amonin = 1000.0
  ELSE
    amonin = -utau2*utau/(batagk*qstar(1))
  ENDIF
!
  hol   = ABS(zi)/amonin
  zol   = ABS(z1)/amonin
  uwsfc = -utau*utau
  vwsfc = -utau*utau
!
! CHECK TSFCC IS GREATER THAN T00 FOR ISFC = 0,1
  tsfcc(1) = 265
!
  IF (l_root) THEN
      WRITE (6,80)
      WRITE (6,2) wtsfc(1),utau,amonin,dtdzf(1),zody,zo,cdbtm,ugcont
  ENDIF
!
  IF (l_debug) THEN
    WRITE (nprt,80)
    WRITE (nprt,2) wtsfc(1),utau,amonin,dtdzf(1),zody,zo,cdbtm,ugcont
  ENDIF
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
    2 FORMAT(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/, &
            10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,     &
            ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4,                &
            ',  UG = ',e12.4)
   80 FORMAT(///,' ***** SCRATCH RUN ***** ',//)
 4040 FORMAT(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043 FORMAT(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 ! 7676 FORMAT(' in init qw_tot_aw = ',e15.6)
 8002 FORMAT(' zwstrt = ',e12.4)
 8003 FORMAT(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
 9901 FORMAT(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,' zody = ',e15.6)
!------------------------------------------------------------------------------!
!
END SUBROUTINE init
