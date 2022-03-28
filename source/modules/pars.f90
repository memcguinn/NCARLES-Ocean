! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
MODULE pars
!
  INTEGER, PARAMETER :: flg_stokes = 1      ! Stokes drift (0/1)
  INTEGER, PARAMETER :: flg_wavebreak = 0   ! Wave breaking (0/1)
  INTEGER, PARAMETER :: flg_diurnal = 1     ! Diurnal forcing (0/1)
  INTEGER, PARAMETER :: flg_reaction = 0    ! Carbonate chemistry reactions (0/1)
  INTEGER, PARAMETER :: flg_asflux = 0      ! Air-sea flux (0/1)
!
  INTEGER, PARAMETER :: itmax=100000,     & ! Maximum number of iterations
                        imean=1,          & ! Time increment ???
                        ihst=1,           & ! Frequency of history files
                        itape=500,        & ! Frequency for data outputs
                        itstr=1,          & ! Frequency of iterations
                        it_his=9e6,       & ! Iteration history (keep large, do not change for restart)
                        i_viz=9e6
!
  INTEGER, PARAMETER :: nscl = 8,         & ! Number of scalars (1 iff flg_reaction = 0, 8 if else)
                        nvar = (4+nscl)     ! Total number of variables
  INTEGER, PARAMETER :: nxg1  = 64,       & ! Number of cells in x-direction
                        nyg1  = 64,       & ! Number of cells in y-direction
                        nzg1  = 64          ! Number of cells in z-direction
  INTEGER, PARAMETER :: maxnx = 256,      & ! Maximum number of points in x-direction
                        maxny = 256,      & ! Maximum number of points in y-direction
                        maxnz = 256         ! Maximum number of points in z-direction
  INTEGER, PARAMETER :: maxnz1 = maxnz + 1, maxnz2 = maxnz + 2, &
                        maxnx2 = maxnx + 2, maxny2 = maxny + 2
!
  INTEGER, PARAMETER :: noalis=1,         & ! No aliasing
                        ismlt=0,          & ! Businger similarity constants (0/1)
                        ifree=0,          & ! Lower boundary condition
                        isfc=0,           & ! Temperature boundary condition
                        iupwnd=1,         & ! Skew symmetric advection form for vertical flux (0/1)
                        ibuoy=1,          & ! Buoyancy effects
                        itcut=1,          & ! Iteration cutoff ???
                        method=3,         & ! See pbltop for description
                        idebug=0,         & ! Write debug file (0/1)
                        iz_space=0,       & ! Varied vertical spacing (0/1)
                        ivis0=1,          & ! THIS CAN PROB BE REMOVED
                        ifix_dt=0,        & ! Fixed time step (0/1)
                        new_vis=-1,       & ! THIS CAN PROB BE REMOVED
                        i_dear = 0          ! Deardorff stability (0/1)
!
  INTEGER, PARAMETER :: j_recl=4            ! Record length in "bytes" for history file
  INTEGER, PARAMETER :: k8=8                ! Kind parameter for integers in mpi_io routines
!
!
! ------------------------------- DEFINITIONS -------------------------------- !
  INTEGER ::    nnx, nny, nnz, nxy, ncx, nnxp1, nnyp1, ncy,         &
                nnxp2, nnyp2, nnzp1, ivis, nnzm1, isize, krec,      &
                izs, ize, ixs, ixe, jxs, jxe, kxs, kxe,             &
                mxs, mxe, iss, ise, iys, iye, jys, jye
!
  INTEGER ::    it_counter, itn
!
!------------------------------------------------------------------------------!
  CHARACTER case*3
!
!------------------------------------------------------------------------------!
  INTEGER  ::   nvel, npre, nhis1, nprt,                            &
                nhisp, nvelc,                                       &
                nviz_z, nviz_y,                                     &
                nviz_x, nviz_s,                                     &
                kfile, jfile, ibcl, ibcu,                           &
                igrdr, imach, it_nxt
!
  LOGICAL ::    mnout, micut, mtape, mhis, msave,                   &
                l_root, l_debug, msave_v, mviz
!
!------------------------------------------------------------------------------!
  REAL    ::    windm,u1xy,v1xy,t1xy(nscl),                         &
                t10xy(nscl),au13m,au23m,aut3m(nscl),tsfcm(nscl),    &
                thstar(nscl), eavg(maxnz), tr_tau(0:maxnz),         &
                zi_min, dir_x, dir_y
!
  INTEGER ::    izi, iz_min
  REAL, ALLOCATABLE ::                                              &
                wind(:,:), tau13m(:,:), tau23m(:,:),                &
                taut3m(:,:,:), t_grnd(:,:,:)
!
!------------------------------------------------------------------------------!
  REAL ::       u_mn(0:maxnz1), v_mn(0:maxnz1),                     &
                w_mn(0:maxnz1), t_mn(0:maxnz1,nscl)
!
!------------------------------------------------------------------------------!
  REAL ::       dzw(0:maxnz2), dzu(0:maxnz2),                       &
                dzw_i(0:maxnz2), dzu_i(0:maxnz2)
!
!------------------------------------------------------------------------------!
  REAL ::       t_factor, t_ref, t_di
!
!------------------------------------------------------------------------------!
  REAL ::       dfac(maxnz), dsl_z(0:maxnz1),                       &
                xksurf, viscon, vise, almin_c,stabmin,              &
                ck,ceps,csmag,stab_c,vis_mean(0:maxnz)
  INTEGER ::    nmatch
!
!------------------------------------------------------------------------------!
  REAL ::       zetas(3), gama(3), etas(4), dt_new,                 &
                umax,vmax,wmax, wabs, emax, vismax,                 &
                cfl, tzero,                                         &
                ucfl, vcfl, wcfl
!
!------------------------------------------------------------------------------!
  REAL ::       wave_height, R_hw, k_nonbreak,                      &
                q0_cool, qd_heat, t_heat
  REAL, DIMENSION(2) ::                                             &
                cos_arr
!
!------------------------------------------------------------------------------!
  CHARACTER*80  path_res, path_sav, path_his, path_prt,             &
                path_hp, path_sav_hp,                               &
                path_v, path_c, path_p, path_h,                     &
                path_sav_v, path_sav_c,                             &
                path_sav_p, path_sav_h,                             &
                bad_news
  CHARACTER case_inp*3
  CHARACTER*80 path_viz_xy, path_viz_xz, path_viz_yz, path_stuf
!
!------------------------------------------------------------------------------!
  INTEGER ::    myid, numprocs, i_root,                             &
                ziloc, myid_newvis, ncpu_s, ncpu_z, maxp
  INTEGER, ALLOCATABLE, DIMENSION(:) ::                             &
                ix_s, ix_e, jx_s, jx_e,                             &
                kx_s, kx_e, mx_s, mx_e,                             &
                iy_s, iy_e, jy_s, jy_e,                             &
                is_s, is_e, iz_s, iz_e
!------------------------------------------------------------------------------!
!
END MODULE pars
