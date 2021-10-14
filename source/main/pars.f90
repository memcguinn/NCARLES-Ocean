! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
MODULE pars
!
  INTEGER, PARAMETER :: flg_stokes = 1     ! Stokes drift (0/1)
  INTEGER, PARAMETER :: flg_reaction = 1   ! Carbonate chemistry reactions (0/1)
  INTEGER, PARAMETER :: iti=0,           & ! Start iteration for restart (default 0)
                        itmax=100000,    & ! Maximum number of iterations
                        imean=1,         & ! Time increment ???
                        ihst=1,          & ! Frequency of history files
                        itape=100,       & ! Frequency for data outputs
                        itstr=1,         & ! Frequency of iterations
                        it_his = 99999     ! Iteration history (keep large, do not change for restart)
!
  INTEGER, PARAMETER :: nscl = 8,        & ! Number of scalars (1 iff flg_reaction = 0, 8 if else)
                        nvar = (4+nscl)    ! Total number of variables
  INTEGER, PARAMETER :: nxg1  = 128,     & ! Number of cells in x-direction
                        nyg1  = 128,     & ! Number of cells in y-direction
                        nzg1  = 128        ! Number of cells in z-direction
  INTEGER, PARAMETER :: maxnx = 128,     & ! Maximum number of points in x-direction
                        maxny = 128,     & ! Maximum number of points in y-direction
                        maxnz = 128        ! Maximum number of points in z-direction
  INTEGER, PARAMETER :: maxnz1 = maxnz + 1, maxnz2 = maxnz + 2, &
                        maxnx2 = maxnx + 2, maxny2 = maxny + 2
!
  INTEGER, PARAMETER :: noalis=1,        & ! No aliasing
                        ismlt=0,         & ! Businger similarity constants (0/1)
                        ifree=0,         & ! Lower boundary condition
                        isfc=1,          & ! Temperature boundary condition
                        ibcu=0,          & ! Gradient boundary condition (0/1, on = 1)
                        ibcl=0,          & ! Lower boundary condition
                        iupwnd=1,        & ! Skew symmetric advection form for vertical flux (0/1)
                        ibuoy=1,         & ! Buoyancy effects
                        itcut=1,         & ! Iteration cutoff ???
                        method=3,        & ! See pbltop for description
                        idebug=0,        & ! Write debug file (0/1)
                        iz_space=0,      & ! Varied vertical spacing (0/1)
                        ifix_dt=0,       & ! Fixed time step (0/1)
                        i_dear=0           ! Deardorff stability (0/1)
!
  INTEGER, PARAMETER :: j_recl=4            ! record length in "bytes" for history file
  INTEGER, PARAMETER :: k8=8                ! kind parameter for integers in mpi_io routines
!
!
! ------------------------------- DEFINITIONS -------------------------------- !
  INTEGER ::   nnx, nny, nnz, nxy, ncx, nnxp1, nnyp1, ncy,          &
               nnxp2, nnyp2, nnzp1, nnzm1, isize, krec,             &
               izs, ize, ixs, ixe, jxs, jxe, kxs, kxe,              &
               mxs, mxe, iss, ise, iys, iye, jys, jye
!
  INTEGER ::   it_counter
!
!------------------------------------------------------------------------------!
  CHARACTER case*3
!
!------------------------------------------------------------------------------!
  INTEGER ::   nvel, npre, nhis1, nprt,nhisp, nvelc, kfile, jfile,   &
               imach, itn, it_nxt
!
  LOGICAL ::   mnout, micut, mtape, mhis, msave, l_root, l_debug
!
!------------------------------------------------------------------------------!
  REAL    ::   windm, u1xy, v1xy, t1xy(nscl), t10xy(nscl),          &
               au13m, au23m, aut3m(nscl), tsfcm(nscl),              &
               thstar(nscl), eavg(maxnz), tr_tau(0:maxnz),          &
               zi_min, dir_x, dir_y
!
  INTEGER ::   izi, iz_min
!
  REAL, ALLOCATABLE ::                                              &
               wind(:,:), tau13m(:,:), tau23m(:,:), taut3m(:,:,:), t_grnd(:,:,:)
!
!------------------------------------------------------------------------------!
  REAL    ::   u_mn(0:maxnz1), v_mn(0:maxnz1),                      &
               w_mn(0:maxnz1), t_mn(0:maxnz1,nscl)
!
!------------------------------------------------------------------------------!
  REAL    ::   dzw(0:maxnz2), dzu(0:maxnz2),                        &
               dzw_i(0:maxnz2), dzu_i(0:maxnz2)
!
!------------------------------------------------------------------------------!
  REAL    ::   t_factor, t_ref
!
!------------------------------------------------------------------------------!
  REAL    ::   dfac(maxnz), dsl_z(0:maxnz1), xksurf, viscon,        &
               vise, almin_c, stabmin, ck, ceps, csmag, stab_c,     &
               vis_mean(0:maxnz)
  INTEGER ::   nmatch
!
!------------------------------------------------------------------------------!
  REAL    ::   zetas(3), gama(3), etas(4), dt_new,                  &
               umax,vmax,wmax, wabs, emax, vismax,                  &
               cfl, tzero, ucfl, vcfl, wcfl
!
!------------------------------------------------------------------------------!
  CHARACTER*80  path_res, path_sav, path_his, path_prt, path_hp,    &
                path_sav_hp, path_v, path_c, path_p, path_h,        &
                path_sav_v, path_sav_c, path_sav_p, path_sav_h,     &
                bad_news
  CHARACTER case_inp*3
!
!------------------------------------------------------------------------------!
  INTEGER ::   myid, numprocs, i_root, ziloc, myid_newvis, ncpu_s, ncpu_z, maxp
  INTEGER, ALLOCATABLE, DIMENSION(:) ::                             &
                    ix_s, ix_e, jx_s, jx_e,                         &
                    kx_s, kx_e, mx_s, mx_e,                         &
                    iy_s, iy_e, jy_s, jy_e,                         &
                    is_s, is_e, iz_s, iz_e
!------------------------------------------------------------------------------!
!
END MODULE pars
