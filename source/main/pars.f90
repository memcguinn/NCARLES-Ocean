! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
MODULE pars
!
  INTEGER, PARAMETER :: flg_stokes = 1
  INTEGER, PARAMETER :: flg_reaction = 1
  INTEGER, PARAMETER :: iti=0,            &
                        itmax=9e8,        &
                        imean=1,          &
                        ihst=01,          &
                        itape=1,          &
                        itstr=1,          &
                        it_his=99999,     &
                        i_viz=99999
!
  INTEGER, PARAMETER :: nscl = 8,         &
                        nvar = (4+nscl)
  INTEGER, PARAMETER :: nxg1  = 128,      &
                        nyg1  = 128,      &
                        nzg1  = 128
  INTEGER, PARAMETER :: maxnx = 128,      &
                        maxny = 128,      &
                        maxnz = 128
  INTEGER, PARAMETER :: maxnz1 = maxnz + 1, maxnz2 = maxnz + 2, &
                        maxnx2 = maxnx + 2, maxny2 = maxny + 2
  REAL, PARAMETER    :: airseaflux_time = 60.0, airsea_val = 8.57
!
  INTEGER, PARAMETER :: noalis=1,         &
                        ismlt=0,          &
                        ifree=0,          &
                        isfc=0,           &
                        iradup=0,         &
                        iupwnd=1,         &
                        ibuoy=1,          &
                        ifilt=0,          &
                        itcut=1,          &
                        isubs=0,          &
                        ibrcl=0,          &
                        iocean=1,         &
                        method=3,         &
                        idebug=0,         &
                        iz_space=0,       &
                        ihurr=0,          &
                        ivis0=1,          &
                        ifix_dt=1,        &
                        new_vis=-1,       &
                        i_dear = 0,       &
                        force_tracer=0
!
  INTEGER, PARAMETER :: j_recl=4
  INTEGER, PARAMETER :: k8=8 
!
!
! ------------------------------- DEFINITIONS -------------------------------- !
  INTEGER ::    nnx, nny, nnz, nxy, ncx, nnxp1, nnyp1, ncy,         &
                nnxp2, nnyp2, nnzp1, ivis, nnzm1, isize, krec,      &
                izs, ize, ixs, ixe, jxs, jxe, kxs, kxe,             &
                mxs, mxe, iss, ise, iys, iye, jys, jye
!
  INTEGER ::    it_counter
!
! ----------------------------------------------------------------------
  CHARACTER case*3
!
! ----------------------------------------------------------------------
  INTEGER  ::   nvel, npre, nhis1, nprt,                            &
                nhisp, nvelc,                                       &
                nviz_z, nviz_y,                                     &
                nviz_x, nviz_s,                                     &
                kfile, jfile, ibcl, ibcu,                           &
                igrdr, imach, itn, it_nxt
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
  REAL ::       t_factor, t_ref
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
