      module pars
c ----------------------------------------------------------------------

cccc
cccc  PARAMETERS ARE DESCRIBED IN PARAMS.pdf
cccc

c     nslab=nproc/ncpu_s=256/24=6
c     nys=ny/nslab=256/6=24
c     integer, parameter :: ncpu_s=24 !set number of x-y slab cpus
c     character*3, parameter :: case_inp = 'yel'
c
c      integer, parameter :: flg_strang = 1 ! strang splitting for reaction
      integer, parameter :: flg_stokes = 1 ! stokes on or off
      integer, parameter :: flg_reaction = 1 ! 3.1536e8 reaction model on or off
      integer, parameter :: iti=0, itmax=20, imean=1, ihst=01,
     +            itape=20, itstr=1, it_his=99999, i_viz=99999
c
      integer, parameter :: nscl = 8, nvar = (4+nscl) !number of scalars and vars
      integer, parameter :: nxg1  = 64, nyg1  = 64, nzg1  = 64 !size of problem
      integer, parameter :: maxnx = 128, maxny = 128, maxnz = 128 !max size
c     integer, parameter :: wid_x=0, wid_y=maxny/4 !max processors
      integer, parameter :: maxnz1 = maxnz + 1, maxnz2 = maxnz + 2,
     +                      maxnx2 = maxnx + 2, maxny2 = maxny + 2
c
c
c ------------- definition of internal flags
c
c       igrdr   =  3; data comes from restart file
c               =  2; data comes from initialization (random)
c               =  1; data comes from coarser grid (or otherwise)
c
c       ibcu    = 1 ; upper boundary condition set by radiation bc
c               = 0 ; fixed value = 0.
c               = -1; value defined by coarser mesh for all variables
c
c       ibcl    = 0 ; lower boundary condition set by similarity theory (sr. setup)
c               = -1; value defined by coarser mesh for all variables
c
c       ifix_dt = 0 ; variable time step with fixed cfl number in setcon
c               = 1 ; fixed time step set in sr. get_dt
c
c       ifree   = 0 ; use spatially averaged surface conditions for MO (call lower)
c               = 1 ; use point-by-point conditions for MO free convection (call lower_free)
c
c       ihst    = nn; frequency at which global variables are output in history file
c               < 0 ; no history files
c
c       it_his  = time step where history files start, incremented by itape
c
c       ismlt   = 1 ; use businger formulas in MO
c                 0 ; use large and everyone elses formulas in MO
c
c       iupwnd  = 0;  use skew symmetric formulas for all derivatives
c                     in scalar equations
c               = 1;  use hybrid upwind scheme for all derivatives
c                     in scalar equations
c
c       ivis0   = 0; old eddy viscosity model
c               = 1; new eddy viscosity model
c
c       new_vis = step; the iteration step for which the new model
c                       is turned on when ivis0=1
c               < 0; new model is on at all steps for ivis0=1
c
c       nscl  .ge. 1   number of scalars to be followed set in parameter statements
c                      change entries in sr. init, and sr. suft for surface bc's
c
      integer, parameter :: noalis=1, ismlt=0, ifree=0, isfc=0,
     +          iradup=0, iupwnd=1, ibuoy=1, ifilt=0, itcut=1, isubs=0,
     +          ibrcl=0, method=3, idebug=0, iz_space=0,
     +          ivis0=1, ifix_dt=1, new_vis=-1, i_dear = 0,
     +          force_tracer=0
      integer, parameter :: j_recl=4 !record length in "bytes" for history file
      integer, parameter :: k8=8 !kind parameter for integers in mpi_io routines
c
c -------- number of vars, size of problem in (x,y,z), max size, max processors
c
c
c ------------ leave following definitions as is
c
c ----------------------------------------------------------------------
      integer ::    nnx, nny, nnz, nxy, ncx, nnxp1, nnyp1, ncy,
     +              nnxp2, nnyp2, nnzp1, ivis, nnzm1, isize, krec,
     +              izs, ize, ixs, ixe, jxs, jxe, kxs, kxe,
     +              mxs, mxe, iss, ise, iys, iye, jys, jye
      integer ::    it_counter
c ----------------------------------------------------------------------
      character case*3
c ----------------------------------------------------------------------
      integer  ::   nvel, npre, nhis1, nprt,
     +              nhisp, nvelc,
     +              nviz_z, nviz_y,
     +              nviz_x, nviz_s,
     +              kfile, jfile, ibcl, ibcu,
     +              igrdr, imach, itn, it_nxt
      logical ::    mnout, micut, mtape, mhis, msave,
     +              l_root, l_debug, msave_v, mviz
c ----------------------------------------------------------------------
      real    ::    windm,u1xy,v1xy,t1xy(nscl),
     +              t10xy(nscl),au13m,au23m,aut3m(nscl),tsfcm(nscl),
     +              thstar(nscl), eavg(maxnz), tr_tau(0:maxnz),
     +              zi_min, dir_x, dir_y
      integer ::    izi, iz_min
      real, allocatable ::
     +              wind(:,:), tau13m(:,:), tau23m(:,:),
     +              taut3m(:,:,:), t_grnd(:,:,:)
c ----------------------------------------------------------------------
      real ::       u_mn(0:maxnz1), v_mn(0:maxnz1),
     +              w_mn(0:maxnz1), t_mn(0:maxnz1,nscl)
c ----------------------------------------------------------------------
      real ::       dzw(0:maxnz2), dzu(0:maxnz2),
     +              dzw_i(0:maxnz2), dzu_i(0:maxnz2)
c ----------------------------------------------------------------------
      real ::       t_factor, t_ref
c , c_rate, t_surf_i
c ----------------------------------------------------------------------
      real ::       dfac(maxnz), dsl_z(0:maxnz1),
     +              xksurf, viscon, vise, almin_c,stabmin,
     +              ck,ceps,csmag,stab_c,vis_mean(0:maxnz)
      integer ::    nmatch
c ----------------------------------------------------------------------
      real ::       zetas(3), gama(3), etas(4), dt_new,
     +              umax,vmax,wmax, wabs, emax, vismax,
     +              cfl, tzero,
     +              ucfl, vcfl, wcfl
c ----------------------------------------------------------------------
      character*80  path_res, path_sav, path_his, path_prt,
     +              path_hp, path_sav_hp,
     +              path_v, path_c, path_p, path_h,
     +              path_sav_v, path_sav_c,
     +              path_sav_p, path_sav_h,
     +              bad_news
      character case_inp*3
      character*80 path_viz_xy, path_viz_xz, path_viz_yz, path_stuf
c ----------------------------------------------------------------------
      integer ::    myid, numprocs, i_root,
     +              ziloc, myid_newvis, ncpu_s, ncpu_z, maxp
      integer, allocatable, dimension(:) ::
     +              ix_s, ix_e, jx_s, jx_e,
     +              kx_s, kx_e, mx_s, mx_e,
     +              iy_s, iy_e, jy_s, jy_e,
     +              is_s, is_e, iz_s, iz_e
      end module pars
