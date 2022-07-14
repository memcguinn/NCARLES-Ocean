      program les_mpi_2d
c
      use pars
      use fields
      use con_data
      use con_stats
      use tracerbc, only: applytracerbc
      include 'mpif.h'
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
c -------------------------------------------------------------------------------
c
c ---------- initialize MPI, get myid, numprocs,
c            test if on root process
c
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)
c
      i_root = 0
      l_root = .false.
      if(myid .eq. i_root) l_root = .true.
c
      l_debug = .false.
      if(idebug .eq. 1) l_debug = .true.
c
      ts_mpi = mpi_wtime()
c
c     -------- set number of x-y slab cpus c
c     ncpu_s = 1
c     ncpu_s = 2
c      ncpu_s = 8
      ncpu_s = 2
c
      itn      = 0
      case_inp = '30L'
c
      call get_units
      call gridd
      call setcon
      call set_paths
      istop = 1
c
c -------------- scratch run
c
      if (iti.eq.0)  then
         igrdr = 2
         case = case_inp
         call init
c         call init_hurr
         call setup(it)
c
c ---------- choose routine for getting initial guess
c
         if(iocean .eq. 1) then
            call randoc
c           call get_fields
         else
c           call random_f
c           call get_fields
            call random
         endif
         call get_max
      else
         igrdr = 3
         call restart
         call get_max
         call setup(it)
      endif
c
c --------------- time loop ------------
c
      tzero = time
      call get_dt(it,iti)
 9000 continue
      call set_sav(it,iti)
c
c ------- update position of vortex
c
c      call hurr_move
c      call stokesv
c
      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif
c
c ---------------- 3 stage runge-kutta time stepping
c
      do  8999 istage=1,3
c
      dtzeta = dt*zetas(istage)
      dtgama = dt*gama(istage)
c
c ---------- compute derivatives of (u,v,w)
c
      call exchange
c
      call get_derv
c
c --------- new eddy viscosity, and bcs
c
      if(iss .eq. 0 .and. ifree .eq. 0) then
         call lower(it)
      elseif(ifree .eq. 1) then
         call lower_free(it)
      endif
      if(ise .eq. numprocs-1) then
         call upper
      endif
      call applytracerbc(it)
      call bcast_pbc
      call get_means(istage)
      if(ivis .eq. 1) then
         call iso(it)
         call surfvis(it)
      endif
      if(istage .eq. 1)then
        call xy_stats
        call tke_budget
        call pbltop(itop)
      endif
c
c ------------ save velocity field
c
c      if(msave .and. istage .eq. 1) then
c         call save_v(it)
c      endif
      
c      if(msave_v .and. istage .eq. 1) then
c         call save_viz(it)
c      endif
c
c ------------ save pressure field
c
c      if(msave .and. istage .eq. 1) then
c         call save_p
c      endif
c
c --------- get rhs for all equations
c
      if(istage.eq.1 .and. flg_reaction.eq.1)then
         call strang1(it)
      endif
      call comp1(istage,it)
c      if(istage .eq. 1) then
c         if(msave .and. l_root) call save_c(it)
c      endif
c
c --------- solve for pressure
c
      call comp_p
c
c --------- add pressure gradient and dealias
c
      call comp2
      if(istage.eq.3 .and. flg_reaction.eq.1)then
         call strang1(it)
      endif
      
      if(msave .and. istage .eq. 3) then
         call save_v(it)
      endif
      if(istage .eq. 3) then
         if(msave .and. l_root) call save_c(it)
      endif

      if(micut) then
         call dealias
      endif
      if(mnout .and. istage .eq. 1)  then
          if(l_debug) then
             call print(nprt,it,izs,ize)
          endif
          if(l_root) call print(6,it,1,nnz)
      endif
      if(l_root) then
         if(mhis  .and. istage .eq. 1)  call write_his(itop)
         if(mhis  .and. istage .eq. 1 .and. mtape) call close_his
      endif
c
 8999 continue
      call get_max
      call get_dt(it,iti)
      if (it.ge.itmax) go to 99000
      go to 9000
c
99000 continue
      te_mpi = mpi_wtime()
      write(6,9997) (te_mpi - ts_mpi)
 9997 format(' Job Execution Time = ',e15.6)
c
 9998 continue
      call mpi_finalize(ierr)
c
      stop
      end
