PROGRAM les_mpi

  USE pars
  USE fields
  USE con_data
  USE con_stats
  USE tracerbc, ONLY: applytracerbc

  INCLUDE 'mpif.h'

  ! INITIALIZE MPI, GET MYID, NUMPROCS
  ! TEST IF ON ROOT PROCESS

  CALL mpi_init(ierr)
  CALL mpi_comm_rank(mpi_comm_world,myid,ierr)
  CALL mpi_comm_size(mpi_comm_world,numprocs,ierr)

  i_root = 0
  l_root = .FALSE.
  IF(myid .EQ. i_root) l_root = .TRUE.

  l_debug = .FALSE.
  IF(idebug .EQ. 1) l_debug = .TRUE.

  ts_mpi = mpi_wtime()

  ! SET NUMBER OF X-Y SLAB CPUS
  ncpu_s   = 2
  itn      = 0
  case_inp = '30L'

  CALL get_units
  CALL gridd
  CALL setcon
  CALL set_paths
  istop = 1

  ! SCRATCH RUN
  IF (iti.EQ.0)  THEN
    igrdr = 2
    case = case_inp
    CALL init
    CALL setup(it)

    ! CHOOSE ROUTINE FOR GETTING INITIAL GUESS
    CALL randoc
    CALL get_max
  ELSE
    igrdr = 3
    CALL restart
    CALL get_max
    CALL setup(it)
  ENDIF

  ! TIME LOOP
  tzero = time
  CALL get_dt(it,iti)
  9000 CONTINUE
  CALL set_sav(it,iti)

  ! UPDATE POSITION OF VORTEX
  IF(it .GE. new_vis .AND. ivis0 .EQ. 1) THEN
    ivis = 1
  ELSE
    ivis = 0
  ENDIF

  ! 3 STAGE RUNGE-KUTTA TIME STEPPING
  DO 8999 istage=1,3
  dtzeta = dt*zetas(istage)
  dtgama = dt*gama(istage)

  ! COMPUTE DERIVATIVES OF (U,V,W)
  CALL exchange
  CALL get_derv

  ! NEW EDDY VISCOSITY, AND BCS
  IF(iss .EQ. 0 .AND. ifree .EQ. 0) THEN
    CALL lower(it)
  ELSEIF(ifree .EQ. 1) THEN
    CALL lower_free(it)
  ENDIF

  IF(ise .EQ. numprocs-1) THEN
    CALL upper
  ENDIF

  CALL applytracerbc(it)
  CALL bcast_pbc
  CALL get_means(istage)

  IF(ivis .EQ. 1) THEN
    CALL iso(it)
    CALL surfvis(it)
  ENDIF

  IF(istage .EQ. 1)THEN
    CALL xy_stats
    CALL tke_budget
    CALL pbltop(itop)
  ENDIF

  ! GET RHS FOR ALL EQUATIONS
  IF(istage.EQ.1 .AND. flg_reaction.EQ.1)THEN
    CALL strang1(it)
  ENDIF

  CALL comp1(istage,it)

  ! SOLVE FOR PRESSURE
  CALL comp_p

  ! ADD PRESSURE GRADIENT AND DEALIAS
  CALL comp2

  IF(istage.EQ.3 .AND. flg_reaction.EQ.1)THEN
    CALL strang1(it)
  ENDIF

  IF(msave .AND. istage .EQ. 3) THEN
    CALL save_v(it)
  ENDIF

  IF(istage .EQ. 3) THEN
    IF(msave .AND. l_root) CALL save_c(it)
  ENDIF

  IF(micut) THEN
    CALL dealias
  ENDIF

  IF(mnout .AND. istage .EQ. 1)  THEN
    IF(l_debug) THEN
      CALL print(nprt,it,izs,ize)
    ENDIF
    IF(l_root) CALL print(6,it,1,nnz)
  ENDIF

  IF(l_root) THEN
    IF(mhis  .AND. istage .EQ. 1)  CALL write_his(itop)
    IF(mhis  .AND. istage .EQ. 1 .AND. mtape) CALL close_his
  ENDIF

  8999 CONTINUE
  CALL get_max
  CALL get_dt(it,iti)

  IF (it.GE.itmax) GO TO 99000
  GO TO 9000

  99000 CONTINUE
  te_mpi = mpi_wtime()

  WRITE(6,9997) (te_mpi - ts_mpi)
  9997 FORMAT(' Job Execution Time = ',e15.6)

  9998 CONTINUE
  CALL mpi_finalize(ierr)

  STOP
END PROGRAM
