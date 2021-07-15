! ============================================================================ !
! ABOUT:                                                                       !
!         les_mpi is the primary program which runs NCARLES-Ocean.             !
! ============================================================================ !
!
PROGRAM les_mpi
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
    USE tracerbc, ONLY: applytracerbc
!
    INCLUDE 'mpif.h'
!
! --------------------------------------------------------------------------- !
!
! INITIALIZE MPI, GET MYID, NUMPROCS, TEST IF ON ROOT PROCESS
  CALL mpi_init(ierr)
  CALL mpi_comm_rank(mpi_comm_world,myid,ierr)
  CALL mpi_comm_size(mpi_comm_world,numprocs,ierr)
!
  i_root = 0
  l_root = .FALSE.
  IF (myid .EQ. i_root) l_root = .TRUE.
!
  l_debug = .FALSE.
  IF (idebug .EQ. 1) l_debug = .TRUE.
!
  ts_mpi   = mpi_wtime()
  ncpu_s   = 2                         ! Set number of x-y slab cpus
  itn      = 20                        ! Input file number for restart (0, default)
  case_inp = '30L'
!
  CALL get_units
  CALL gridd
  CALL setcon
  CALL set_paths
!
  istop = 1
!
! ------------------------------- SCRATCH RUN -------------------------------- !
!
! CHOOSE ROUTINE FOR GETTING INITIAL GUESS
  IF (iti .EQ. 0) THEN
    case = case_inp
    CALL init
    CALL setup(it)
    CALL randoc
    CALL get_max
  ELSE
    CALL restart
    CALL get_max
    CALL setup(it)
  END IF
!
! --------------------------------- TIME LOOP -------------------------------- !
  tzero = time
!
  CALL get_dt(it,iti)
!
  DO WHILE (it .LE. itmax)
    CALL set_sav(it,iti)
!
!   3 STAGE RUNGE-KUTTA TIME-STEPPING
    DO istage=1,3
!
!     R-K: FIRST STAGE
      IF (istage .EQ. 1) THEN
        dtzeta = dt*zetas(istage)
        dtgama = dt*gama(istage)
!
!       COMPUTE DERIVATIES OF (U,V,W)
        CALL exchange
        CALL get_derv
!
!       NEW EDDY VISCOSITY, BOUNDARY CONDITIONS
        CALL lower(it)
!
        IF (ise .EQ. numprocs-1) THEN
          CALL upper
        END IF
!
        CALL applytracerbc(it)
        CALL bcast_pbc
        CALL get_means(istage)
        CALL iso(it)
        CALL surfvis(it)
        CALL xy_stats
        CALL tke_budget
        CALL pbltop(itop)
!
        IF (flg_reaction .EQ. 1) THEN
          CALL strang1(it)
        END IF
!
        CALL comp1(istage,it)
!
! ----------------------------- PRESSURE SOLVER ------------------------------ !
        CALL comp_p
!
!       ADD PRESSURE GRADIENT AND DEALIAS
        CALL comp2
!
        IF (micut) THEN
          CALL dealias
        END IF
!
        IF (mnout) THEN
          IF (l_debug) THEN
            CALL print(nprt,it,izs,ize)
          END IF
          IF (l_root) CALL print(6,it,1,nnz)
        END IF
!
        IF (l_root) THEN
          IF (mhis) CALL write_his(itop)
          IF (mhis .and. mtape) THEN
            CALL close_his
          END IF
        END IF
!
!     R-K: SECOND STAGE
      ELSE IF (istage .EQ. 2) THEN
!
        dtzeta = dt*zetas(istage)
        dtgama = dt*gama(istage)
!
!       COMPUTE DERIVATIVE OF (U,V,W)
        CALL exchange
        CALL get_derv
!
!       NEW EDDY VISCOSITY AND BOUNDARY CONDITIONS
        CALL lower(it)
!
        IF (ise .EQ. numprocs-1) THEN
          CALL upper
        END IF
!
        CALL applytracerbc(it)
        CALL bcast_pbc
        CALL get_means(istage)
        CALL iso(it)
        CALL surfvis(it)
        CALL comp1(istage,it)
        CALL comp_p
        CALL comp2
!
        IF (micut) THEN
          CALL dealias
        END IF
!
!     R-K: THIRD STAGE
      ELSE
!
        dtzeta = dt*zetas(istage)
        dtgama = dt*gama(istage)
!
!       COMPUTE DERIVATIVES OF (U,V,W)
        CALL exchange
        CALL get_derv
!
!       NEW EDDY VISCOSITY AND BCS
        CALL lower(it)
!
        IF (ise .EQ. numprocs-1) THEN
          CALL upper
        END IF
!
        CALL applytracerbc(it)
        CALL bcast_pbc
        CALL get_means(istage)
        CALL iso(it)
        CALL surfvis(it)
        CALL comp1(istage,it)
        CALL comp_p
        CALL comp2
!
        IF (flg_reaction .EQ. 1) THEN
          CALL strang1(it)
        END IF
!
        IF (msave) THEN
          CALL save_v(it)
        END IF
!
        IF (msave .and. l_root) CALL save_c(it)
        IF (micut) THEN
          CALL dealias
        END IF
      END IF
    END DO
!
    CALL get_max
    CALL get_dt(it,iti)
  END DO
!
  te_mpi = mpi_wtime()
!
  WRITE (6,9997) (te_mpi - ts_mpi)
  CALL mpi_finalize(ierr)
!
STOP
!
! ---------------------------------- FORMAT ---------------------------------- !
 9997 FORMAT(' Job Execution Time = ',e15.6)
!------------------------------------------------------------------------------!
!
END PROGRAM les_mpi
