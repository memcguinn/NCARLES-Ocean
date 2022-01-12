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
  USE inputs
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
  ts_mpi = mpi_wtime()
!
  ncpu_s = 32
  itn=0
!
  case_inp = '30L'
!
  CALL get_units
  CALL gridd
  CALL setcon
  CALL set_paths
  istop = 1
!
! ------------------------------- SCRATCH RUN -------------------------------- !
!
! CHOOSE ROUTINE FOR GETTING INITIAL GUESS
  IF (iti.EQ.0)  THEN
     igrdr = 2
     case = case_inp
     CALL init
     CALL setup(it)
!
!    CHOOSE ROUTINE FOR GETTING INITIAL GUESS
     CALL randoc
!
     CALL get_max
  ELSE
     igrdr = 3
     CALL restart
     CALL get_max
     CALL setup(it)
  END IF
!
! --------------------------------- TIME LOOP -------------------------------- !
  tzero = time
  CALL get_dt(it,iti)
  9000 CONTINUE
  CALL set_sav(it,iti)
!
  IF (it .GE. new_vis .AND. ivis0 .EQ. 1) THEN
      ivis = 1
  ELSE
      ivis = 0
  END IF
!
! 3 STAGE RUNGE-KUTTA TIME-STEPPING
  DO 8999 istage=1,3
!
  dtzeta = dt*zetas(istage)
  dtgama = dt*gama(istage)
!
! COMPUTE DERIVATIES OF (U,V,W)
  CALL exchange
  CALL get_derv
!
! NEW EDDY VISCOSITY, BOUNDARY CONDITIONS
  IF (iss .EQ. 0 .AND. ifree .EQ. 0) THEN
     CALL lower(it)
  ELSEIF (ifree .EQ. 1) THEN
     CALL lower_free(it)
  END IF
!
  IF (ise .EQ. numprocs-1) THEN
     CALL upper
  END IF
!
  CALL applytracerbc(it)
  CALL bcast_pbc
  CALL get_means(istage)
!
  IF (ivis .EQ. 1) THEN
     CALL iso(it)
     CALL surfvis(it)
  END IF
!
  IF (istage .EQ. 1) THEN
    CALL xy_stats
    CALL tke_budget
    CALL pbltop(itop)
  END IF
!
  IF (istage.EQ.1 .AND. flg_reaction.EQ.1) THEN
     CALL strang1(it)
  END IF
!
  CALL comp1(istage,it)
!
! ----------------------------- PRESSURE SOLVER ------------------------------ !
  CALL comp_p
!
! ADD PRESSURE GRADIENT AND DEALIAS
  CALL comp2
!
  IF (istage.EQ.3 .AND. flg_reaction.EQ.1) THEN
     CALL strang1(it)
  END IF
!
  IF (msave .AND. istage .EQ. 3) THEN
     CALL save_v(it)
  END IF
!
  IF (istage .EQ. 3) THEN
     IF (msave .AND. l_root) CALL save_c(it)
  END IF
!
  IF (micut) THEN
     CALL dealias
  END IF
!
  IF (mnout .AND. istage .EQ. 1) THEN
      IF (l_debug) THEN
         CALL PRINT(nprt,it,izs,ize)
      END IF
      IF (l_root) CALL PRINT(6,it,1,nnz)
  END IF
!
  IF (l_root) THEN
     IF (mhis .AND. istage .EQ. 1) CALL write_his(itop)
     IF (mhis .AND. istage .EQ. 1 .AND. mtape) CALL close_his
  END IF
!
  8999 CONTINUE
!
  CALL get_max
  CALL get_dt(it,iti)
!
  IF (it.GE.itmax) GO TO 99000
  GO TO 9000
!
  99000 CONTINUE
!
  te_mpi = mpi_wtime()
!
! ---------------------------------- FORMAT ---------------------------------- !
  WRITE(6,9997) (te_mpi - ts_mpi)
  9997 FORMAT(' Job Execution Time = ',e15.6)
!------------------------------------------------------------------------------!
!
  9998 CONTINUE
!
  CALL mpi_finalize(ierr)
!
STOP
!
END PROGRAM les_mpi
