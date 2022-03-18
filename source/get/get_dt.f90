! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine computes max time step for a given CFL number,       !
!         previously computed.                                                 !
! ============================================================================ !
!
SUBROUTINE get_dt(it,istart)
!
    USE inputs, ONLY: iti
    USE pars
    USE con_data
    USE con_stats
!
    DATA dt_max /30.0/
    SAVE dt_max
!
! ---------------------------------------------------------------------------- !
!
  ucfl = umax
  vcfl = vmax
  wcfl = wmax
  vel_max = MAX(ucfl,vcfl,wcfl)
!
  IF (it .EQ. istart) THEN
    vel_max = cfl/dt_max
  END IF
!
  IF (vel_max .LE. 0.0) THEN
    WRITE (6,6000) ucfl, vcfl,wcfl, vel_max, cfl
    STOP
  END IF
!
! CHOOSE FIXED/VARIABLE TIME STEP
  IF (ifix_dt .NE. 0) THEN
    dt_new = 10                       ! WARNING: IF USED, CHANGE TO FIT PROBLEM
  ELSE
!
! NEW ESTIMATE OF BEST TIME STEP FROM CFL CONSTRAIN
    dt_new = AMIN1(cfl/vel_max, dt_max)
  END IF
!
! COMPARE AGAINST VISCOUS STABILITY LIMIT
  IF (vismax*dt_new .GT. 0.5) THEN
    dt_cfl = dt_new
    dt_new = 0.5/vismax
    IF (l_root) THEN
      WRITE (6,6200) dt_new, dt_cfl, vismax
    END IF
  END IF
!
! CHECK: IF RESTART, SET TIMESTEP = SAVED TIMESTEP IN RESTART FILE
  IF (it .EQ. iti .AND. iti .NE. 0) THEN
    dt_new = dt
  END IF

RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
6000 FORMAT ('6000, sr. get_dt bad news, umax = ',e15.6,/, &
             ' vmax = ',e15.6,' wmax = ',e15.6,/,          &
             ' vel_max = ',e15.5,/,                        &
             ' it = ', e15.6,/,                            &
             ' infinite time step !!!')
6200 FORMAT (' 6200 get_dt: cfl time step too large',/,    &
             '   viscous time step = ',e15.6,              &
             ' cfl time step = ',e15.6,' vismax = ',e15.6)
! ---------------------------------------------------------------------------- !
!
END SUBROUTINE get_dt
