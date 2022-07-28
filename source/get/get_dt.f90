SUBROUTINE get_dt(it,istart)
! ROUTINE COMPUTES MAX TIME STEP FOR GIVEN CFL NUMBER

  USE pars
  USE con_data
  USE con_stats

  DATA dt_max /30.0/
  SAVE dt_max

  ucfl = umax
  vcfl = vmax
  wcfl = wmax
  ucflm = ucfl
  vcflm = vcfl
  wcflm = wcfl
  vel_max = wcflm
  vel_max = AMAX1(ucflm,vel_max)
  vel_max = AMAX1(vcflm,vel_max)

  IF(it == istart) THEN
    vel_max = cfl/dt_max
  ENDIF

  IF(vel_max <= 0.0) THEN
    WRITE(6,6000) ucflm, vcflm,wcflm, vel_max, cfl
    STOP
  ENDIF

  ! CHOOSE FIXED OR VARIABLE TIME STEP
  IF(ifix_dt /= 0) THEN
    ! IF USED, CHANGE TO FIT PROBLEM
    dt_new = 10.0
  ELSE
    ! NEW ESTIMATES OF BEST TIME STEP FROM CFL CONSTRAINT
    dt_new = AMIN1(cfl/vel_max, dt_max)
  ENDIF

  ! COMPARE AGAINST VISCOUS STABILITY LIMIT
  IF(vismax*dt_new > 0.5) THEN
    dt_cfl = dt_new
    dt_new = 0.5/vismax
    IF(l_root) THEN
      WRITE(6,6200) dt_new, dt_cfl, vismax
    ENDIF
  ENDIF

  ! FOR SAFETY IF RESTART SET TIMESTEP = SAVED TIMESTEP IN RESTART FILE
  IF(it == iti .AND. iti /= 0) THEN
    dt_new = dt
  ENDIF

  RETURN

! FORMAT
6000  FORMAT('6000, sr. get_dt bad news, umax = ',e15.6,/,' vmax = ',       &
            e15.6,' wmax = ',e15.6,/,' vel_max = ',e15.5,/,' it = ',        &
             e15.6,/, ' infinite time step !!!')
6200  FORMAT(' 6200 get_dt: cfl time step too large',/,                     &
            '   viscous time step = ',e15.6,' cfl time step = ',e15.6,'     &
            vismax = ',e15.6)

END SUBROUTINE
