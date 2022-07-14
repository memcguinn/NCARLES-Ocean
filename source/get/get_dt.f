      subroutine get_dt(it,istart)
c
c ---------- routine computes max time step for given cfl number
c            from max's found previously
c
      use pars
      use con_data
      use con_stats
c
      data dt_max /30.0/
      save dt_max
c
      ucfl = umax
      vcfl = vmax
      wcfl = wmax
      ucflm = ucfl
      vcflm = vcfl
      wcflm = wcfl
      vel_max = wcflm
      vel_max = amax1(ucflm,vel_max)
      vel_max = amax1(vcflm,vel_max)
      if(it .eq. istart) then
         vel_max = cfl/dt_max
      endif
      if(vel_max .le. 0.0) then
          write(6,6000) ucflm, vcflm,wcflm, vel_max, cfl
 6000     format('6000, sr. get_dt bad news, umax = ',e15.6,/,
     +           ' vmax = ',e15.6,' wmax = ',e15.6,/,
     +           ' vel_max = ',e15.5,/,
     +           ' it = ', e15.6,/,
     +           ' infinite time step !!!')
          stop
      endif
c
c ---------------- choose fixed or variable time step
c
      if(ifix_dt .ne. 0) then
c
c ------------- if used, change to fit your problem
c
        dt_new = 10.0
      else
c
c ------------------- new estimate of best time step
c                     from cfl constraint
c
        dt_new = amin1(cfl/vel_max, dt_max)
      endif
c
c ---------------- compare against viscous stability limit
c
      if(vismax*dt_new .gt. 0.5) then
         dt_cfl = dt_new
         dt_new = 0.5/vismax
         if(l_root) then
            write(6,6200) dt_new, dt_cfl, vismax
 6200       format(' 6200 get_dt: cfl time step too large',/,
     +      '   viscous time step = ',e15.6,
     +      ' cfl time step = ',e15.6,' vismax = ',e15.6)
         endif
      endif
c
c -------- for safety if restart set timestep = saved timestep in
c          restart file
c
      if(it .eq. iti .and. iti .ne. 0) then
        dt_new = dt
      endif
c	if(dt_new < 0.25) then
c		print *, 'dt too small,dt=',dt_new
c		stop
c	endif

c
      return
      end
