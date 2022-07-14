      subroutine init_hurr
c
c ------ initialization information for water les
c        driven by hurricane vortex
c
      use pars
      use fields
      use con_data
      use con_stats
c
c -------- initial location of vortex
c
      x_hurr    = 0.0
      y_hurr_i  = 0.0
      y_hurr    = y_hurr_i
c
c -------- location of les domain, relative to vortex location
c
      xpt_les = 55.0e03
      ypt_les = 700.0e03
c
c -------- vertical translation speed of storm
c
      h_speed = 5.5
c
c -------- start moving storm at time
c
      t_move = 10.0
c
      return
      end
