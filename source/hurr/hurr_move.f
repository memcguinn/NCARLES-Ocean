      subroutine hurr_move
c
c ----------- update movement of vortex location
c             in the vertical direction
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      if(time .gt. t_move) then
        y_hurr = y_hurr_i + h_speed*(time - t_move)
      else
        y_hurr = y_hurr_i
      endif
c
      return
      end
