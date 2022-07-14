      subroutine stokes_hurr
c
c ----------- update stokes drift for a hurricane
c             assuming simple rules, constant wave age,
c             constant langmuir number, and exponential
c             decay profile
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
c --------- need to get the water stress
c           for use in turbulent langmuir number
c
      call f_vortex(x_hurr,y_hurr,xpt_les,ypt_les,
     +              u_10,v_10,cd_10,tau_x,tau_y)
      tau_unif = rho_a*sqrt(tau_x**2 + tau_y**2)
      utau_st  = sqrt(tau_unif/rho_w)
c
c -------- get amplitude of the stokes drift
c          and wavenumber of the peak
c
      stokess = utau_st/(turb_la**2)
      wind_st = sqrt(u_10**2 + v_10**2)
      stokesw = grav/((cpou10*wind_st)**2)
c
c ------- direction of stokes drift
c
      dir_x = u_10/wind_st
      dir_y = v_10/wind_st
c
c ------- every cpu has stokes profile, no stokes here
c
c     do iz=1,nnzp1
c        stokes(iz) = 0.0
c     enddo
      do iz=1,nnzp1
         stokes(iz) = stokess*exp(2.0*stokesw*zz(iz))
      enddo
c     if(l_root) then
c       write(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
c6000   format(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))
c     endif
c
      return
      end
