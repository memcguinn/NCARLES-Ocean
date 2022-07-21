      subroutine speed2stress(u_10,v_10,cd_10,tau_x,tau_y)

c
c ----------- Specify speed here KEEP in 5-10 m/s range
c
      u_10 = 5.75
      v_10 = 0.0
      s_10 = sqrt(u_10**2 + v_10**2)
c
c -------- fancy drag rule with limited cd
c
      cd_fac = 0.7
      arg1   = (s_10 - 25.0)/5.0
      arg2   = -25.0/5.0
      a1     = tanh(arg1)
      a2     = tanh(arg2)
      cd_10  = cd_fac*(a1*0.5 - a2*0.5 + 1.2)*0.001
      tau_x  = cd_10*s_10*u_10
      tau_y  = cd_10*s_10*v_10
c
      return
      end
