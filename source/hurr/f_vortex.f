      subroutine f_vortex(xo,yo,x,y,u_10,v_10,cd_10,tau_x,tau_y)
c
c ------------- simple routine to generate wind and stress
c               fields for hurricane frances according to the
c               zedler thesis on page 112. wind speeds in m/s
c               and radius in meters
c
      vr_rule(r) = 53.86*exp(-5.4e-06*r) + 2.86
      ur_rule(r) = 13.27*exp(-3.5e-06*r) - 3.33
c
c ------ get wind speeds and stress vectors
c
      rpos   = sqrt((x-xo)**2 + (y-yo)**2)
      cos_th = (x-xo)/rpos
      sin_th = (y-yo)/rpos
c
c ---------- radial and azmuthal formula
c
      rmax   = 40000.0
      ur_max = ur_rule(rmax)
      vr_max = vr_rule(rmax)
      ur_20  = ur_rule(rmax*20.0)
      vr_20  = vr_rule(rmax*20.0)
c
      if(rpos .le. rmax) then
c
c ----------- inside radius of max winds
c
         alph = alog(ur_max + 1.0)/rmax
         ur   = exp(alph*rpos) - 1.0
         alph = alog(vr_max + 1.0)/rmax
         vr   = exp(alph*rpos) - 1.0
      elseif(rpos .gt. rmax .and. rpos .le. 20.0*rmax) then
c
c ---------- decaying region
c
         ur = ur_rule(rpos)
         vr = vr_rule(rpos)
      elseif(rpos .gt. 20.0*rmax .and. rpos .le. 22.5*rmax) then
c
c ---------- decaying tail, zero at 22.5*rmax
c
         weit = 1.0 - (rpos - 20.0*rmax)/(2.5*rmax)
         ur   = ur_20*weit
         vr   = vr_20*weit
      else
         ur   = 0.0
         vr   = 0.0
      endif
c
c ---------- convert these into u_10,v_10
c
      u_10 = -(ur*cos_th + vr*sin_th)
      v_10 = -ur*sin_th + vr*cos_th
      s_10 = sqrt(u_10**2 + v_10**2)
c
c -------- fancy drag rule with limited cd
c
c     cd_fac = 1.1
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
