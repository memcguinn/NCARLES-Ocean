      subroutine init
c
      use pars
      use fields
      use con_data
      use con_stats
c
      pi     = 4.0*atan(1.0)
      pi2    = 2.0*pi
      d_to_r = pi2/360.0
      grav   = 9.81
      bfac   = 1.0
      if(ibuoy.eq.0) bfac = 0.
c
c -------------------- case specific data
c
         rho_a   = 1.0
         rho_w   = 1000.0
         t00     = 283.
         t00b    = 5000.0
         cp_a    = 1.0057e03
         cp_w    = 4.20e03
         gcp     = grav/cp_w
         batag   = bfac*grav/t00b
c
c --------- specify stokes drift parameters
c
         cpou10  = 0.6
         turb_la = 0.3
c
         rlat    = 30

         fcor    = 2.0*pi2*sin(rlat*d_to_r)/(24.0*3600.0)
         fcor_h  = 0.0
         ugcont  = 0.0
         vgcont  = 0.
        wtsfc(1) = 0.0 !       5.0e-7
         qstar(1) = wtsfc(1)
c
c ------- other thermodynamic variables for hurricane
c         time varying heat fluxes
c
         c_e    = 0.0015         ! exchange coefficient for moisture
         c_h    = 0.0015         ! exchange coefficient for temperature
         rlv    = 2.45e06        ! latent heat of vaporization J/Kg
         p_surf = 1000.0         ! surface pressure in millibars
         t_surf = 273.15 + 26.0  ! surface temperature in deg K
         q_rel  = 0.8            ! relative humidity in % at 10m height
c
         t_10m  = t_surf - 2.5    ! 10m temperature
c
c -------- saturation vapor pressure from Gill, p. 606, A.4.5
c
         e_val = (0.7859 + 0.03477*(t_surf-273.15))/
     +           (1.0 + 0.00412*(t_surf-273.15))
         e_sat = 10.0**e_val
         e_rat = e_sat/p_surf
c
         r_sat = e_rat*0.62197/(1.0 - e_rat)
         q_sat = r_sat/(1.0 + r_sat)
c
c -------- find mixing ratio at 10m height
c
         r_10m = q_rel*r_sat
         q_10m = r_10m/(1.0 + r_10m)
c
c -------- get atmosphere and ocean heat fluxes without wind
c
         q_lat_a   = rho_a*c_e*rlv*(q_sat - q_10m)
         q_sen_a   = rho_a*c_h*(t_surf - t_10m)*cp_a
         qw_tot_aw = (q_lat_a + q_sen_a)/(rho_w*cp_w)
         if(l_root) write(6,7676) qw_tot_aw
 7676    format(' in init qw_tot_aw = ',e15.6)
c
         dtdzf(1)=0.010
         dtjump  = 0.
         divgls  = 0.
         zo      = 0.0001
         zi      = -30
         xl      = 320.0
         yl      = 320.0
         zl      = -96.0
         izi     = nint((zi/zl)*nnz)
c
c ---------- if stretched grid specify location of first point
c
         zw1 = -0.5

c
c ----------- donelan parameters
c
         ann      = 0.00615
         bnn      = 1.0
c        f2w      = 0.123
         f2w      = 0.13
         f_p      = f2w*grav/u_10
         npm      = 4
         sigma_p  = pi2*f_p
         grav_w   = grav

c
c -------- ratio of k_1/k_p ala Phillips and Donelan
c
         r_kp     = grav*(f2w*pi2/u_10)**2
         r_k1     = r_fac*grav/(cd_10*u_10*u_10)
         rk_ratio = r_k1/r_kp

c
      time  = 0.0
c
c ---------- outermost coarse grid  indicies are bounds of grid
c
      izlow = 1
      izup  = nnz
      dz    = zl/nnz
      dzg   = abs(dz)
      if(l_root) write(6,4040) zl,nnz,dzg
c
c --------------- generate z grids for particular mesh from
c                 iz = 0,1,...,nnz+1; this allows indexing
c                 to array elements z(0), etc.
c
      zwstrt = 0.0
c
c ------------ if uniform vertical spacing then
c
      if(iz_space .eq. 0) then
c
c ------------ build z grid for w points
c
         do iz=0,nnz+1
            z(iz) = dz*float(iz) + zwstrt
         enddo

      else
        call vgrid(zw1,zi,zl,nnz,z(0),l_root,l_debug)
      endif
c
      call get_dz
c
      if(l_root) then
         write(6,8002) zwstrt
         write(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
      endif
c
      nnzm1 = nnz-1
      dx    = xl/nnx
      dy    = yl/nny
      fnxy  = 1./float(nxy)
      dzdz  = dzw(1)*dzw(1)
      z1    = zz(1)
c
      c23  = 2.0/3.0
      dsl  = (dx*1.5*dy*1.5*abs(dzw(1)))**(1./3.)
      dslg = dsl
      cs   = 0.2
c
      vk     = 0.4
      batagk = batag*vk
      vkin   = 1./vk
      ttmean = 0.
      zody   = alog(abs(z1/zo))
      write(nprt, 9901) z1,zo,zody
 9901 format(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,
     +       ' zody = ',e15.6)
      zodyin = 1./zody
      wstar  = abs(batag*zi*wtsfc(1))**(1./3.)
      if(ismlt .eq. 1) then
c
c ---- set constants for businger similarity functions
c
         vk74   = vk*0.74
         vk74in = 0.74/vk
         zody74 = zody*0.74
      else
c
c ---- set constants for large similarity functions
c
        vk74    = vk
        vk74in  = 1.0/vk
        zody74  = zody
      endif
      ugal   = ugcont*0.5
      cdbtm  = vk*vk/zody/zody
c ----------- set surface friction velocity here and in sr. sufto
         utau = sqrt(rho_a*(8.5e-4)*5.75*5.75/rho_w)

      utau2    = utau*utau
      if(ibuoy .eq. 0 .or. qstar(1) .eq. 0.) then
        amonin = 1000.0
      else
        amonin = -utau2*utau/(batagk*qstar(1))
      endif
      hol   = abs(zi)/amonin
      zol   = abs(z1)/amonin
      uwsfc = -utau*utau
      vwsfc = -utau*utau
c ------- make sure tsfcc is gt than t00 for both isfc=0 or 1
      tsfcc(1) = 265.00
c
      if(l_root) then
         write(6,80)
         write(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      if(l_debug) then
         write(nprt,80)
         write(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      return
c ------------------------
   2  format(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/,
     +       10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,
     +       ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4,
     +       ',  UG = ',e12.4)
  80  format(///,' ***** SCRATCH RUN ***** ',//)
 4040 format(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043 format(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 8002 format(' zwstrt = ',e12.4)
 8003 format(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
      end
