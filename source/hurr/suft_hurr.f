      subroutine suft_hurr(it)
c
c --------- surface layer routine for high wind
c           hurricane driven obls
c
      use pars
      use fields
      use con_data
      use con_stats
      real buf(3+nscl)
c
c ------- version of similarity theory adpated for ocean flows
c      option to use businger or large version of similarity theory
c
      iz      = 1
      izm1    = iz - 1
      izp1    = iz + 1
      z1_a    = abs(z1)
      buf(1)  = 0.0
      buf(2)  = 0.0
      buf(3)  = 0.0
      tol     = 0.01
      do iy=iys,iye
      do ix=1,nnx
         buf(1) = buf(1) + u(ix,iy,iz)
         buf(2) = buf(2) + v(ix,iy,iz)
         wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +                    +v(ix,iy,iz)*v(ix,iy,iz))
         buf(3) = buf(3) + wind(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         buf(3+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
         enddo
         enddo
      enddo
c
c -------- get x-y slab sums
c
      call mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
      u1xy  = buf(1)*fnxy + ugal
      v1xy  = buf(2)*fnxy
      windm = buf(3)*fnxy
      do iscl=1,nscl
         t1xy(iscl) = buf(3+iscl)*fnxy
      enddo
      vsfc  = sqrt(u1xy*u1xy+v1xy*v1xy)
      windm = amax1(windm,ufree)
      vsfc  = amax1(vsfc,ufree)
c
c      call f_vortex(x_hurr,y_hurr,xpt_les,ypt_les,
c     +              u_10,v_10,cd_10,tau_x,tau_y)
c

      call speed2stress(u_10,v_10,cd_10,tau_x,tau_y)

c -------- get time varying surface heat flux
c          based on atmospheric winds!!
c

      qstar(1) = qw_tot_aw*sqrt(u_10**2 + v_10**2)

c ----------- get water u_* assuming rho_a/rho_w = 1/1000
c             and assign a fraction to unresolved breaking
c             and viscous piece (DBRK)
c
c     tau_unif = rho_a*cd_10*u_10*u_10*frac_v + tau_frac
      tau_unif = rho_a*sqrt(tau_x**2 + tau_y**2)
      utau     = sqrt(tau_unif/rho_w)
c
      t10xy(1)=-qstar(1)/utau*zody*vk74in

c
c ---- check for temperature boundary condition
c
      if(isfc .eq. 0 ) then
         tsfcc(1) = t1xy(1) - t10xy(1)
      endif
c
c ---- save old utau
c
      utausv = utau
      utau2  = utau*utau
      if (ibuoy.eq.0 .or. qstar(1) .eq. 0.) then
          amonin    = 1000.
          zeta      = 0.
          thstar(1) = 0.0
          t10xy(1)  = 0.0
      else
          amonin = -utau2*utau/(batagk*qstar(1))
          zeta   = z1_a/amonin
      endif
      if (t10xy(1).lt.0. .and. qstar(1) .lt. 0.) then
         write(6,1234)u1xy,v1xy,t1xy(1),tsfcc(1),amonin,utau,it
 1234    format(' ** check sfc u=',e12.3,' v=',e12.3,' t,ts=',2f10.3,
     +     ' l=',e12.3,' u*=',e12.3,' at it=',i5)
         go to 9999
      endif
c
c -------- for stable,neutral and unstable pbl get drift velocity
c
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      udrift = windm + stokes(1) - stokess + utau*(zody-psim)*vkin
      vdrift = 0.0
      dnom      = (zody-psis)*vk74in
      if (isfc.eq.1) then
         thstar(1) = (t1xy(1) - tsfcc(1))/dnom
         t10xy(1)  = thstar(1)*dnom
         qstar(1)  = - utau*thstar(1)
      else
         thstar(1)  = -qstar(1)/utau
         tsfcc(1)   = t1xy(1) - thstar(1)*dnom
         t10xy(1)   = thstar(1)*dnom
      endif
      zol = zeta
      hol = zol*zi/z1
c
c ------ pick correct sign for stress, account for density
c
      au13m    = rho_a*tau_x/rho_w
      au23m    = rho_a*tau_y/rho_w
      aut3m(1) = wtsfc(1)
c qstar(1)
c
      return
c
c --------- trouble in sl routine
c
 9999 continue
c
      write(nprt,9000)
 9000 format(' Trouble in SR. suft_hurr')
      call mpi_finalize(ierr)
      stop
      end
