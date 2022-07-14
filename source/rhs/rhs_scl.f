      subroutine rhs_scl(istep,iscl)
c
c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize)
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
c
      real fnt1(nnx,iys:iye,izs:ize)
      real tx(nnx,iys:iye), ty(nnx,iys:iye,izs:ize)
      real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
      real taut3_u(nnx,iys:iye,nscl), taut3_l(nnx,iys:iye,nscl) 
      real Sc, tscal, kconst
c      integer avalue_on = 0
c
c --------- set sign for ocean simulations that use monotone
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      upwn = 2.0
      if(iupwnd .ne. 1) upwn = 1.0
	
c
c --------- outer loop over z
c
      do iz=izs,ize
c
      izm2 = iz - 2
      izm1 = iz - 1
      izp1 = iz + 1
      izp2 = iz + 2
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
      weit3 = dzw(izm1)/(dzw(iz) + dzw(izm1))
      weit4 = 1.0 - weit3
      dzw2_i = 1.0/(dzw(iz) + dzw(izp1))
      dzw3_i = 2.0*dzw2_i
c
      do iy=iys,iye
      do ix=1,nnx
         tx(ix,iy) = t(ix,iy,iscl,iz)
      enddo
      enddo
      call xderivp(tx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
c
c --------- compute tau_t3 at iz-1
c
      if (iz.ne.1 .or. ibcl.ne.0) then
         do iy=iys,iye
         do ix=1,nnx
            taut3_l(ix,iy,iscl) = -vis_sv(ix,iy,izm1)*
     +              (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz)
         enddo
         enddo
      else
         do iy=iys,iye
         do ix=1,nnx
            taut3_l(ix,iy,iscl) = taut3m(ix,iy,iscl)
         enddo
         enddo
      endif
c
c ---------- SGS tau_t1, tau_t3 and resolved u*theta scalar fluxes
c            skew symmetric advective term 0.5(udt/dx + d/dx(ut))
c
      do iy=iys,iye
      do ix=1,nnx
         taut3_u(ix,iy,iscl) = -vis_sv(ix,iy,iz)*(t(ix,iy,iscl,izp1) -
     +                      t(ix,iy,iscl,iz))*dzu_i(izp1)
         fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iz)+vis_s(ix,iy,izm1))*
     +                    tx(ix,iy) - upwn*t(ix,iy,iscl,iz)
     +                      *(u(ix,iy,iz)+stokes(iz)*dir_x))
      enddo
      enddo

      call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)

      do iy=iys,iye
      do ix=1,nnx
         r4(ix,iy,iscl,iz) = -fnt1(ix,iy,iz)
     +           -(taut3_u(ix,iy,iscl)-taut3_l(ix,iy,iscl))*dzw_i(iz)
      enddo
      enddo
c
      if(iupwnd .ne. 1) then
c
c --------- skew symmetric advective form for
c           vertical flux = 0.5(wdt/dz + d/dz(wt))
c
      do iy=iys,iye
      do ix=1,nnx
         theta_u = weit1*t(ix,iy,iscl,iz) +
     +                weit*t(ix,iy,iscl,izp1)
         theta_l = weit3*t(ix,iy,iscl,iz) +
     +                weit4*t(ix,iy,iscl,izm1)
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +     -0.5*(u(ix,iy,iz)+stokes(iz)*dir_x)*tx(ix,iy)
     +     -0.5*(w(ix,iy,iz)*theta_u - w(ix,iy,izm1)*theta_l)*dzw_i(iz)
c
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +     -0.25*(w(ix,iy,iz)*
     +       (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1) +
     +            w(ix,iy,izm1)*
     +       (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz))
c
      enddo
      enddo
c
      else
c
c ----------- z-direction special
c
         if(iz .eq. 1) then
              do iy=iys,iye
              do ix=1,nnx 
                if(iscl.eq.2)then ! air-sea flux bc
                   Sc = 0.0d0
                   kconst = 0.0d0
                   tscal = 0.0d0
                   tscal = t(ix,iy,1,iz) - 273.15d0
                   Sc = 2073.1d0 - 125.62d0*tscal +
     +                  3.6276d0*tscal*tscal - 
     +                  0.043219d0*tscal*tscal*tscal ! calculate schmidt number (Table A1 in Wanninkof, 1992 for co2)
                   kconst = (2.77778d-6)*
     +                  0.31d0*5.75d0*5.75d0*sqrt(660.0d0/Sc) ! calculate piston velocity (Eq. 3 in Wanninkof, 1992)
                   flux_l(ix,iy) = (kconst)*
     +                  (8.325d0-t(ix,iy,iscl,iz)) ! calculate surface flux rate
                else
                   flux_l(ix,iy) = sgn*0.5*w(ix,iy,izm1)*
     +                  (t(ix,iy,iscl,izm1)+t(ix,iy,iscl,iz))
                endif
c
                flux_u(ix,iy) =
     +               amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +               rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +               t(ix,iy,iscl,izm1))) +
     +               amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +               rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +               t(ix,iy,iscl,izp2)))
             enddo
             enddo
         else if(iz .eq. nnz) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) = sgn*0.5*w(ix,iy,iz)*
     +                        (t(ix,iy,iscl,izp1)+t(ix,iy,iscl,iz))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         else
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         endif
c
c ---------- sum vertical monotone flux
c
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +           - sgn*(flux_u(ix,iy) - flux_l(ix,iy))*dzw_i(iz)
         enddo
         enddo
c
c -------- end monotone if block
c
      endif
c
c -------- save SGS fluxes for printout, gather sums on exit
c
      if(istep .eq. 1) then
         utsb(iz,iscl) = 0.0
         wtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            wtsb(iz,iscl) = wtsb(iz,iscl) + taut3_u(ix,iy,iscl)
            utsb(iz,iscl) = utsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iz)+vis_s(ix,iy,izm1))*tx(ix,iy)
         enddo
         enddo
         utsb(iz,iscl) = utsb(iz,iscl)*fnxy
         wtsb(iz,iscl) = wtsb(iz,iscl)*fnxy
      endif
c
c ---------- end z loop
c
      enddo
c
c --------- outer loop over z for y-depenence
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         ty(ix,iy,iz)  = t(ix,iy,iscl,iz)
      enddo
      enddo
      enddo
c
c --------- y-derivative of t for [izs:ize]
c
      call yd_mpi(ty(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c ------------- add skew symmetric advective flux and SGS flux
c               to y-derivative computation. check for monotone
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iz)+vis_s(ix,iy,izm1))*
     +                 ty(ix,iy,iz) - upwn*t(ix,iy,iscl,iz)
     +                      *(v(ix,iy,iz)+stokes(iz)*dir_y))
         enddo
         enddo
         if(iupwnd .ne. 1) then
           do iy=iys,iye
           do ix=1,nnx
              r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) -
     +                   0.5*(v(ix,iy,iz)+stokes(iz)*dir_y)*ty(ix,iy,iz)
           enddo
           enddo
         endif
      enddo
c
c --------- y-derivatives of scalar fluxes for [izs:ize]
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - fnt1(ix,iy,iz)
         enddo
         enddo
      enddo
c
c --------- reaction sources for for [izs:ize], for slow reactions (tau>1000)
c           see SOLVE/STRANG.F for fast reaction sources
c
!      if((flg_reaction.eq.1).and.(iscl.ge.2)) then
!        do iz=izs,ize
!           do iy=iys,iye
!              do ix=1,nnx
!                 r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) + 
!     +                react_src(ix,iy,iscl,iz)
!              enddo
!           enddo
!        enddo
!      end if
c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
      do iz=izs,ize
         vtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vtsb(iz,iscl) = vtsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iz)+vis_s(ix,iy,izm1))*ty(ix,iy,iz)
         enddo
         enddo
         vtsb(iz,iscl) = vtsb(iz,iscl)*fnxy
      enddo
      endif
c
      return
      end
