      subroutine xy_stats
c
c ------------ get statistics
c
      use pars
      use fields
      use con_data
      use con_stats
c
c ------- indices for indexing array stat(.,.)
c         js = number of non-scalar stats
c         ns = number of scalar stats
c
      parameter(js = 10, ns = 5, nstat = js + ns*nscl)
      real stat(1:nnz,nstat)
c
c -------- stat(.,1) = u*u = ups
c          stat(.,2) = v*v = vps
c          stat(.,3) = w*w = wps
c          stat(.,4) = w**3 = wcube
c          stat(.,5) = w**4 = wfour
c          stat(.,6) = resolved tke at zw = englez
c          stat(.,7) = sgs e at zu = engsbz
c          stat(.,8) = sgs e at zw = eavg
c          stat(.,9) = resolved uw at zw = uwle
c          stat(.,10) = resolved vw at zw = vwle
c          stat(.,m1) = resolved scalar flux wt at zw = wtle
c          stat(.,m2) = resolved scalar flux ut at zw = utle
c          stat(.,m3) = resolved scalar flux vt at zw = vtle
c          stat(.,m4) = scalar t*t at zu = tps
c          stat(.,m5) = scalar t*t*t at zu = tcube
c
c --------- use a trick with mpi reduce over all z to get averages
c           by setting stat array = 0 for all z on each process
c
      do i=1,nstat
      do iz=1,nnz
         stat(iz,i) = 0.0
      enddo
      enddo
c
c -------- indices for scalars
c
      m1 = js
      m2 = js + nscl
      m3 = js + 2*nscl
      m4 = js + 3*nscl
      m5 = js + 4*nscl
c
      sgn = 1.0
      if(iocean .eq. 1 .and. iupwnd .eq. 1) sgn = -1.0
c
      do iz=izs,ize
c
      izp2 = iz + 2
      izp1 = iz + 1
      izm1 = iz - 1
c
      do iy=iys,iye
      do ix=1,nnx
         stat(iz,1) = stat(iz,1) + (u(ix,iy,iz) - uxym(iz))**2
         stat(iz,2) = stat(iz,2) + (v(ix,iy,iz) - vxym(iz))**2
         stat(iz,3) = stat(iz,3) + (w(ix,iy,iz) - wxym(iz))**2
         stat(iz,4) = stat(iz,4) + (w(ix,iy,iz) - wxym(iz))**3
         stat(iz,5) = stat(iz,5) + (w(ix,iy,iz) - wxym(iz))**4
         stat(iz,6) = stat(iz,6) +
     +                ((w(ix,iy,iz)-wxym(iz))**2 +
     +                (0.5*(u(ix,iy,iz)-uxym(iz) +
     +                      u(ix,iy,izp1)-uxym(izp1)))**2 +
     +                (0.5*(v(ix,iy,iz)-vxym(iz) +
     +                      v(ix,iy,izp1)-vxym(izp1)))**2)*0.5
c
         stat(iz,7) = stat(iz,7) + 0.5*(e(ix,iy,iz)+e(ix,iy,izm1))
         stat(iz,8) = stat(iz,8) + e(ix,iy,iz)
         stat(iz,9) = stat(iz,9) + (w(ix,iy,iz)-wxym(iz))*
     +              0.5*((u(ix,iy,iz)-uxym(iz))+
     +                   (u(ix,iy,izp1)-uxym(izp1)))
         stat(iz,10) = stat(iz,10) + (w(ix,iy,iz)-wxym(iz))*
     +              0.5*((v(ix,iy,iz)-vxym(iz))+
     +                   (v(ix,iy,izp1)-vxym(izp1)))
      enddo
      enddo
c
c ------------ get scalar resolved fluxes and variances
c
      do l=1,nscl
         if(iupwnd .ne. 1 .or. iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
               stat(iz,m1+l)=stat(iz,m1+l) +
     +               (w(ix,iy,iz)-wxym(iz))*
     +               0.5*(t(ix,iy,l,iz)-txym(iz,l) +
     +                    t(ix,iy,l,izp1)-txym(izp1,l))
            enddo
            enddo
         else
c
c ------------------- monotone fluxes
c
           do iy=iys,iye
           do ix=1,nnx
              stat(iz,m1+l) = stat(iz,m1+l) +
     +    amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,l,iz) +
     + rlim(t(ix,iy,l,izp1),t(ix,iy,l,iz),t(ix,iy,l,izm1))) +
     +    amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,l,izp1) +
     + rlim(t(ix,iy,l,iz),t(ix,iy,l,izp1),t(ix,iy,l,izp2)))
           enddo
           enddo
         endif
         stat(iz,m1+l)= sgn*stat(iz,m1+l)
c
c ------------ get horizontal scalar resolved fluxes
c
         do iy=iys,iye
         do ix=1,nnx
            stat(iz,m2+l) = stat(iz,m2+l)+
     +               (u(ix,iy,iz)-uxym(iz))*
     +               (t(ix,iy,l,iz)-txym(iz,l))
            stat(iz,m3+l) = stat(iz,m3+l)+
     +               (v(ix,iy,iz)-vxym(iz))*
     +               (t(ix,iy,l,iz)-txym(iz,l))
         enddo
         enddo
c
c ------------------- scalar variances & higher moments
c
         do iy=iys,iye
         do ix=1,nnx
            stat(iz,m4+l) = stat(iz,m4+l) +
     +                (t(ix,iy,l,iz) - txym(iz,l))**2
            stat(iz,m5+l) = stat(iz,m5+l) +
     +                (t(ix,iy,l,iz) - txym(iz,l))**3
         enddo
         enddo
c
c ------ end scalar loop
c
      enddo
c
c ------ end z loop
c
      enddo
c
c -------- add partial sums and send it to all
c
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*nstat,1)
c
c ------ fill arrays for printout and constant file
c
      do iz=1,nnz
c
      ups(iz)    = stat(iz,1)*fnxy
      vps(iz)    = stat(iz,2)*fnxy
      wps(iz)    = stat(iz,3)*fnxy
      wcube(iz)  = stat(iz,4)*fnxy
      wfour(iz)  = stat(iz,5)*fnxy
      englez(iz) = stat(iz,6)*fnxy
      engsbz(iz) = stat(iz,7)*fnxy
      eavg(iz)   = stat(iz,8)*fnxy
      uwle(iz)   = stat(iz,9)*fnxy
      vwle(iz)   = stat(iz,10)*fnxy
      uw_tot(iz) = uwle(iz) + uwsb(iz)
      vw_tot(iz) = vwle(iz) + vwsb(iz)
c
c ------------ get scalar resolved fluxes and variances
c
      do l=1,nscl
         wtle(iz,l)   = stat(iz,m1+l)*fnxy
         utle(iz,l)   = stat(iz,m2+l)*fnxy
         vtle(iz,l)   = stat(iz,m3+l)*fnxy
         tps(iz,l)    = stat(iz,m4+l)*fnxy
         tcube(iz,l)  = stat(iz,m5+l)*fnxy
         wt_tot(iz,l) = wtle(iz,l) + wtsb(iz,l)
      enddo
      enddo
c
      return
      end
