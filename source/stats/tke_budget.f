      subroutine tke_budget
c
c -------- get terms in resolved scale tke budget
c          as in gabls writeup at w-points
c          at istage = 1.
c          t_diss and tr_tau are in comp1
c
      use pars
      use fields
      use con_data
      use con_stats
c
      real pxym(1:nnz), stat(1:nnz,2)
c
c -------- stat(.,1) = tke transport  = wq
c          stat(.,2) = pressure transport  = wp
c
      do iz=1,nnz
         pxym(iz)   = 0.0
         stat(iz,1) = 0.0
         stat(iz,2) = 0.0
      enddo
c
      do iz=izs,ize
c
c ----------------- get "mean p_star pressure"
c
         izm1 = iz - 1
c
         pxym(iz) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            pxym(iz) = pxym(iz) + p(ix,iy,iz)
     +           -(e(ix,iy,iz)+e(ix,iy,izm1))/3.0
     +           -0.5*((u(ix,iy,iz)+stokes(iz)*dir_x)**2 +
     +                 (v(ix,iy,iz)+stokes(iz)*dir_y)**2 +
     +      0.5*(w(ix,iy,iz)*w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
         enddo
         enddo
         pxym(iz) = pxym(iz)*fnxy
      enddo
      call mpi_sum_z(pxym(1),i_root,myid,nnz,1)
c
c --------------- get transport terms as vertical arrays
c
      do iz=izs,ize
c
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
c
c --------- get estimate of turbulent transport term
c
            ufluc   = u(ix,iy,iz) - uxym(iz)
            vfluc   = v(ix,iy,iz) - vxym(iz)
            wfluc   = w(ix,iy,iz) - wxym(iz)
            wfluc_l = w(ix,iy,izm1) - wxym(izm1)
            stat(iz,1)  = stat(iz,1) + 0.25*(wfluc + wfluc_l)*
     +             (ufluc**2 + vfluc**2 + 0.5*(wfluc_l**2 + wfluc**2))
c
c --------- get estimate of pressure transport term
c
            pfluc = p(ix,iy,iz) - pxym(iz)
     +           -(e(ix,iy,iz)+e(ix,iy,izm1))/3.0
     +           -0.5*((u(ix,iy,iz)+stokes(iz)*dir_x)**2 +
     +                 (v(ix,iy,iz)+stokes(iz)*dir_y)**2 +
     +      0.5*(w(ix,iy,iz)*w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
            stat(iz,2) = stat(iz,2) + pfluc*0.5*(wfluc_l + wfluc)
         enddo
         enddo
         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
      enddo
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*2,1)
c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
c
         izp1 = iz + 1
         izm1 = iz - 1
c
c ------- treat tr_tau at bottom special, tr_tau = 0.0
c
         if(iz .eq. 1) tr_tau(izm1) = 0.0
         if(iz .eq. nnz) then
            t_tau(iz) = 0.0
            t_wp(iz)  = 0.0
            t_wq(iz)  = 0.0
         else
            t_tau_u   = 0.5*(tr_tau(izp1) + tr_tau(iz))
            t_tau_l   = 0.5*(tr_tau(izm1) + tr_tau(iz))
            t_tau(iz) = -(t_tau_u - t_tau_l)*dzu_i(izp1)
            t_wq(iz)  = -(stat(izp1,1) - stat(iz,1))*dzu_i(izp1)
            t_wp(iz)  = -(stat(izp1,2) - stat(iz,2))*dzu_i(izp1)
         endif
         dudz = (uxym(izp1) - uxym(iz))*dzu_i(izp1)
         dvdz = (vxym(izp1) - vxym(iz))*dzu_i(izp1)
c
c ------------- gather all the budget terms
c
         t_tran(iz)  = t_wq(iz) + t_wp(iz) + t_tau(iz)
         t_rprod(iz) = -(dudz*uwle(iz) + dvdz*vwle(iz))
         t_sprod(iz) =  (dudz*uwsb(iz) + dvdz*vwsb(iz))
         t_buoy(iz)  =  batag*wtle(iz,1)
c
      enddo
c
      return
      end
