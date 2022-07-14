      subroutine get_means(istage)
c
c ------------ get means for all variables
c              for use in iso, surfvis, comp1, compmn.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      do iz=0,nnz+1
         u_mn(iz)   = 0.0
         v_mn(iz)   = 0.0
         w_mn(iz)   = 0.0
         engz(iz)   = 0.0
         engsbz(iz) = 0.0
         divz(iz)   = 0.0
      enddo
      do iscl=1,nscl
         do iz=0,nnz+1
            t_mn(iz,iscl) = 0.0
         enddo
      enddo
      iz_ee = ize
      if(ize .eq. nnz) iz_ee = nnzp1
      do iz=izs,iz_ee
         do iy=iys,iye
         do ix=1,nnx
            u_mn(iz) = u_mn(iz) + u(ix,iy,iz)
            v_mn(iz) = v_mn(iz) + v(ix,iy,iz)
            w_mn(iz) = w_mn(iz) + w(ix,iy,iz)
         enddo
         enddo
         u_mn(iz) = u_mn(iz)*fnxy
         v_mn(iz) = v_mn(iz)*fnxy
         w_mn(iz) = w_mn(iz)*fnxy
         do iscl=1,nscl
            t_mn(iz,iscl) = 0.0
            do iy=iys,iye
            do ix=1,nnx
               t_mn(iz,iscl) = t_mn(iz,iscl) + t(ix,iy,iscl,iz)
            enddo
            enddo
            t_mn(iz,iscl) = t_mn(iz,iscl)*fnxy
         enddo
      enddo
      call mpi_sum_z(u_mn(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(v_mn(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(w_mn(1),i_root,myid,nnzp1,1)
      do iscl=1,nscl
         call mpi_sum_z(t_mn(1,iscl),i_root,myid,nnzp1,1)
      enddo
c
c -------- set e to minimum value
c
      do iz=izs-1,ize+1
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz )=amax1(e(ix,iy,iz ),sml_eg)
         enddo
         enddo
      enddo
c
c -------------- get terms which contribute to mean pressure
c                careful with the sum
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            engz(iz)  = engz(iz) +
     +   .5*((u(ix,iy,iz) + stokes(iz)*dir_x)**2 +
     +       (v(ix,iy,iz) + stokes(iz)*dir_y)**2 +
     +    .5*(w(ix,iy,iz)*w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
            engsbz(iz) = engsbz(iz) +
     +                    0.5*(e(ix,iy,iz)+e(ix,iy,izm1))
         enddo
         enddo
         engz(iz)   = engz(iz)*fnxy
         engsbz(iz) = engsbz(iz)*fnxy
      enddo
      call mpi_sum_z(engz(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(engsbz(1),i_root,myid,nnzp1,1)
c
c ------------ save means and divergence for printout and compmn
c              all cpus have means over all z
c
      if(istage .eq. 1) then
        do iz=izs,ize
           izm1 = iz - 1
           do iy=iys,iye
           do ix=1,nnx
              divz(iz) = divz(iz) +
     +                  (ux(ix,iy,iz)+vy(ix,iy,iz)+
     +                  (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz))**2
           enddo
           enddo
           divz(iz) = divz(iz)*fnxy
        enddo
        call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
        do iz=0,nnz+1
           uxym(iz) = u_mn(iz)
           vxym(iz) = v_mn(iz)
           wxym(iz) = w_mn(iz)
        enddo
        do iscl=1,nscl
           do iz=0,nnz+1
              txym(iz,iscl) = t_mn(iz,iscl)
           enddo
        enddo
      endif
c
      return
      end
