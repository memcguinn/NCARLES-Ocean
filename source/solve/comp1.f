      subroutine comp1(istep,it)
c
c ----- 3-order runge-kutta time stepping and monotone scalar fluxes in x,y,z.
c       designed to use mpi in x & y directions.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      parameter(js = 7, ns = 3, nstat = js + ns*nscl)
      real stat(1:nnz,nstat)
c
c ------ temp arrays to hold rhs from step n-1 and
c        field variables from step n
c
      real urhs(nnx,iys:iye,izs:ize),
     +     vrhs(nnx,iys:iye,izs:ize),
     +     wrhs(nnx,iys:iye,izs:ize),
     +     erhs(nnx,iys:iye,izs:ize),
     +     trhs(nnx,iys:iye,nscl,izs:ize)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            urhs(ix,iy,iz) = u(ix,iy,iz) + dtzeta*r1(ix,iy,iz)
            vrhs(ix,iy,iz) = v(ix,iy,iz) + dtzeta*r2(ix,iy,iz)
            wrhs(ix,iy,iz) = w(ix,iy,iz) + dtzeta*r3(ix,iy,iz)
            erhs(ix,iy,iz) = e(ix,iy,iz) + dtzeta*r5(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx

               trhs(ix,iy,l,iz) = t(ix,iy,l,iz) + dtzeta*r4(ix,iy,l,iz)

         enddo
         enddo
         enddo
      enddo
c
c --------- get viscosity and rhs of (e,u,v,w)-equations
c           at next step
c
      call tke_vis(istep)
      call rhs_uvw(istep)

c -------- evaluate rhs of scalar equations

c
      do l=1,nscl
         call rhs_scl(istep,l,it)
      enddo

c ---------- gather stat sums on root processor
c            using mpi_reduction over all processors
c
      if(istep .eq. 1) then
c
        do j=1,nstat
        do iz=1,nnz
           stat(iz,j) = 0.0
        enddo
        enddo
        do iz=izs,ize
           stat(iz,1) = uwsb(iz)
           stat(iz,2) = vwsb(iz)
           stat(iz,3) = wwsb(iz)
           stat(iz,4) = tr_tau(iz)
           stat(iz,5) = triz(iz)
           stat(iz,6) = shrz(iz)
           stat(iz,7) = t_diss(iz)
        enddo
        m1 = js
        m2 = js + nscl
        m3 = js + 2*nscl
        do l=1,nscl
           do iz=izs,ize
              stat(iz,m1+l) = utsb(iz,l)
              stat(iz,m2+l) = vtsb(iz,l)
              stat(iz,m3+l) = wtsb(iz,l)
           enddo
        enddo
        call mpi_sum_z(stat(1,1),i_root,myid,nstat*nnz,1)
        do iz=1,nnz
           uwsb(iz)   = stat(iz,1)
           vwsb(iz)   = stat(iz,2)
           wwsb(iz)   = stat(iz,3)
           tr_tau(iz) = stat(iz,4)
           triz(iz)   = stat(iz,5)
           shrz(iz)   = stat(iz,6)
           t_diss(iz) = stat(iz,7)
        enddo
        do l=1,nscl
           do iz=1,nnz
              utsb(iz,l) = stat(iz,m1+l)
              vtsb(iz,l) = stat(iz,m2+l)
              wtsb(iz,l) = stat(iz,m3+l)
           enddo
        enddo
        do iz=1,nnz
           buyz(iz) = batag*wtsb(iz,1)
        enddo
c
c -------- end if block
c
      endif
c
c ------- save old rhs in field variables for RK-advancement
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = urhs(ix,iy,iz)
            v(ix,iy,iz) = vrhs(ix,iy,iz)
            w(ix,iy,iz) = wrhs(ix,iy,iz)
            e(ix,iy,iz) = erhs(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx

            t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
