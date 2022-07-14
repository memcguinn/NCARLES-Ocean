      subroutine dear_vis(alk)
c
c --------- get deardorff stability corrected length
c           scales, possible new vis model
c
      use pars
      use fields
      use con_data
      use con_stats
      real alk(nnx,iys:iye,izs-1:ize+1)
c
      do iz=izs-1,ize+1
c
      izp1 = iz + 1
      dslk  = dsl_z(iz)
      if(iz .gt. 0) dslk  = amin1(dsl_z(iz),vk*abs(z(iz))/csmag)
      almin = almin_c*dsl_z(iz)
      if(iz .eq. 0 .or. iz .eq. nnz+1) then
         dfack = 1.0
      else
         dfack = dfac(iz)
      endif
c
      if(ivis .eq. 1 .and. iz .le. nmatch) then
c
c --------------- no stability corrected length scales
c
         do j=iys,iye
         do i=1,nnx
            alk(i,j,iz) = dslk
         end do
         end do
      else
         do j=iys,iye
         do i=1,nnx
            alk(i,j,iz) = dslk
            stab = batag*(t(i,j,1,izp1) - t(i,j,1,iz))*dzu_i(izp1)
            if(stab.gt.stabmin) then
              als = stab_c*sqrt(e(i,j,iz)/stab)
              alk(i,j,iz) = amin1(dslk,als)
            endif
            alk(i,j,iz)  = amax1(almin,alk(i,j,iz))
         enddo
         enddo
      endif
c
      do j=iys,iye
      do i=1,nnx
         vis_m(i,j,iz)  = ck*alk(i,j,iz)*sqrt(e(i,j,iz))*dfack
         vis_s(i,j,iz)  = (1.+2.*alk(i,j,iz)/dslk)*vis_m(i,j,iz)
         vis_sv(i,j,iz) = vis_s(i,j,iz)
      enddo
      enddo
c
c
c -------------- special case for iz = 1
c
      if(iz.eq.1 .and. ibcl .eq. 0) then
         do iy=iys,iye
         do ix=1,nnx
            vis_m(ix,iy,iz-1)  = vis_m(ix,iy,iz)
            vis_s(ix,iy,iz-1)  = vis_s(ix,iy,iz)
            vis_sv(ix,iy,iz-1) = vis_sv(ix,iy,iz)
         enddo
         enddo
      endif
c
c -------------- end z loop
c
      enddo
c
      return
      end
