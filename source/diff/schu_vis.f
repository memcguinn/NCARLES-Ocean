      subroutine schu_vis(alk)
c
c --------- get schumann stability corrected length
c           scales, possible new vis model.
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
c --------------- no stability corrected length scales
c
      do j=iys,iye
      do i=1,nnx
         alk(i,j,iz) = dslk
      end do
      end do
c
c -------------- stability correction for vertical scalar flux
c
      do j=iys,iye
      do i=1,nnx
         vis_m(i,j,iz)  = ck*dslk*sqrt(e(i,j,iz))*dfack
         vis_s(i,j,iz)  = 3.0*vis_m(i,j,iz)
c
         stab = amax1(batag*(t(i,j,1,izp1) - t(i,j,1,iz))*dzu_i(izp1),
     +                0.0)
         vis_sv(i,j,iz) = vis_s(i,j,iz)*(e(i,j,iz)/
     +               (e(i,j,iz) + 0.3*stab*dslk**2))
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
