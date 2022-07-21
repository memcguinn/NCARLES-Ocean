      subroutine get_derv
c
c ------- get ux,uy,vx,vy at all z for this node
c         using parallel fft. can be improved (?)
c         by using exchange to send derivatives
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      iz_ss = izs-1
      iz_ee = ize+1
      if(iss .eq. 0) then
         iz_ss = izs
      endif
      if(ise .eq. numprocs-1) then
         iz_ee = ize
      endif
c
c ------- make sure <w> = 0
c
      do iz=izs-1,ize+1
         w_sum = 0.0
         do iy=iys,iye
         do ix=1,nnx
            w_sum = w_sum + w(ix,iy,iz)
         enddo
         enddo
         w_sum = w_sum*fnxy
         call mpi_sum_xy(w_sum,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz) = w(ix,iy,iz) - w_sum
         enddo
         enddo
      enddo
c
      do iz=izs-1,ize+1

         do iy=iys,iye
         do ix=1,nnx
            ux(ix,iy,iz) = u(ix,iy,iz)
            vx(ix,iy,iz) = v(ix,iy,iz)
            wx(ix,iy,iz) = w(ix,iy,iz)
            uy(ix,iy,iz) = u(ix,iy,iz)
            vy(ix,iy,iz) = v(ix,iy,iz)
            wy(ix,iy,iz) = w(ix,iy,iz)
         enddo
         enddo
c        endif
         call xderivp(ux(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
         call xderivp(vx(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
         call xderivp(wx(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
      enddo
c
c ---------- get y derivatives for (u,v,w)

      call yd_mpi(uy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(vy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(wy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
      return
      end
