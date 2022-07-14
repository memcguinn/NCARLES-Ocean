      subroutine gridd
c
c ----------- allocate space and pass arrays using modules
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
c ------------- establish association between pointers
c               and data structures, hurricane data also
c
      call fill_cc
      call fill_ch
      call fill_cs
c
c      write(6,5001) isize
c 5001 format(' size of stats array = ',i8)
c
c ---------------- debug for arrays
c
      big = -99.0e+300
c
c ---------------- setup grid
c
      nnx = nxg1
      nny = nyg1
      nnz = nzg1
c     izs = 1
c     ize = nnz
c
c
c ----------- make sure problem and cpu's match
c
      maxp   = numprocs-1
      ncpu_z = numprocs/ncpu_s
      if(mod(numprocs,ncpu_s) .ne. 0 .or.
     +   ncpu_z .gt. nnz) then
         go to 999
      endif
      if(l_root) write(6, 1100) ncpu_s, ncpu_z, numprocs,
     +                          maxp
      write(nprt,1100) ncpu_s, ncpu_z, numprocs, maxp
 1100 format(' Number of x-y slab cpus = ',i5,/,
     +       ' Number of z-level cpus  = ',i5,/,
     +       ' Total number of cpus    = ',i5,/,
     +       ' Max-p for index arrays  = ',i5)
c
c ---------------- allocate arrays for (i,j,k)-indexing on
c                  each processor (see set_range)
c
      allocate(ix_s(0:maxp), ix_e(0:maxp),
     +         jx_s(0:maxp), jx_e(0:maxp),
     +         kx_s(0:maxp), kx_e(0:maxp),
     +         mx_s(0:maxp), mx_e(0:maxp),
     +         iy_s(0:maxp), iy_e(0:maxp),
     +         jy_s(0:maxp), jy_e(0:maxp),
     +         is_s(0:maxp), is_e(0:maxp),
     +         iz_s(0:maxp), iz_e(0:maxp))
c
c ---------------- setup array sizes and variable dimensions
c
      nxy   = nnx*nny
      ncx   = nnx/2 + 1
      ncy   = nny/2 + 1
      nnxp1 = nnx + 1
      nnyp1 = nny + 1
      nnxp2 = nnx + 2
      nnyp2 = nny + 2
      nnzp1 = nnz + 1
      nnzm1 = nnz - 1
      ivis = ivis0
      fnxy  = 1.0/float(nnx*nny)
c
      write(nprt,7001) nnx,nny,nnz
 7001 format(' 7001 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      call set_range
c
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      num_y = iye + 1 - iys
c
c ------------- allocate solution arrays
c               account for nnxp2 for fields but not in rhs
c               and possible monotone for scalars
c
      allocate(u(nnxp2,iys:iye,izs-1:ize+1),
     +         v(nnxp2,iys:iye,izs-1:ize+1),
     +         w(nnxp2,iys:iye,izs-1:ize+1),
     +         t(nnxp2,iys:iye,nscl,izs-2:ize+2),
     +         e(nnxp2,iys:iye,izs-1:ize+1),
     +         r1(nnx,iys:iye,izs-1:ize+1),
     +         r2(nnx,iys:iye,izs-1:ize+1),
     +         r3(nnx,iys:iye,izs-1:ize+1),
     +         r4(nnx,iys:iye,nscl,izs-1:ize+1),
     +         r5(nnx,iys:iye,izs-1:ize+1))
c
c ------------- allocate space for boundary condition arrays
c               on top and bottom of domain
c
      allocate(ubc(nnx,iys:iye,2),
     +         vbc(nnx,iys:iye,2),
     +         wbc(nnx,iys:iye,2),
     +         tbc(nnx,iys:iye,nscl,2),
     +         ebc(nnx,iys:iye,2),
     +         pbc(nnx,iys:iye,2),
     +         pbc2(nnx,iys:iye,2))
c
c ------------ allocate space for wind and surface arrays
c
      allocate(wind(nnx,iys:iye),
     +         tau13m(nnx,iys:iye),
     +         tau23m(nnx,iys:iye),
     +         taut3m(nnx,iys:iye,nscl),
     +         t_grnd(nnx,iys:iye,nscl))
c
c ------------------- allocate space for derivative arrays
c
      allocate(ux(nnx,iys:iye,izs-1:ize+1),
     +         uy(nnx,iys:iye,izs-1:ize+1),
     +         vx(nnx,iys:iye,izs-1:ize+1),
     +         vy(nnx,iys:iye,izs-1:ize+1),
     +         wx(nnx,iys:iye,izs-1:ize+1),
     +         wy(nnx,iys:iye,izs-1:ize+1))
c
c ------------- allocate space for pressure, pressure bcs
c
      allocate(p(nnxp2,iys:iye,izs-1:ize+1),
     +         ptop(nnxp2,iys:iye,2))
c
c ------------- allocate space for viscosity and diffusivity
c
      allocate(vis_m(nnx,iys:iye,izs-1:ize+1),
     +         vis_s(nnx,iys:iye,izs-1:ize+1),
     +         vis_sv(nnx,iys:iye,izs-1:ize+1))
c
c ------------- allocate space for fft trig factors
c
      nq_trig = max(nnx,nny)
      allocate(trigx(2*nq_trig+15,2),
     +         trigc(4*nq_trig+15))
c
      return
  999 continue
c
      if(l_root) write(6,1000) numprocs, ncpu_s, mmz
      write(nprt,1000) numprocs, ncpu_s, nnz
 1000 format(' Gridd Trouble number of processors and grid',
     +          ' partitioning do not match!',/,
     +          ' Total num of cpus   = ',i5,
     +          ' Num cpu on x-y slab = ',i5,/,
     +          ' Num of z-levels     = ',i5)
      call mpi_finalize(ierr)
      end
