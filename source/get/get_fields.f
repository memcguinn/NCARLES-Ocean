      subroutine get_fields
c
c ----------- special routine to read just 3d fields
c             as an initial guess, easy to customize
c             if missing data, etc..
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c#if defined(SWAP)
c      use module_byteswap
c#endif
      include 'mpif.h'
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
      real, allocatable, dimension(:,:,:) :: temp
      real, allocatable, dimension(:) ::     ui_u, ui_v
      real :: max_TKE
      parameter(max_TKE=5E-4)
      character*80 path_ran, path_inertial
      logical there
c
      allocate(temp(nvar,nnx,iys:iye))
      allocate(ui_u(nnz),ui_v(nnz))
c
c ---------- input file to read from
c
      path_ran = '/lustre/janus_scratch/peha5625/restart/bluefire_restart'
c      path_inertial = '/les/hurr/big/inertial.xy'
c
c --------------------- get restart file from local directory
c                       reuse unit number
c
      close(nvel)
c
      inquire(file=path_ran,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_ran
      else
         if(l_root) write(6,6005) path_ran
         stop
      endif
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_ran,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ------------ read 3d fields, make rhs*8
c
      nsize  = int(nvar,k8)*nnx*nny
      nsize2 = int(nvar,k8)*nnx*(iys-1)
      n_read = nvar*nnx*(iye+1-iys)
c
      do k=izs,ize
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_read_at_all(nvel,offset,temp,n_read,
     +                              mpi_real8,status,ierr)
         if (ierr /= 0) goto 9992
c#if defined(SWAP)
c         call byteswap(temp)
c#endif
         do j=iys,iye
         do i=1,nnx
            u(i,j,k) = temp(1,i,j)
            v(i,j,k) = temp(2,i,j)
            w(i,j,k) = temp(3,i,j)
            e(i,j,k) = min(max_TKE,temp(nvar,i,j))

         enddo
         enddo
         do is = 1,nscl
            do j = iys,iye
            do i = 1,nnx
               t(i,j,is,k) = temp(3+is,i,j)+273.15
            enddo
            enddo
         enddo
c
      enddo
c
c ---- close file
c
      call mpi_file_close(nvel, ierr)
c
      deallocate(temp)
c
      do iz=izs,ize
c
         ug(iz) = 0.0
         vg(iz) = 0.0
c
c ---------------- initial guess for pressure
c
         do iy=iys,iye
         do ix=1,nnx
            p(ix,iy,iz) = 0.0
         enddo
         enddo
      enddo
c
c ------------ fix currents for inertial oscillation
c
c     if(l_root) then
c        close(99)
c        open(99,file=path_inertial,form='formatted')
c        read(99,*) (ui_u(i),ui_v(i),i=1,nnz)
c        close(99)
c     endif
c     call mpi_bcast(ui_u,nnz,mpi_real8,
c    +              i_root,mpi_comm_world,ierr)
c     call mpi_bcast(ui_v,nnz,mpi_real8,
c    +              i_root,mpi_comm_world,ierr)
c
c     write(nprt,1910) (i, ui_u(i), ui_v(i),i=1,nnz)
c1910 format(' Inertial oscillation correction',/,
c    +       ' i ',5x,' ui u ',5x,' ui v ',/,
c    +      (i5,2e15.6))
c
c     do k=izs,ize
c     do j=iys,iye
c     do i=1,nnx
c        u(i,j,k) = u(i,j,k) + ui_u(k)
c        v(i,j,k) = v(i,j,k) + ui_v(k)
c     enddo
c     enddo
c     enddo
c
      deallocate(ui_u,ui_v)
c
      return
c ---------------------------- process errors
  100 continue
      write(6,9000) path_ran, nvel
      call mpi_finalize(ierr)
      stop
c
 9992 continue
      write(6,6100) nvel,iz
      call mpi_finalize(ierr)
      stop
c ---------------------
 6001 format(' SR. GET_FIELDS: FILE READ FOR INITIALIZATION = ',a80)
 6005 format(' 6005, SR. GET_FIELDS: cannot find restart file = ',a80)
 6100 format(' SR. GET_FIELDS: file read error on unit number = ',i2,/,
     +       '               at iz = ',i4)
 9000 format(' 9000, SR. GET_FIELDS: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
      end
