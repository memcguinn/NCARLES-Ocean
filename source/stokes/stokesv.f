      subroutine stokesv
c
c ----------- get stokes drift velocity using donelan spectrum
c             matched to wave data see alves et al jpo, 2003, vol. 33 .
c             compute all pdf variables in sr. pdf_ndot
c
c      include 'par.f'
c      include 'breaker_par.f'
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
c
c ----------- donelan shape (constants set in init)
c             recompute peak if we change u_10
c
	 call speed2stress(u_10,v_10,cd_10,tau_x,tau_y)

	 speedval = sqrt(u_10**2 + v_10**2)
         f_p     = f2w*grav/speedval
         sigma_p = pi2*f_p
c
c ----------- set parameters though not used here
c
         stokesw = 0.0
         stokesa = 0.0
         stokess = 0.0
c
         range_min = 0.1
         range_max = 5000.0
         do iz=1,nnzp1
            z_pt = zz(iz)
            call s_int(range_min,range_max,value)
            stokes(iz) = value
c
c ---------- for no stokes
c
            if(flg_stokes /= 1) then
                stokes(iz) = 0.0
            endif
         enddo
         if(l_root) then
            write(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
 6000       format(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))

         endif
c

         dir_x = 1.0
	 dir_y = 0.0


      return
 2212 format('#k ',/,
     +       '#m 4',/,
     +       '#lw 1.0',/,
     +       (2e15.6))
      end
