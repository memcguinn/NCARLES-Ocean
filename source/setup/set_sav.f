      subroutine set_sav(it,istart)
c
      use pars
      use fields
      use con_data
      use con_stats
c
      data ionce /0/
      save ionce
c
      if(it .ne. istart) then
c
c ------------------- increment time if not first time through
c
         time=time+dt
      endif
c
      it=it+1
      it_counter = it - iti
c
      dt    = dt_new
      dtg   = dt
      mnout = (mod(it_counter,imean).eq.0).or. (it.eq.1)
      mtape = (mod(it_counter,itape).eq.0)
      micut = (mod(it_counter,itcut).eq.0)
      mviz  = (mod(it_counter,i_viz).eq.0)
      if(ihst .lt. 0) then
         mhis = .false.
      else
         mhis = (mod(it_counter,ihst).eq.0 .and. it .ge. it_his)
      endif
c
c ---------- decide whether velocity fields are saved
c
      msave = .false.
      if(it_counter .ge. itstr .and. mtape) then
         itn=itn+1
         msave = .true.
         call get_output_filenames
      endif
c
c ---------- decide whether viz fields are saved
c
      msave_v = .false.
      if(it_counter .ge. itstr .and. mviz .and. i_viz .gt. 0) then
         msave_v = .true.
         if(ionce .eq. 0) then
            ionce = 1
            call open_viz
         endif
      endif
c
c --------- decide whether history files are to be saved
c
      if((ihst .gt. 0) .and. (it_counter .ge. it_nxt)) then
         call open_his(it)
         it_nxt = it_nxt + itape
      endif
c
      return
      end
