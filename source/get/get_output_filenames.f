      subroutine get_output_filenames
c
c ----------- build file names for velocity, pressure, and constants
c
      use pars
      include 'mpif.h'
      character cgrid*10, num*5
c
c --------------- build character strings for file name
c
      cgrid = '.mp.'
      write(num,'(i5.5)') itn
      iblnk = index(path_sav,' ')
      call blnk(path_sav_v)
      call blnk(path_sav_p)
      call blnk(path_sav_c)
      path_sav_v = path_sav(1:iblnk-1)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:5)
      path_sav_p = path_sav(1:iblnk-1)//'/p'//
     +                 cgrid(1:4)//case(1:3)//num(1:5)
      path_sav_c = path_sav(1:iblnk-1)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:5)//'.con'
c
      return
      end
