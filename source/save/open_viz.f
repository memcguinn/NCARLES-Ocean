      subroutine open_viz
c
c ------------------- set visualization files,
c                     leaves files in scratch directory
c
      use pars
      include 'mpif.h'
      character iblks*16
c
c --------------- build character strings for file names
c                 with time step
c
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:7),'(i6.6)') iti
      iblks(8:8) = '_'
      write(iblks(9:14),'(i6.6)') itmax
c
      iloc = index(path_viz_xy,' ')
      path_viz_xy = path_viz_xy(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:14)//'.xy.data'
      iloc = index(path_viz_xz,' ')
      path_viz_xz = path_viz_xz(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:14)//'.xz.data'
      iloc = index(path_viz_yz,' ')
      path_viz_yz = path_viz_yz(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:14)//'.yz.data'
      path_stuf = path_stuf(1:iloc-1)
     +      //'/stuff.'//case(1:3)//iblks(1:14)//'.data'
      close(nviz_z)
      close(nviz_y)
      close(nviz_x)
      close(nviz_s)
c
c ----------- do not actually open the files here since
c             not all processors may have been picked and
c             its unknown how many variables are selected.
c             customized in sr. save_viz
c
      return
      end
