      subroutine set_paths
c
c ------------- set file path for RESTART, and
c               directories for saving, history, and viz files
c
c     path_res  --- the path and file name of the velocity restart file
c     path_sav  --- the path where the new 3d volumes are to be saved
c     path_his  --- the path where the new history files are to be saved
c     path_viz  --- the path where xy, xz, or yz planes of data will be stored
c     path_stuf --- the path where fun facts about the viz planes of
c                   data will be stored
c
      use pars
c
c ----------------- restart at step 05000
c
      path_res='./data/u.mp.30L00005'

c
      path_sav    = './data'
      path_his    = './data'
      path_viz_xy = './data'
      path_viz_xz = './data'
      path_viz_yz = './data'
      path_stuf   = './data'
c
c -------------- typical habu names
c
c     path_res    = '/scr/pps/data/les/2dmpi/xxx/u.mp.xxx004'
c     path_sav    = '/scr/pps/data/les/2dmpi/xxx'
c     path_his    = '/scr/pps/data/les/2dmpi/xxx'
c     path_viz_xy = '/scr/pps/data/les/2dmpi/xxx'
c     path_viz_xz = '/scr/pps/data/les/2dmpi/xxx'
c     path_viz_yz = '/scr/pps/data/les/2dmpi/xxx'
c     path_stuf   = '/scr/pps/data/les/2dmpi/xxx'
c
      return
      end
