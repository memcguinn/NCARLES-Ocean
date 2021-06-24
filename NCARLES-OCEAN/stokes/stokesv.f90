! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine computes the Stokes drift velocity using the         !
!         Donelan spectral wave parameters (Alves et al (2003), JPO).          !
! ============================================================================ !
!
SUBROUTINE stokesv
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
! ---------------------------------------------------------------------------- !
!
! USING DONELAN SHAPE (SET IN INIT), PEAK RECOMPUTES WITH VARYING U10
  CALL speed2stress(u_10,v_10,cd_10,tau_x,tau_y)
!
  speedval = SQRT(u_10**2 + v_10**2)
  f_p     = f2w*grav/speedval
  sigma_p = pi2*f_p
!
! SET PARAMETERS
  range_min = 0.1
  range_max = 5000.0
  DO iz=1,nnzp1
    z_pt = zz(iz)
    CALL s_int(range_min,range_max,value)
    stokes(iz) = value
!
! STOKES DRIFT OFF
    IF (flg_stokes .NE. 1) THEN
      stokes(iz) = 0.0
    END IF
  END DO
  IF (l_root) THEN
    WRITE (6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
  END IF
  dir_x = 1.0
  dir_y = 0.0
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
 6000 FORMAT(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))
 2212 FORMAT('#k ',/,'#m 4',/,'#lw 1.0',/,(2e15.6))
!------------------------------------------------------------------------------!
!
END SUBROUTINE stokesv
