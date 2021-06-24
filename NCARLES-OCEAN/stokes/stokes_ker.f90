! ============================================================================ !
! ABOUT:                                                                       !
!         This function evaluates the kernel of the Stokes integral for the    !
!         Donelan spectral shape of the wave. This is then converted into the  !
!         "sigma" space using factors of 2pi.                                  !
! ============================================================================ !
!
FUNCTION stokes_ker(sigma)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
! ---------------------------------------------------------------------------- !
!
  wave_spec =  (ann*grav_w*grav_w/(sigma_p*sigma**4))*       &
               EXP(-bnn*(sigma_p/sigma)**4)
  stokes_ker = 2.0*(wave_spec*sigma**3)*                     &
               EXP(2.0*sigma*sigma*z_pt/grav_w)/grav_w
!
RETURN
END FUNCTION stokes_ker
