SUBROUTINE stokes_ideal
! GET STOKES DRIFT VELOCITY FOR ASSUMED WAVELENGTH STOKESW AND WAVE AMPLITUDE
! STOKESA.
! CHANGED SIGN FOR Z

  USE pars
  USE con_data
  USE con_stats
  INCLUDE 'mpif.h'

  ! COMPUTE STOKES VELOCITY FOR OCEAN PBLS
  stokesw = pi2/76.5
  ak      = 0.00
  stokesa = ak/stokesw
  sigma = SQRT(ABS(grav)*stokesw)
  stokess = sigma*stokesw*stokesa**2

  DO iz=1,nnzp1
    stokes(iz) = stokess*EXP(2.0*stokesw*zz(iz))
  ENDDO

  IF(l_root) THEN
    WRITE(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
    6000 FORMAT(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))
  ENDIF

  RETURN
END SUBROUTINE
