! ============================================================================ !
! ABOUT:                                                                       !
!         Runge-Kutta coefficients are set in this subroutine. Currently, the  !
!         scheme is of Or(3) which is the default for NCARLES.                 !
! ============================================================================ !
!
SUBROUTINE setcon
!
    USE pars
    USE fftwk
!    USE fft
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    EXTERNAL get_zi
!
! --------------------------------------------------------------------------- !
!
! GET MACHINE TYPE, CAN ERAD DATADIR
! INITIALIZE FFT
  CALL rffti(nnx,trigx(1,1))
  CALL rffti(nny,trigx(1,2))
  CALL cffti(nny,trigc(1))
!
! HISTORY FILE INITIAL STEP
  it_nxt = it_his
!
! DEFINE MINIMUM VALUE OF E
  smal_e = 1.0e-12
!
! SET CONSTANTS IN EDDY VISCOSITY MODEL
  ck       = 0.1
  ceps     = 0.93
  csmag    = SQRT(ck*SQRT(ck/ceps))
  stab_c   = 0.76
!
! SET STABILITY CONSTANT
  stabmin = 1.0e-12
!
! MINIMUM DSL LENGTH CONSTANT
  almin_c = 0.0001
!
! CREATE MPI OPERATION TO FIND MAX LOCATION USING LOCAL GRADIENT
  CALL mpi_op_create(get_zi,.true.,ziloc,ierror)
!
! DEFINE COEFFICIENTS FOR OR(3) R-K TIME-STEP SCHEME - SPALART, MOSER, AND ROGERS (1990)
! NOTE: THIS IS A SIMPLER VERSION SINCE ALL TERMS ARE LUMPED IN NON-LINEAR TERMS.
!       CFL NUMBER IS FOR AN ENTIRE R-K STEP IN THIS CASE THREE STAGES.
  zetas(1) = 0.0
  zetas(2) = -17.0/60.0
  zetas(3) = -5.0/12.0
  gama(1)  = 8.0/15.0
  gama(2)  = 5.0/12.0
  gama(3)  = 3.0/4.0
  etas(1)  = -1.0
  etas(2)  = -1.0 + 8.0/15.0
  etas(3)  = -1.0 + 2.0/3.0                   ! Full step at new time
  etas(4)  =  0.0
!
  cfl = 0.40                                  ! CFL = MAX(u) * dt/dx
!
RETURN
END SUBROUTINE setcon
