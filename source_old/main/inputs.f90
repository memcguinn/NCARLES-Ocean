! ============================================================================ !
! ABOUT:                                                                       !
!         inputs contains the parameters most likely to be changed during the  !
!         initialization and runs of NCARLES-Ocean.                            !
! ============================================================================ !
!
!
MODULE inputs
!
! ---------------------------- RESTART PARAMETERS ---------------------------- !
INTEGER, PARAMETER :: iti = 0,              & ! Initial iteration
                      itn_restart = 0         ! Initial iteration save value (default 0)
CHARACTER*80       :: its_path = './data/u.mp.30L00000'
                                              ! Path to save file for initial iteration
! ---------------------------------------------------------------------------- !
REAL, PARAMETER ::   wtsfc_l1 = 0.0           ! Define surface cooling
!
! ----------------------------- STOKES PARAMETERS ---------------------------- !
REAL, PARAMETER ::  u10_stokes = 7.0,       & ! Wind strength to change Stokes velocity and La_t
                    u10_windspeed = 10.0      ! 10m wind speed
! ---------------------------------------------------------------------------- !
!
! equil_val = (/7.56903, 1.67006e03, 3.14655e02, 2.96936e02, 1.18909e02, 6.30928e-03, 9.60492/)
END MODULE inputs
