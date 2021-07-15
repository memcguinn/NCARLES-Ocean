! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine writes history file with global parameters with      !
!         tsfcc to preserve digits.                                            !
! ============================================================================ !
!
SUBROUTINE write_his(iloc)
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
!
! --------------------------------------------------------------------------- !
!
  divgmax = 0.0
!
  DO iz=1,nnz
    divgmax = AMAX1(divgmax, divz(iz))
  END DO
!
  ziavg = zi
  holtop = hol
  wt_min = wtsb(iloc,1)
  wt_le  = wtle(iloc,1)
  krec = krec + 1
  mid = nnz/4
!
  WRITE (nhis1,6000) time,dt,utau,ziavg,amonin,holtop,         &
               (tsfcc(1)-t_ref),uwsfc,vwsfc,divgmax, wt_min,   &
               wt_le, ucfl, vcfl, wcfl, wtsfc(1),ups(mid),     &
               vps(mid),wps(mid),tps(mid,1),uwle(mid),         &
               uwsb(mid),uw_tot(mid),vwle(mid),vwsb(mid),      &
               vw_tot(mid),wtle(mid,1),wtsb(mid,1),            &
               wt_tot(mid,1),englez(mid),eavg(mid), wabs,      &
               u_10, v_10, cd_10, tau_x, tau_y, stokess,       &
               cpou10
!
  CALL write_prof(nhisp,krec,isize,c_s%wwsb)
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
6000 FORMAT (5e17.8)
!------------------------------------------------------------------------------!
!
END SUBROUTINE write_his
