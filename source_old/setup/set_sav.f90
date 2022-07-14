! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE set_sav(it,istart)
!
    USE inputs, ONLY: iti
    USE pars
    USE fields
    USE con_data
    USE con_stats
!
    DATA ionce /0/
    SAVE ionce
!
! --------------------------------------------------------------------------- !
!
! INCREMENT TIME IF NOT FIRST TIME THROUGH
  IF (it .NE. istart) THEN
    time=time+dt
  END IF
!
  it=it+1
  it_counter = it - iti
  dt    = dt_new
  dtg   = dt
!
  mnout = (MOD(it_counter,imean) .EQ. 0) .OR. (it .EQ. 1)
  mtape = (MOD(it_counter,itape) .EQ. 0)
  micut = (MOD(it_counter,itcut) .EQ. 0)
!
  IF (ihst .LT. 0) THEN
    mhis = .FALSE.
  ELSE
    mhis = (MOD(it_counter,ihst).eq.0 .and. it .ge. it_his)
  END IF
!
! BOOLEAN OPERATOR: VELOCITY FIELDS SAVED
  msave = .FALSE.
  IF (it_counter .GE. itstr .AND. mtape) THEN
    itn=itn+1
    msave = .TRUE.
    CALL get_output_filenames
  END IF
!
! BOOLEAN OPERATOR: HISTORY FIELDS SAVED
  IF ((ihst .GT. 0) .AND. (it_counter .GE. it_nxt)) THEN
    CALL open_his(it)
    it_nxt = it_nxt + itape
  END IF
!
RETURN
END SUBROUTINE set_sav
