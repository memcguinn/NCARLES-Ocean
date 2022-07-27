SUBROUTINE set_sav(it,istart)

  USE pars
  USE fields
  USE con_DATA
  USE con_stats

  DATA ionce /0/
  SAVE ionce

  IF(it .NE. istart) THEN
    ! INCREMENT TIME IF NOT FIRST TIME THROUGH
    time=time+dt
  ENDIF

  it=it+1
  it_counter = it - iti

  dt    = dt_new
  dtg   = dt
  mnout = (MOD(it_counter,imean).EQ.0) .OR. (it.EQ.1)
  mtape = (MOD(it_counter,itape).EQ.0)
  micut = (MOD(it_counter,itcut).EQ.0)
  mviz  = (MOD(it_counter,i_viz).EQ.0)

  IF(ihst .LT. 0) THEN
    mhis = .FALSE.
  ELSE
    mhis = (MOD(it_counter,ihst).EQ.0 .AND. it .GE. it_his)
  ENDIF

  ! DECIDE WHETHER VELOCITY FIELDS ARE SAVED
  msave = .FALSE.
  IF(it_counter .GE. itstr .AND. mtape) THEN
    itn=itn+1
    msave = .TRUE.
    CALL get_output_filenames
  ENDIF

  ! DECIDE WHETHER VIZ FIELDS ARE SAVED
  msave_v = .FALSE.
  IF(it_counter .GE. itstr .AND. mviz .AND. i_viz .GT. 0) THEN
    msave_v = .TRUE.
    IF(ionce .EQ. 0) THEN
      ionce = 1
      CALL open_viz
    ENDIF
  ENDIF

  ! DECIDE WHETHER HISTORY FILES ARE TO BE SAVED
  IF((ihst .GT. 0) .AND. (it_counter .GE. it_nxt)) THEN
    CALL open_his(it)
    it_nxt = it_nxt + itape
  ENDIF

  RETURN
END SUBROUTINE
