SUBROUTINE diurnal
!
  USE pars
  USE inputs
  USE con_data
!
  t_les = time
  day = t_les/86400
  pi     = 4.0*ATAN(1.0)            ! Pi = 3.1417
  pi2    = 2.0*pi                   !
  d_to_r = pi2/360.0

  IF (day .LE. 1) THEN
      t_di = t_les
  ELSE
      t_di = t_les - (FLOOR(day)*86400)
  END IF

  q0_cool = .0002
  qd_heat = .2
  t_heat  = 86400/3                 ! (s)
  cos_arr = (/(COS((2*pi*t_di/(3600*24))*d_to_r)), (COS((pi*t_heat)*d_to_r))/)
  qstar(1) = q0_cool + qd_heat * (MAXVAL(cos_arr) - COS(pi*t_heat))
  wtsfc(1) = qstar(1)
!
  RETURN
END SUBROUTINE diurnal
