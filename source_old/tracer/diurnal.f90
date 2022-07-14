SUBROUTINE diurnal
!
  USE pars
  USE inputs
  USE con_data
!
  t_les = time
  day   = t_les/86400
  pi    = 4.0*ATAN(1.0)            ! Pi = 3.1417

  t_di = day - FLOOR(day)

  q0_cool = 2e-5
  qd_heat = -2e-5
  t_heat  = 1              ! (d)
  cos_arr = (/(COS(2*pi*t_di)), (COS(pi*t_heat))/)
  qstar(1) = q0_cool + qd_heat * (MAXVAL(cos_arr) - COS(pi*t_heat))
  wtsfc(1) = qstar(1)
!
  RETURN
END SUBROUTINE diurnal
