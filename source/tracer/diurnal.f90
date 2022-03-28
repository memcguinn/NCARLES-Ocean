SUBROUTINE diurnal
!
  USE pars
  USE inputs
  USE con_data
!
  t_di = time
  pi     = 4.0*ATAN(1.0)            ! Pi = 3.1417
  q0_cool = 0.200
  qd_heat = -1.834
  t_heat  = 86400/3                 ! (s)
  cos_arr = (/(COS(2*pi*t_di/(3600*24))), (COS(pi*t_heat))/)
  qstar(1) = q0_cool + qd_heat * (MAXVAL(cos_arr) - COS(pi*t_heat))
  wtsfc(1) = qstar(1)
!
  RETURN
END SUBROUTINE diurnal
