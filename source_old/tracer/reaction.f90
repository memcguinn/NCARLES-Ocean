! ============================================================================ !
! ABOUT:                                                                       !
!         Stabilized Runge-Kutta Chebyshev solver of Or(2). Tolerances,        !
!         uround, etc. set here.
! ============================================================================ !
!
MODULE reaction
!
  USE fields, ONLY: t
  USE con_data, ONLY: time,dt
  USE pars, ONLY: nscl
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: flg_debug = 0       ! Write debug file (0/1)
!
CONTAINS
!
! --------------------------------------------------------------------------- !
!
! CALCULATE THE SCALAR REACTION SOURCE TERM FOR GIVEN SCALAR AND POINT (rhs_scl)
  FUNCTION react_src(ix,iy,iscl,iz)
!
    REAL, DIMENSION(nscl-1) :: react_src
    INTEGER, INTENT(in) :: ix,            & ! Location index of interest in x
                           iy,            & ! Location index of interest in y
                           iscl,          & ! Scalar index of interest
                           iz               ! Location index of interest in z
!
    INTEGER :: zzi, i
    REAL, DIMENSION(0:nscl-2) :: co2, co2tmp
    REAL :: dtst, temper, K1star, K2star, Kwater, Kboron, Rgas, salt
    REAL :: h, t_rkc, t_next, t_end, task
    INTEGER :: steps
!
    t_rkc  = time
    t_end  = time + dt*0.5
    task   = 1                              ! Integrate to end (0/1), single if else
!
    co2(0) = t(ix,iy,2,iz)                  ! Carbon Dioxide, [CO2]
    co2(1) = t(ix,iy,3,iz)                  ! Bicarbonate, [HCO3-]
    co2(2) = t(ix,iy,4,iz)                  ! Carbonate, [CO32-]
    co2(3) = t(ix,iy,5,iz)                  ! Boric Acid, [B(OH)3]
    co2(4) = t(ix,iy,6,iz)                  ! Tetrahydroxyborate, [B(OH)4-]
    co2(5) = t(ix,iy,7,iz)                  ! Hydrogen Ion, [H+]
    co2(6) = t(ix,iy,8,iz)                  ! Hydroxide, [OH-]
    temper = t(ix,iy,1,iz)
!
    co2tmp = intDriver(t_rkc, t_end, co2, temper)
!
    DO i = 0,nscl-2
       co2(i) = co2tmp(i)
    ENDDO
!
    DO i = 0,nscl-2
       react_src(i+1) = co2(i)
    ENDDO
!
  END FUNCTION react_src
!
!
! --------------------------------------------------------------------------- !
!
  FUNCTION intDriver(t_rkc, t_end, yGlobal, temper)
!
    REAL, INTENT(in) :: t_rkc                 ! Chemistry time
    REAL, INTENT(in) :: t_end, temper         ! t_end = time + dt/2
    REAL, INTENT(in), DIMENSION(0:nscl-2) :: yGlobal
    REAL, DIMENSION(0:nscl-2) :: intDriver
    REAL, DIMENSION(0:nscl-2) :: yLocal, yLocal2
    REAL, DIMENSION(0:4+nscl-1) :: workLocal
!
    INTEGER i
!
    workLocal(:) = 0.0
!
    DO i = 0,nscl-2
       yLocal(i)    = yGlobal(i)
    ENDDO
!
    yLocal2 = rkc_driver(t_rkc, t_end, workLocal, yLocal, temper)
!
    DO i = 0,nscl-2
       yLocal(i)      = yLocal2(i)
       intDriver(i)   = yLocal2(i)
    ENDDO
!
  END FUNCTION intDriver
!
!
! --------------------------------------------------------------------------- !
!
! DRIVER FUNCTION FOR RUNGE-KUTTA CHEBYSHEV INTEGRATOR
  FUNCTION rkc_driver(t_rkc2, t_end, work, yLocal, temper)
!
    REAL, INTENT(in) :: t_rkc2                ! Next/new time
    REAL, INTENT(in) :: t_end,              & !
                        temper
    REAL, INTENT(inout), DIMENSION(0:nscl-2) :: &
                        yLocal                ! Dependent variable array, integrated values replace initial conditions
    REAL, DIMENSION(0:nscl-2) :: rkc_driver   !
    REAL, INTENT(inout), DIMENSION(0:4+nscl-1) :: &
                        work                  ! Real work array, size 3
    REAL, DIMENSION(0:nscl-2) :: y_n, F_n, temp_arr, temp_arr2
!
    INTEGER nstep, m_max, i, m
!
    REAL abs_tol, rel_tol, uround, hmax, hmin, err, est
    REAL fac, temp1, temp2, t_rkc
!
!
    t_rkc = t_rkc2
    nstep   = 0
    rel_tol = 1.0e-7
    abs_tol = 1.0e-10
    uround  = 2.22e-16
    m_max   = NINT(SQRT(rel_tol / (10.0 * uround)))
    hmax    = ABS(t_end - t_rkc)
    hmin    = 10.0 * uround * MAX(ABS(t_rkc), hmax)
!
    IF(m_max < 2)THEN
       m_max = 2
    ENDIF
!
    DO i = 0,nscl-2
       y_n(i) = yLocal(i)
    ENDDO
!
!   CALCULATE F_N FOR INITIAL Y
    F_n = dydt(t_rkc, y_n, temper)
!
!   LOAD INITIAL ESTIMATE FOR EIGENVECTOR
    IF(work(2) < uround) THEN
       DO i = 0,nscl-2
          work(4+i) = F_n(i)
       ENDDO
    ENDIF
!
!   ESTIMATE JACOBIAN SPECTRAL RADIUS, USE TIME STEP FROM WORK(3)
!   SPEC_RAD = WORK(4) IFF 25 STEPS PASS
    DO WHILE (t_rkc < t_end)
!
       temp_arr(:)  = 0.0
       temp_arr2(:) = 0.0
       err          = 0.0
!
       IF(MOD(nstep,25) == 0)THEN
          work(3) = rkc_spec_rad(t_rkc, hmax, y_n, F_n, work(4), temp_arr2, temper)
       ENDIF
!
!      FIRST STEP - ESTIMATE STEP SIZE
       IF(work(2) < uround)THEN
          work(2) = hmax
!
          IF((work(3) * work(2)) > 1.0)THEN
             work(2) = 1.0/work(3)
          ENDIF
!
          work(2) = MAX(work(2), hmin)
!
          DO i = 0,nscl-2
             temp_arr(i) = y_n(i) + (work(2) * F_n(i))
          ENDDO
!
          temp_arr2 = dydt(t_rkc + work(2), temp_arr, temper)
          err = 0.0
!
          DO i = 0,nscl-2
             est = (temp_arr2(i) - F_n(i)) / (abs_tol + rel_tol * ABS(y_n(i)))
             err = err + est*est
          ENDDO
!
          err = work(2) * SQRT(err/REAL(nscl-2))
          IF((0.1 * work(2)) < (hmax * SQRT(err)))THEN
             work(2) = MAX((0.1 * work(2)) / SQRT(err), hmin)
          ELSE
             work(2) = hmax
          ENDIF
       ENDIF
!
!      CHECK: LAST STEP
       IF((1.1 * work(2)) .GE. ABS(t_end - t_rkc))THEN
          work(2) = ABS(t_end - t_rkc)
       ENDIF
!
!      CALCULATE NUMBER OF STEPS
       m = 1 + NINT(SQRT(1.54 * work(2) * work(3) + 1.0))
!
       IF(m > m_max)THEN
          m = m_max
          work(2) = REAL((m*m - 1) / (1.54*work(3)))
       ENDIF
!
       hmin = 10.0 * uround * MAX(ABS(t_rkc), ABS(t_rkc + work(2)))
!
!      PERFORM TENTATIVE TIME STEP
       yLocal = rkc_step(t_rkc, work(2), y_n, F_n, m, temper)
!
!      CALCULATE F_NP1 WITH TENTATIVE Y_NP1
       temp_arr = dydt(t_rkc + work(2), yLocal, temper)
!
!      ESTIMATE ERROR
       err = 0.0
       DO i = 0,nscl-2
          est = 0.0
          est = 0.8 * (y_n(i) - yLocal(i)) + 0.4 * work(2) * (F_n(i) + temp_arr(i))
          est = est / (abs_tol + rel_tol * MAX(ABS(yLocal(i)), ABS(y_n(i))))
          err = err + est*est
       ENDDO
!
       err = SQRT(err / 7.0)
       IF (err > 1.0) THEN
!
!         REJECT STEP IF ERROR TOO LARGE AND SELECT SMALLER STEP SIZE
          work(2) = 0.8 * work(2) / (err**(1.0/3.0))
!
!         REEVALUATE SPECTRAL RADIUS
          work(3) = rkc_spec_rad(t_rkc, hmax, y_n, F_n, work(4), temp_arr2, temper)
!
       ELSE
!
!         ACCEPT STEP IF ERROR SMALL
          t_rkc = t_rkc + work(2)
          nstep = nstep + 1
!
          fac   = 10.0
          temp1 = 0.0
          temp2 = 0.0
!
          IF(work(1) < uround)THEN
             temp2 = err**(1.0/3.0)
             IF(0.8 < (fac * temp2))THEN
                fac = 0.8 /  temp2
             ENDIF
!
          ELSE
!
             temp1 = 0.8 * work(2) * (work(0)**(1.0/3.0))
             temp2 = work(1) * (err**(2.0/3.0))
!
             IF(temp1 < (fac * temp2))THEN
                fac = temp1 / temp2
             ENDIF
          ENDIF
!
!         SET PREVIOUS VALUES FOR CURRENT TIME STEP
          work(0) = err
          work(1) = work(2)
!
          DO i = 0,nscl-2
             y_n(i) = yLocal(i)
             F_n(i) = temp_arr(i)
          ENDDO
!
!         STORE NEXT TIME STEP
          work(2) = work(2) * MAX(0.1, fac)
          work(2) = MAX(hmin, min(hmax, work(2)))
!
       ENDIF
    ENDDO
!
    DO i = 0,nscl-2
       rkc_driver(i) = yLocal(i)
    ENDDO
!
  END FUNCTION rkc_driver
!
!
! --------------------------------------------------------------------------- !
!
! ESTIMATE SPECTRAL RADIUS
  REAL FUNCTION rkc_spec_rad(t_rkc, hmax, yLocal, F, v, Fv, temper)
!
    REAL, INTENT(in) :: t_rkc
    REAL, INTENT(in) :: hmax,                 & ! Max time step size
                        temper
    REAL, INTENT(inout), DIMENSION(0:nscl-2) :: &
                        v,                    & !
                        Fv,                   & !
                        F                       ! Derivative evaluated at current state
    REAL, INTENT(in), DIMENSION(0:nscl-2) ::  &
                        yLocal                  ! Array of dependent variables
!
    INTEGER itmax, i, iter, ind
!
    REAL uround, small, nrm1, nrm2, dynrm, sigma
!
!
    uround  = 2.22e-16
    itmax   = 50
    small   = 1.0 / hmax
    nrm1    = 0.0
    nrm2    = 0.0
    sigma   = 0.0
!
    DO i = 0,nscl-2
       nrm1 = nrm1 + yLocal(i) * yLocal(i)
       nrm2 = nrm2 + v(i) * v(i)
    ENDDO
!
    nrm1 = SQRT(nrm1)
    nrm2 = SQRT(nrm2)
!
    IF((nrm1 .NE. 0.0) .AND. (nrm2 .NE. 0.0))THEN
       dynrm = nrm1 * SQRT(uround)
       DO i = 0,nscl-2
          v(i) = yLocal(i) + v(i) * (dynrm / nrm2)
       ENDDO
!
    ELSEIF(nrm1 .NE. 0.0)THEN
!
       dynrm = nrm1 * SQRT(uround)
       DO i = 0,nscl-2
          v(i) = yLocal(i) * (1.0 + SQRT(uround))
       ENDDO
!
    ELSEIF(nrm2 .NE. 0.0)THEN
!
       dynrm = uround
       DO i = 0,nscl-2
          v(i) = v(i) * (dynrm / nrm2)
       ENDDO
!
    ELSE
!
       dynrm = uround
       DO i = 0,nscl-2
          v(i) = uround
       ENDDO
    ENDIF
!
!   ITERATE USING NONLINEAR POWER METHOD
    sigma = 0.0
    DO iter = 1,itmax
       Fv = dydt(t_rkc, v, temper)
!
       nrm1 = 0.0
!
       DO i = 0,nscl-2
          nrm1 = nrm1 + ((Fv(i) - F(i)) * (Fv(i) - F(i)))
       ENDDO
!
       nrm1  = SQRT(nrm1)
       nrm2  = sigma
       sigma = nrm1 / dynrm
!
       IF((iter .GE. 2) .AND. (ABS(sigma - nrm2) .le. (MAX(sigma, small) * 0.01)))THEN
          DO i = 0,nscl-2
             v(i) = v(i) - yLocal(i)
          ENDDO
          rkc_spec_rad = 1.2 * sigma
       ENDIF
!
       IF(nrm1 .NE. 0.0)THEN
          DO i = 0,nscl-2
             v(i) = yLocal(i) + ((Fv(i) - F(i)) * (dynrm / nrm1))
          ENDDO
!
       ELSE
!
          ind = MOD(iter, INT(nscl-1))
          v(ind) = yLocal(ind) - (v(ind) - yLocal(ind))
       ENDIF
    ENDDO
!
    rkc_spec_rad = 1.2 * sigma
!
  END FUNCTION rkc_spec_rad
!
!
! --------------------------------------------------------------------------- !
!
! TAKE A SINGLE RKC INTEGRATION STEP
  FUNCTION rkc_step(t_rkc, h, y_0, F_0, s, temper)
!
    REAL, INTENT(in) :: t_rkc                   ! Start time
    REAL, INTENT(in) :: h,                    & ! Time-step size
                        temper                  !
    REAL, INTENT(inout), DIMENSION(0:nscl-2) :: &
                        y_0,                  & ! Initial conditions
                        F_0                     ! Derivative function at IC
!
    INTEGER, INTENT(in) :: s                    ! Number of steps
!
    REAL, DIMENSION(0:nscl-2) :: rkc_step       ! Integrated variables
    REAL, DIMENSION(0:nscl-2) :: y_j
    REAL w0, temp1, temp2, arg, w1, b_jm1, b_jm2, mu_t
    REAL c_jm2, c_jm1, zjm1, zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2
    REAL zj, dzj, d2zj, b_j, gamma_t, nu, mu, c_j
    REAL, DIMENSION(0:nscl-2) :: y_jm1, y_jm2
!
    INTEGER i, j
!
!
    w0    = 1.0 + 2.0 / (13.0 * REAL(s * s))
    temp1 = (w0 * w0) - 1.0
    temp2 = SQRT(temp1)
    arg   = REAL(s) * LOG(w0 + temp2)
    w1    = SINH(arg) * temp1 / (COSH(arg) * REAL(s) * temp2 - w0 * SINH(arg))
!
    b_jm1 = 1.0 / (4.0 * (w0 * w0))
    b_jm2 = b_jm1
!
!   CALCULATE Y_1
    mu_t = w1 * b_jm1
    DO i = 0,nscl-2
       y_jm2(i) = y_0(i)
       y_jm1(i) = y_0(i) + (mu_t * h * F_0(i))
    ENDDO
!
    c_jm2 = 0.0
    c_jm1 = mu_t
    zjm1 = w0
    zjm2 = 1.0
    dzjm1 = 1.0
    dzjm2 = 0.0
    d2zjm1 = 0.0
    d2zjm2 = 0.0
!
    DO j = 2,s
!
       zj = 2.0 * w0 * zjm1 - zjm2
       dzj = 2.0 * w0 * dzjm1 - dzjm2 + 2.0 * zjm1
       d2zj = 2.0 * w0 * d2zjm1 - d2zjm2 + 4.0 * dzjm1
       b_j = d2zj / (dzj * dzj)
       gamma_t = 1.0 - (zjm1 * b_jm1)
!
       nu = -b_j / b_jm2
       mu = 2.0 * b_j * w0 / b_jm1
       mu_t = mu * w1 / w0
!
!      CALCULATE DERIVATIVE, USE Y ARRAY FOR TEMPORARY STORAGE
       y_j = dydt(t_rkc + (h * c_jm1), y_jm1, temper)
!
       DO i = 0,nscl-2
          y_j(i) = (1.0 - mu - nu) * y_0(i) + (mu * y_jm1(i)) + (nu * y_jm2(i)) &
               + h * mu_t * (y_j(i) - (gamma_t * F_0(i)))
       ENDDO
!
       c_j = (mu * c_jm1) + (nu * c_jm2) + mu_t * (1.0 - gamma_t)
!
       IF(j < s)THEN
          DO i = 0,nscl-2
             y_jm2(i) = y_jm1(i)
             y_jm1(i) = y_j(i)
          ENDDO
       ENDIF
!
       c_jm2  = c_jm1
       c_jm1  = c_j
       b_jm2  = b_jm1
       b_jm1  = b_j
       zjm2   = zjm1
       zjm1   = zj
       dzjm2  = dzjm1
       dzjm1  = dzj
       d2zjm2 = d2zjm1
       d2zjm1 = d2zj
    ENDDO
!
    DO i = 0,nscl-2
       rkc_step(i) = y_j(i)
    ENDDO
!
  END FUNCTION rkc_step
!
!
! --------------------------------------------------------------------------- !
!
  FUNCTION dydt(t_rkc, y, temper)
    REAL, INTENT(in),  DIMENSION(0:nscl-2) :: y
    REAL, DIMENSION(0:nscl-2) :: dydt, dy
    REAL, DIMENSION(nscl-1) :: c
    REAL, INTENT(in) :: t_rkc, temper
    REAL K1s, K2s, Kw, Kb, Rgas, salt
!
    INTEGER i
!
    LOGICAL reduced
!
    REAL  forward_coeff1, forward_coeff2, forward_coeff3, forward_coeff4,   &
          forward_coeff5, forward_coeff6, forward_coeff7
    REAL  reverse_coeff1, reverse_coeff2, reverse_coeff3, reverse_coeff4,   &
          reverse_coeff5, reverse_coeff6, reverse_coeff7
!
    reduced = .false.
    salt   = 35.0
    Rgas = 0.0083143
!
    DO i = 0,nscl-2
       c(i+1) = y(i)
    ENDDO
!
    K1s = EXP(-2307.1266/temper + 2.83655 - 1.5529413*LOG(temper) +         &
         (-4.0484/temper - 0.20760841)*(salt**0.5) + 0.08468345*salt -      &
         0.00654208*(salt**1.5) + LOG(1.0-0.001005*salt))*(1.0e6)
    K2s = EXP(-3351.6106/temper - 9.226508 - 0.2005743*LOG(temper) +        &
         (-23.9722/temper - 0.106901773)*(salt**0.5) + 0.1130822*salt -     &
         0.00846934*(salt**1.5) + LOG(1.0-0.001005*salt))*(1.0e6)
!
!   SOLVE KW USING DOE (1994)
    Kw = EXP(148.96502 - 13847.26/temper - 23.65218*LOG(temper) +           &
         (118.67/temper - 5.977 + 1.0495*LOG(temper))*(salt**0.5) -         &
         0.01615*salt)*(1.0e6)
!
!   SOLVE KB USING DICKSON (1990)
    Kb = EXP((-8966.9 - 2890.53*(salt**0.5) - 77.942*salt +                 &
         1.728*(salt**1.5) - 0.0996*(salt**2))/temper                       &
         + 148.0248 + 137.1942*(salt**0.5) + 1.62142*salt -                 &
         (24.4344 + 25.085*(salt**0.5) + 0.2474*salt)*LOG(temper) +         &
         0.053105*(salt**0.5)*temper)*(1.0e6)
!
    forward_coeff1 = EXP(1246.98-6.19*(10.0**4)/temper - 183.0*LOG(temper))
    forward_coeff2 = (4.7e7)*EXP(-23.3/(Rgas*temper))/(1.0e6)
    forward_coeff3 = (5.0e10)/(1.0e6)
    forward_coeff4 = (6.0e9)/(1.0e6)
    forward_coeff5 = (1.4e-3)*(1.0e6)
    forward_coeff6 = (4.58e10)*EXP(-(20.8/(Rgas*temper)))/(1.0e6)
    forward_coeff7 = (3.05e10)*EXP(-(20.8/(Rgas*temper)))/(1.0e6)
    reverse_coeff1 = forward_coeff1/K1s
    reverse_coeff2 = (Kw*forward_coeff2/K1s)*(1.0e6)
    reverse_coeff3 = forward_coeff3*K2s
    reverse_coeff4 = (forward_coeff4*Kw/K2s)*(1.0e6)
    reverse_coeff5 = (forward_coeff5/Kw)/(1.0e6)
    reverse_coeff6 = (forward_coeff6*Kw/Kb)*(1.0e6)
    reverse_coeff7 = forward_coeff7*K2s/Kb
!
!   REDUCED MODEL, QSS SPECIES H+
    c(6) = (forward_coeff1*c(1) + reverse_coeff3*c(2) + forward_coeff5)/    &
            (reverse_coeff1*c(2) + forward_coeff3*c(3) + reverse_coeff5*c(7))
!
    dy(0) = reverse_coeff1*c(2)*c(6)+reverse_coeff2*c(2)-forward_coeff1*    &
            c(1)-forward_coeff2*c(1)*c(7)
    dy(1) = forward_coeff1*c(1)+forward_coeff2*c(1)*c(7)-reverse_coeff1*    &
            c(2)*c(6)-reverse_coeff2*c(2)+forward_coeff3*c(3)*c(6)-         &
            reverse_coeff3*c(2)-forward_coeff4*c(2)*c(7)+reverse_coeff4*    &
            c(3)+forward_coeff7*c(3)*c(4)-reverse_coeff7*c(5)*c(2)
    dy(2) = -forward_coeff3*c(3)*c(6)+ reverse_coeff3*c(2)+forward_coeff4*  &
            c(2)*c(7)-reverse_coeff4*c(3)-forward_coeff7*c(3)*c(4)+         &
            reverse_coeff7*c(5)*c(2)
    dy(3) = -forward_coeff6*c(4)*c(7)+ reverse_coeff6*c(5)-forward_coeff7*  &
            c(3)*c(4)+reverse_coeff7*c(5)*c(2)
    dy(4) = forward_coeff6*c(4)*c(7)- reverse_coeff6*c(5)+forward_coeff7*   &
            c(3)*c(4)-reverse_coeff7*c(5)*c(2)
    dy(5) = 0
!    dy(5) = forward_coeff1*c(1) + c(2)*(reverse_coeff3 -reverse_coeff1*c(6) &
!            ) - c(3)*forward_coeff3*c(6)+ forward_coeff5 - reverse_coeff5*  &
!            c(6)*c(7)
    dy(6) = reverse_coeff2*c(2)-forward_coeff2*c(1)*c(7)-forward_coeff4*    &
            c(2)*c(7)+reverse_coeff4*c(3)+forward_coeff5-reverse_coeff5*    &
            c(6)*c(7)-forward_coeff6*c(4)*c(7)+reverse_coeff6*c(5)
!
    DO i = 0,nscl-2
       dydt(i) = dy(i)
    ENDDO
!
  END FUNCTION dydt
!
!
END MODULE reaction
