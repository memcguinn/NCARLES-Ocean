FUNCTION ran1(idum)
! STOLEN FROM NUMERICAL RECIPES, P. 271

  INTEGER :: idum, ia, im, iq, ir, ntab, ndiv
  REAL :: ran1, am, eps, rnmx
  PARAMETER :: ia=16807,im=2147483647,am=1.0/im,iq=127773,ir=2836.0,        &
        ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-07,rnmx=1.0-eps)
  INTEGER :: j, k, iv(ntab), iy

  SAVE iv, iy
  DATA iv /ntab*0/, iy /0/

  IF(idum .LE. 0 .OR. iy .EQ. 0) THEN
    idum = MAX(-idum,1)
    DO j=ntab+8,1,-1
      k = idum/iq
      idum = ia*(idum - k*iq) - ir*k
      IF(idum .LT. 0) idum = idum + im
      IF(j .LE. ntab) iv(j) = idum
    ENDDO
    iy = iv(1)
  ENDIF

  k     = idum/iq
  idum  = ia*(idum - k*iq) - ir*k

  IF(idum .LT. 0) idum = idum + im

  j     = 1 + iy/ndiv
  iy    = iv(j)
  iv(j) = idum
  ran1  = MIN(am*iy, rnmx)

  RETURN
END

! ============================================================================ !
FUNCTION ranf()

  DATA inc /1/
  SAVE inc, ix, ia, m, fm

  IF(inc.EQ.1) THEN
    inc = 2
    m = 2**20
    fm = FLOAT(m)
    ix = 566387
    ia = 2**10 + 3
  ENDIF

  ix = MOD(ia*ix,m)
  fx = FLOAT(ix)
  ranf = fx/fm

  RETURN
END
