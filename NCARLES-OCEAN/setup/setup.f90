! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE setup(it)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
! --------------------------------------------------------------------------- !
!
  it = iti
  it_counter = it - iti
!
! TURN ON NEW SGS MODEL AT DESIGNATED STEP 
  IF (it .NE. 0) THEN
    IF (l_root) THEN
      WRITE (6,6) iti, utau, tsfcc(1), qstar(1)
      WRITE (6,510)
      WRITE (6,520) wtsfc(1), utau, amonin, dtdzf(1), zody, zo, cdbtm, ugcont
      CALL print(6,it,1,nnz)
    END IF
!
    IF (l_debug) THEN
      WRITE (nprt,6) iti, utau, tsfcc(1), qstar(1)
      WRITE (nprt,510)
      WRITE (nprt,520) wtsfc(1), utau, amonin, dtdzf(1), zody, zo, cdbtm, ugcont
      CALL print(nprt,it,izs,ize)
    END IF
  END IF
!
  IF (l_root) THEN
    WRITE (6,1) nnx,nny,nnz,ismlt,iti,itmax,iupwnd,ibuoy, &
                noalis,itcut,dt,zo,tsfcc(1),method
  END IF
!
  IF (l_debug) THEN
    WRITE (nprt,1) nnx,nny,nnz,ismlt,iti,itmax,iupwnd,ibuoy, &
                   noalis,itcut,dt,zo,tsfcc(1),method
  END IF
!
  DO i=1,nnx
    xkn(i) = float(i-1)*pi2/xl
    IF(i .GT. ncx) xkn(i) = -float(nnx-i+1)*pi2/xl
  END DO
!
  fn = 1.0/float(nnx)
  DO i=1,nnx
    xk(i) = xkn(i)*fn
  END DO
!
  DO i=1,nny
    ykn(i) = float(i-1)*pi2/yl
    IF (i .GT. ncy) ykn(i) = -float(nny-i+1)*pi2/yl
  END DO
!
  fn = 1.0/float(nny)
  DO i=1,nny
    yk(i) = ykn(i)*fn
  END DO
!
  ii = -1
  DO i=1,ncx
    ii = ii + 2
    temp = xkn(i)**2
    DO j=1,nny
      temp1       = temp + ykn(j)**2
      xks(ii,j)   = temp1
      xks(ii+1,j) = temp1
    END DO
  END DO
!
  xnn = ABS(batag*dtdzf(1))
!
! ENSURE GRAVITY WAVES PROPAGATE OUT OF DOMAIN
  sgn = -1.0
  IF (ibcu .EQ. 1) THEN
    DO iy=1,nny
      DO ix=1,nnxp2
        IF (xks(ix,iy) .LE. 0.) THEN
          wavexy(ix,iy) = 0.0
        ELSE
          wavexy(ix,iy) = sgn*SQRT(xnn/xks(ix,iy))
        END IF
      END DO
    END DO
  END IF
!
! SET LENGTH SCALE FOR SGS MODEL
  IF (iz_space .EQ. 0) THEN
!
!   UNIFORM VERTICAL SPACING, NEED FOR REACTION MODULES
    dx32 = dx*3./2.
    dy32 = dy*3./2.
    dsl  = (ABS(dx32*dy32*dzw(1)))**(1./3.)
    dslg = dsl
!
    IF (l_root) WRITE (6,2000) dsl
    IF (l_debug) WRITE (nprt,2000) dsl
!
!   CREATE DSL ARRAY FOR EASY INDEXING, COMP1
    DO iz=0,nnzp1
      dsl_z(iz) = dslg
    END DO
!
!   VARIABLE VERTICAL SPACING
  ELSE
!
!   ESTIMATE DSL ARRAY, USE AVERAGE SPACING
    dx32 = dx*3./2.
    dy32 = dy*3./2.
    dsl_max = (ABS(dx32*dy32*dzw(0)))**(1./3.)
!
    DO iz=0,nnzp1
      dsl_z(iz) = (ABS(dx32*dy32*dzw(iz)))**(1./3.)
      IF (dsl_z(iz) .GT. dsl_max) dsl_max = dsl_z(iz)
    END DO
!
    dsl  = dsl_max
    dslg = dsl
  END IF
!
  gridr = 1.0
  sml_eg = smal_e*gridr
!
! SET STOKES VELOCITY FOR OCEANIC FLOW
  CALL stokesv
!
! CAN ADD TIME FACTOR TO SKIP ANY PART OF SPECIFIFIED GEOSTROPHIC ARRAYS [SEC]
  t_factor = 7200.0
!
! GET MORE DIGITS FOR PRINT OUT
  t_ref = 290.16
!
! SET LOWER BOUND AT ZI_MIN FOR ZI
  zi_min = -5.0
  iz_min = 1
!
  DO iz=1,nnz-1
    IF (zz(iz) .LT. zi_min .AND. zz(iz+1) .GE. zi_min) iz_min = iz
  END DO
!
  IF (l_root) THEN
    WRITE (6,9000) zi_min, iz_min
  END IF
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
    6 FORMAT (///,' DATA FROM RESTART FILE AT STEP =',I5,           &
             ' U_* = ',e15.6,' TS = ',e15.6,' Q_* = ',e15.6,///)
  510 FORMAT (' RESTART ***** CASE WITH : ******',/)
  520 FORMAT (' WT = ',e12.4,', U_* = ',e12.4,', L = ',e12.4,       &
              ', DTDZ FREE = ',e12.4,', ZODY = ',e12.4,/,10x,       &
              '  ZO(BTM) = ',e12.4,', CDBTM = ',e12.4,', UG = ',e12.4)
    1 FORMAT (10x,' NNX = ',i5,',  NNY = ',i5,',  NNZ = ',i5,/,     &
              10x,' SFC SMLT = ',i1,',  FILTER = ',i1,',  ITI = ',  &
              i6,',  ITMAX = ',i6,/,10x,' IUPWIND = ',i1,           &
              ',  BUYNCY = ',i1,',  NO ALISNG = ',i1,',  ITCUT = ', &
              i1,/,10x,' DT = ',e15.6,',  ZO = ',e15.6,',  TS = ',  &
              e15.6,',  SUBSD = ',i1,/,10x,' BRCLICITY = ',i1,',  METHOD = ',i1)
 2000 FORMAT (10x,' DSL = ',e15.6)
 9000 FORMAT (' Search for zi above the height = ',e15.6,/,' iz_min = ',i5)
!------------------------------------------------------------------------------!
!
END SUBROUTINE setup
