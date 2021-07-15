SUBROUTINE sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
! sharp cutoff filter for field wve stored in 2d-fft order

    USE pars

    REAL wve(nny,jxs:jxe,izs:ize)

  DO iz=izs,ize
    DO ix=jxs,jxe
      DO iy=iy_cut_l,iy_cut_u
        wve(iy,ix,iz) = 0.0
      END DO
    END DO
  END DO
  IF (jxe .LT. ix_cut) THEN
    CONTINUE
  ELSE
    DO iz=izs,ize
      DO ix=MAX(jxs,ix_cut),jxe
        DO iy=1,nny
          wve(iy,ix,iz) = 0.0
        END DO
      END DO
    END DO
  END IF

RETURN
END SUBROUTINE sharp
