SUBROUTINE cfftf1 (N,C,CH,WA,IFAC)

    DIMENSION CH(*), C(*), WA(*), IFAC(*)

  NF = IFAC(2)
  NA = 0
  L1 = 1
  IW = 1
  DO K1=1,NF
    IP = IFAC(K1+2)
    L2 = IP*L1
    IDO = N/L2
    IDOT = IDO+IDO
    IDL1 = IDOT*L1
    IF (IP .NE. 4) THEN
      IF (IP .NE. 2) THEN
        IF (IP .NE. 3) THEN
          IF (IP .NE. 5) THEN
            IF (NA .NE. 0) THEN
              CALL passf (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
            ELSE
              CALL passf (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
            END IF
            IF (NAC .NE. 0) NA = 1-NA
          ELSE
            IX2 = IW+IDOT
            IX3 = IX2+IDOT
            IX4 = IX3+IDOT
            IF (NA .NE. 0) THEN
              CALL passf5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            ELSE
              CALL passf5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            END IF
            NA = 1-NA
          END IF
        ELSE
          IX2 = IW+IDOT
          IF (NA .NE. 0) THEN
            CALL passf3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
          ELSE
            CALL passf3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
          END IF
          NA = 1-NA
        END IF
      ELSE
        IF (NA .NE. 0) THEN
          CALL passf2 (IDOT,L1,CH,C,WA(IW))
        ELSE
          CALL passf2 (IDOT,L1,C,CH,WA(IW))
        END IF
        NA = 1-NA
      END IF
    ELSE
      IX2 = IW+IDOT
      IX3 = IX2+IDOT
      IF (NA .NE. 0) THEN
        CALL passf4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
      ELSE
        CALL passf4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
      END IF
      NA = 1-NA
    END IF
    L1 = L2
    IW = IW+(IP-1)*IDOT
  END DO
  IF (NA .EQ. 0) THEN
    CONTINUE
  ELSE
    N2 = N+N
    DO I=1,N2
      C(I) = CH(I)
    END DO
  END IF

RETURN
END SUBROUTINE cfftf1
