SUBROUTINE rfftb1 (N,C,CH,WA,IFAC)

    DIMENSION CH(*), C(*), WA(*), IFAC(*)

  NF = IFAC(2)
  NA = 0
  L1 = 1
  IW = 1
  DO K1=1,NF
    IP = IFAC(K1+2)
    L2 = IP*L1
    IDO = N/L2
    IDL1 = IDO*L1
    IF (IP .NE. 4) THEN
      IF (IP .NE. 2) THEN
        IF (IP .NE. 3) THEN
          IF (IP .NE. 5) THEN
            IF (NA .NE. 0) THEN
              CALL radbg (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
            ELSE
              CALL radbg (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
            END IF
            IF (IDO .EQ. 1) NA = 1-NA
          ELSE
            IX2 = IW+IDO
            IX3 = IX2+IDO
            IX4 = IX3+IDO
            IF (NA .NE. 0) THEN
              CALL radb5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            ELSE
              CALL radb5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            END IF
            NA = 1-NA
          END IF
        ELSE
          IX2 = IW+IDO
          IF (NA .NE. 0) THEN
            CALL radb3 (IDO,L1,CH,C,WA(IW),WA(IX2))
          ELSE
            CALL radb3 (IDO,L1,C,CH,WA(IW),WA(IX2))
          END IF
          NA = 1-NA
        END IF
      ELSE
        IF (NA .NE. 0) THEN
          CALL radb2 (IDO,L1,CH,C,WA(IW))
        ELSE
          CALL radb2 (IDO,L1,C,CH,WA(IW))
        END IF
        NA = 1-NA
      END IF
    ELSE
      IX2 = IW+IDO
      IX3 = IX2+IDO
      IF (NA .NE. 0) THEN
        CALL radb4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
      ELSE
        CALL radb4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
      END IF
      NA = 1-NA
    END IF
    L1 = L2
    IW = IW+(IP-1)*IDO
  END DO
  IF (NA .EQ. 0) THEN
    CONTINUE
  ELSE
    DO I=1,N
      C(I) = CH(I)
    END DO
  END IF

RETURN
END SUBROUTINE rfftb1
