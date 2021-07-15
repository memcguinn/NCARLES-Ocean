SUBROUTINE rfftf1 (N,C,CH,WA,IFAC)

    DIMENSION CH(*), C(*), WA(*), IFAC(*)

  NF = IFAC(2)
  NA = 1
  L2 = N
  IW = N
  DO K1=1,NF
    KH = NF-K1
    IP = IFAC(KH+3)
    L1 = L2/IP
    IDO = N/L2
    IDL1 = IDO*L1
    IW = IW-(IP-1)*IDO
    NA = 1-NA
    IF (IP .NE. 4) THEN
      IF (IP .NE. 2) THEN
        IF (IP .NE. 3) THEN
          IF (IP .NE. 5) THEN
            IF (IDO .EQ. 1) NA = 1-NA
            IF (NA .NE. 0) THEN
              CALL radfg (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
              NA = 0
            ELSE
              CALL radfg (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
              NA = 1
            END IF
          ELSE
            IX2 = IW+IDO
            IX3 = IX2+IDO
            IX4 = IX3+IDO
            IF (NA .NE. 0) THEN
              CALL radf5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            ELSE
              CALL radf5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            END IF
          END IF
        ELSE
          IX2 = IW+IDO
          IF (NA .NE. 0) THEN
            CALL radf3 (IDO,L1,CH,C,WA(IW),WA(IX2))
          ELSE
            CALL radf3 (IDO,L1,C,CH,WA(IW),WA(IX2))
          END IF
        END IF ! END OF 104
      ELSE
        IF (NA .NE. 0) THEN
          CALL radf2 (IDO,L1,CH,C,WA(IW))
        ELSE
          CALL radf2 (IDO,L1,C,CH,WA(IW))
        END IF
      END IF ! END OF 102
    ELSE
      IX2 = IW+IDO
      IX3 = IX2+IDO
      IF (NA .NE. 0) THEN
        CALL radf4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
      ELSE
        CALL radf4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
      END IF
    END IF
    L2 = L1
  END DO
  IF (NA .EQ. 1) THEN
    CONTINUE
  ELSE
    DO I=1,N
      C(I) = CH(I)
    END DO
  END IF

RETURN
END SUBROUTINE rfftf1
