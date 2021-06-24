SUBROUTINE ezfft1 (N,WA,IFAC)

    DIMENSION WA(*), IFAC(*), NTRYH(4)

    DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/

  TPI = 2.0*PIMACH(DUM)
  NL = N
  NF = 0
  J = 0
  J = J+1
  IF (J-4 .LE. 0) THEN
    NTRY = NTRYH(J)
  ELSE !(J-4 .GT. 0)
    NTRY = NTRY+2
  END IF
  DO WHILE (NL .NE. 1)
    NQ = NL/NTRY
    NR = NL-NTRY*NQ
    IF (NR .LT. 0 .OR. NR .GT. 0) THEN
      J = J+1
      IF (J-4 .LE. 0) THEN
        NTRY = NTRYH(J)
      ELSE !(J-4 .GT. 0)
        NTRY = NTRY+2
      END IF
      NQ = NL/NTRY
      NR = NL-NTRY*NQ
    ELSE !(NR .EQ. 0)
      CONTINUE
    END IF
    NF = NF+1
    IFAC(NF+2) = NTRY
    NL = NQ
    IF (NTRY .NE. 2) THEN
      CONTINUE
    ELSE
      IF (NF .EQ. 1) THEN
        CONTINUE
      ELSE
        DO I=2,NF
          IB = NF-I+2
          IFAC(IB+2) = IFAC(IB+1)
        END DO
      END IF
    END IF
    IFAC(3) = 2
  END DO
  IFAC(1) = N
  IFAC(2) = NF
  ARGH = TPI/FLOAT(N)
  IS = 0
  NFM1 = NF-1
  L1 = 1
  IF (NFM1 .EQ. 0) THEN
    CONTINUE
  ELSE
    DO K1=1,NFM1
      IP = IFAC(K1+2)
      L2 = L1*IP
      IDO = N/L2
      IPM = IP-1
      ARG1 = FLOAT(L1)*ARGH
      CH1 = 1.
      SH1 = 0.
      DCH1 = COS(ARG1)
      DSH1 = SIN(ARG1)
      DO J=1,IPM
        CH1H = DCH1*CH1-DSH1*SH1
        SH1 = DCH1*SH1+DSH1*CH1
        CH1 = CH1H
        I = IS+2
        WA(I-1) = CH1
        WA(I) = SH1
        IF (IDO .LT. 5) THEN
          CONTINUE
        ELSE
          DO II=5,IDO,2
            I = I+2
            WA(I-1) = CH1*WA(I-3)-SH1*WA(I-2)
            WA(I) = CH1*WA(I-2)+SH1*WA(I-3)
          END DO
        END IF
        IS = IS+IDO
      END DO
      L1 = L2
    END DO
  END IF

RETURN
END SUBROUTINE ezfft1
