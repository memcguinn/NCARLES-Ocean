SUBROUTINE cffti1 (N,WA,IFAC)

    DIMENSION WA(*), IFAC(*), NTRYH(4)

    DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/

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
        IFAC(3) = 2
      END IF
    END IF
  END DO
  IFAC(1) = N
  IFAC(2) = NF
  TPI = 2.*PIMACH(DUM)
  ARGH = TPI/FLOAT(N)
  I = 2
  L1 = 1
  DO K1=1,NF
    IP = IFAC(K1+2)
    LD = 0
    L2 = L1*IP
    IDO = N/L2
    IDOT = IDO+IDO+2
    IPM = IP-1
    DO J=1,IPM
      I1 = I
      WA(I-1) = 1.
      WA(I) = 0.
      LD = LD+L1
      FI = 0.
      ARGLD = FLOAT(LD)*ARGH
      DO II=4,IDOT,2
        I = I+2
        FI = FI+1.
        ARG = FI*ARGLD
        WA(I-1) = COS(ARG)
        WA(I) = SIN(ARG)
      END DO
      IF (IP .LE. 5) THEN
        CONTINUE
      ELSE
        WA(I1-1) = WA(I-1)
        WA(I1) = WA(I)
      END IF
    END DO
    L1 = L2
  END DO

RETURN
END SUBROUTINE cffti1
