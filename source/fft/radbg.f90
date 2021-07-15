SUBROUTINE radbg (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)

    DIMENSION CH(IDO,L1,IP), CC(IDO,IP,L1), C1(IDO,L1,IP), C2(IDL1,IP),  &
               CH2(IDL1,IP), WA(*)

  TPI = 2.0*PIMACH(DUM)
  ARG = TPI/FLOAT(IP)
  DCP = COS(ARG)
  DSP = SIN(ARG)
  IDP2 = IDO+2
  NBD = (IDO-1)/2
  IPP2 = IP+2
  IPPH = (IP+1)/2
  IF (IDO .LT. L1) THEN
    DO I=1,IDO
      DO K=1,L1
        CH(I,K,1) = CC(I,1,K)
      END DO
    END DO
  ELSE
    DO K=1,L1
      DO I=1,IDO
        CH(I,K,1) = CC(I,1,K)
      END DO
    END DO
  END IF
  DO J=2,IPPH
    JC = IPP2-J
    J2 = J+J
    DO K=1,L1
      CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
      CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
    END DO
  END DO
  IF (IDO .EQ. 1) THEN
    CONTINUE
  ELSE
    IF (NBD .LT. L1) THEN
      DO J=2,IPPH
        JC = IPP2-J
        DO I=3,IDO,2
          IC = IDP2-I
          DO K=1,L1
            CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
            CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
            CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
            CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
          END DO
        END DO
      END DO ! END OF 112
    ELSE
      DO J=2,IPPH
        JC = IPP2-J
        DO K=1,L1
          DO I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
            CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
            CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
            CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
          END DO
        END DO
      END DO
    END IF
  END IF
  AR1 = 1.
  AI1 = 0.
  DO L=2,IPPH
    LC = IPP2-L
    AR1H = DCP*AR1-DSP*AI1
    AI1 = DCP*AI1+DSP*AR1
    AR1 = AR1H
    DO IK=1,IDL1
      C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
      C2(IK,LC) = AI1*CH2(IK,IP)
    END DO
    DC2 = AR1
    DS2 = AI1
    AR2 = AR1
    AI2 = AI1
    DO J=3,IPPH
      JC = IPP2-J
      AR2H = DC2*AR2-DS2*AI2
      AI2 = DC2*AI2+DS2*AR2
      AR2 = AR2H
      DO IK=1,IDL1
        C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
        C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
      END DO
    END DO
  END DO
  DO J=2,IPPH
    DO IK=1,IDL1
      CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
    END DO
  END DO
  DO J=2,IPPH
    JC = IPP2-J
    DO K=1,L1
      CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
      CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
    END DO
  END DO
  IF (IDO .EQ. 1) THEN
    CONTINUE
  ELSE
    IF (NBD .LT. L1) THEN
      DO J=2,IPPH
        JC = IPP2-J
        DO I=3,IDO,2
          DO K=1,L1
            CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
            CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
            CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
            CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
          END DO
        END DO
      END DO ! END OF 128
    ELSE
      DO J=2,IPPH
        JC = IPP2-J
        DO K=1,L1
          DO I=3,IDO,2
            CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
            CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
            CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
            CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
          END DO
        END DO
      END DO
    END IF
  END IF
  IF (IDO .EQ. 1) THEN
    CONTINUE
  ELSE
    DO IK=1,IDL1
      C2(IK,1) = CH2(IK,1)
    END DO
    DO J=2,IP
      DO K=1,L1
        C1(1,K,J) = CH(1,K,J)
      END DO
    END DO
    IF (NBD .GT. L1) THEN
      IS = -IDO
      DO J=2,IP
        IS = IS+IDO
        DO K=1,L1
          IDIJ = IS
          DO I=3,IDO,2
            IDIJ = IDIJ+2
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
          END DO
        END DO
      END DO ! END OF 139
    ELSE
      IS = -IDO
      DO J=2,IP
        IS = IS+IDO
        IDIJ = IS
        DO I=3,IDO,2
          IDIJ = IDIJ+2
          DO K=1,L1
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
          END DO
        END DO
      END DO
    END IF
  END IF

RETURN
END SUBROUTINE radbg
