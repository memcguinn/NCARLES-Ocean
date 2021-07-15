SUBROUTINE passb (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)

    DIMENSION CH(IDO,L1,IP), CC(IDO,IP,L1), C1(IDO,L1,IP), WA(*), &
              C2(IDL1,IP), CH2(IDL1,IP)

  IDOT = IDO/2
  NT = IP*IDL1
  IPP2 = IP+2
  IPPH = (IP+1)/2
  IDP = IP*IDO
  IF (IDO .LT. L1) THEN
    DO J=2,IPPH
      JC = IPP2-J
      DO I=1,IDO
        DO K=1,L1
          CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
          CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
        END DO
      END DO
    END DO
    DO I=1,IDO
      DO K=1,L1
        CH(I,K,1) = CC(I,1,K)
      END DO
    END DO
  ELSE
    DO J=2,IPPH
      JC = IPP2-J
      DO K=1,L1
        DO I=1,IDO
          CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
          CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
        END DO
      END DO
    END DO
    DO K=1,L1
      DO I=1,IDO
        CH(I,K,1) = CC(I,1,K)
      END DO
    END DO
  END IF
  IDL = 2-IDO
  INC = 0
  DO L=2,IPPH
    LC = IPP2-L
    IDL = IDL+IDO
    DO IK=1,IDL1
      C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
      C2(IK,LC) = WA(IDL)*CH2(IK,IP)
    END DO
    IDLJ = IDL
    INC = INC+IDO
    DO J=3,IPPH
      JC = IPP2-J
      IDLJ = IDLJ+INC
      IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
      WAR = WA(IDLJ-1)
      WAI = WA(IDLJ)
      DO IK=1,IDL1
        C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
        C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
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
    DO IK=2,IDL1,2
      CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
      CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
      CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
      CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
    END DO
  END DO
  NAC = 1
  IF (IDO .EQ. 2) THEN
    CONTINUE
  ELSE
    NAC = 0
    DO IK=1,IDL1
      C2(IK,1) = CH2(IK,1)
    END DO
    DO J=2,IP
      DO K=1,L1
        C1(1,K,J) = CH(1,K,J)
        C1(2,K,J) = CH(2,K,J)
      END DO
    END DO
    IF (IDOT .GT. L1) THEN
      IDJ = 2-IDO
      DO J=2,IP
        IDJ = IDJ+IDO
        DO K=1,L1
          IDIJ = IDJ
          DO I=4,IDO,2
            IDIJ = IDIJ+2
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
          END DO
        END DO
      END DO
    ELSE
      IDIJ = 0
      DO J=2,IP
        IDIJ = IDIJ+2
        DO I=4,IDO,2
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
END SUBROUTINE passb
