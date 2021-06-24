SUBROUTINE radf4 (IDO,L1,CC,CH,WA1,WA2,WA3)

    DIMENSION CC(IDO,L1,4), CH(IDO,4,L1), WA1(*), WA2(*), WA3(*)

    DATA HSQT2 /.7071067811865475/

  DO K=1,L1
    TR1 = CC(1,K,2)+CC(1,K,4)
    TR2 = CC(1,K,1)+CC(1,K,3)
    CH(1,1,K) = TR1+TR2
    CH(IDO,4,K) = TR2-TR1
    CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
    CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  END DO
  IF (IDO-2 .LT. 0) THEN
    CONTINUE
  ELSE IF (IDO-2 .EQ. 0) THEN
    DO K=1,L1
      TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
      TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
      CH(IDO,1,K) = TR1+CC(IDO,K,1)
      CH(IDO,3,K) = CC(IDO,K,1)-TR1
      CH(1,2,K) = TI1-CC(IDO,K,3)
      CH(1,4,K) = TI1+CC(IDO,K,3)
    END DO
  ELSE !(ID O .GT. 0)
    IDP2 =  IDO+2
    DO K=1,L1
      DO I=3 ,IDO,2
        IC  = IDP2-I
        CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        TR1 = CR2+CR4
        TR4 = CR4-CR2
        TI1 = CI2+CI4
        TI4 = CI2-CI4
        TI2 = CC(I,K,1)+CI3
        TI3 = CC(I,K,1)-CI3
        TR2 = CC(I-1,K,1)+CR3
        TR3 = CC(I-1,K,1)-CR3
        CH(I-1,1,K) = TR1+TR2
        CH(IC-1,4,K) = TR2-TR1
        CH(I,1,K) = TI1+TI2
        CH(IC,4,K) = TI1-TI2
        CH(I-1,3,K) = TI4+TR3
        CH(IC-1,2,K) = TR3-TI4
        CH(I,3,K) = TR4+TI3
        CH(IC,2,K) = TR4-TI3
      END DO
    END DO
    IF (MOD(IDO,2) .EQ. 1) THEN
      CONTINUE
    ELSE
      DO K=1,L1
        TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
        TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
        CH(IDO,1,K) = TR1+CC(IDO,K,1)
        CH(IDO,3,K) = CC(IDO,K,1)-TR1
        CH(1,2,K) = TI1-CC(IDO,K,3)
        CH(1,4,K) = TI1+CC(IDO,K,3)
      END DO
    END IF
  END IF

RETURN
END SUBROUTINE radf4
