SUBROUTINE radf3 (IDO,L1,CC,CH,WA1,WA2)

    DIMENSION CH(IDO,3,L1), CC(IDO,L1,3), WA1(*), WA2(*)

    DATA TAUR,TAUI /-.5,.866025403784439/

  DO K=1,L1
    CR2 = CC(1,K,2)+CC(1,K,3)
    CH(1,1,K) = CC(1,K,1)+CR2
    CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
    CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  END DO
  IF (IDO .EQ. 1) THEN
    CONTINUE
  ELSE
    IDP2 = IDO+2
    DO K=1,L1
      DO I=3,IDO,2
        IC = IDP2-I
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR2 = DR2+DR3
        CI2 = DI2+DI3
        CH(I-1,1,K) = CC(I-1,K,1)+CR2
        CH(I,1,K) = CC(I,K,1)+CI2
        TR2 = CC(I-1,K,1)+TAUR*CR2
        TI2 = CC(I,K,1)+TAUR*CI2
        TR3 = TAUI*(DI2-DI3)
        TI3 = TAUI*(DR3-DR2)
        CH(I-1,3,K) = TR2+TR3
        CH(IC-1,2,K) = TR2-TR3
        CH(I,3,K) = TI2+TI3
        CH(IC,2,K) = TI3-TI2
      END DO
    END DO
  END IF

RETURN
END SUBROUTINE radf3
