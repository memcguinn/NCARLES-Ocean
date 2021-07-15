SUBROUTINE passb2 (IDO,L1,CC,CH,WA1)

    DIMENSION CC(IDO,2,L1), CH(IDO,L1,2), WA1(1)

  IF (IDO .GT. 2) THEN
    DO K=1,L1
      DO I=2,IDO,2
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
      END DO
    END DO
  ELSE
    DO K=1,L1
      CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
      CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
      CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
      CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
    END DO
  END IF

RETURN
END SUBROUTINE passb2
