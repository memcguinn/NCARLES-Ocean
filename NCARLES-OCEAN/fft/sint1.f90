SUBROUTINE sint1 (N,WAR,WAS,XH,X,IFAC)

    DIMENSION WAR(*),WAS(*),X(*),XH(*),IFAC(*)

    DATA SQRT3 /1.73205080756888/

  DO I=1,N
    XH(I) = WAR(I)
    WAR(I) = X(I)
  END DO
  IF (N-2 .LT. 0) THEN
    XH(1) = XH(1)+XH(1)
  ELSE IF (N-2 .EQ. 0) THEN
    XHOLD = SQRT3*(XH(1)+XH(2))
    XH(2) = SQRT3*(XH(1)-XH(2))
    XH(1) = XHOLD
  ELSE !(N-2 .GT. 0)
    NP1 = N+1
    NS2 = N/2
    X(1) = 0.
    DO K=1,NS2
      KC = NP1-K
      T1 = XH(K)-XH(KC)
      T2 = WAS(K)*(XH(K)+XH(KC))
      X(K+1) = T1+T2
      X(KC+1) = T2-T1
    END DO
    MODN = MOD(N,2)
    IF (MODN .NE. 0) X(NS2+2) = 4.*XH(NS2+1)
    CALL rfftf1 (NP1,X,XH,WAR,IFAC)
    XH(1) = .5*X(1)
    DO I=3,N,2
      XH(I-1) = -X(I)
      XH(I) = XH(I-2)+X(I-1)
    END DO
    IF (MODN .NE. 0) THEN
      CONTINUE
    ELSE
      XH(N) = -X(N+1)
    END IF
  END IF
  DO I=1,N
    X(I) = WAR(I)
    WAR(I) = XH(I)
  END DO

RETURN
END
