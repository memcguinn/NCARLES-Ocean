SUBROUTINE cosqf1 (N,X,W,XH)

    DIMENSION X(*), W(*), XH(*)

  NS2 = (N+1)/2
  NP2 = N+2
  DO K=2,NS2
    KC = NP2-K
    XH(K) = X(K)+X(KC)
    XH(KC) = X(K)-X(KC)
  END DO
  MODN = MOD(N,2)
  IF (MODN .EQ. 0) XH(NS2+1) = X(NS2+1)+X(NS2+1)
  DO K=2,NS2
    KC = NP2-K
    X(K) = W(K-1)*XH(KC)+W(KC-1)*XH(K)
    X(KC) = W(K-1)*XH(K)-W(KC-1)*XH(KC)
  END DO
  IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*XH(NS2+1)
  CALL rfftf (N,X,XH)
  DO I=3,N,2
    XIM1 = X(I-1)-X(I)
    X(I) = X(I-1)+X(I)
    X(I-1) = XIM1
  END DO

RETURN
END SUBROUTINE cosqf1
