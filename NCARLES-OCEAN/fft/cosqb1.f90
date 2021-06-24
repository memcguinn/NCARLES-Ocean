SUBROUTINE cosqb1 (N,X,W,XH)

    DIMENSION X(*), W(*), XH(*)

  NS2 = (N+1)/2
  NP2 = N+2
  DO I=3,N,2
    XIM1 = X(I-1)+X(I)
    X(I) = X(I)-X(I-1)
    X(I-1) = XIM1
  END DO
  X(1) = X(1)+X(1)
  MODN = MOD(N,2)
  IF (MODN .EQ. 0) X(N) = X(N)+X(N)
  CALL rfftb (N,X,XH)
  DO K=2,NS2
    KC = NP2-K
    XH(K) = W(K-1)*X(KC)+W(KC-1)*X(K)
    XH(KC) = W(K-1)*X(K)-W(KC-1)*X(KC)
  END DO
  IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
  DO K=2,NS2
    KC = NP2-K
    X(K) = XH(K)+XH(KC)
    X(KC) = XH(K)-XH(KC)
  END DO
  X(1) = X(1)+X(1)

RETURN
END SUBROUTINE cosqb1
