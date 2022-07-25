SUBROUTINE nblnk(word)

  PARAMETER :: nmax=304
  CHARACTER wordt*304, word*(*)

  nchar = len(word)

  IF(nchar .GT. nmax) THEN
    WRITE(6,6000) nchar,nmax
    6000  FORMAT(' TROUBLE, IN SR. NBLNK : NCHAR = ',i6,' EXCEEDS NMAX = ',i6)
    STOP
  ENDIF

  jj = 0
  DO j=1,nchar
    IF(word(j:j) .NE. ' ') THEN
      jj = jj + 1
      wordt(jj:jj) = word(j:j)
    ENDIF
    word(j:j) = ' '
  ENDDO

  DO j=1,jj
    word(j:j) = wordt(j:j)
  ENDDO

  RETURN
END SUBROUTINE
