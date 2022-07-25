SUBROUTINE range(n1,n2,nprocs,irank,ista,iEND)
! THE IBM RANGE FINDER TO BALANCE LOAD

  iwork1 = (n2 - n1 + 1)/nprocs
  iwork2 = MOD(n2 - n1 +1, nprocs)
  ista = irank*iwork1 + n1 + MIN(irank,iwork2)
  iEND = ista + iwork1 - 1

  IF(iwork2 .gt. irank) iEND = iEND + 1

  RETURN
END SUBROUTINE
