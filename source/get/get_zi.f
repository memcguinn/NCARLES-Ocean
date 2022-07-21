      subroutine get_zi(gradmax,gradout,len,itype)
c
      use pars
      real gradmax(*), gradout(*)
c
      do i=1,len,2
         if(gradmax(i) .gt. gradout(i)) then
              gradout(i)   = gradmax(i)
              gradout(i+1) = gradmax(i+1)
         endif
      enddo
c
      return
      end
