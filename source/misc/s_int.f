      subroutine s_int(r_min,r_max,value)
c
c ---------- get integral using a mid-point rule
c            stolen from numerical recipes
c
c
      iter = 10
      value = 0
      do j=1,iter
         call midpnt(r_min,r_max,value,j)

 1000 format(' j = ',i5,' value = ',e15.6)
      enddo
c
      return
      end
