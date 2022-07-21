      function rlim(d1,d2,d3)
c
c ------------- Cees's kappa=1/3 scheme
c
      r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
      rlim = (d2-d3)*amax1(0.,amin1(r,amin1(1./6.+1./3.*r,1.)))
c
      return
      end
