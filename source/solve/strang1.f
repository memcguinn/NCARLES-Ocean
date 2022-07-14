      subroutine strang1(it)
c
c ----- strang splitting of scalar reaction term - 0.5*react, advect, 0.5*react
c       for fast reactions (tau<=1000), see RHS/RHS_SCL.F for slow reaction sources
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use reaction, only: react_src
      include 'mpif.h'
c
c ------ temp scalar array to hold rhs from step n-1 and
c        field variables from step n
c
      real trhs(nnx,iys:iye,nscl,izs:ize)
      real, dimension(nscl-1) :: tmp
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         tmp = react_src(ix,iy,1,iz)
         do l=2,nscl
            trhs(ix,iy,l,iz) = tmp(l-1)
                  
            if(trhs(ix,iy,l,iz).le.1.0e-20)then
               trhs(ix,iy,l,iz) = 1.0e-20
            endif
         enddo
      enddo
      enddo
      enddo
c
      do iz=izs,ize
         do l=2,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo   
c
      return
      end
