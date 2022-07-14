!     this module contains the interface to a generic byteswapping
!     routine called "byteswap", that will switch the argument "data"
!     between big and little endian format.
!     usage:
!     use module_byteswap          
!     ...
!     call byteswap(data)
!     ...      
!     data can be real or integer, kind=4 or kind=8. the routine handles
!     scalar and 1d-5d real arrays, scalar and 1d-3d integer arrays.
!     the whole implementation rests on vendor extensions, but since
!     endianness is a vendor thing anyway... seems to work for IEEE
!     unformatted data.
!     
!     written by Helge Avlesen <avle@ii.uib.no>, para//ab
!     
module module_byteswap

  interface byteswap
     module procedure &
          byteswap_real4,&
          byteswap_real8,&
          byteswap_1dreal4,&
          byteswap_1dreal8,&
          byteswap_2dreal4,&
          byteswap_2dreal8,&
          byteswap_3dreal4,&
          byteswap_3dreal8,&
          byteswap_4dreal4,&
          byteswap_4dreal8,&
          byteswap_5dreal4,&
          byteswap_5dreal8,&
          byteswap_integer4,&
          byteswap_integer8,&
          byteswap_1dinteger4,&
          byteswap_1dinteger8,&
          byteswap_2dinteger4,&
          byteswap_2dinteger8,&
          byteswap_3dinteger4,&
          byteswap_3dinteger8
  end interface

  integer, parameter :: kind4=4, kind8=8

contains

!     
!     slightly ugly hack for converting the reals in data(im,jm,km)
!     between little and big endian. the equivalence trick is legal
!     Fortran 77, but not 9X. we also assume kind to return the number
!     of bytes used in the representation of data. this may not always
!     be the case but all compilers I know (ifort,pgf,xlf,sgi
!     f90,g95)use this scheme. the routine of course also assume the
!     same floating point representation (except for endianess), but
!     this is also true for most compilers as variations on the IEEE
!     floating point format is more and more used.
!     
  subroutine byteswap_real4array(data,im,jm,km,lm,mm)
    implicit none
    integer im,jm,km,lm,mm,i,j,k,l,m,n
    real(kind4) data(im,jm,km,lm,mm), target
    integer, parameter :: ul=5, hl=2
    character tmpchar, swaparr(4)
    equivalence (swaparr(1),target)

    do m=1,mm
        do l=1,lm         
            do k=1,km
                do j=1,jm
                    do i=1,im
                        target = data(i,j,k,l,m)
                        do n=1,hl
                            tmpchar = swaparr(n)
                            swaparr(n) = swaparr(ul-n)
                            swaparr(ul-n) = tmpchar
                        end do
                        data(i,j,k,l,m) = target
                    end do
                end do
            end do
        end do
    end do
  end subroutine byteswap_real4array

  subroutine byteswap_real8array(data,im,jm,km,lm,mm)
    implicit none
    integer im,jm,km,lm,mm,i,j,k,l,m,n
    real(kind8) data(im,jm,km,lm,mm), target
    integer, parameter :: ul=9, hl=4
    character tmpchar, swaparr(8)
    equivalence (swaparr(1),target)
    do m=1,mm
        do l=1,lm
            do k=1,km
                do j=1,jm
                    do i=1,im
                        target = data(i,j,k,l,m)
                        do n=1,hl
                            tmpchar = swaparr(n)
                            swaparr(n) = swaparr(ul-n)
                            swaparr(ul-n) = tmpchar
                        end do
                        data(i,j,k,l,m) = target
                    end do
                end do
            end do
        end do
    end do
  end subroutine byteswap_real8array

! the same as above for integers.

  subroutine byteswap_integer4array(data,im,jm,km)
    implicit none
    integer im,jm,km,i,j,k,n
    integer(kind4) data(im,jm,km), target
    integer, parameter :: ul=5, hl=2
    character tmpchar, swaparr(4)
    equivalence (swaparr(1),target)

    do k=1,km
        do j=1,jm
            do i=1,im
                target = data(i,j,k)
                do n=1,hl
                    tmpchar = swaparr(n)
                    swaparr(n) = swaparr(ul-n)
                    swaparr(ul-n) = tmpchar
                end do
                data(i,j,k) = target
            end do
        end do
    end do
  end subroutine byteswap_integer4array

  subroutine byteswap_integer8array(data,im,jm,km)
    implicit none
    integer im,jm,km,i,j,k,n
    integer(kind8) data(im,jm,km), target
    integer, parameter :: ul=9, hl=4
    character tmpchar, swaparr(8)
    equivalence (swaparr(1),target)

    do k=1,km
        do j=1,jm
            do i=1,im
                target = data(i,j,k)
                do n=1,hl
                    tmpchar = swaparr(n)
                    swaparr(n) = swaparr(ul-n)
                    swaparr(ul-n) = tmpchar
                end do
                data(i,j,k) = target
            end do
        end do
    end do
  end subroutine byteswap_integer8array

! so. let's cut&paste some macros to fill in the generic interface:

  subroutine byteswap_real4(data)
    real(kind4) data, tmpdata(1,1,1,1,1)
    tmpdata(1,1,1,1,1)=data
    call byteswap_real4array(tmpdata,1,1,1,1,1)
    data=tmpdata(1,1,1,1,1)
  end subroutine byteswap_real4

  subroutine byteswap_real8(data)
    real(kind8) data, tmpdata(1,1,1,1,1)
    tmpdata(1,1,1,1,1)=data
    call byteswap_real8array(tmpdata,1,1,1,1,1)
    data=tmpdata(1,1,1,1,1)
  end subroutine byteswap_real8

  subroutine byteswap_1dreal4(data)
    real(kind4) data(:)
    call byteswap_real4array(data,size(data,1),1,1,1,1)
  end subroutine byteswap_1dreal4

  subroutine byteswap_1dreal8(data)
    real(kind8) data(:)
    call byteswap_real8array(data,size(data,1),1,1,1,1)
  end subroutine byteswap_1dreal8

  subroutine byteswap_2dreal4(data)
    real(kind4) data(:,:)
    call byteswap_real4array(data,size(data,1),size(data,2),1,1,1)
  end subroutine byteswap_2dreal4

  subroutine byteswap_2dreal8(data)
    real(kind8) data(:,:)
    call byteswap_real8array(data,size(data,1),size(data,2),1,1,1)
  end subroutine byteswap_2dreal8

  subroutine byteswap_3dreal4(data)
    real(kind4) data(:,:,:)
    call byteswap_real4array(&
         data,size(data,1),size(data,2),size(data,3),1,1)
  end subroutine byteswap_3dreal4

  subroutine byteswap_3dreal8(data)
    real(kind8) data(:,:,:)
    call byteswap_real8array(&
         data,size(data,1),size(data,2),size(data,3),1,1)
  end subroutine byteswap_3dreal8

  subroutine byteswap_4dreal4(data)
    real(kind4) data(:,:,:,:)
    call byteswap_real4array(&
         data,size(data,1),size(data,2),size(data,3),size(data,4),1)
  end subroutine byteswap_4dreal4

  subroutine byteswap_4dreal8(data)
    real(kind8) data(:,:,:,:)
    call byteswap_real8array(&
         data,size(data,1),size(data,2),size(data,3),size(data,4),1)
  end subroutine byteswap_4dreal8

  subroutine byteswap_5dreal4(data)
    real(kind4) data(:,:,:,:,:)
    call byteswap_real4array(data,size(data,1),size(data,2),&
         size(data,3),size(data,4),size(data,5))
  end subroutine byteswap_5dreal4

  subroutine byteswap_5dreal8(data)
    real(kind8) data(:,:,:,:,:)
    call byteswap_real8array(data,size(data,1),size(data,2),&
         size(data,3),size(data,4),size(data,5))
  end subroutine byteswap_5dreal8

  subroutine byteswap_integer4(data)
    integer(kind4) data,tmpdata(1,1,1)
    tmpdata(1,1,1)=data
    call byteswap_integer4array(tmpdata,1,1,1)
    data=tmpdata(1,1,1)
  end subroutine byteswap_integer4

  subroutine byteswap_1dinteger4(data)
    integer(kind4) data(:)
    call byteswap_integer4array(data,size(data,1),1,1)
  end subroutine byteswap_1dinteger4

  subroutine byteswap_2dinteger4(data)
    integer(kind4) data(:,:)
    call byteswap_integer4array(data,size(data,1),size(data,2),1)
  end subroutine byteswap_2dinteger4

  subroutine byteswap_3dinteger4(data)
    integer(kind4) data(:,:,:)
    call byteswap_integer4array(&
         data,size(data,1),size(data,2),size(data,3))
  end subroutine byteswap_3dinteger4

  subroutine byteswap_integer8(data)
    integer(kind8) data,tmpdata(1,1,1)
    tmpdata(1,1,1)=data
    call byteswap_integer8array(tmpdata,1,1,1)
    data=tmpdata(1,1,1)
  end subroutine byteswap_integer8

  subroutine byteswap_1dinteger8(data)
    integer(kind8) data(:)
    call byteswap_integer8array(data,size(data,1),1,1)
  end subroutine byteswap_1dinteger8

  subroutine byteswap_2dinteger8(data)
    integer(kind8) data(:,:)
    call byteswap_integer8array(data,size(data,1),size(data,2),1)
  end subroutine byteswap_2dinteger8

  subroutine byteswap_3dinteger8(data)
    integer(kind8) data(:,:,:)
    call byteswap_integer8array(&
         data,size(data,1),size(data,2),size(data,3))
  end subroutine byteswap_3dinteger8

end module module_byteswap

