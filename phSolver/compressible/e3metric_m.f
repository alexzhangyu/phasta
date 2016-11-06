      module e3metric_m
c
        use workfc_m
        use global_const_m
        use propar_m
        use number_def_m
        implicit none
c
      contains
c
      subroutine e3metric(shg,shgl,xl,WdetJ,qwt)
c
        implicit none
c
        real*8, dimension(:,:,:), pointer, intent(out) :: shg
        real*8, dimension(:,:,:), pointer, intent(in)  :: xl, shgl
        real*8, dimension(:), pointer, optional, intent(out)  :: WdetJ
        real*8, optional, intent(in) :: qwt
c
        real*8, dimension(:,:,:), allocatable  :: dxdxi, dxidx
        real*8, dimension(npro)  :: tmp
c
        integer :: n, nshl, err0,err1
c
        allocate(dxdxi(npro,nsd,nsd),dxidx(npro,nsd,nsd))
c
        call calc_deform_grad(dxdxi,xl,shgl)
c
c.... compute the inverse of deformation gradient
c
       dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3) 
     &                - dxdxi(:,3,2) * dxdxi(:,2,3)
       dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3) 
     &                - dxdxi(:,1,2) * dxdxi(:,3,3)
       dxidx(:,1,3) =  dxdxi(:,1,2) * dxdxi(:,2,3) 
     &                - dxdxi(:,1,3) * dxdxi(:,2,2)
       tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1) 
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)  
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
       dxidx(:,1,1) = dxidx(:,1,1) * tmp
       dxidx(:,1,2) = dxidx(:,1,2) * tmp
       dxidx(:,1,3) = dxidx(:,1,3) * tmp
       dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1) 
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
       dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3) 
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
       dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3) 
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
       dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2) 
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
       dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2) 
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
       dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2) 
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
       if (present(WdetJ) .or. present(qwt))  then
         if (present(WdetJ) .and. present(qwt)) then
           WdetJ = qwt / tmp
         else
           call error ('e3metric_m:  ', 'both WdetJ and qwt should be present', qwt)
         endif
       endif
c
c.... compute the global gradient of shape-functions
c
       nshl = size(shgl,3)
c
       do n = 1, nshl
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) + 
     &                 shgl(:,2,n) * dxidx(:,2,1) +
     &                 shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) + 
     &                 shgl(:,2,n) * dxidx(:,2,2) +
     &                 shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) + 
     &                 shgl(:,2,n) * dxidx(:,2,3) +
     &                 shgl(:,3,n) * dxidx(:,3,3) 
       enddo
c
       deallocate(dxdxi,stat=err0)
       deallocate(dxidx,stat=err1)
c      print*, myrank,npro,nsd,err0,err1
c
       return
c
      end subroutine e3metric      
c
      subroutine calc_deform_grad(dxdxi,xl,shgl)
c
        real*8, dimension(npro,nsd,nsd), intent(out) :: dxdxi
        real*8, dimension(:,:,:), intent(in)  :: xl, shgl
c
        integer :: nenl,n,i,j
c
        nenl = size(xl,2)
c
        dxdxi = zero
c
        do n = 1, nenl
          do i = 1, nsd
            do j = 1, nsd
              dxdxi(:,i,j) = dxdxi(:,i,j) + xl(:,n,i)*shgl(:,j,n)
            enddo
          enddo
        enddo
c
      end subroutine calc_deform_grad
c
      end module e3metric_m
