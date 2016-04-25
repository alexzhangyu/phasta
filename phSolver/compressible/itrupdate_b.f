      subroutine itrupdate_b
c....Update the left Cauchy-green tensor at the end of each time step
c       
      
c      implicit none
      use pointer_data
      use mattype_m
      include "common.h"
c
      real*8, dimension(npro,ngauss,6) :: b_temp,b_dot_temp
c    
      material_loop: do imattype = 1, nummat
c
         blocks_loop: do iblk = 1, nelblk
c
!for solid block only
           if (mat_eos(imattype,1).eq.ieos_solid_1)then
c.....save the temp
              b_temp(:,:,:) = b(iblk)%p(:,:,:)
              b_dot_temp(:,:,:) = b_dot(iblk)%p(:,:,:)
c
c.....Do the update
              b(iblk)%p(:,:,:) = (one/alfBi) * (b_af(iblk)%p(:,:,:) - b_temp(:,:,:) )
     &+                          b_temp(:,:,:)
              b_dot(iblk)%p(:,:,:)  =  -one/(alfBi * gamBi * Delt(1))* b_temp(:,:,:)
     &+                           (one - one/gamBi) * b_dot_temp(:,:,:)
     &+                           one/(alfBi * gamBi * Delt(1)) * b_af(iblk)%p(:,:,:)
           endif
c..
        enddo blocks_loop
      enddo material_loop 
c
      return
      end subroutine itrupdate_b
