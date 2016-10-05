      subroutine itrCorrectElasSolid (old_x, y)
c....Update the displacement within solid at end of each time step
c       
      
c      implicit none
      use pointer_data
      use mattype_m
      include "common.h"
c
      real*8  old_x(numnp,nsd)
      real*8  y(nshg,ndof)
      real*8  temp_count(nshg,nsd)
      integer:: mater_s, npro_s, nshl_s
      integer:: mater_sb
c.....initialize the count array      
      temp_count = zero
c..
c    
c.....for interior blocks
         blocks_loop: do iblk = 1, nelblk
           mater_s = lcblk(7,iblk)
c
!for solid block only
           if (mat_eos(mater_s,1).eq.ieos_solid_1)then
c.....save the temp
              npro_s = SIZE(mien(iblk)%p,1)
              nshl_s = SIZE(mien(iblk)%p,2)
c.....set the global number of solid node within this block to 1.0
              do ipro = 1, npro_s !loop over all elements within solid
                do ishl =1, nshl_s !loop over all local dof
                  temp_count( mien(iblk)%p( ipro, ishl),:) = 1.0
                
                enddo
              enddo

           endif
c..
        enddo blocks_loop 
c
c......adding the displacment from solution field
              do ith =1,nsd
              temp_count(:,ith) = temp_count(:,ith) *
     &                          y( :,ith)
              enddo
             old_x = old_x + temp_count * Delt(1)

      return
      end subroutine itrCorrectElasSolid
