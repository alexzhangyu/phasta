      subroutine itrupdate_b(elmb1,elmb2)
c....Update the left Cauchy-green tensor at the end of each time step
c       
      
c      implicit none
      use pointer_data
      use mattype_m
      include "common.h"
c
      real*8, dimension(:,:,:), allocatable :: b_temp,b_dot_temp
      real*8, dimension(:,:,:), allocatable :: bdy_b_temp, bdy_b_dot_temp
c      real*8, dimension(npro,ngaussb,6) :: bdy_b_temp,bdy_b_dot_temp
      integer:: mater_s
      integer:: mater_sb
      integer:: npro_temp, npro_b_temp
c
      real*8 elmb1(numel,3) !store elm-wise solid array
      real*8 elmb2(numel,3) !store elm-wise solid array
      integer::iel_temp, lcsyst_temp, ngauss_temp

c    
c.....for interior blocks
         blocks_loop: do iblk = 1, nelblk
           mater_s = lcblk(7,iblk)
c
!for solid block only
           if (mat_eos(mater_s,1).eq.ieos_solid_1)then
c.....save the temp
              npro_temp = SIZE(b(iblk)%p,1)
              iel_temp = lcblk(1,iblk)
              lcsyst_temp = lcblk(3,iblk)
              ngauss_temp = nint(lcsyst_temp)
c
              allocate(b_temp(npro_temp,ngauss,6))
              allocate(b_dot_temp(npro_temp,ngauss,6))
              b_temp(:,:,:) = b(iblk)%p(:,:,:)
              b_dot_temp(:,:,:) = b_dot(iblk)%p(:,:,:)
c
c.....Do the update
              b(iblk)%p(:,:,:) = (one/alfBi) * (b_af(iblk)%p(:,:,:) - b_temp(:,:,:) )
     &+                          b_temp(:,:,:)
              b_dot(iblk)%p(:,:,:)  =  -one/(alfBi * gamBi * Delt(1))* b_temp(:,:,:)
     &+                           (one - one/gamBi) * b_dot_temp(:,:,:)
     &+                           one/(alfBi * gamBi * Delt(1)) * b_af(iblk)%p(:,:,:)
              deallocate(b_temp,b_dot_temp)
c
c.....fill the elm-wise solid arrays
              call fillelmb( elmb1, elmb2, iel_temp,    npro_temp,
     &                       lcsyst_temp,  ngauss_temp, b(iblk)%p )
c.....end of fill elm-wise solid arrays
c
c.....end of updates
           endif
c..
        enddo blocks_loop 
c
c.....for boundary blocks
         boundary_blocks_loop: do iblk = 1, nelblb
           mater_sb = lcblkb(7,iblk)
c
!for solid block only
           if (mat_eos(mater_sb,1).eq.ieos_solid_1)then
c.....save the temp
              npro_b_temp = SIZE(bdy_b(iblk)%p,1)
              allocate(bdy_b_temp(npro_b_temp,ngaussb,6))
              allocate(bdy_b_dot_temp(npro_b_temp,ngaussb,6)) 
              bdy_b_temp(:,:,:) = bdy_b(iblk)%p(:,:,:)
              bdy_b_dot_temp(:,:,:) = bdy_b_dot(iblk)%p(:,:,:)
c
c.....Do the update
              bdy_b(iblk)%p(:,:,:) = (one/alfBi) * (bdy_b_af(iblk)%p(:,:,:) - bdy_b_temp(:,:,:) )
     &+                          bdy_b_temp(:,:,:)
              bdy_b_dot(iblk)%p(:,:,:)  =  -one/(alfBi * gamBi * Delt(1))* bdy_b_temp(:,:,:)
     &+                           (one - one/gamBi) * bdy_b_dot_temp(:,:,:)
     &+                           one/(alfBi * gamBi * Delt(1)) * bdy_b_af(iblk)%p(:,:,:)
              deallocate(bdy_b_temp, bdy_b_dot_temp)
c.....end of updates
            endif
c..
        enddo boundary_blocks_loop 
c
      return
      end subroutine itrupdate_b
