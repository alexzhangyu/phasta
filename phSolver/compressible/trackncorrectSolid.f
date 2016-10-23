      subroutine trackncorrectSolid( iblk, npro_temp, nshl_temp, dt, old_x, y)
c      
      use conpar_m
      use global_const_m
      use pointer_data, only: mien
      use solid_m, only:is_solid
      implicit none
c      
      integer :: npro_temp, nshl_temp
      integer :: iblk, ipro, ishl, ith
      real*8  dt
      real*8  old_x(numnp,nsd)
      real*8  y(nshg,ndof)
c
      allocate( is_solid(nshg) ) 
c.....set the global number of solid node within this block to 1.0
      do ipro = 1, npro_temp !loop over all elements within solid
        do ishl =1, nshl_temp !loop over all local dof
          is_solid( mien(iblk)%p( ipro, ishl)) = 1
        enddo
      enddo
c......adding the displacment from solution field
      do ith =1,nsd
        old_x(:,ith) = old_x(:,ith)+ is_solid(:)
     &*                y( :,ith) * dt
      enddo
      end subroutine trackncorrectSolid 
