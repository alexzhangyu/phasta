      subroutine trackncorrectSolid(dt, old_x, y, old_um)
c      
      use conpar_m
      use global_const_m
      use solid_m, only:is_solid,disp_solid_temp
      implicit none
c      
      integer :: ith
      real*8  dt
      real*8  old_x(numnp,nsd)
c.....solid debug
      real*8  old_um(nshg,nsd)
c.....solid debug
      real*8  y(nshg,ndof)
c
c......adding the displacment from solution field
      do ith =1,nsd
c......solid debug
c        old_um(:,ith) = is_solid(:) * y( :,ith)
c......solid debug
c        old_x(:,ith) = old_x(:,ith)+ is_solid(:)
c     &*                y( :,ith) * dt
        if (iSOLID == 1)then
          disp_solid_temp(:,ith) =  disp_solid_temp(:,ith)+ is_solid(:)
     &*                y( :,ith) * dt
        endif
      enddo
c
      end subroutine trackncorrectSolid 
