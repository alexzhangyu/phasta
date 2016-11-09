      subroutine e3ivar_solid(g1yi_,g2yi_,g3yi_,almBi_,alfBi_,gamBi_)
        use e3_solid_func_m
        implicit none
        real*8, dimension(npro,nflow),target, intent(in) :: g1yi_,g2yi_,g3yi_
        real*8, intent(in) :: almBi_, alfBi_, gamBi_
c
c        almBi = almBi_
c        alfBi = alfBi_
c        gamBi = gamBi_  
cc
c        allocate(dudx(npro,3))
c        allocate(dudy(npro,3))
c        allocate(dudz(npro,3))
c
c        dudx = g1yi_(:,2:4)
c        dudy = g2yi_(:,2:4)
c        dudz = g3yi_(:,2:4)
c
c        call calc_solid
cc       
c        deallocate(dudx)
c        deallocate(dudy)
c        deallocate(dudz)

      end subroutine e3ivar_solid
c
      subroutine e3bvar_solid(g1yb_,g2yb_,g3yb_,almBi_,alfBi_,gamBi_)
        use e3_solid_func_m
        implicit none
        real*8, dimension(npro,nflow),target, intent(in) :: g1yb_,g2yb_,g3yb_
        real*8, intent(in) :: almBi_, alfBi_, gamBi_
c
c        almBi = almBi_
c        alfBi = alfBi_
c        gamBi = gamBi_  
cc
cc        dudx => g1yi_(:,2:4)
cc        dudy => g2yi_(:,2:4)
cc        dudz => g3yi_(:,2:4)
c        allocate(dudx(npro,3))
c        allocate(dudy(npro,3))
c        allocate(dudz(npro,3))
cc
c        dudx = g1yb_(:,2:4)
c        dudy = g2yb_(:,2:4)
c        dudz = g3yb_(:,2:4)
cc
c        call calc_solid_bdy
cc
c        deallocate(dudx)
c        deallocate(dudy)
c        deallocate(dudz)
c
      end subroutine e3bvar_solid
