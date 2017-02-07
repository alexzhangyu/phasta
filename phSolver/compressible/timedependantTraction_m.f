      module timedependantTraction_m
c.....used to test the time dependant traction condition for solid
c
        implicit none
c
c
        contains
c
        subroutine settimedepandantflux(x)
c.....apply the time dependant flux BC for solid, test only
c
        use number_def_m
        use pointer_data
        use propar_m
        use elmpar_m
        use shpdat_m
        use conpar_m
        use global_const_m
        use inpdat_m
        use timdat_m
        use blkdat_m
        implicit none

c
        real*8, dimension(numnp,nsd),intent(in)::x 
        integer, dimension(27) :: lnode
        integer :: n,ipro,iblk,temp,nodlc,temp1,npro_temp
        real*8, parameter :: r = 1.0d-5, tt = 2.0d3,f =1.264d3,
     &                       tp =1.0d3
c
        call getbnodes(lnode)
c
        do iblk = 1, nelblb
          npro_temp = lcblkb(1,iblk+1) -lcblkb(1,iblk)
          do ipro = 1, npro_temp
            do n = 1, nshlb
               nodlc = lnode(n)
               temp = mienb(iblk)%p(ipro,nodlc)
               temp1 = miBCB(iblk)%p(ipro,1)
               if ( (abs(x(temp,1)) .lt. r) .and.
     &              (btest(temp1,2)) .and.
     &              (btest(temp1,1))     )then
               mBCB(iblk)%p(ipro,n,2)= 1.0d5 
     &                               + tp*dsin(two*pi*f*lstep * Delt(1))
               mBCB(iblk)%p(ipro,n,3)= tt*dsin(two*pi*f*lstep * Delt(1))
               endif
            enddo
          enddo
        enddo
c
        end subroutine settimedepandantflux
c
      end module timedependantTraction_m
