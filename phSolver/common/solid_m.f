      module solid_data_m
        use pointer_data
        use outpar_m
        implicit none
c
        type (r3d), dimension(MAXBLK2) ::  b !for solid,left Cauchy_green tensor,added
        type (r3d), dimension(MAXBLK2) ::  b_dot!time derivative of b,added
        type (r3d), dimension(MAXBLK2) ::  b_af! b at time step n+af,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b !left Cauchy_green tensor on the boundary,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b_dot!time derivative of b on the boudary,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b_af! b at time step n+af on the boudary,added
c
        integer, pointer :: is_solid(:)
        integer, parameter :: b_size = 6
        real*8, dimension(:,:),pointer::disp_solid_temp !track the solid displacement field
        real*8, dimension(:),pointer :: elm_b1, elm_b2, elm_b3,
     &                          elm_b4, elm_b5, elm_b6 ! element averaged b array
c
        type solid_t
          logical :: is_active, restart
          integer :: nel,nelb
        end type solid_t
c
        integer :: b_array_size
        real*8, dimension(:), pointer :: temp_b, temp_b_dot, temp_b_af  ! temp arrays to read b, b_dot and b_af
        real*8, dimension(:), pointer :: bdy_temp_b, bdy_temp_b_dot, bdy_temp_b_af  ! on the boundary...
c
        type(solid_t) :: solid_p
c
      end module solid_data_m
c
      module solid_m
c
        use solid_data_m
        implicit none
c
        contains
c
        subroutine init_block(b,bdot,baf,tempb,tempbdot,tempbaf)
          use intpt_m
          use propar_m
          implicit none
          real*8, dimension(:,:,:), pointer, intent(out) :: b,bdot,baf
          real*8, dimension(:), pointer, intent(in) :: tempb,tempbdot,tempbaf
          integer :: ib,ifirst,ilast,incr
          do intp = 1,ngauss
            do ib = 1,b_size
              ifirst = ib + (intp-1)*ngauss
              incr   = ngauss*b_size
              ilast  = ifirst + (npro-1)*incr
              b(1:npro,intp,ib)    = tempb(ifirst:ilast:incr)
              bdot(1:npro,intp,ib) = tempbdot(ifirst:ilast:incr)
              baf(1:npro,intp,ib)  = tempbaf(ifirst:ilast:incr)
            enddo
          enddo
        end subroutine init_block
c
        subroutine dump_block(tempb,tempbdot,tempbaf,b,bdot,baf)
          use intpt_m
          use propar_m
          implicit none
          real*8, dimension(:), pointer, intent(out) :: tempb,tempbdot,tempbaf
          real*8, dimension(:,:,:), pointer, intent(in) :: b,bdot,baf
          integer :: ib,ifirst,ilast,incr
          do intp = 1,ngauss
            do ib = 1,b_size
              ifirst = ib + (intp-1)*ngauss
              incr   = ngauss*b_size
              ilast  = ifirst + (npro-1)*incr
              tempb(ifirst:ilast:incr)    = b(1:npro,intp,ib)
              tempbdot(ifirst:ilast:incr) = bdot(1:npro,intp,ib)
              tempbaf(ifirst:ilast:incr)  = baf(1:npro,intp,ib)
            enddo
          enddo
        end subroutine dump_block
c
        subroutine malloc_solid
          use matdat_def_m
          use number_def_m
          use elmpar_m
          use blkdat_m
          use intpt_m
          implicit none
          integer :: mattype, iblk, npro
c
c.... allocate space for solid arrays 
c
          solid_p%nel = 0
          b_array_size = 0  ! keep this size for later write_field...
c
          blocks_loop: do iblk = 1, nelblk
c
            mattype = lcblk(i_mattype,iblk)
            lcsyst = lcblk(3,iblk)
            ngauss = nint(lcsyst)
            npro = lcblk(1,iblk+1) - lcblk(1,iblk)
c
            if (mat_eos(mattype,1).eq.ieos_solid_1)then
c
              allocate (b(iblk)%p(npro,ngauss,b_size))
              allocate (b_dot(iblk)%p(npro,ngauss,b_size))
              allocate (b_af(iblk)%p(npro,ngauss,b_size))
c
              solid_p%nel = solid_p%nel + npro
              b_array_size = b_array_size + npro*ngauss*b_size
c
              if (solid_p%restart) then
c
                call init_block(b(iblk)%p,b_dot(iblk)%p,b_af(iblk)%p,temp_b,temp_b_dot,temp_b_af)
c
              else
c
                b(iblk)%p(:,:,:)= one
                b(iblk)%p(:,:,4)= zero
                b(iblk)%p(:,:,5)= zero
                b(iblk)%p(:,:,6)= zero
                b_af(iblk)%p(:,:,:) = one
                b_af(iblk)%p(:,:,4) = zero
                b_af(iblk)%p(:,:,5) = zero
                b_af(iblk)%p(:,:,6) = zero
c
                b_dot(iblk)%p(:,:,:) = zero
c
              endif
            else
              cycle   
            endif
c
          enddo blocks_loop
c
          if (solid_p%restart) then
            deallocate(temp_b)
            deallocate(temp_b_dot)
            deallocate(temp_b_af)
          endif
c
          solid_p%nelb = 0
c
          boundary_blocks_loop: do iblk = 1, nelblb
c
            mattype = lcblkb(i_mattype,iblk)
            lcsyst = lcblkb(3,iblk)
            ngaussb = nintb(lcsyst)
            npro = lcblkb(1,iblk+1) - lcblkb(1,iblk)
c
            if (mat_eos(mattype,1).eq.ieos_solid_1)then
c
              allocate (bdy_b(iblk)%p(npro,ngaussb,b_size))
              allocate (bdy_b_dot(iblk)%p(npro,ngaussb,b_size))
              allocate (bdy_b_af(iblk)%p(npro,ngaussb,b_size))
c
              solid_p%nelb = solid_p%nelb + npro
c
              if (solid_p%restart) then
c
                call init_block(bdy_b(iblk)%p,bdy_b_dot(iblk)%p,bdy_b_af(iblk)%p,bdy_temp_b,bdy_temp_b_dot,bdy_temp_b_af)
c
              else
c
                bdy_b(iblk)%p(:,:,:)= one
                bdy_b(iblk)%p(:,:,4)= zero
                bdy_b(iblk)%p(:,:,5)= zero
                bdy_b(iblk)%p(:,:,6)= zero
                bdy_b_af(iblk)%p(:,:,:) = one
                bdy_b_af(iblk)%p(:,:,4) = zero
                bdy_b_af(iblk)%p(:,:,5) = zero
                bdy_b_af(iblk)%p(:,:,6) = zero
c
                bdy_b_dot(iblk)%p(:,:,:) = zero
c
              endif
            else
              cycle
            endif
c..
            if (solid_p%restart) then
              deallocate(bdy_temp_b)
              deallocate(bdy_temp_b_dot)
              deallocate(bdy_temp_b_af)
            endif
c
          enddo boundary_blocks_loop
c
c
        end subroutine malloc_solid
c
        subroutine free_solid
c
        end subroutine free_solid
c
        subroutine read_field(fieldtag,b)
          use iso_c_binding 
          use phio
          implicit none
          character(len=*), intent(in) :: fieldtag
          real*8, pointer, intent(out) :: b(:)
          character(len=1024) :: dataInt, dataDbl
          integer, target :: intfromfile(50) ! integers read from headers
          intfromfile=0
c          call phio_readheader(fhandle,
c     &     c_char_'trim(fieldtag)' // char(0),
c     &     c_loc(intfromfile), 1, dataInt, iotype)
c          if (intfromfile(1) > 0) then
c            allocate(b(intfromfile(1)))
c            call phio_readdatablock(fhandle,
c     &       c_char_'trim(fieldtag)' // char(0),
c     &       c_loc(b), intfromfile(1), dataDbl, iotype)
c          endif
        end subroutine read_field
c
        subroutine read_restart_solid
          use iso_c_binding 
          use phio
          implicit none
          call read_field('solid b',temp_b)
          call read_field('solid b_dot',temp_b_dot)
          call read_field('solid b_af',temp_b_af)
          if (associated(temp_b)) then
            solid_p%restart   = .true.
          endif
        end subroutine read_restart_solid
c
        subroutine write_restart_solid
          use matdat_def_m
          use number_def_m
          use elmpar_m
          use blkdat_m
          use intpt_m
          use workfc_m
          use timdat_m
          implicit none
          integer :: mattype, iblk, npro, b_len
c
          allocate(temp_b(b_array_size))
          allocate(temp_b_dot(b_array_size))
          allocate(temp_b_af(b_array_size))
c
          b_len = 0
          blocks_loop: do iblk = 1, nelblk
c
            mattype = lcblk(i_mattype,iblk)
            lcsyst = lcblk(3,iblk)
            ngauss = nint(lcsyst)
            npro = lcblk(1,iblk+1) - lcblk(1,iblk)
c
            if (mat_eos(mattype,1).eq.ieos_solid_1)then
              call dump_block(temp_b,temp_b_dot,temp_b_af,b(iblk)%p,b_dot(iblk)%p,b_af(iblk)%p)
            endif
c
          enddo blocks_loop
c
          call write_field(myrank,'a'//char(0),'solid b'//char(0),7,
     &     temp_b,'d'//char(0),b_array_size,1,lstep)
c
          deallocate(temp_b)
          deallocate(temp_b_dot)
          deallocate(temp_b_af)
c
        end subroutine write_restart_solid
c
c
       subroutine itrSetupSolid
c
       use inpdat_m
       use timdat_m
       use genpar_m
       use number_def_m
       implicit none
c
c  Setting the generalized alpha method parameters for solid
c
      if( rhoinf_B(itseq).lt.0.or.rhoinf_B(itseq).gt.1) then ! backward Euler
         almBi   = one
         alfBi   = one
         gamBi   = one
         ipred  = 1
      else           !second order family
         almBi   = (three-rhoinf_B(itseq))/(one+rhoinf_B(itseq))/two
         alfBi   = one/(one+rhoinf_B(itseq))
         gamBi   = pt5+almBi-alfBi
      endif
c    
        end subroutine itrSetupSolid
c
c

       subroutine fillelmb( iel,      npro_s,
     &                      lcsyst_s,   ngauss_s,
     &                      blk_b)
c....fill the elm-wise solid array
c       
      use solid_data_m, only: elm_b1, elm_b2, elm_b3,
     &                   elm_b4, elm_b5, elm_b6
      use intpt_m, only: intp, Qwt
      use number_def_m  
      implicit none
c
      real*8 blk_b(npro_s,ngauss_s,6)
      real*8 sum_temp(npro_s)
      integer::npro_s, iel 
      integer::lcsyst_s, ngauss_s
      integer::ilast, ith
c
      sum_temp = zero !initialization
c
      ilast = iel+npro_s-1
c    
c.....loop over all quadrature point
         quad_loop: do intp = 1, ngauss_s
c
                   elm_b1(iel:ilast) = elm_b1(iel:ilast)
     &                                 + blk_b(:,intp,1)* Qwt(lcsyst_s,intp)
                   elm_b2(iel:ilast) = elm_b2(iel:ilast)
     &                                 + blk_b(:,intp,2)* Qwt(lcsyst_s,intp)
                   elm_b3(iel:ilast) = elm_b3(iel:ilast)
     &                                 + blk_b(:,intp,3)* Qwt(lcsyst_s,intp)
                   elm_b4(iel:ilast) = elm_b4(iel:ilast)
     &                                 + blk_b(:,intp,4)* Qwt(lcsyst_s,intp)
                   elm_b5(iel:ilast) = elm_b5(iel:ilast)
     &                                 + blk_b(:,intp,5)* Qwt(lcsyst_s,intp)
                   elm_b6(iel:ilast) = elm_b6(iel:ilast)
     &                                 + blk_b(:,intp,6)* Qwt(lcsyst_s,intp)
                   sum_temp(:) = sum_temp(:) + Qwt(lcsyst_s,intp)
c..
        enddo quad_loop 
c.. Normalize the elm-wise field
        elm_b1(iel:ilast) = elm_b1(iel:ilast)/sum_temp(:)
        elm_b2(iel:ilast) = elm_b2(iel:ilast)/sum_temp(:)
        elm_b3(iel:ilast) = elm_b3(iel:ilast)/sum_temp(:)
        elm_b4(iel:ilast) = elm_b4(iel:ilast)/sum_temp(:)
        elm_b5(iel:ilast) = elm_b5(iel:ilast)/sum_temp(:)
        elm_b6(iel:ilast) = elm_b6(iel:ilast)/sum_temp(:)
c
c
       return
       end subroutine fillelmb
c
      end module solid_m
