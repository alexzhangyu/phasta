       subroutine allocate_initial_solid(i_iniSolid)
c...   calculate the left Cauchy-green tensor at time step n+af
C       implicit none 
       use pointer_data
c       use mattype_m
       include "common.h"
c
       integer:: i_iniSolid
       integer:: mattyp_s, lcsyst_s, npro_s, ngauss_s
       integer:: mattyp_sb, lcsyst_sb, npro_sb, ngauss_sb
c..................
c.... allocate space for solid arrays 
 
       blocks_loop: do iblk = 1, nelblk
         mattyp_s = lcblk(7,iblk)
         lcsyst_s = lcblk(3,iblk)
         ngauss_s = nint(lcsyst_s)
         npro_s = SIZE(mien(iblk)%p,1)
!for solid block only
         if (mat_eos(mattyp_s,1).eq.ieos_solid_1)then
           allocate (b(iblk)%p(npro_s,ngauss_s,6))!for solid
           allocate (b_dot(iblk)%p(npro_s,ngauss_s,6)) !for solid,added
           allocate (b_af(iblk)%p(npro_s,ngauss_s,6)) !for solid,added
c......... these arrays need initialization
           if (i_iniSolid .eq. 1)then
              b(iblk)%p    = zero 
              b_dot(iblk)%p = zero
              b(iblk)%p(:,:,1)= one
              b(iblk)%p(:,:,2)= one
              b(iblk)%p(:,:,3)= one
              b_af(iblk)%p = b(iblk)%p
           endif
c
          endif
c..
       enddo blocks_loop
c
c
c
c 
      boundary_ blocks_loop: do iblk = 1, nelblb
         mattyp_sb = lcblkb(7,iblk)
         lcsyst_sb = lcblkb(3,iblk)
         ngauss_sb = nintb(lcsyst_sb)
         npro_sb = SIZE(mienb(iblk)%p,1)
!for solid block only
         if (mat_eos(mattyp_sb,1).eq.ieos_solid_1)then
           allocate (bdy_b(iblk)%p(npro_sb,ngauss_sb,6))!for solid
           allocate (bdy_b_dot(iblk)%p(npro_sb,ngauss_sb,6)) !for solid,added
           allocate (bdy_b_af(iblk)%p(npro_sb,ngauss_sb,6)) !for solid,added
c......... these arrays need initialization
           if (i_iniSolid .eq. 1)then
             bdy_b(iblk)%p    = zero
             bdy_b_dot(iblk)%p = zero
             bdy_b(iblk)%p(:,:,1)= one
             bdy_b(iblk)%p(:,:,2)= one
             bdy_b(iblk)%p(:,:,3)= one
             bdy_b_af(iblk)%p = bdy_b(iblk)%p
           endif
c
          endif
c..
      enddo boundary_blocks_loop

      end
