      subroutine update_solid_blocks(dt, old_x, y)
c
        use conpar_m
        use timdat_m
        use solid_m
        use matdat_def_m
        use number_def_m
        use blkdat_m
        use elmpar_m, only: nelblk, nelblb
        use pointer_data, only: mien
c
        implicit none
c
       integer :: iblk
       integer :: mater_s
       integer :: mater_sb
       integer :: npro_temp, npro_b_temp
       integer :: nshl_temp
       real*8 :: dt
c
       real*8  old_x(numnp,nsd)
       real*8  y(nshg,ndof)
c
c.....for interior blocks
         do iblk = 1, nelblk
           mater_s = lcblk(7,iblk)
!for solid block only
           if (mat_eos(mater_s,1).eq.ieos_solid_1)then
c.....save the temp
              npro_temp = SIZE(b(iblk)%p,1)
              nshl_temp = SIZE(mien(iblk)%p,2)
c              iel_temp = lcblk(1,iblk)
c              lcsyst_temp = lcblk(3,iblk)
c              ngauss_temp = nint(lcsyst_temp)
c
c.....Do the update
              b_dot(iblk)%p(:,:,:)  =  -one/(alfBi * gamBi * dt)* b(iblk)%p(:,:,:)
     &+                           (one - one/gamBi) * b_dot(iblk)%p(:,:,:)
     &+                           one/(alfBi * gamBi * dt) * b_af(iblk)%p(:,:,:)
              b(iblk)%p(:,:,:) = (one/alfBi) * (b_af(iblk)%p(:,:,:) - b(iblk)%p(:,:,:) )
     &+                          b(iblk)%p(:,:,:)  
c....Track the displacment for solid
              call trackncorrectSolid( iblk, npro_temp, nshl_temp, dt, old_x, y)
c
cc.....fill the elm-wise solid arrays

c              call fillelmb( elmb1, elmb2, iel_temp,    npro_temp,
c     &                       lcsyst_temp,  ngauss_temp, b(iblk)%p )
cc.....end of fill elm-wise solid arrays
c
c.....end of updates
           endif
c..
        enddo 

c.....for boundary blocks
         do iblk = 1, nelblb
           mater_sb = lcblkb(7,iblk)
c
!for solid block only
           if (mat_eos(mater_sb,1).eq.ieos_solid_1)then
c.....save the temp
              npro_b_temp = SIZE(bdy_b(iblk)%p,1)
c.....Do the update
              bdy_b_dot(iblk)%p(:,:,:)  =  -one/(alfBi * gamBi * dt)*bdy_b(iblk)%p(:,:,:)
     &+                           (one - one/gamBi) * bdy_b_dot(iblk)%p(:,:,:)
     &+                           one/(alfBi * gamBi * dt) * bdy_b_af(iblk)%p(:,:,:)

              bdy_b(iblk)%p(:,:,:) = (one/alfBi) * (bdy_b_af(iblk)%p(:,:,:) - bdy_b(iblk)%p(:,:,:) )
     &+                              bdy_b(iblk)%p(:,:,:)

c.....end of updates
            endif
        enddo  
c
      return
      end subroutine update_solid_blocks
