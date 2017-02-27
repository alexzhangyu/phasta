      subroutine update_solid_blocks(old_x, y, old_um)
c
        use conpar_m
        use timdat_m
        use solid_m
        use matdat_def_m
        use number_def_m
        use blkdat_m
        use elmpar_m, only: nelblk, nelblb, nelblif
        use pointer_data, only: mien
        use inpdat_m
        use intpt_m, only: nint
c
        implicit none
c
       integer :: iblk,ipro,ishl
       integer :: mater_s
       integer :: mater_sb
       integer :: npro_temp, npro_b_temp
       integer :: nshl_temp
       integer :: iel_temp, lcsyst_temp, ngauss_temp
       real*8 :: dt
! local array for interface elements operations
       integer :: mater_if0, mater_if1
       integer :: npro_if_temp      
c
       real*8  old_x(numnp,nsd)
       real*8  y(nshg,ndof)
c..... solid debug
       real*8  old_um(nshg,nsd)
c..... solid debug
c
c......passing or initialize parameters
       dt = Delt(1)
       allocate( is_solid(nshg))
       is_solid = zero
c
       allocate(elm_b1(numel) )  
       allocate(elm_b2(numel) )  
       allocate(elm_b3(numel) )  
       allocate(elm_b4(numel) )  
       allocate(elm_b5(numel) )  
       allocate(elm_b6(numel) )
       elm_b1 = zero  
       elm_b2 = zero
       elm_b3 = zero
       elm_b4 = zero
       elm_b5 = zero
       elm_b6 = zero
      
c-----------------------for interior blocks-----------------------------------
         do iblk = 1, nelblk
           mater_s = lcblk(7,iblk)
!for solid block only
           if (mat_eos(mater_s,1).eq.ieos_solid_1)then
c.....save the temp
              npro_temp = SIZE(b(iblk)%p,1)
              nshl_temp = SIZE(mien(iblk)%p,2)
              iel_temp = lcblk(1,iblk)
              lcsyst_temp = lcblk(3,iblk)
              ngauss_temp = nint(lcsyst_temp)
c
c.....Do the update
              b_dot(iblk)%p(:,:,:)  =  -one/(alfBi * gamBi * dt)* b(iblk)%p(:,:,:)
     &+                           (one - one/gamBi) * b_dot(iblk)%p(:,:,:)
     &+                           one/(alfBi * gamBi * dt) * b_af(iblk)%p(:,:,:)
              b(iblk)%p(:,:,:) = (one/alfBi) * (b_af(iblk)%p(:,:,:) - b(iblk)%p(:,:,:) )
     &+                          b(iblk)%p(:,:,:)
c
c.....set the global number of solid node within this block to 1.0     
                                      
              do ipro = 1, npro_temp !loop over all elements within solid      
                do ishl =1, nshl_temp !loop over all local dof                 
                  is_solid( mien(iblk)%p( ipro, ishl)) = 1             
                enddo
              enddo                                               
c
c.....fill the elm-wise solid arrays

              call fillelmb( iel_temp,     npro_temp,
     &                       lcsyst_temp,  ngauss_temp,
     &                       b(iblk)%p )
c.....end of fill elm-wise solid arrays

c.....end of updates
           endif
c..
        enddo
c....Track the displacment for solid
c       call trackncorrectSolid( dt, old_x, y)
c........solid debug
       call trackncorrectSolid( dt, old_x, y, old_um)
c........solid debug
c
        deallocate(is_solid) 
c
c--------------------for boundary blocks---------------------------------------
         do iblk = 1, nelblb
           mater_sb = lcblkb(7,iblk)
c
!for solid block only
           if (mat_eos(mater_sb,1).eq.ieos_solid_1)then
c.....save the temp
              npro_b_temp = SIZE(bdy_b(iblk)%p,1)
c.....Do the update
              bdy_b_dot(iblk)%p(:,:,:)  = -one/(alfBi * gamBi * dt)
     &                                  * bdy_b(iblk)%p(:,:,:)
     &                                  + (one - one/gamBi) 
     &                                  * bdy_b_dot(iblk)%p(:,:,:)
     &                                  + one/(alfBi * gamBi * dt) 
     &                                  * bdy_b_af(iblk)%p(:,:,:)

              bdy_b(iblk)%p(:,:,:) = (one/alfBi) 
     &                             * ( bdy_b_af(iblk)%p(:,:,:) 
     &                             - bdy_b(iblk)%p(:,:,:) )
     &                             + bdy_b(iblk)%p(:,:,:)

c.....end of updates
            endif
        enddo  
c
c-----------------------for interface blocks-----------------------------------
c.....for interface blocks
         do iblk = 1, nelblif
           mater_if0 = lcblkif(9, iblk)
           mater_if1 = lcblkif(10,iblk)
           npro_if_temp = lcblkif(1,iblk+1) - lcblkif(1,iblk)
c
! for phase_0
           if (mat_eos(mater_if0,1).eq.ieos_solid_1)then
c.....Do the update
              if_b0_dot(iblk)%p(:,:,:)  = -one/(alfBi * gamBi * dt)
     &                                  * if_b0(iblk)%p(:,:,:)
     &                                  + (one - one/gamBi)
     &                                  * if_b0_dot(iblk)%p(:,:,:)
     &                                  + one/(alfBi * gamBi * dt) 
     &                                  * if_b0_af(iblk)%p(:,:,:)

              if_b0(iblk)%p(:,:,:) = (one/alfBi) 
     &                             * ( if_b0_af(iblk)%p(:,:,:)
     &                             - if_b0(iblk)%p(:,:,:) )
     &                             + if_b0(iblk)%p(:,:,:)

c.....end of updates
           endif
! for phase_1           
           if (mat_eos(mater_if1,1).eq.ieos_solid_1)then
c.....Do the update
              if_b1_dot(iblk)%p(:,:,:)  = -one/(alfBi * gamBi * dt)
     &                                  * if_b1(iblk)%p(:,:,:)
     &                                  + (one - one/gamBi)
     &                                  * if_b1_dot(iblk)%p(:,:,:)
     &                                  + one/(alfBi * gamBi * dt) 
     &                                  * if_b1_af(iblk)%p(:,:,:)

              if_b1(iblk)%p(:,:,:) = (one/alfBi) 
     &                             * ( if_b1_af(iblk)%p(:,:,:)
     &                             - if_b1(iblk)%p(:,:,:) )
     &                             + if_b1(iblk)%p(:,:,:)

c.....end of updates
            endif
c            
        enddo
c--------------------------------------
c          
      return
      end subroutine update_solid_blocks
