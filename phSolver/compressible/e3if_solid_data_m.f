      module e3if_solid_data_m
c-------------------------------------------------------------------------------
! This module includes the data(parameters, arrays) which would be used in the
! the calculation in the interface element blocks for solid.
!
! Fuctions included:
! e3if0_malloc_solid : memory allocation of solid data for phase_0
! e3if1_malloc_solid : memory allocation of solid data for phase_1 
! 
! e3if0_mfree_solid : memory deallocation
! e3if1_mfree_solid : memory deallocation
c-------------------------------------------------------------------------------
        implicit none
c
        integer :: iblkif_solid
c
        real*8, dimension(:,:), pointer :: if_d0, if_d1
        real*8, dimension(:), pointer :: if_det_d0, if_det_d1
        real*8, dimension(:), pointer :: if_det_baf0, if_det_baf1
        real*8, dimension(:), pointer :: if_Ja_def0, if_Ja_def1
        real*8, dimension(:), pointer :: bulkMod0, shearMod0
        real*8, dimension(:), pointer :: bulkMod1, shearMod1
        real*8, dimension(:,:), pointer :: dudx_0, dudy_0, dudz_0
        real*8, dimension(:,:), pointer :: dudx_1, dudy_1, dudz_1
! stress -temperature modulus
        real*8, dimension(:), pointer :: stress_T_Mod0, stress_T_Mod1
! temperature gradient
        real*8, dimension(:), pointer :: dtpdx_0, dtpdy_0, dtpdz_0   
        real*8, dimension(:), pointer :: dtpdx_1, dtpdy_1, dtpdz_1        
c
        contains
c............................................................................
c
        subroutine e3if0_malloc_solid
          use global_const_m
          use propar_m
          use solid_data_m,only: b_size
          implicit none
c
          allocate(if_d0(npro,b_size))
          allocate(if_det_d0(npro),if_det_baf0(npro))
          allocate(if_Ja_def0(npro))
          allocate(bulkMod0(npro),shearMod0(npro))
          allocate(dudx_0(npro,nsd))
          allocate(dudy_0(npro,nsd))
          allocate(dudz_0(npro,nsd))
          allocate( stress_T_Mod0(npro) )
          allocate( dtpdx_0(npro), dtpdy_0(npro), dtpdz_0(npro) )
c
        end subroutine e3if0_malloc_solid
c............................................................................
c
        subroutine e3if1_malloc_solid
          use global_const_m
          use propar_m
          use solid_data_m,only: b_size
          implicit none
c
          allocate(if_d1(npro,b_size))
          allocate(if_det_d1(npro),if_det_baf1(npro))
          allocate(if_Ja_def1(npro))
          allocate(bulkMod1(npro),shearMod1(npro))
          allocate(dudx_1(npro,nsd))
          allocate(dudy_1(npro,nsd))
          allocate(dudz_1(npro,nsd))
          allocate( stress_T_Mod1(npro) )
          allocate( dtpdx_1(npro), dtpdy_1(npro), dtpdz_1(npro) )
c
        end subroutine e3if1_malloc_solid

c............................................................................
c
        subroutine e3if0_mfree_solid
c
          implicit none
c
          deallocate(if_d0)
          deallocate(if_det_d0, if_det_baf0)
          deallocate(if_Ja_def0)
          deallocate(bulkMod0,shearMod0)
          deallocate(dudx_0,dudy_0,dudz_0)
          deallocate( stress_T_Mod0 )
          deallocate( dtpdx_0, dtpdy_0, dtpdz_0)         
c
        end subroutine e3if0_mfree_solid
c.............................................................................
c
        subroutine e3if1_mfree_solid
c
          implicit none        
c
          deallocate(if_d1)
          deallocate(if_det_d1, if_det_baf1)
          deallocate(if_Ja_def1)
          deallocate(bulkMod1,shearMod1)
          deallocate(dudx_1,dudy_1,dudz_1)
          deallocate( stress_T_Mod1 )
          deallocate( dtpdx_1, dtpdy_1, dtpdz_1)         
c
        end subroutine e3if1_mfree_solid
c
c
      end module e3if_solid_data_m
