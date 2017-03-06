      module e3_solid_m
c
        use e3_param_m
        use solid_data_m
c
        implicit none
c
        integer :: iblk_solid
        integer :: iblk_first
        integer :: iblk_last
c
        real*8, dimension(:,:), pointer :: d
        real*8, dimension(:), pointer :: det_d, det_baf
        real*8, dimension(:), pointer :: bulkMod, shearMod, Ja_def
        real*8, dimension(:), pointer :: stress_T_Mod
        real*8, dimension(:,:), pointer :: dudx, dudy, dudz
! temperature gradient
        real*8, dimension(:), pointer :: dtpdx, dtpdy, dtpdz        
c
        contains
c
        subroutine e3_malloc_solid
          use global_const_m
c
          call e3_malloc
c
          allocate(d(npro,b_size))
          allocate(det_d(npro),det_baf(npro))
          allocate(Ja_def(npro),bulkMod(npro),shearMod(npro))
          allocate( stress_T_Mod(npro) )
          allocate(dudx(npro,nsd))
          allocate(dudy(npro,nsd))
          allocate(dudz(npro,nsd))
          allocate( dtpdx(npro), dtpdy(npro), dtpdz(npro) )
c
        end subroutine e3_malloc_solid
c
        subroutine e3_mfree_solid
c
          call e3_mfree
c
          deallocate(d)
          deallocate(det_d,det_baf)
          deallocate(Ja_def,bulkMod,shearMod)
          deallocate(stress_T_Mod)
          deallocate(dudx,dudy,dudz)
          deallocate( dtpdx, dtpdy, dtpdz)
c
        end subroutine e3_mfree_solid
c
      end module e3_solid_m
c
      module e3_solid_func_m
c
        use e3_solid_m
        use number_def_m
        use global_const_m
        use intpt_m
        use inpdat_m
        use conpar_m
        use timdat_m
c
        implicit none
c
        real*8, dimension(:,:,:), pointer :: AS
c
      contains
c
      subroutine calc_solid
c
        implicit none
c
        allocate(AS(npro,6,6))
c
        call calc_as_matrix
        d(:,:) = almBi * b(iblk_solid)%p(:,intp,:)
     &+      alfBi * Delt(1) * (almBi - gamBi) 
     &*      b_dot(iblk_solid)%p(:,intp,:)
        call setB_af
        call get_det(b_af(iblk_solid)%p(:,intp,:),det_baf,npro)
        Ja_def= (det_baf)**0.5
        call get_det(d,det_d,npro)
c
        deallocate(AS)
c
      end subroutine calc_solid
c
      subroutine calc_solid_bdy
c
        implicit none
c
        allocate(AS(npro,6,6))
c
        call calc_as_matrix
c.......debug
c        AS(:,5,3) = zero
c        AS(:,6,2) = zero
c.......debug
        d(:,:) = almBi * bdy_b(iblk_solid)%p(:,intp,:)
     &+      alfBi * Delt(1) * (almBi - gamBi) 
     &*      bdy_b_dot(iblk_solid)%p(:,intp,:)
        call setB_af_bdy
        call get_det(bdy_b_af(iblk_solid)%p(:,intp,:),det_baf,npro)
        Ja_def= (det_baf)**0.5
        call get_det(d,det_d,npro)
c
        deallocate(AS)
c
      end subroutine calc_solid_bdy
c
      subroutine calc_as_matrix
c
        implicit none
c
c initialize the AS matrix
         AS = zero !add
c
         AS(:,1,1) = 2 * dudx(:,1)
         AS(:,1,5) = 2 * dudz(:,1)
         AS(:,1,6) = 2 * dudy(:,1)
c
         AS(:,2,2) = 2 * dudy(:,2)
         AS(:,2,4) = 2 * dudz(:,2)
         AS(:,2,6) = 2 * dudx(:,2)
c
         AS(:,3,3) = 2 * dudz(:,3)
         AS(:,3,4) = 2 * dudy(:,3)
         AS(:,3,5) = 2 * dudx(:,3)
c
         AS(:,4,2) = dudy(:,3)
         AS(:,4,3) = dudz(:,2)
         AS(:,4,4) = dudy(:,2) + dudz(:,3)
         AS(:,4,5) = dudx(:,2)
         AS(:,4,6) = dudx(:,3)
c
         AS(:,5,1) = dudx(:,3)
         AS(:,5,3) = dudz(:,1)
         AS(:,5,4) = dudy(:,1)
         AS(:,5,5) = dudx(:,1) + dudz(:,3)
         AS(:,5,6) = dudy(:,3)
c         
         AS(:,6,1) = dudx(:,2)
         AS(:,6,2) = dudy(:,1)
         AS(:,6,4) = dudz(:,1)
         AS(:,6,5) = dudz(:,2)
         AS(:,6,6) = dudx(:,1) + dudy(:,2)
c
c
      end subroutine calc_as_matrix
c
      subroutine setB_af
c
c... calculate the left Cauchy-green tensor at time step n+af
c
c         use function_and_complex
         implicit none 
c
         integer,parameter :: nsize = 6 
c
         real*8, dimension(6,6) :: ident
         real*8, dimension(6,6) :: temp_matrix
         real*8, dimension(6) :: temp_baf
         integer :: i
c
         ident = zero
         do i = 1, nsize
           ident(i,i) = one
         enddo
c
         do i = 1, npro
          temp_matrix(:,:) = (one/almBi) * ident(:,:) + ( gamBi * Delt(1)
     &                     *alfBi/(almBi)**2 ) * AS(i,:,:) !check here
          b_af(iblk_solid)%p(i,intp,:) = matmul(temp_matrix(:,:) , d(i,:))
c           temp_matrix(:,:) = (almBi) * ident(:,:) - ( gamBi * Delt(1)
c     &                     *alfBi) * AS(i,:,:) !check here
c
c           call lu_solve(temp_matrix, d(i,:), b_af(iblk_solid)%p(i,intp,:), 6)
         enddo
c
       end subroutine setB_af 
c
      subroutine setB_af_bdy
c
c... calculate the left Cauchy-green tensor at time step n+af
c
c         use function_and_complex
         implicit none 
c
         integer,parameter :: nsize = 6 
c
         real*8, dimension(6,6) :: ident
         real*8, dimension(6,6) :: temp_matrix
         integer :: i

c
         ident = zero
         do i = 1, nsize
           ident(i,i) = one
         enddo
c
         do i = 1, npro
          temp_matrix(:,:) = (one/almBi) * ident(:,:) + ( gamBi * Delt(1)
     &                     *alfBi/(almBi)**2 ) * AS(i,:,:) !check here
          bdy_b_af(iblk_solid)%p(i,intp,:) = matmul(temp_matrix(:,:) , d(i,:))
c           temp_matrix(:,:) = (almBi) * ident(:,:) - ( gamBi * Delt(1)
c     &                     *alfBi) * AS(i,:,:) !check here
c           call lu_solve(temp_matrix, d(i,:), bdy_b_af(iblk_solid)%p(i,intp,:), 6)
         enddo
c
       end subroutine setB_af_bdy
c
       subroutine get_det(matr,det_matr,npro)
c
        implicit none 
c
        real*8, dimension(npro,6),intent(in) :: matr
        real*8, dimension(npro) :: det_matr
        integer, intent(in) :: npro
c
        det_matr(:) = matr(:,1) * matr(:,2) * matr(:,3)
     &+             matr(:,6) * matr(:,4) * matr(:,5)
     &+             matr(:,5) * matr(:,6) * matr(:,4) 
     &-             matr(:,1) * matr(:,4) * matr(:,4)
     &-             matr(:,5) * matr(:,2) * matr(:,5)
     &-             matr(:,6) * matr(:,6) * matr(:,3)
c
      end subroutine get_det
c
c
      subroutine modify_solid_pres(mater_s)
c
        use matdat_def_m
        use e3_param_m, only:pres
        implicit none
c
        integer, intent(in) :: mater_s
        real*8 :: pres_ref_s
c
        pres_ref_s   = mat_prop(mater_s,iprop_solid_1_p_ref,  1)
        pres(:) = pres(:) - pres_ref_s
      end subroutine modify_solid_pres
c
c
      subroutine set_solid_kij_interior ( stiff, u1, u2, u3,
     &                                     con)
c-----------------------------------------------------------------
! setting the K_ij matrix for solid interior blocks so that
! F^{diff} = K_ij Y,{j}
! input:
! 
! output:
! stiff  (npro,nsd*nflow,nsd*nflow) : stiffness matrix,K_ij
c-----------------------------------------------------------------
c
        implicit none
c
        real*8, dimension(npro,nsd*nflow,nsd*nflow),intent(inout)::stiff
        real*8, dimension(npro),intent(in) :: u1, u2, u3, con
! local arrays
        real*8, dimension(npro) :: d_temp1,d_temp2,d_temp3
        real*8, dimension(npro,6) :: bq_af
c.................................................................
        d_temp1(:) = 2.0/3.0*d(:,1) - 1.0/3.0*d(:,2) - 1.0/3.0*d(:,3)
        d_temp2(:) = -1.0/3.0*d(:,1) + 2.0/3.0*d(:,2) - 1.0/3.0*d(:,3)
        d_temp3(:) = -1.0/3.0*d(:,1) - 1.0/3.0*d(:,2) + 2.0/3.0*d(:,3)
c
c.... K11
c
         stiff(:, 2, 2) = ShearMod * (det_d)**(-5.0/6.0) 
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d(:,1)
     &                  - (5.0/3.0)* (almBi)**(-2.0) * d_temp1 )  
         stiff(:, 2, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (-2.0/3.0)*d(:,6)
         stiff(:, 2, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (-2.0/3.0)*d(:,5)
         stiff(:, 3, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d(:,6)- (5.0/3.0)* (almBi)**(-2) * d(:,6))
         stiff(:, 3, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,1) )
         stiff(:, 4, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d(:,5)- (5.0/3.0)* (almBi)**(-2) * d(:,5))
         stiff(:, 4, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,1))
         stiff(:, 5, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d(:,1) * u1 + d(:,6) * u2 + d(:,5) * u3
     &                  - (5.0/3.0)* (almBi)**(-2)
     &                  *( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3) )
         stiff(:, 5, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d(:,6) * u1 + d(:,1) * u2 )
         stiff(:, 5, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d(:,5) * u1 + d(:,1) * u3 )
         stiff(:, 5, 5) = con ! notice the + or -
c
c.... K12                 
c     
         stiff(:, 2, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d(:,6) )
         stiff(:, 2, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d(:,2) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_temp1 )
         stiff(:, 2, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d(:,4) )
         stiff(:, 3, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d(:,2))
         stiff(:, 3, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,6) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d(:,6) )
         stiff(:, 4, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,4) )
         stiff(:, 4, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( 0.0 
     &                  - (5.0/3.0)* (almBi)**(-2) * d(:,5) )    
         stiff(:, 4, 9) = ShearMod * (det_d)**(-5.0/6.0) 
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,6) )
         stiff(:, 5, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d(:,6) * u1 + d(:,2) * u2 + d(:,4) * u3 )
         stiff(:, 5, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d(:,2) * u1 + d(:,6) * u2 
     &                  - (5.0/3.0)* (almBi)**(-2)
     &                  * ( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3) )
         stiff(:, 5, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d(:,4) * u1 + d(:,6) * u3 )
c
c.... K13
c     
         stiff(:, 2,12) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0 /3.0)*d(:,5) )
         stiff(:, 2,13) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0 /3.0)*d(:,4) )
         stiff(:, 2,14) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0 /3.0)*d(:,3) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_temp1)
         stiff(:, 3,12) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d(:,4))
         stiff(:, 3,13) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,5) )
         stiff(:, 3,14) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( - (5.0/3.0)* (almBi)**(-2) * d(:,6) )
         stiff(:, 4,12) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,3) )
         stiff(:, 4,14) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,5)- (5.0/3.0)* (almBi)**(-2) * d(:,5) )
         stiff(:, 5,12) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0 /3.0)*d(:,5) * u1 + d(:,4) * u2 + d(:,3) * u3 )
         stiff(:, 5,13) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0 /3.0)*d(:,4) * u1 + d(:,5) * u2  )
         stiff(:, 5,14) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0 /3.0)*d(:,3) * u1 + d(:,5) * u3
     &                  - (5.0/3.0)* (almBi)**(-2)
     &                  * ( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3))
c     
c.... K21
c     
         stiff(:, 7, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d(:,6)- (5.0/3.0)* (almBi)**(-2) * d(:,6))
         stiff(:, 7, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,1) )
         stiff(:, 8, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d(:,1) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_temp2 )  
         stiff(:, 8, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (4.0/3.0)*d(:,6)
         stiff(:, 8, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (-2.0/3.0)*d(:,5)
         stiff(:, 9, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (- (5.0/3.0)* (almBi)**(-2) * d(:,4))
         stiff(:, 9, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  *  d(:,5)
         stiff(:, 9, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * d(:,6)
         stiff(:, 10, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,6) * u1 + (-2.0/3.0)*d(:,1) * u2 
     &                   - (5.0/3.0)* (almBi)**(-2)
     &                   * ( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3) )
         stiff(:, 10, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * (d(:,1) * u1 
     &                   + (4.0/3.0)*d(:,6) * u2 +d(:,5) * u3)
         stiff(:, 10, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d(:,5) * u2 +d(:,6) * u3 )
c
c.... K22
c     
         stiff(:, 7, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,2) )
         stiff(:, 7, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,6) -(5.0/3.0)* (almBi)**(-2) * d(:,6) )
         stiff(:, 8, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( -(2.0/3.0)*d(:,6) )
         stiff(:, 8, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d(:,2)
     &                  -(5.0/3.0)* (almBi)**(-2) * d_temp2 )
         stiff(:, 8, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d(:,4) )
         stiff(:, 9, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
         stiff(:, 9, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d(:,2) )
         stiff(:, 10, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,2) * u1 + (-2.0/3.0)*d(:,6) * u2 )
         stiff(:, 10, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * (d(:,6) * u1 + (4.0/3.0)*d(:,2) * u2 +d(:,4) * u3
     &                   - (5.0/3.0)* (almBi)**(-2)
     &                   *( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3))
         stiff(:, 10, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d(:,4) * u2 +d(:,2) * u3 )
         stiff(:, 10, 10) = con ! notice the + or -
c
c.... K23
c     
         stiff(:, 7, 12) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,4) )
         stiff(:, 7, 13) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,5) )
         stiff(:, 7, 14) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( 0.0 -(5.0/3.0)* (almBi)**(-2) * d(:,6) )
         stiff(:, 8, 12) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( -(2.0/3.0)*d(:,5) )
         stiff(:, 8, 13) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (4.0/3.0)*d(:,4) )
         stiff(:, 8, 14) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d(:,3)
     &                   - (5.0/3.0)* (almBi)**(-2) * d_temp2 )   
         stiff(:, 9, 13) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,3) )
         stiff(:, 9, 14) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
         stiff(:, 10, 12) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d(:,4) * u1 + (-2.0/3.0)*d(:,5) * u2 )
         stiff(:, 10, 13) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d(:,5) * u1 + (4.0/3.0)*d(:,4) * u2 
     &                    +d(:,3) * u3 )
         stiff(:, 10, 14) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( (-2.0/3.0)*d(:,3) * u2 +d(:,4) * u3 
     &                    - (5.0/3.0)* (almBi)**(-2) 
     &                    * ( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3) )
c
c.... K31
c
         stiff(:, 12, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,5) - (5.0/3.0)* (almBi)**(-2) *d(:,5) )
         stiff(:, 12, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,1) )
         stiff(:, 13, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( 0.0 - (5.0/3.0)* (almBi)**(-2) *d(:,4))
         stiff(:, 13, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,5) )
         stiff(:, 13, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,6) )
         stiff(:, 14, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d(:,1) 
     &                   - (5.0/3.0)* (almBi)**(-2) *d_temp3 )   
         stiff(:, 14, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d(:,6) )
         stiff(:, 14, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (4.0/3.0)*d(:,5) )
         stiff(:, 15, 2) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,5) * u1 + (-2.0/3.0)*d(:,1) * u3 
     &                   - (5.0/3.0)* (almBi)**(-2) 
     &                   *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3 ) )
         stiff(:, 15, 3) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,5) * u2 +(-2.0/3.0)*d(:,6) * u3 )
         stiff(:, 15, 4) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,1) * u1 + (4.0/3.0)*d(:,5) * u3 
     &                   + d(:,6) * u2 )
c
c     
c.... K32
c     
         stiff(:, 12, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,4) )
         stiff(:, 12, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( 0.0 - (5.0/3.0)* (almBi)**(-2) *d(:,5) )
         stiff(:, 12, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                     * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                     * ( d(:,6) )
         stiff(:, 13, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,4) - (5.0/3.0)* (almBi)**(-2) *d(:,4) )
         stiff(:, 13, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,2) )
         stiff(:, 14, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( -(2.0/3.0)*d(:,6) )
         stiff(:, 14, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( -(2.0/3.0)*d(:,2) 
     &                   - (5.0/3.0)* (almBi)**(-2) *d_temp3)
         stiff(:, 14, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (4.0/3.0)*d(:,4) )
         stiff(:, 15, 7) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,4) * u1 + (-2.0/3.0)*d(:,6) * u3 )
         stiff(:, 15, 8) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,4) * u2 + (-2.0/3.0)*d(:,2) * u3 
     &                   - (5.0/3.0)* (almBi)**(-2) 
     &                   *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3 ) )
         stiff(:, 15, 9) = ShearMod * (det_d)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d(:,6) * u1 + d(:,2) * u2 
     &                   + (4.0/3.0)*d(:,4) * u3 )
c
c     
c.... K33
c     
         stiff(:, 12, 12) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d(:,3) )
         stiff(:, 12, 14) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d(:,5)-(5.0/3.0)* (almBi)**(-2) * d(:,5) )
         stiff(:, 13, 13) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d(:,3) )
         stiff(:, 13, 14) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
         stiff(:, 14, 12) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( -(2.0/3.0)*d(:,5) )
         stiff(:, 14, 13) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( -(2.0/3.0)*d(:,4) )
         stiff(:, 14, 14) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( (4.0/3.0)*d(:,3)
     &                    -(5.0/3.0)* (almBi)**(-2) * d_temp3 )   
         stiff(:, 15, 12) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d(:,3) * u1 + (-2.0/3.0)*d(:,5) * u3 )
         stiff(:, 15, 13) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d(:,3) * u2 + (-2.0/3.0)*d(:,4) * u3 )
         stiff(:, 15, 14) = ShearMod * (det_d)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * (d(:,5) * u1 + d(:,4) * u2 + (4.0/3.0)*d(:,3) * u3
     &                    - (5.0/3.0)* (almBi)**(-2)
     &                    *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3) )
         stiff(:, 15, 15) = con ! notice the + or -
c
c     
      end subroutine set_solid_kij_interior
c
c
      subroutine set_solid_difflux_interior(rmi, ri,  u1,   u2,
     &                                      u3,  con)
c------------------------------------------------------------------------------
! setting the diffusive flux for solid interior blocks( not from K_{ij} Y,j )
!
!input:
! u1 : x component of velocity
! u2 : y component of velocity
! u3 : z component of velocity
! con : thermal conductivity
!output:
! rmi : modified residual
! ri : residual

c------------------------------------------------------------------------------
        implicit none
c
        real*8, dimension(npro,nflow*(nsd+1)),intent(inout) :: rmi
        real*8, dimension(npro,nflow*(nsd+1)),intent(inout) :: ri
        real*8, dimension(npro),intent(in) :: u1, u2, u3, con
! local array
        real*8, dimension(npro,6) :: bq_af
c...........................................................................
        bq_af(:,:)= b_af(iblk_solid)%p(:,intp,:)
c.... diffusive flux in x1-direction

c         rmi(:,1) = zero ! already initialized
         rmi(:,2) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * (1.0/3.0)
     &            * ( 2 * bq_af(:,1) - bq_af(:,2) - bq_af(:,3) )
         rmi(:,3) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * bq_af(:,6)
         rmi(:,4) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * bq_af(:,5)
         rmi(:,5) =  (det_baf)**(-5.0/6.0) * ShearMod
     &            * ( (1.0/3.0) 
     &            * ( 2 * bq_af(:,1) - bq_af(:,2) - bq_af(:,3) ) *u1
     &            + bq_af(:,6) * u2 + bq_af(:,5) * u3 )
     &            + con * dtpdx
c                
        ri (:,2:5) = ri (:,2:5) + rmi(:,2:5)
c       rmi(:,2:5) = rmi(:,2:5) + qdi(:,2:5)
c
c!      flops = flops + 74*npro
c
c.... diffusive flux in x2-direction
c
c       rmi(:, 6) = zero
        rmi(:, 7) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * bq_af(:,6)
        rmi(:, 8) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * (1.0/3.0)
     &            *(-bq_af(:,1) + 2.0* bq_af(:,2) - bq_af(:,3) )
        rmi(:, 9) =  (det_baf)**(-5.0/6.0) * ShearMod
     &            * bq_af(:,4)
        rmi(:,10) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * ( bq_af(:,6) * u1 
     &            + (1.0/3.0) 
     &            * ( -bq_af(:,1) + 2.0* bq_af(:,2) - bq_af(:,3) ) *u2 
     &            + bq_af(:,4) * u3 )
     &            + con * dtpdy
c
      ri (:,7:10) = ri (:,7:10) + rmi(:,7:10)
c     rmi(:,7:10) = rmi(:,7:10) + qdi(:,2:5)
c
c!      flops = flops + 74*npro
c
c.... diffusive flux in x3-direction
c
c       rmi(:,11) = zero
        rmi(:,12) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * bq_af(:,5)
        rmi(:,13) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * bq_af(:,4)
        rmi(:,14) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * (1.0/3.0) 
     &            *( -bq_af(:,1) - bq_af(:,2) + 2.0* bq_af(:,3) )
        rmi(:,15) = (det_baf)**(-5.0/6.0) * ShearMod
     &            * ( bq_af(:,5) * u1 + bq_af(:,4) * u2
     &            + (1.0/3.0) 
     &            * ( -bq_af(:,1) - bq_af(:,2) + 2.0* bq_af(:,3) ) *u3 )
     &            + con * dtpdz
c
       ri (:,12:15) = ri (:,12:15) + rmi(:,12:15)
c!      flops = flops + 74*npro
c
      end subroutine set_solid_difflux_interior
c
c
      end module e3_solid_func_m
