      module e3if_solid_func_m
c-------------------------------------------------------------------------------
!
c-------------------------------------------------------------------------------
!        use e3if_solid_m
!        use number_def_m
!        use global_const_m
!        use intpt_m
        use solid_data_m
        use e3if_solid_data_m
        use propar_m, only: npro
        use global_const_m, only: nsd
        use timdat_m, only: almBi, alfBi, gamBi
        use number_def_m
        use inpdat_m, only: Delt
        implicit none
c
! arrays used commonly in inetegration calucation for solid
! matrix used in calculation the strain field at time n+af
        real*8, dimension(:,:,:), allocatable :: AS_0
        real*8, dimension(:,:,:), allocatable :: AS_1
c        
        contains
c
c..............................................................................
c~         subroutine set_solid_kij_if ()
c~ c-----------------------------------------------------------------
c~ ! setting the K_ij matrix for solid interface blocks so that
c~ ! F^{diff} = K_ij Y,{j}
c~ ! input:
c~ ! 
c~ ! output:
c~ ! stiff  (npro,nsd*nflow,nsd*nflow) : stiffness matrix,K_ij
c~ c-----------------------------------------------------------------
c~ c
c~          implicit none
c~ c
c~ !        real*8, dimension(npro,nsd*nflow,nsd*nflow),intent(inout)::stiff
c~ !        real*8, dimension(npro),intent(in) :: u1, u2, u3, con
c~ ! local arrays
c~          real*8, dimension(npro) :: d_temp1,d_temp2,d_temp3
c~          real*8, dimension(npro,b_size) :: bq_af
c~ c.................................................................
c~          d_temp1(:) = 2.0/3.0*d(:,1) - 1.0/3.0*d(:,2) - 1.0/3.0*d(:,3)
c~          d_temp2(:) = -1.0/3.0*d(:,1) + 2.0/3.0*d(:,2) - 1.0/3.0*d(:,3)
c~          d_temp3(:) = -1.0/3.0*d(:,1) - 1.0/3.0*d(:,2) + 2.0/3.0*d(:,3)
c~ c
c~ c.... K11
c~ c
c~          stiff(:, 2, 2) = ShearMod * (det_d)**(-5.0/6.0) 
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (4.0/3.0)*d(:,1)
c~      &                  - (5.0/3.0)* (almBi)**(-2.0) * d_temp1 )  
c~          stiff(:, 2, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (-2.0/3.0)*d(:,6)
c~          stiff(:, 2, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (-2.0/3.0)*d(:,5)
c~          stiff(:, 3, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (d(:,6)- (5.0/3.0)* (almBi)**(-2) * d(:,6))
c~          stiff(:, 3, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,1) )
c~          stiff(:, 4, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (d(:,5)- (5.0/3.0)* (almBi)**(-2) * d(:,5))
c~          stiff(:, 4, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,1))
c~          stiff(:, 5, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (4.0/3.0)*d(:,1) * u1 + d(:,6) * u2 + d(:,5) * u3
c~      &                  - (5.0/3.0)* (almBi)**(-2)
c~      &                  *( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3) )
c~          stiff(:, 5, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0/3.0)*d(:,6) * u1 + d(:,1) * u2 )
c~          stiff(:, 5, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0/3.0)*d(:,5) * u1 + d(:,1) * u3 )
c~          stiff(:, 5, 5) = con ! notice the + or -
c~ c
c~ c.... K12                 
c~ c     
c~          stiff(:, 2, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (4.0/3.0)*d(:,6) )
c~          stiff(:, 2, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0/3.0)*d(:,2) 
c~      &                  - (5.0/3.0)* (almBi)**(-2) * d_temp1 )
c~          stiff(:, 2, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0/3.0)*d(:,4) )
c~          stiff(:, 3, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (d(:,2))
c~          stiff(:, 3, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,6) 
c~      &                  - (5.0/3.0)* (almBi)**(-2) * d(:,6) )
c~          stiff(:, 4, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,4) )
c~          stiff(:, 4, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( 0.0 
c~      &                  - (5.0/3.0)* (almBi)**(-2) * d(:,5) )    
c~          stiff(:, 4, 9) = ShearMod * (det_d)**(-5.0/6.0) 
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,6) )
c~          stiff(:, 5, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (4.0/3.0)*d(:,6) * u1 + d(:,2) * u2 + d(:,4) * u3 )
c~          stiff(:, 5, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0/3.0)*d(:,2) * u1 + d(:,6) * u2 
c~      &                  - (5.0/3.0)* (almBi)**(-2)
c~      &                  * ( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3) )
c~          stiff(:, 5, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0/3.0)*d(:,4) * u1 + d(:,6) * u3 )
c~ c
c~ c.... K13
c~ c     
c~          stiff(:, 2,12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (4.0 /3.0)*d(:,5) )
c~          stiff(:, 2,13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0 /3.0)*d(:,4) )
c~          stiff(:, 2,14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0 /3.0)*d(:,3) 
c~      &                  - (5.0/3.0)* (almBi)**(-2) * d_temp1)
c~          stiff(:, 3,12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (d(:,4))
c~          stiff(:, 3,13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,5) )
c~          stiff(:, 3,14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( - (5.0/3.0)* (almBi)**(-2) * d(:,6) )
c~          stiff(:, 4,12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,3) )
c~          stiff(:, 4,14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,5)- (5.0/3.0)* (almBi)**(-2) * d(:,5) )
c~          stiff(:, 5,12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (4.0 /3.0)*d(:,5) * u1 + d(:,4) * u2 + d(:,3) * u3 )
c~          stiff(:, 5,13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0 /3.0)*d(:,4) * u1 + d(:,5) * u2  )
c~          stiff(:, 5,14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0 /3.0)*d(:,3) * u1 + d(:,5) * u3
c~      &                  - (5.0/3.0)* (almBi)**(-2)
c~      &                  * ( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3))
c~ c     
c~ c.... K21
c~ c     
c~          stiff(:, 7, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (d(:,6)- (5.0/3.0)* (almBi)**(-2) * d(:,6))
c~          stiff(:, 7, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,1) )
c~          stiff(:, 8, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0/3.0)*d(:,1) 
c~      &                  - (5.0/3.0)* (almBi)**(-2) * d_temp2 )  
c~          stiff(:, 8, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (4.0/3.0)*d(:,6)
c~          stiff(:, 8, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (-2.0/3.0)*d(:,5)
c~          stiff(:, 9, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * (- (5.0/3.0)* (almBi)**(-2) * d(:,4))
c~          stiff(:, 9, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  *  d(:,5)
c~          stiff(:, 9, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * d(:,6)
c~          stiff(:, 10, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,6) * u1 + (-2.0/3.0)*d(:,1) * u2 
c~      &                   - (5.0/3.0)* (almBi)**(-2)
c~      &                   * ( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3) )
c~          stiff(:, 10, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * (d(:,1) * u1 
c~      &                   + (4.0/3.0)*d(:,6) * u2 +d(:,5) * u3)
c~          stiff(:, 10, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( (-2.0/3.0)*d(:,5) * u2 +d(:,6) * u3 )
c~ c
c~ c.... K22
c~ c     
c~          stiff(:, 7, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,2) )
c~          stiff(:, 7, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,6) -(5.0/3.0)* (almBi)**(-2) * d(:,6) )
c~          stiff(:, 8, 7) = ShearMod * (almBi * det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( -(2.0/3.0)*d(:,6) )
c~          stiff(:, 8, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (4.0/3.0)*d(:,2)
c~      &                  -(5.0/3.0)* (almBi)**(-2) * d_temp2 )
c~          stiff(:, 8, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( (-2.0/3.0)*d(:,4) )
c~          stiff(:, 9, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
c~          stiff(:, 9, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                  * ( d(:,2) )
c~          stiff(:, 10, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,2) * u1 + (-2.0/3.0)*d(:,6) * u2 )
c~          stiff(:, 10, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * (d(:,6) * u1 + (4.0/3.0)*d(:,2) * u2 +d(:,4) * u3
c~      &                   - (5.0/3.0)* (almBi)**(-2)
c~      &                   *( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3))
c~          stiff(:, 10, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( (-2.0/3.0)*d(:,4) * u2 +d(:,2) * u3 )
c~          stiff(:, 10, 10) = con ! notice the + or -
c~ c
c~ c.... K23
c~ c     
c~          stiff(:, 7, 12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,4) )
c~          stiff(:, 7, 13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,5) )
c~          stiff(:, 7, 14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( 0.0 -(5.0/3.0)* (almBi)**(-2) * d(:,6) )
c~          stiff(:, 8, 12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( -(2.0/3.0)*d(:,5) )
c~          stiff(:, 8, 13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( (4.0/3.0)*d(:,4) )
c~          stiff(:, 8, 14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( (-2.0/3.0)*d(:,3)
c~      &                   - (5.0/3.0)* (almBi)**(-2) * d_temp2 )   
c~          stiff(:, 9, 13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,3) )
c~          stiff(:, 9, 14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
c~          stiff(:, 10, 12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( d(:,4) * u1 + (-2.0/3.0)*d(:,5) * u2 )
c~          stiff(:, 10, 13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( d(:,5) * u1 + (4.0/3.0)*d(:,4) * u2 
c~      &                    +d(:,3) * u3 )
c~          stiff(:, 10, 14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( (-2.0/3.0)*d(:,3) * u2 +d(:,4) * u3 
c~      &                    - (5.0/3.0)* (almBi)**(-2) 
c~      &                    * ( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3) )
c~ c
c~ c.... K31
c~ c
c~          stiff(:, 12, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,5) - (5.0/3.0)* (almBi)**(-2) *d(:,5) )
c~          stiff(:, 12, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,1) )
c~          stiff(:, 13, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( 0.0 - (5.0/3.0)* (almBi)**(-2) *d(:,4))
c~          stiff(:, 13, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,5) )
c~          stiff(:, 13, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,6) )
c~          stiff(:, 14, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( (-2.0/3.0)*d(:,1) 
c~      &                   - (5.0/3.0)* (almBi)**(-2) *d_temp3 )   
c~          stiff(:, 14, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( (-2.0/3.0)*d(:,6) )
c~          stiff(:, 14, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( (4.0/3.0)*d(:,5) )
c~          stiff(:, 15, 2) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,5) * u1 + (-2.0/3.0)*d(:,1) * u3 
c~      &                   - (5.0/3.0)* (almBi)**(-2) 
c~      &                   *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3 ) )
c~          stiff(:, 15, 3) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,5) * u2 +(-2.0/3.0)*d(:,6) * u3 )
c~          stiff(:, 15, 4) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,1) * u1 + (4.0/3.0)*d(:,5) * u3 
c~      &                   + d(:,6) * u2 )
c~ c
c~ c     
c~ c.... K32
c~ c     
c~          stiff(:, 12, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,4) )
c~          stiff(:, 12, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( 0.0 - (5.0/3.0)* (almBi)**(-2) *d(:,5) )
c~          stiff(:, 12, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                     * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                     * ( d(:,6) )
c~          stiff(:, 13, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,4) - (5.0/3.0)* (almBi)**(-2) *d(:,4) )
c~          stiff(:, 13, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,2) )
c~          stiff(:, 14, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( -(2.0/3.0)*d(:,6) )
c~          stiff(:, 14, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( -(2.0/3.0)*d(:,2) 
c~      &                   - (5.0/3.0)* (almBi)**(-2) *d_temp3)
c~          stiff(:, 14, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( (4.0/3.0)*d(:,4) )
c~          stiff(:, 15, 7) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,4) * u1 + (-2.0/3.0)*d(:,6) * u3 )
c~          stiff(:, 15, 8) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,4) * u2 + (-2.0/3.0)*d(:,2) * u3 
c~      &                   - (5.0/3.0)* (almBi)**(-2) 
c~      &                   *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3 ) )
c~          stiff(:, 15, 9) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                   * ( d(:,6) * u1 + d(:,2) * u2 
c~      &                   + (4.0/3.0)*d(:,4) * u3 )
c~ c
c~ c     
c~ c.... K33
c~ c     
c~          stiff(:, 12, 12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( d(:,3) )
c~          stiff(:, 12, 14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( d(:,5)-(5.0/3.0)* (almBi)**(-2) * d(:,5) )
c~          stiff(:, 13, 13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( d(:,3) )
c~          stiff(:, 13, 14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
c~          stiff(:, 14, 12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( -(2.0/3.0)*d(:,5) )
c~          stiff(:, 14, 13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( -(2.0/3.0)*d(:,4) )
c~          stiff(:, 14, 14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( (4.0/3.0)*d(:,3)
c~      &                    -(5.0/3.0)* (almBi)**(-2) * d_temp3 )   
c~          stiff(:, 15, 12) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( d(:,3) * u1 + (-2.0/3.0)*d(:,5) * u3 )
c~          stiff(:, 15, 13) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * ( d(:,3) * u2 + (-2.0/3.0)*d(:,4) * u3 )
c~          stiff(:, 15, 14) = ShearMod * (det_d)**(-5.0/6.0)
c~      &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
c~      &                    * (d(:,5) * u1 + d(:,4) * u2 + (4.0/3.0)*d(:,3) * u3
c~      &                    - (5.0/3.0)* (almBi)**(-2)
c~      &                    *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3) )
c~          stiff(:, 15, 15) = con ! notice the + or -
c~ c
c~ c     
c~         end subroutine set_solid_kij_if
c~ c
c~ c..............................................................................
c~         subroutine set_solid_difflux_interior()
c~ c------------------------------------------------------------------------------
c~ ! setting the diffusive flux for solid interface blocks( not from K_{ij} Y,j )
c~ !
c~ !input:
c~ ! u1 : x component of velocity
c~ ! u2 : y component of velocity
c~ ! u3 : z component of velocity
c~ ! con : thermal conductivity
c~ !output:
c~ ! rmi : modified residual
c~ ! ri : residual
c~ 
c~ c------------------------------------------------------------------------------
c~          implicit none
c~ c
c~ !        real*8, dimension(npro,nflow*(nsd+1)),intent(inout) :: rmi
c~ !        real*8, dimension(npro,nflow*(nsd+1)),intent(inout) :: ri
c~ !        real*8, dimension(npro),intent(in) :: u1, u2, u3, con
c~ ! local array
c~          real*8, dimension(npro,6) :: bq_af
c~ c...........................................................................
c~          bq_af(:,:)= b_af(iblk_solid)%p(:,intp,:)
c~ c.... diffusive flux in x1-direction
c~ 
c~ c         rmi(:,1) = zero ! already initialized
c~           rmi(:,2) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * (1.0/3.0)
c~      &             * ( 2 * bq_af(:,1) - bq_af(:,2) - bq_af(:,3) )
c~           rmi(:,3) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * bq_af(:,6)
c~           rmi(:,4) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * bq_af(:,5)
c~           rmi(:,5) =  (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * ( (1.0/3.0) 
c~      &             * ( 2 * bq_af(:,1) - bq_af(:,2) - bq_af(:,3) ) *u1
c~      &             + bq_af(:,6) * u2 + bq_af(:,5) * u3 )
c~      &             + con * dtpdx
c~ c                
c~           ri (:,2:5) = ri (:,2:5) + rmi(:,2:5)
c~ c       rmi(:,2:5) = rmi(:,2:5) + qdi(:,2:5)
c~ c
c~ c!      flops = flops + 74*npro
c~ c
c~ c.... diffusive flux in x2-direction
c~ c
c~ c       rmi(:, 6) = zero
c~          rmi(:, 7) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * bq_af(:,6)
c~          rmi(:, 8) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * (1.0/3.0)
c~      &             *(-bq_af(:,1) + 2.0* bq_af(:,2) - bq_af(:,3) )
c~          rmi(:, 9) =  (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * bq_af(:,4)
c~          rmi(:,10) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * ( bq_af(:,6) * u1 
c~      &             + (1.0/3.0) 
c~      &             * ( -bq_af(:,1) + 2.0* bq_af(:,2) - bq_af(:,3) ) *u2 
c~      &             + bq_af(:,4) * u3 )
c~      &             + con * dtpdy
c~ c
c~          ri (:,7:10) = ri (:,7:10) + rmi(:,7:10)
c~ c     rmi(:,7:10) = rmi(:,7:10) + qdi(:,2:5)
c~ c
c~ c!      flops = flops + 74*npro
c~ c
c~ c.... diffusive flux in x3-direction
c~ c
c~ c       rmi(:,11) = zero
c~          rmi(:,12) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * bq_af(:,5)
c~          rmi(:,13) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * bq_af(:,4)
c~          rmi(:,14) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * (1.0/3.0) 
c~      &             *( -bq_af(:,1) - bq_af(:,2) + 2.0* bq_af(:,3) )
c~          rmi(:,15) = (det_baf)**(-5.0/6.0) * ShearMod
c~      &             * ( bq_af(:,5) * u1 + bq_af(:,4) * u2
c~      &             + (1.0/3.0) 
c~      &             * ( -bq_af(:,1) - bq_af(:,2) + 2.0* bq_af(:,3) ) *u3 )
c~      &             + con * dtpdz
c~ c
c~          ri (:,12:15) = ri (:,12:15) + rmi(:,12:15)
c~ c!      flops = flops + 74*npro
c~ c
c~         end subroutine set_solid_difflux_interior  
c~ c
c.............................................................................
        subroutine e3if0_var_solid(intp, var)
c-----------------------------------------------------------------------------
! for phase_0
! getting the varibles needed for solid calculation at each quadrature points
! for interface elements
c-----------------------------------------------------------------------------
c       
        use e3if_param_m
        implicit none
c
        integer :: intp ! could be improved
        type(var_t), dimension(:),pointer, intent(in) :: var
! local
        integer :: ipro
c        
        allocate(AS_0(npro,b_size,b_size))     
! pass the array
        do ipro = 1, npro
          dudx_0(ipro,:) = var(ipro)%grad_y(1,2:4)
          dudy_0(ipro,:) = var(ipro)%grad_y(2,2:4)
          dudz_0(ipro,:) = var(ipro)%grad_y(3,2:4)
  !        
          dtpdx_0(ipro) = var(ipro)%grad_y(1,5)
          dtpdy_0(ipro) = var(ipro)%grad_y(2,5)
          dtpdz_0(ipro) = var(ipro)%grad_y(3,5)
        enddo 
c               
        call calc_as_matrix_if( AS_0, dudx_0, dudy_0, dudz_0 )
        if_d0(:,:) = almBi * if_b0(iblkif_solid)%p(:,intp,:)
     &             + alfBi * Delt(1) * (almBi - gamBi) 
     &             * if_b0_dot(iblkif_solid)%p(:,intp,:)
        call setB_af_if( AS_0, if_d0, 
     &                   if_b0_af(iblkif_solid)%p(:,intp,:) )
        call get_det_if(if_b0_af(iblkif_solid)%p(:,intp,:),if_det_baf0)
        if_Ja_def0 = (if_det_baf0)**0.5
        call get_det_if(if_d0,if_det_d0)
c
        deallocate(AS_0)        
c        
        end subroutine e3if0_var_solid
c
c..............................................................................
        subroutine e3if1_var_solid(intp,var)
c-----------------------------------------------------------------------------
! for phase_1
! getting the varibles needed for solid calculation at each quadrature points
! for interface elements
c-----------------------------------------------------------------------------
c
        use e3if_param_m
        implicit none
c
        integer :: intp ! could be improved
        type(var_t), dimension(:), pointer, intent(in) :: var
! local
         integer :: ipro
c                  
        allocate(AS_1(npro,b_size,b_size))
! pass the array
        do ipro = 1, npro
          dudx_1(ipro,:) = var(ipro)%grad_y(1,2:4)
          dudy_1(ipro,:) = var(ipro)%grad_y(2,2:4)
          dudz_1(ipro,:) = var(ipro)%grad_y(3,2:4)
  !        
          dtpdx_1(ipro) = var(ipro)%grad_y(1,5)
          dtpdy_1(ipro) = var(ipro)%grad_y(2,5)
          dtpdz_1(ipro) = var(ipro)%grad_y(3,5)
        enddo
c               
        call calc_as_matrix_if( AS_1, dudx_1, dudy_1, dudz_1 )
        if_d1(:,:) = almBi * if_b1(iblkif_solid)%p(:,intp,:)
     &             + alfBi * Delt(1) * (almBi - gamBi) 
     &             * if_b1_dot(iblkif_solid)%p(:,intp,:)
        call setB_af_if( AS_1, if_d1, 
     &                   if_b1_af(iblkif_solid)%p(:,intp,:) )
        call get_det_if( if_b1_af(iblkif_solid)%p(:,intp,:),if_det_baf1 )
        if_Ja_def1 = (if_det_baf1)**0.5
        call get_det_if( if_d1,if_det_d1 )
c
        deallocate(AS_1)        
c  
        
        end subroutine e3if1_var_solid    
c
c..............................................................................
        subroutine calc_as_matrix_if(A_f,  dudx_f,  dudy_f,
     &                               dudz_f)
c------------------------------------------------------------------------------
! calculating the A matrix used in calculation the strain field at time n+af
c------------------------------------------------------------------------------     
c
        implicit none
c
! pass the arrays
        real*8, dimension(npro,b_size,b_size), intent(inout)::A_f
        real*8, dimension(npro,nsd), intent(in) :: dudx_f, dudy_f, 
     &                                             dudz_f
! initialize the AS matrix
         A_f = zero !add
c
c
         A_f(:,1,1) = 2 * dudx_f(:,1)
         A_f(:,1,5) = 2 * dudz_f(:,1)
         A_f(:,1,6) = 2 * dudy_f(:,1)
c
         A_f(:,2,2) = 2 * dudy_f(:,2)
         A_f(:,2,4) = 2 * dudz_f(:,2)
         A_f(:,2,6) = 2 * dudx_f(:,2)
c
         A_f(:,3,3) = 2 * dudz_f(:,3)
         A_f(:,3,4) = 2 * dudy_f(:,3)
         A_f(:,3,5) = 2 * dudx_f(:,3)
c
         A_f(:,4,2) = dudy_f(:,3)
         A_f(:,4,3) = dudz_f(:,2)
         A_f(:,4,4) = dudy_f(:,2) + dudz_f(:,3)
         A_f(:,4,5) = dudx_f(:,2)
         A_f(:,4,6) = dudx_f(:,3)
c
         A_f(:,5,1) = dudx_f(:,3)
         A_f(:,5,3) = dudz_f(:,1)
         A_f(:,5,4) = dudy_f(:,1)
         A_f(:,5,5) = dudx_f(:,1) + dudz_f(:,3)
         A_f(:,5,6) = dudy_f(:,3)
c         
         A_f(:,6,1) = dudx_f(:,2)
         A_f(:,6,2) = dudy_f(:,1)
         A_f(:,6,4) = dudz_f(:,1)
         A_f(:,6,5) = dudz_f(:,2)
         A_f(:,6,6) = dudx_f(:,1) + dudy_f(:,2)
c
        end subroutine calc_as_matrix_if
c
c..............................................................................
      subroutine setB_af_if(A_f, d_f, bf_af_qt)
c------------------------------------------------------------------------------
! calculating the left Cauchy-green tensor at time step n+af
c------------------------------------------------------------------------------
c         use function_and_complex
         implicit none 
c
! passing arrays
         real*8, dimension(npro,b_size,b_size),intent(in) :: A_f
         real*8, dimension(npro,b_size), intent(in) :: d_f
         real*8, dimension(npro,b_size), intent(inout) :: bf_af_qt
! local varibles         
         real*8, dimension(b_size,b_size) :: ident
         real*8, dimension(b_size,b_size) :: temp_matrix
         integer :: i
c
         ident = zero
         do i = 1, b_size
           ident(i,i) = one
         enddo
c
         do i = 1, npro
          temp_matrix(:,:) = (one/almBi) * ident(:,:) + ( gamBi * Delt(1)
     &                     *alfBi/(almBi)**two ) * A_f(i,:,:) !check here
          bf_af_qt(i,:) = matmul(temp_matrix(:,:) , d_f(i,:))
c           temp_matrix(:,:) = (almBi) * ident(:,:) - ( gamBi * Delt(1)
c     &                     *alfBi) * AS(i,:,:) !check here
c
c           call lu_solve(temp_matrix, d(i,:), b_af(iblk_solid)%p(i,intp,:), b_size)
         enddo
c
       end subroutine setB_af_if
c
c...............................................................................
       subroutine get_det_if(matr,det_matr)
c
        implicit none 
c
        real*8, dimension(npro,b_size),intent(in) :: matr
        real*8, dimension(npro) :: det_matr
c
        det_matr(:) = matr(:,1) * matr(:,2) * matr(:,3)
     &+             matr(:,6) * matr(:,4) * matr(:,5)
     &+             matr(:,5) * matr(:,6) * matr(:,4) 
     &-             matr(:,1) * matr(:,4) * matr(:,4)
     &-             matr(:,5) * matr(:,2) * matr(:,5)
     &-             matr(:,6) * matr(:,6) * matr(:,3)
c
      end subroutine get_det_if
c
c
      end module e3if_solid_func_m  
