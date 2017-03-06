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
        subroutine set_solid_kij_if (Kij_f, con_f,   u_f, 
     &                               d_f,   det_d_f, shearMod_f)
c-----------------------------------------------------------------
! setting the K_ij matrix for solid interface blocks so that
! F^{diff} = K_ij Y,{j}
! input:
! 
! output:
! Kij_f  (npro,nsd*nflow,nsd*nflow) : stiffness matrix,K_ij
c-----------------------------------------------------------------
c
         implicit none
c
         real*8, dimension(:,:,:,:,:), pointer, intent(out) :: Kij_f
         real*8, dimension(npro,nsd), intent(in) :: u_f
         real*8, dimension(npro), intent(in) :: con_f
         real*8, dimension(npro,b_size) :: d_f
         real*8, dimension(npro) :: det_d_f
         real*8, dimension(npro) :: shearMod_f
         
! local arrays
         real*8, dimension(npro) :: u1, u2, u3
         real*8, dimension(npro) :: d_temp1,d_temp2,d_temp3
         real*8, dimension(npro,b_size) :: bq_af
c.................................................................
         d_temp1(:) = 2.0/3.0*d_f(:,1) - 1.0/3.0*d_f(:,2) - 1.0/3.0*d_f(:,3)
         d_temp2(:) = -1.0/3.0*d_f(:,1) + 2.0/3.0*d_f(:,2) - 1.0/3.0*d_f(:,3)
         d_temp3(:) = -1.0/3.0*d_f(:,1) - 1.0/3.0*d_f(:,2) + 2.0/3.0*d_f(:,3)
c
c
         u1 = u_f(:,1)
         u2 = u_f(:,2)
         u3 = u_f(:,3)
c
         Kij_f = zero
c
c.... K11
c
         Kij_f(:,1,1,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0) 
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d_f(:,1)
     &                  - (5.0/3.0)* (almBi)**(-2.0) * d_temp1 )  
         Kij_f(:,1,1,2,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (-2.0/3.0)*d_f(:,6)
         Kij_f(:,1,1,2,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (-2.0/3.0)*d_f(:,5)
         Kij_f(:,1,1,3,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d_f(:,6)- (5.0/3.0)* (almBi)**(-2) * d_f(:,6))
         Kij_f(:,1,1,3,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,1) )
         Kij_f(:,1,1,4,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d_f(:,5)- (5.0/3.0)* (almBi)**(-2) * d_f(:,5))
         Kij_f(:,1,1,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,1))
         Kij_f(:,1,1,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d_f(:,1) * u1 + d_f(:,6) * u2 + d_f(:,5) * u3
     &                  - (5.0/3.0)* (almBi)**(-2)
     &                  *( d_temp1 * u1 + d_f(:,6) * u2 + d_f(:,5) * u3) )
         Kij_f(:,1,1,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d_f(:,6) * u1 + d_f(:,1) * u2 )
         Kij_f(:,1,1,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d_f(:,5) * u1 + d_f(:,1) * u3 )
         Kij_f(:,1,1,5,5) = con_f ! notice the + or -
c
c.... K12                 
c     
         Kij_f(:,1,2,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d_f(:,6) )
         Kij_f(:,1,2,2,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d_f(:,2) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_temp1 )
         Kij_f(:,1,2,2,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d_f(:,4) )
         Kij_f(:,1,2,3,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d_f(:,2))
         Kij_f(:,1,2,3,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,6) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_f(:,6) )
         Kij_f(:,1,2,4,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,4) )
         Kij_f(:,1,2,4,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( 0.0 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_f(:,5) )    
         Kij_f(:,1,2,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0) 
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,6) )
         Kij_f(:,1,2,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d_f(:,6) * u1 + d_f(:,2) * u2 + d_f(:,4) * u3 )
         Kij_f(:,1,2,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d_f(:,2) * u1 + d_f(:,6) * u2 
     &                  - (5.0/3.0)* (almBi)**(-2)
     &                  * ( d_temp1 * u1 + d_f(:,6) * u2 + d_f(:,5) * u3) )
         Kij_f(:,1,2,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d_f(:,4) * u1 + d_f(:,6) * u3 )
c
c.... K13
c     
         Kij_f(:,1,3,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0 /3.0)*d_f(:,5) )
         Kij_f(:,1,3,2,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0 /3.0)*d_f(:,4) )
         Kij_f(:,1,3,2,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0 /3.0)*d_f(:,3) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_temp1)
         Kij_f(:,1,3,3,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d_f(:,4))
         Kij_f(:,1,3,3,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,5) )
         Kij_f(:,1,3,3,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( - (5.0/3.0)* (almBi)**(-2) * d_f(:,6) )
         Kij_f(:,1,3,4,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,3) )
         Kij_f(:,1,3,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,5)- (5.0/3.0)* (almBi)**(-2) * d_f(:,5) )
         Kij_f(:,1,3,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0 /3.0)*d_f(:,5) * u1 + d_f(:,4) * u2 + d_f(:,3) * u3 )
         Kij_f(:,1,3,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0 /3.0)*d_f(:,4) * u1 + d_f(:,5) * u2  )
         Kij_f(:,1,3,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0 /3.0)*d_f(:,3) * u1 + d_f(:,5) * u3
     &                  - (5.0/3.0)* (almBi)**(-2)
     &                  * ( d_temp1 * u1 + d_f(:,6) * u2 + d_f(:,5) * u3))
c     
c.... K21
c     
         Kij_f(:,2,1,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (d_f(:,6)- (5.0/3.0)* (almBi)**(-2) * d_f(:,6))
         Kij_f(:,2,1,2,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,1) )
         Kij_f(:,2,1,3,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d_f(:,1) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_temp2 )  
         Kij_f(:,2,1,3,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (4.0/3.0)*d_f(:,6)
         Kij_f(:,2,1,3,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (-2.0/3.0)*d_f(:,5)
         Kij_f(:,2,1,4,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * (- (5.0/3.0)* (almBi)**(-2) * d_f(:,4))
         Kij_f(:,2,1,4,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  *  d_f(:,5)
         Kij_f(:,2,1,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * d_f(:,6)
         Kij_f(:,2,1,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,6) * u1 + (-2.0/3.0)*d_f(:,1) * u2 
     &                   - (5.0/3.0)* (almBi)**(-2)
     &                   * ( d_f(:,6) * u1 + d_temp2 * u2 + d_f(:,4) * u3) )
         Kij_f(:,2,1,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * (d_f(:,1) * u1 
     &                   + (4.0/3.0)*d_f(:,6) * u2 +d_f(:,5) * u3)
         Kij_f(:,2,1,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d_f(:,5) * u2 +d_f(:,6) * u3 )
c
c.... K22
c     
         Kij_f(:,2,2,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,2) )
         Kij_f(:,2,2,2,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,6) 
     &                  - (5.0/3.0)* (almBi)**(-2) * d_f(:,6) )
         Kij_f(:,2,2,3,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( -(2.0/3.0)*d_f(:,6) )
         Kij_f(:,2,2,3,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (4.0/3.0)*d_f(:,2)
     &                  -(5.0/3.0)* (almBi)**(-2) * d_temp2 )
         Kij_f(:,2,2,3,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( (-2.0/3.0)*d_f(:,4) )
         Kij_f(:,2,2,4,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,4) -(5.0/3.0)* (almBi)**(-2) * d_f(:,4))
         Kij_f(:,2,2,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                  * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                  * ( d_f(:,2) )
         Kij_f(:,2,2,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,2) * u1 + (-2.0/3.0)*d_f(:,6) * u2 )
         Kij_f(:,2,2,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * (d_f(:,6) * u1 + (4.0/3.0)*d_f(:,2) * u2 +d_f(:,4) * u3
     &                   - (5.0/3.0)* (almBi)**(-2)
     &                   *( d_f(:,6) * u1 + d_temp2 * u2 + d_f(:,4) * u3))
         Kij_f(:,2,2,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d_f(:,4) * u2 +d_f(:,2) * u3 )
         Kij_f(:,2,2,5,5) = con_f ! notice the + or -
c
c.... K23
c     
         Kij_f(:,2,3,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,4) )
         Kij_f(:,2,3,2,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,5) )
         Kij_f(:,2,3,2,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( 0.0 -(5.0/3.0)* (almBi)**(-2) * d_f(:,6) )
         Kij_f(:,2,3,3,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( -(2.0/3.0)*d_f(:,5) )
         Kij_f(:,2,3,3,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (4.0/3.0)*d_f(:,4) )
         Kij_f(:,2,3,3,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d_f(:,3)
     &                   - (5.0/3.0)* (almBi)**(-2) * d_temp2 )   
         Kij_f(:,2,3,4,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,3) )
         Kij_f(:,2,3,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,4) -(5.0/3.0)* (almBi)**(-2) * d_f(:,4))
         Kij_f(:,2,3,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d_f(:,4) * u1 + (-2.0/3.0)*d_f(:,5) * u2 )
         Kij_f(:,2,3,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d_f(:,5) * u1 + (4.0/3.0)*d_f(:,4) * u2 
     &                    +d_f(:,3) * u3 )
         Kij_f(:,2,3,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( (-2.0/3.0)*d_f(:,3) * u2 +d_f(:,4) * u3 
     &                    - (5.0/3.0)* (almBi)**(-2) 
     &                    * ( d_f(:,6) * u1 + d_temp2 * u2 + d_f(:,4) * u3) )
c
c.... K31
c
         Kij_f(:,3,1,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,5) - (5.0/3.0)* (almBi)**(-2) *d_f(:,5) )
         Kij_f(:,3,1,2,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,1) )
         Kij_f(:,3,1,3,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( 0.0 - (5.0/3.0)* (almBi)**(-2) *d_f(:,4))
         Kij_f(:,3,1,3,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,5) )
         Kij_f(:,3,1,3,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,6) )
         Kij_f(:,3,1,4,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d_f(:,1) 
     &                   - (5.0/3.0)* (almBi)**(-2) *d_temp3 )   
         Kij_f(:,3,1,4,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (-2.0/3.0)*d_f(:,6) )
         Kij_f(:,3,1,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (4.0/3.0)*d_f(:,5) )
         Kij_f(:,3,1,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,5) * u1 + (-2.0/3.0)*d_f(:,1) * u3 
     &                   - (5.0/3.0)* (almBi)**(-2) 
     &                   *( d_f(:,5) * u1 + d_f(:,4) * u2 + d_temp3 * u3 ) )
         Kij_f(:,3,1,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,5) * u2 +(-2.0/3.0)*d_f(:,6) * u3 )
         Kij_f(:,3,1,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,1) * u1 + (4.0/3.0)*d_f(:,5) * u3 
     &                   + d_f(:,6) * u2 )
c
c     
c.... K32
c     
         Kij_f(:,3,2,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,4) )
         Kij_f(:,3,2,2,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( 0.0 - (5.0/3.0)* (almBi)**(-2) *d_f(:,5) )
         Kij_f(:,3,2,2,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                     * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                     * ( d_f(:,6) )
         Kij_f(:,3,2,3,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,4) - (5.0/3.0)* (almBi)**(-2) *d_f(:,4) )
         Kij_f(:,3,2,3,5) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,2) )
         Kij_f(:,3,2,4,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( -(2.0/3.0)*d_f(:,6) )
         Kij_f(:,3,2,4,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( -(2.0/3.0)*d_f(:,2) 
     &                   - (5.0/3.0)* (almBi)**(-2) *d_temp3)
         Kij_f(:,3,2,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( (4.0/3.0)*d_f(:,4) )
         Kij_f(:,3,2,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,4) * u1 + (-2.0/3.0)*d_f(:,6) * u3 )
         Kij_f(:,3,2,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,4) * u2 + (-2.0/3.0)*d_f(:,2) * u3 
     &                   - (5.0/3.0)* (almBi)**(-2) 
     &                   *( d_f(:,5) * u1 + d_f(:,4) * u2 + d_temp3 * u3 ) )
         Kij_f(:,3,2,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                   * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                   * ( d_f(:,6) * u1 + d_f(:,2) * u2 
     &                   + (4.0/3.0)*d_f(:,4) * u3 )
c
c     
c.... K33
c     
         Kij_f(:,3,3,2,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d_f(:,3) )
         Kij_f(:,3,3,2,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d_f(:,5)-(5.0/3.0)* (almBi)**(-2) * d_f(:,5) )
         Kij_f(:,3,3,3,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d_f(:,3) )
         Kij_f(:,3,3,3,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d_f(:,4) -(5.0/3.0)* (almBi)**(-2) * d_f(:,4))
         Kij_f(:,3,3,4,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( -(2.0/3.0)*d_f(:,5) )
         Kij_f(:,3,3,4,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( -(2.0/3.0)*d_f(:,4) )
         Kij_f(:,3,3,4,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( (4.0/3.0)*d_f(:,3)
     &                    -(5.0/3.0)* (almBi)**(-2) * d_temp3 )   
         Kij_f(:,3,3,5,2) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d_f(:,3) * u1 + (-2.0/3.0)*d_f(:,5) * u3 )
         Kij_f(:,3,3,5,3) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * ( d_f(:,3) * u2 + (-2.0/3.0)*d_f(:,4) * u3 )
         Kij_f(:,3,3,5,4) = shearMod_f * (det_d_f)**(-5.0/6.0)
     &                    * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)
     &                    * (d_f(:,5) * u1 + d_f(:,4) * u2 + (4.0/3.0)*d_f(:,3) * u3
     &                    - (5.0/3.0)* (almBi)**(-2)
     &                    *( d_f(:,5) * u1 + d_f(:,4) * u2 + d_temp3 * u3) )
         Kij_f(:,3,3,5,5) = con_f ! notice the + or -
c
c     
        end subroutine set_solid_kij_if
c
c..............................................................................
        subroutine set_solid_difflux_if( flux,         shearMod_f, 
     &                                   con_f,         det_baf_f,
     &                                   u_f,           dtpdx_f,
     &                                   dtpdy_f,       dtpdz_f,                      
     &                                   bq_af )
c------------------------------------------------------------------------------
! setting the diffusive flux for solid interface blocks( not from K_{ij} Y,j )
!
!input:
! u_f : velocity field
! con_f : thermal conductivity
!output:
! fdiff : residual

c------------------------------------------------------------------------------
         use conpar_m, only: nflow
         implicit none
c
         real*8, dimension(nsd,nflow), intent(out) :: flux
         real*8, intent(in) :: shearMod_f, con_f, det_baf_f
         real*8, dimension(nsd), intent(in) :: u_f
         real*8, intent(in) :: dtpdx_f, dtpdy_f, dtpdz_f
         real*8, dimension(b_size), intent(in) :: bq_af
! local array
         real*8 :: u1, u2, u3
c...........................................................................
         u1 = u_f(1)
         u2 = u_f(2)
         u3 = u_f(3)
! intialize the flux
         flux = zero
c.... diffusive flux in x1-direction
c
          flux(1,2) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * (1.0/3.0)
     &             * ( 2 * bq_af(1) - bq_af(2) - bq_af(3) )
          flux(1,3) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * bq_af(6)
          flux(1,4) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * bq_af(5)
          flux(1,5) =  (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * ( (1.0/3.0) 
     &             * ( 2 * bq_af(1) - bq_af(2) - bq_af(3) ) *u1
     &             + bq_af(6) * u2 + bq_af(5) * u3 )
     &             + con_f * dtpdx_f
c                
c.... diffusive flux in x2-direction
c
         flux(2, 2) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * bq_af(6)
         flux(2, 3) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * (1.0/3.0)
     &             *(-bq_af(1) + 2.0* bq_af(2) - bq_af(3) )
         flux(2, 4) =  (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * bq_af(4)
         flux(2, 5) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * ( bq_af(6) * u1 
     &             + (1.0/3.0) 
     &             * ( -bq_af(1) + 2.0* bq_af(2) - bq_af(3) ) *u2 
     &             + bq_af(4) * u3 )
     &             + con_f * dtpdy_f
c
c.... diffusive flux in x3-direction
c
         flux(3,2) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * bq_af(5)
         flux(3,3) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * bq_af(4)
         flux(3,4) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * (1.0/3.0) 
     &             *( -bq_af(1) - bq_af(2) + 2.0* bq_af(3) )
         flux(3,5) = (det_baf_f)**(-5.0/6.0) * shearMod_f
     &             * ( bq_af(5) * u1 + bq_af(4) * u2
     &             + (1.0/3.0) 
     &             * ( -bq_af(1) - bq_af(2) + 2.0* bq_af(3) ) *u3 )
     &             + con_f * dtpdz_f
c
c
        end subroutine set_solid_difflux_if  
c
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
      subroutine set_kinematic_mu_solid( mu_f,      con_f,    shearMod_f,
     &                                   bulkMod_f, mater_f)
c-------------------------------------------------------------------------------
!
c-------------------------------------------------------------------------------
c
        use conpar_m, only:nflow
        use number_def_m
        use eqn_state_m
        use e3if_param_m, only: length_h
        implicit none
! pass arrays
        real*8, dimension(npro,nflow), intent(out) :: mu_f
        real*8, dimension(npro), intent(in) :: con_f
        real*8, dimension(npro), intent(in) :: shearMod_f, bulkMod_f
        integer :: mater_f
! local varibles
        real*8, dimension(npro) :: c_s
        real*8, dimension(npro) :: nu_p !poisson's ratio
        real*8, dimension(npro) :: temp_solid        
        
c....................
        nu_p = ( three*bulkMod_f - two*shearMod_f)
     &       / (two * (three*bulkMod_f + shearMod_f) )
        rho_ref_s = mat_prop(mater_f,iprop_solid_1_rho_ref,1)
        c_s = ( two * shearMod_f * (one - nu_p)
     &      / (rho_ref_s * (one - two*nu_p )) )**pt5
        temp_solid = shearMod_f / c_s * length_h
c         
        mu_f = zero
        mu_f(:,2) = temp_solid
        mu_f(:,3) = temp_solid
        mu_f(:,4) = temp_solid
        mu_f(:,5) = con_f   
c
c
      end subroutine set_kinematic_mu_solid
c
c
      end module e3if_solid_func_m  
