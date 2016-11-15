       subroutine getLHS_bs(shape,         shgb,   d,      det_d,
     &                      bdy_shearMod,  u1,     u2,     u3,
     &                      WdetJb,        con,    bnorm,  lnode,
     &                      EGmass_bs)
c
       include "common.h"
c      
       real*8, dimension(npro,nshl)::shape
       real*8, dimension(npro,nshl,nsd)::shgb
       real*8, dimension(npro,6)::d
       real*8, dimension(npro)::det_d,bdy_shearMod
       real*8, dimension(npro)::u1, u2, u3
       real*8, dimension(npro):: WdetJb
       real*8, dimension(npro):: con
       real*8, dimension(npro,nsd)::bnorm
       integer, dimension(27):: lnode
c
       real*8, dimension(npro,nedof,nedof)::EGmass_bs
c    
       real*8, dimension(npro,3*nflow,3*nflow)::stiff
       real*8, dimension(npro) :: d_temp1,d_temp2,d_temp3
       real*8, dimension(npro) :: Wshg1,Wshg2,Wshg3
       real*8, dimension(npro,nflow,nflow) :: stif1,stif2,stif3
       integer :: j,j0,jdof,jdof2, jdof3, ib
       integer :: i,i0,idof,idof2, idof3, ia

c
c
        stiff = zero!initialize
c
        d_temp1(:) = 2.0/3.0*d(:,1) - 1.0/3.0*d(:,2) - 1.0/3.0*d(:,3)
        d_temp2(:) = -1.0/3.0*d(:,1) + 2.0/3.0*d(:,2) - 1.0/3.0*d(:,3)
        d_temp3(:) = -1.0/3.0*d(:,1) - 1.0/3.0*d(:,2) + 2.0/3.0*d(:,3)
c        
c.... K11
c
         stiff(:, 2, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0/3.0)*d(:,1) - (5.0/3.0)* (almBi)**(-2) * d_temp1 )  
         stiff(:, 2, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    (-2.0/3.0)*d(:,6)
         stiff(:, 2, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    (-2.0/3.0)*d(:,5)
         stiff(:, 3, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    (d(:,6)- (5.0/3.0)* (almBi)**(-2) * d(:,6))
         stiff(:, 3, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,1) )
         stiff(:, 4, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    (d(:,5)- (5.0/3.0)* (almBi)**(-2) * d(:,5))
         stiff(:, 4, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,1))
         stiff(:, 5, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0/3.0)*d(:,1) * u1 + d(:,6) * u2 + d(:,5) * u3 - 
     &                    (5.0/3.0)* (almBi)**(-2) *( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3) )
         stiff(:, 5, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0/3.0)*d(:,6) * u1 + d(:,1) * u2 )
         stiff(:, 5, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( (-2.0/3.0)*d(:,5) * u1 + d(:,1) * u3 )
         stiff(:, 5, 5) = con ! notice the + or -
c     
c.... K12
c     
         stiff(:, 2, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0/3.0)*d(:,6) )  
         stiff(:, 2, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0/3.0)*d(:,2) - (5.0/3.0)* (almBi)**(-2) * d_temp1 )
         stiff(:, 2, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0/3.0)*d(:,4) )
         stiff(:, 3, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    (d(:,2))
         stiff(:, 3, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,6) - (5.0/3.0)* (almBi)**(-2) * d(:,6) )
         stiff(:, 4, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,4) )
         stiff(:, 4, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( 0.0 - (5.0/3.0)* (almBi)**(-2) * d(:,5) )    
         stiff(:, 4, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,6) )
          stiff(:, 5, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0/3.0)*d(:,6) * u1 + d(:,2) * u2 + d(:,4) * u3 )
          stiff(:, 5, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0/3.0)*d(:,2) * u1 + d(:,6) * u2 - (5.0/3.0)* (almBi)**(-2) *
     &                    ( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3) )
         stiff(:, 5, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0/3.0)*d(:,4) * u1 + d(:,6) * u3 )
c         stiff(:, 5, 10) = con ! notice the + or -
c     
c.... K13
c     
         stiff(:, 2,12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0 /3.0)*d(:,5) )  
         stiff(:, 2,13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0 /3.0)*d(:,4) )
         stiff(:, 2,14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0 /3.0)*d(:,3) - (5.0/3.0)* (almBi)**(-2) * d_temp1)
         stiff(:, 3,12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    (d(:,4))
         stiff(:, 3,13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,5) )
         stiff(:, 3,14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( - (5.0/3.0)* (almBi)**(-2) * d(:,6) )   
         stiff(:, 4,12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,3) )  
         stiff(:, 4,14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,5)- (5.0/3.0)* (almBi)**(-2) * d(:,5) )
         stiff(:, 5,12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0 /3.0)*d(:,5) * u1 + d(:,4) * u2 + d(:,3) * u3 )
         stiff(:, 5,13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0 /3.0)*d(:,4) * u1 + d(:,5) * u2  )
         stiff(:, 5,14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( (-2.0 /3.0)*d(:,3) * u1 + d(:,5) * u3 - (5.0/3.0)* (almBi)**(-2) *
     &                    ( d_temp1 * u1 + d(:,6) * u2 + d(:,5) * u3))
c         stiff(:, 5,15) = con ! notice the + or -
c     
c.... K21
c     
         stiff(:, 7, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    (d(:,6)- (5.0/3.0)* (almBi)**(-2) * d(:,6))
         stiff(:, 7, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    (d(:,1))
         stiff(:, 8, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ((-2.0/3.0)*d(:,1) - (5.0/3.0)* (almBi)**(-2) * d_temp2 )  
         stiff(:, 8, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     (4.0/3.0)*d(:,6)
         stiff(:, 8, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     (-2.0/3.0)*d(:,5)
         stiff(:, 9, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    (- (5.0/3.0)* (almBi)**(-2) * d(:,4))
         stiff(:, 9, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    d(:,5)
         stiff(:, 9, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    d(:,6)
         stiff(:, 10, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,6) * u1 + (-2.0/3.0)*d(:,1) * u2 - (5.0/3.0)* (almBi)**(-2) *
     &                    ( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3) )
         stiff(:, 10, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    (d(:,1) * u1 + (4.0/3.0)*d(:,6) * u2 +d(:,5) * u3)
         stiff(:, 10, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ((-2.0/3.0)*d(:,5) * u2 +d(:,6) * u3)
c         stiff(:, 10, 5) = con ! notice the + or -
c     
c.... K22
c     
         stiff(:, 7, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,2) )
         stiff(:, 7, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,6) -(5.0/3.0)* (almBi)**(-2) * d(:,6) )
         stiff(:, 8, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( -(2.0/3.0)*d(:,6) )
         stiff(:, 8, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0/3.0)*d(:,2) -(5.0/3.0)* (almBi)**(-2) * d_temp2 )
         stiff(:, 8, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0/3.0)*d(:,4) )   
         stiff(:, 9, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
         stiff(:, 9, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     ( d(:,2) )
         stiff(:, 10, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,2) * u1 + (-2.0/3.0)*d(:,6) * u2 )
         stiff(:, 10, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    (d(:,6) * u1 + (4.0/3.0)*d(:,2) * u2 +d(:,4) * u3 -
     &                    (5.0/3.0)* (almBi)**(-2) *( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3))
         stiff(:, 10, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ((-2.0/3.0)*d(:,4) * u2 +d(:,2) * u3)
         stiff(:, 10, 10) = con ! notice the + or -
c     
c.... K23
c     
         stiff(:, 7, 12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,4) )
         stiff(:, 7, 13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,5) )
         stiff(:, 7, 14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( 0.0 -(5.0/3.0)* (almBi)**(-2) * d(:,6) )
         stiff(:, 8, 12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( -(2.0/3.0)*d(:,5) )
         stiff(:, 8, 13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0/3.0)*d(:,4) )
         stiff(:, 8, 14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0/3.0)*d(:,3) -(5.0/3.0)* (almBi)**(-2) * d_temp2 )   
         stiff(:, 9, 13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     ( d(:,3) )
         stiff(:, 9, 14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
         stiff(:, 10, 12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,4) * u1 + (-2.0/3.0)*d(:,5) * u2 )
         stiff(:, 10, 13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    (d(:,5) * u1 + (4.0/3.0)*d(:,4) * u2 +d(:,3) * u3 )
         stiff(:, 10, 14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ((-2.0/3.0)*d(:,3) * u2 +d(:,4) * u3 -
     &                    (5.0/3.0)* (almBi)**(-2) *( d(:,6) * u1 + d_temp2 * u2 + d(:,4) * u3))
c         stiff(:, 10, 15) = con ! notice the + or -
c     
c.... K31
c
         stiff(:, 12, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,5) - (5.0/3.0)* (almBi)**(-2) *d(:,5) )
         stiff(:, 12, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,1) )
         stiff(:, 13, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( 0.0 - (5.0/3.0)* (almBi)**(-2) *d(:,4))
         stiff(:, 13, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,5) )
         stiff(:, 13, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,6) )
         stiff(:, 14, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (-2.0/3.0)*d(:,1) - (5.0/3.0)* (almBi)**(-2) *d_temp3 )   
         stiff(:, 14, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     ( (-2.0/3.0)*d(:,6) )
         stiff(:, 14, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     ( (4.0/3.0)*d(:,5) )
         stiff(:, 15, 2) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,5) * u1 + (-2.0/3.0)*d(:,1) * u3 -
     &                    (5.0/3.0)* (almBi)**(-2) *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3) )
         stiff(:, 15, 3) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,5) * u2 +(-2.0/3.0)*d(:,6) * u3 )     
         stiff(:, 15, 4) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    (d(:,1) * u1 + (4.0/3.0)*d(:,5) * u3 +d(:,6) * u2 )
c         stiff(:, 15, 5) = con ! notice the + or -     
c
c.... K32
c     
         stiff(:, 12, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,4) )
         stiff(:, 12, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( 0.0 - (5.0/3.0)* (almBi)**(-2) *d(:,5) )
         stiff(:, 12, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( d(:,6) )
         stiff(:, 13, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     ( d(:,4) - (5.0/3.0)* (almBi)**(-2) *d(:,4) )
         stiff(:, 13, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                     ( d(:,2) )
         stiff(:, 14, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( -(2.0/3.0)*d(:,6) )
         stiff(:, 14, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( -(2.0/3.0)*d(:,2) - (5.0/3.0)* (almBi)**(-2) *d_temp3)
         stiff(:, 14, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*
     &                    ( (4.0/3.0)*d(:,4) )   
         stiff(:, 15, 7) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,4) * u1 + (-2.0/3.0)*d(:,6) * u3 )
         stiff(:, 15, 8) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,4) * u2 + (-2.0/3.0)*d(:,2) * u3 -
     &                    (5.0/3.0)* (almBi)**(-2) *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3))
         stiff(:, 15, 9) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,6) * u1 + d(:,2) * u2 + (4.0/3.0)*d(:,4) * u3)
c         stiff(:, 15, 10) = con ! notice the + or -
c     
c.... K33
c     
         stiff(:, 12, 12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,3) )
         stiff(:, 12, 14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( d(:,5)-(5.0/3.0)* (almBi)**(-2) * d(:,5)  )
         stiff(:, 13, 13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                     ( d(:,3) )
         stiff(:, 13, 14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                     ( d(:,4) -(5.0/3.0)* (almBi)**(-2) * d(:,4))
         stiff(:, 14, 12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( -(2.0/3.0)*d(:,5) )
         stiff(:, 14, 13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( -(2.0/3.0)*d(:,4) )
         stiff(:, 14, 14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)* 
     &                    ( (4.0/3.0)*d(:,3) -(5.0/3.0)* (almBi)**(-2) * d_temp3 )   
         stiff(:, 15, 12) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*  
     &                    ( d(:,3) * u1 + (-2.0/3.0)*d(:,5) * u3 )
         stiff(:, 15, 13) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*  
     &                    ( d(:,3) * u2 + (-2.0/3.0)*d(:,4) * u3 )
         stiff(:, 15, 14) = bdy_shearMod * (det_d)**(-5.0/6.0) * (gamBi*Delt(1) * alfBi)*(almBi)**(-7.0/6.0)*  
     &                    (d(:,5) * u1 + d(:,4) * u2 + (4.0/3.0)*d(:,3) * u3 -
     &                    (5.0/3.0)* (almBi)**(-2) *( d(:,5) * u1 + d(:,4) * u2 + d_temp3 * u3))
         stiff(:, 15, 15) = con ! notice the + or -
         stiff(:,:,:) = -stiff(:,:,:)
c         stiff = zero 

c..........> end assembly K_ij for solid <.........
c..........> assembly element stiffness matrix <.........   
c..loop through columns (nodes j)
c
          do ia = 1, nshlb
             j =lnode(ia)
             j0 = nflow * (j - 1)
c
c.... set up factors
c
             Wshg1 = WdetJb * shgb(:,j,1)
             Wshg2 = WdetJb * shgb(:,j,2)
             Wshg3 = WdetJb * shgb(:,j,3)
c
c.... loop through d.o.f.'s
c
             do jdof = 1, nflow
                do idof = 1, nflow
                   idof2 = idof  + nflow
                   jdof2 = jdof  + nflow

                   idof3 = idof2 + nflow
                   jdof3 = jdof2 + nflow
c
c.... calculate weighted stiffness matrix (first part)
c
                   stif1(:,idof,jdof) = Wshg1 * stiff(:,idof,jdof)
     &                                + Wshg2 * stiff(:,idof,jdof2)
     &                                + Wshg3 * stiff(:,idof,jdof3)
                   stif2(:,idof,jdof) = Wshg1 * stiff(:,idof2,jdof)
     &                                + Wshg2 * stiff(:,idof2,jdof2)
     &                                + Wshg3 * stiff(:,idof2,jdof3)
                   stif3(:,idof,jdof) = Wshg1 * stiff(:,idof3,jdof)
     &                                + Wshg2 * stiff(:,idof3,jdof2)
     &                                + Wshg3 * stiff(:,idof3,jdof3)
                enddo
             enddo
c
c.... loop through rows (nodes i)
c
             do ib = 1, nshlb
                i = lnode(ib)
                i0 = nflow * (i - 1)
c
c.... add contribution of stiffness to EGmass
c
                do jdof = 1, nflow
                   EGmass_bs(:,i0+1,j0+jdof) = EGmass_bs(:,i0+1,j0+jdof) 
     &                  + shape(:,i) * stif1(:,1,jdof) * bnorm(:,1)
     &                  + shape(:,i) * stif2(:,1,jdof) * bnorm(:,2)
     &                  + shape(:,i) * stif3(:,1,jdof) * bnorm(:,3)
                   EGmass_bs(:,i0+2,j0+jdof) = EGmass_bs(:,i0+2,j0+jdof) 
     &                  + shape(:,i) * stif1(:,2,jdof) * bnorm(:,1)
     &                  + shape(:,i) * stif2(:,2,jdof) * bnorm(:,2)
     &                  + shape(:,i) * stif3(:,2,jdof) * bnorm(:,3)
                   EGmass_bs(:,i0+3,j0+jdof) = EGmass_bs(:,i0+3,j0+jdof) 
     &                  + shape(:,i) * stif1(:,3,jdof) * bnorm(:,1)
     &                  + shape(:,i) * stif2(:,3,jdof) * bnorm(:,2)
     &                  + shape(:,i) * stif3(:,3,jdof) * bnorm(:,3)
                   EGmass_bs(:,i0+4,j0+jdof) = EGmass_bs(:,i0+4,j0+jdof) 
     &                  + shape(:,i) * stif1(:,4,jdof) * bnorm(:,1)
     &                  + shape(:,i) * stif2(:,4,jdof) * bnorm(:,2)
     &                  + shape(:,i) * stif3(:,4,jdof) * bnorm(:,3)
                   EGmass_bs(:,i0+5,j0+jdof) = EGmass_bs(:,i0+5,j0+jdof) 
     &                  + shape(:,i) * stif1(:,5,jdof) * bnorm(:,1)
     &                  + shape(:,i) * stif2(:,5,jdof) * bnorm(:,2)
     &                  + shape(:,i) * stif3(:,5,jdof) * bnorm(:,3)
                enddo
c
c.... end loop on rows
c
             enddo
c
c.... end loop on columns
c
          enddo
c
c.........> end assembly the element stiffness matrix <.....
     

       end subroutine getLHS_bs
