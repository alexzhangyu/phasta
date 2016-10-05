      subroutine fillelmb( elmb1, elmb2, iel,      npro_s,
     &                     lcsyst_s,     ngauss_s, blk_b)
c....fill the elm-wise solid array
c       
      
c      implicit none
      use pointer_data
      include "common.h"
c
      real*8 elmb1(numel,3) !elm-wise solid array
      real*8 elmb2(numel,3) !elm-wise solid array
      real*8 blk_b(npro_s,ngauss_s,6)
      real*8 sum_temp(npro_s)
      integer::iel_s, npro_s 
      integer::lcsyst_s, ngauss_s
c
      sum_temp = zero !initialization
c
      ilast = iel+npro_s-1
c    
c.....loop over all quadrature point
         quad_loop: do intp = 1, ngauss_s
c
                   elmb1(iel:ilast,:) = elmb1(iel:ilast,:)
     &                                 + blk_b(:,intp,1:3)* Qwt(lcsyst_s,intp)
                   elmb2(iel:ilast,:) = elmb2(iel:ilast,:)
     &                                 + blk_b(:,intp,4:6)* Qwt(lcsyst_s,intp)
                   sum_temp(:) = sum_temp(:) + Qwt(lcsyst_s,intp)
c..
        enddo quad_loop 
c.. Normalize the elm-wise field
        do ith = 1,3
        elmb1(iel:ilast,ith) = elmb1(iel:ilast,ith)/sum_temp(:)
        elmb2(iel:ilast,ith) = elmb2(iel:ilast,ith)/sum_temp(:)
        enddo 
c
c
      return
      end subroutine fillelmb
