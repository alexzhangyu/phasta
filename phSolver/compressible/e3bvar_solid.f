      subroutine e3bvar_solid(bdy_Ja_def_,         bdy_det_baf_, 
     &                        g1yi_,     g2yi_,    g3yi_,
     &                        npro_,     nsd_,     almBi_,   alfBi_,
     &                        gamBi_,    intp_,    Delt_)
c
        use solid_func_m
c
        implicit none
c
        real*8, dimension(npro_,nsd_), target, intent(in) :: g1yi_,g2yi_,g3yi_
        real*8, dimension(npro_), target :: bdy_Ja_def_
        real*8, dimension(npro_), target :: bdy_det_baf_
        real*8, intent(in) :: Delt_
        integer, intent(in) :: npro_, nsd_
        integer, intent(in) :: intp_
        real*8, intent(in) :: almBi_, alfBi_, gamBi_
        real*8, dimension(npro_,6) :: bdy_d
c
        npro = npro_
        nsd = nsd_
        Delt = Delt_
        intp_s = intp_
        almBi = almBi_
        alfBi = alfBi_
        gamBi = gamBi_  
c
        Ja_def => bdy_Ja_def_
        det_baf => bdy_det_baf_
        dudx => g1yi_(:,2:4)
        dudy => g2yi_(:,2:4)
        dudz => g3yi_(:,2:4)
c
        allocate(AS(npro,6,6))
        call calc_as_matrix
        bdy_d(:,:) = almBi * bdy_b(iblk_solid)%p(:,intp_s,:)
     &+      alfBi * Delt * (almBi - gamBi) 
     &*      bdy_b_dot(iblk_solid)%p(:,intp_s,:)
        call setB_af(bdy_d, AS, bdy_b_af(iblk_solid)%p(:,intp_s,:))
        call get_det(bdy_b_af( iblk_solid)%p(:,intp_s,:),det_baf )
        Ja_def= (det_baf)**0.5
c
        deallocate(AS)         
c
      end subroutine e3bvar_solid
