      subroutine getthm_solid_1 (rho,    ei
     &,                           p,      T,     npro, mater
     &,                           h,      cv,    cp
     &,                           alphaP, betaT, bulkMod, shearMod
     &,                           Ja_def)
        use number_def_m
        use matdat_def_m
c
        integer, intent(in) :: npro, mater
        real*8, dimension(npro), intent(in)  :: p, T, Ja_def        
        real*8, dimension(npro), intent(out) :: rho,ei,h
        real*8, dimension(npro), intent(out) :: cv,cp,alphaP,betaT
        real*8, dimension(npro), intent(out) :: bulkMod, shearMod
c
c
        real*8 :: rho_ref_s, p_ref_s, T_ref_s, alpha_P_s, cv_s
        real*8 :: bulkMod_s, shearMod_s 
c
        rho_ref_s = mat_prop(mater,iprop_solid_1_rho_ref,1)
        p_ref_s   = mat_prop(mater,iprop_solid_1_p_ref,  1)
        T_ref_s   = mat_prop(mater,iprop_solid_1_T_ref,  1)
        cv_s     = mat_prop(mater,iprop_solid_1_cv,     1)
        alpha_P_s = mat_prop(mater,iprop_solid_1_alphaP, 1)
        bulkMod_s = mat_prop(mater,iprop_solid_1_bulkMod, 1)
        shearMod_s = mat_prop(mater,iprop_solid_1_shearMod, 1)
!        beta_T  = mat_prop(mater,iprop_solid_1_betaT,  1)
c
        alphaP = alpha_P_s
        betaT  = one /(bulkMod_s * Ja_def) ! double check here
        rho = rho_ref_s * (one - alphaP*(T - T_ref_s) 
     &                    + betaT*(p - p_ref_s))
c        ei  = ( cv_s - P * alpha_P_s/rho)* T 
c     &        + (betaT * P - alphaP * T)/rho * p
        ei  = cv_s * T
        h   = ei + p/rho
        cv  = cv_s
        cp  = cv_s
        bulkMod = bulkMod_s
        shearMod = shearMod_s
c        c =  sqrt(one/(rho_ref*betaT))
C         c =  sqrt(one/(rho*betaT))
C         gamb = zero
c
      end subroutine getthm_solid_1
